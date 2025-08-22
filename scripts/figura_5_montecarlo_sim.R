# 0) librerias
library(tidyr) 
library(dplyr)
library(purrr)
 
## ★ paralelismo
library(future)
library(future.apply)

## usa núcleos-físicos 
workers_use <- max(1, parallel::detectCores(logical = FALSE) - 1)
plan(multisession, workers = workers_use)       

# 1) Cargar todas las funciones base
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
source("R/montecarlo_helpers.R")
source("R/montecarlo_function.R")
source("R/oracles.R")

# 2) Cargar la betting function precomputada (Generada por calibrar_kde.R)
kde_bf_fixed <- readRDS("data/kde_bf_fixed.rds")

make_ncm_LR <- function(mu_shift) {
  function(xi, training_set, ...) {
    Non_conformity_LNR(xi, training_set, mu_r = mu_shift)
  }
}

# 3) Definir escenarios y parámetros
theta_vals <- c(100, 200)
mu1_vals   <- c(1, 1.5, 2)
h_vals     <- seq(1, 6, 0.5)
n_sim      <- 1000
m          <- 200
k_par      <- 7

# 4-5) Bucle principal ---------------------------------------------------------
t_total <- Sys.time()
all_results  <- list()
all_taus_tbl <- list()

for (theta_s in theta_vals) {
  for (mu1_s in mu1_vals) {
    t_scen <- Sys.time()
    # ------ NCM y BFs para este salto de media --------------------------------
    ncm_tbl <- tibble(
      ncm_fun  = list(Non_conformity_KNN,
                      make_ncm_LR(mu1_s), Non_conformity_MAD),
      ncm_lbl  = c("KNN", "LR", "MAD"),
      needs_k  = c(TRUE,  FALSE, FALSE)
    )
    
    bet_tbl <- tibble(
      bet_fun   = list(Constant_BF, Mixture_BF, kde_bf_fixed, histogram_betting_function),
      bet_lbl   = c("Constant BF", "Mixture BF", "Precomputed KDE BF", "Histogram BF"),
      params_bf = list(
        list(),                    # Constant no necesita nada
        list(),                    # Mixture tampoco
        list(),                    # Precomputed KDE tampoco
        list(num_bins = 20)        # Histogram sí
      )
    )
    
    icm_res <- expand_grid(ncm_tbl, bet_tbl) |>
      pmap(function(ncm_fun, ncm_lbl, needs_k,
                        bet_fun, bet_lbl, params_bf) {
        k_val <- if (needs_k) k_par else NULL
        
        # ---- Método ICM clásico
        out_icm <- montecarlo_ICM(
          n_sim        = n_sim,
          h_vals       = h_vals,
          theta_stream = theta_s,
          mu1          = mu1_s,
          m            = m,
          ncm_fun      = ncm_fun,
          bet_fun      = bet_fun,
          k            = k_val,
          params_bf    = params_bf,
          n_stream     = 1000 
        )
        sum_icm <- out_icm$summary |>
          mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM"))
        tau_icm <- out_icm$taus |>
          mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM"),
                 ncm = ncm_lbl, bf = bet_lbl)
        
        # ---- Método ICM Cautious
        out_icm_cbf <- montecarlo_ICM_CBF(
          n_sim        = n_sim,
          h_vals       = h_vals,
          theta_stream = theta_s,
          mu1          = mu1_s,
          m            = m,
          ncm_fun      = ncm_fun,
          bet_fun      = bet_fun,
          k            = k_val,
          params_bf    = params_bf,
          n_stream     = 1000
        )
        sum_cbf <- out_icm_cbf$summary |>
          mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM CBF"))
        tau_cbf <- out_icm_cbf$taus |>
          mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM CBF"),
                 ncm = ncm_lbl, bf = bet_lbl)
        
        list(summary = bind_rows(sum_icm, sum_cbf),
             taus    = bind_rows(tau_icm,  tau_cbf))
      })
    
    icm_res_summary <- bind_rows(purrr::map(icm_res, "summary"))
    icm_res_taus    <- bind_rows(purrr::map(icm_res, "taus"))
    
    # ---- Oráculos (sin cambios) ---------------------------------------------
    df_CUSUM <- montecarlo_oraculo(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      n_stream     = 1000,
      detector_fn  = cusum_oracle,
      detector_label = "CUSUM Oracle"
    )
    
    df_SR <- montecarlo_oraculo(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      n_stream     = 1000,
      detector_fn  = sr_oracle,
      detector_label = "S-R Oracle"
    )
    
    df_POST <- montecarlo_oraculo(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      n_stream     = 1000,
      detector_fn  = posterior_oracle,
      detector_label = "Posterior Oracle",
      p = 1 / 100
    )
    or_summary <- bind_rows(df_CUSUM$summary, df_SR$summary, df_POST$summary)
    or_taus    <- bind_rows(
      mutate(df_CUSUM$taus, Method = "CUSUM Oracle"),
      mutate(df_SR$taus,    Method = "S-R Oracle"),
      mutate(df_POST$taus,  Method = "Posterior Oracle")
    )
    
    combined_summary <- bind_rows(icm_res_summary, or_summary) %>%
      mutate(
        theta_stream = factor(theta_stream, levels = theta_vals),
        mu1          = factor(mu1,        levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    combined_taus <- icm_res_taus %>%
      bind_rows(or_taus) %>%
      mutate(
        theta_stream = factor(theta_s, levels = theta_vals),
        mu1          = factor(mu1_s,   levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    # Guardar por escenario
    all_results[[paste0("θ=", theta_s, "_μ1=", mu1_s)]]  <- combined_summary
    all_taus_tbl[[paste0("θ=", theta_s, "_μ1=", mu1_s)]] <- combined_taus
    secs <- as.numeric(difftime(Sys.time(), t_scen, units = "secs"))
    cat(sprintf("⏱️ Escenario θ=%s, μ1=%s listo en %.1f s (%.1f min)\n",
                theta_s, mu1_s, secs, secs/60))
  }
}

df_all_methods <- bind_rows(all_results,  .id = "scenario_key")
df_all_taus    <- bind_rows(all_taus_tbl, .id = "scenario_key")

saveRDS(df_all_methods, file = "data/prueba_fast.rds")
saveRDS(df_all_taus,    file = "data/prueba_fast_taus.rds") 
cat("✅ Simulación con KNN y LR (parametrizado) completada.\n")

secs_total <- as.numeric(difftime(Sys.time(), t_total, units = "secs"))
cat(sprintf("⏰ Tiempo total: %.1f s (%.1f min)\n", secs_total, secs_total/60))

plan(sequential)