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

# 1) Cargar todas las funciones base-----------------------------------
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
source("R/montecarlo_helpers.R")
source("R/montecarlo_function.R")

# 2) Cargar la betting function precomputada -------------------------------
kde_bf_fixed <- readRDS("data/kde_bf_fixed.rds")

make_ncm_LR <- function(mu_shift) {
  function(xi, training_set, ...) {
    Non_conformity_LNR(xi, training_set, mu_r = mu_shift)
  }
}

# === Parámetros de contaminación (solo shift / scale) =====================
CONTAM_CFG <- list(
  train = list(rate = 0.05, model = "scale", params = list(lambda = 6)),
  pre   = list(rate = 0.05, model = "scale", params = list(lambda = 6)),
  post  = list(rate = 0.05, model = "scale", params = list(lambda = 6))
)

# 3) Definir escenarios y parámetros---------------------------------------
theta_vals <- c(100, 200)
mu1_vals   <- c(1, 1.5, 2)
h_vals     <- seq(1, 9, 0.5)
n_sim      <- 1000
m          <- 200
k_par      <- 7

# 4-5) Bucle principal -----------------------------------------------------
t_total <- Sys.time()
all_results  <- list()
all_taus_tbl <- list()

for (theta_s in theta_vals) {
  for (mu1_s in mu1_vals) {
    t_scen <- Sys.time()
    
    # ------ NCM y BFs para este salto de media ----------------------------
    ncm_tbl <- tibble(
      ncm_fun  = list(Non_conformity_KNN,
                      make_ncm_LR(mu1_s), Non_conformity_MAD, Non_conformity_IQR),
      ncm_lbl  = c("KNN", "LR", "MAD", "IQR"),
      needs_k  = c(TRUE,  FALSE, FALSE, FALSE)
    )
    
    bet_tbl <- tibble(
      bet_fun   = list(Constant_BF, Mixture_BF, kde_bf_fixed, histogram_betting_function),
      bet_lbl   = c("Constant BF", "Mixture BF", "Precomputed KDE BF", "Histogram BF"),
      params_bf = list(
        list(),                    # Constant
        list(),                    # Mixture
        list(),                    # Precomputed KDE
        list(num_bins = 2)         # Histogram
      )
    )
    
    icm_res <- tidyr::expand_grid(ncm_tbl, bet_tbl) |>
      purrr::pmap(function(ncm_fun, ncm_lbl, needs_k,
                           bet_fun, bet_lbl, params_bf) {
        k_val <- if (needs_k) k_par else NULL
        
        # ---- Llamada ÚNICA: ICM CONTAMINADO (sin CBF ni oráculos)
        out_icm <- montecarlo_ICM_contaminado(
          n_sim        = n_sim,
          h_vals       = h_vals,
          theta_stream = theta_s,
          mu1          = mu1_s,
          m            = m,
          ncm_fun      = ncm_fun,
          bet_fun      = bet_fun,
          k            = k_val,
          params_bf    = params_bf,
          n_stream     = 1000,
          # parámetros de contaminación
          train_contam_rate   = CONTAM_CFG$train$rate,
          train_contam_model  = CONTAM_CFG$train$model,
          train_contam_params = CONTAM_CFG$train$params,
          pre_contam_rate     = CONTAM_CFG$pre$rate,
          pre_contam_model    = CONTAM_CFG$pre$model,
          pre_contam_params   = CONTAM_CFG$pre$params,
          post_contam_rate    = CONTAM_CFG$post$rate,
          post_contam_model   = CONTAM_CFG$post$model,
          post_contam_params  = CONTAM_CFG$post$params
        )
        
        sum_icm <- out_icm$summary |>
          dplyr::mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM (contaminado)"))
        tau_icm <- out_icm$taus |>
          dplyr::mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM (contaminado)"),
                        ncm = ncm_lbl, bf = bet_lbl)
        
        list(summary = sum_icm, taus = tau_icm)
      })
    
    icm_res_summary <- dplyr::bind_rows(purrr::map(icm_res, "summary"))
    icm_res_taus    <- dplyr::bind_rows(purrr::map(icm_res, "taus"))
    
    combined_summary <- icm_res_summary %>%
      dplyr::mutate(
        theta_stream = factor(theta_stream, levels = theta_vals),
        mu1          = factor(mu1,        levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    combined_taus <- icm_res_taus %>%
      dplyr::mutate(
        theta_stream = factor(theta_s, levels = theta_vals),
        mu1          = factor(mu1_s,   levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    # Guardar por escenario
    all_results[[paste0("θ=", theta_s, "_μ1=", mu1_s)]]  <- combined_summary
    all_taus_tbl[[paste0("θ=", theta_s, "_μ1=", mu1_s)]] <- combined_taus
    secs <- as.numeric(difftime(Sys.time(), t_scen, units = "secs"))
    cat(sprintf(" Escenario θ=%s, μ1=%s listo en %.1f s (%.1f min)\n",
                theta_s, mu1_s, secs, secs/60))
  }
}

df_all_methods <- dplyr::bind_rows(all_results,  .id = "scenario_key")
df_all_taus    <- dplyr::bind_rows(all_taus_tbl, .id = "scenario_key")

saveRDS(df_all_methods, file = "data/prueba_contaminada.rds")
saveRDS(df_all_taus,    file = "data/prueba_taus_contaminada.rds") 
cat(" Simulación contaminada completada.\n")

secs_total <- as.numeric(difftime(Sys.time(), t_total, units = "secs"))
cat(sprintf("Tiempo total: %.1f s (%.1f min)\n", secs_total, secs_total/60))

plan(sequential)
