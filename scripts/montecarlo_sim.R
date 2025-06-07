# 0) librerias
library(tidyr) 
library(dplyr)
library(purrr)
 
# 1) Cargar todas las funciones base
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
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
n_sim      <- 200
m          <- 200
k_par      <- 7

# 4-5) Bucle principal ---------------------------------------------------------
all_results <- list()

for (theta_s in theta_vals) {
  for (mu1_s in mu1_vals) {
    
    # ------ NCM y BFs para este salto de media --------------------------------
    ncm_tbl <- tibble(
      ncm_fun  = list(Non_conformity_KNN,
                      make_ncm_LR(mu1_s)),
      ncm_lbl  = c("KNN", "LR"),
      needs_k  = c(TRUE,  FALSE)
    )
    
    bet_tbl <- tibble(
      bet_fun = list(Constant_BF, Mixture_BF, kde_bf_fixed),
      bet_lbl = c("Constant BF", "Mixture BF", "Precomputed KDE BF")
    )
    
    icm_res <- expand_grid(ncm_tbl, bet_tbl) |>
      pmap_dfr(function(ncm_fun, ncm_lbl, needs_k,
                        bet_fun, bet_lbl) {
        montecarlo_ICM(
          n_sim        = n_sim,
          h_vals       = h_vals,
          theta_stream = theta_s,
          mu1          = mu1_s,
          m            = m,
          ncm_fun      = ncm_fun,
          bet_fun      = bet_fun,
          k            = if (needs_k) k_par else NULL
        ) |>
          mutate(Method = paste(ncm_lbl, bet_lbl, sep = " + "))
      })
    
    # ---- Oráculos (sin cambios) ---------------------------------------------
    df_CUSUM <- montecarlo_oraculo(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      detector_fn  = cusum_oracle,
      detector_label = "CUSUM Oracle"
    ) %>% mutate(Method = "CUSUM Oracle")
    
    df_SR <- montecarlo_oraculo(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      detector_fn  = sr_oracle,
      detector_label = "S-R Oracle"
    ) %>% mutate(Method = "S-R Oracle")
    
    df_POST <- montecarlo_oraculo(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      detector_fn  = posterior_oracle,
      detector_label = "Posterior Oracle",
      p = 1 / 100
    ) %>% mutate(Method = "Posterior Oracle")
    combined <- bind_rows(icm_res, df_CUSUM, df_SR, df_POST) %>%
      mutate(
        theta_stream = factor(theta_stream, levels = theta_vals),
        mu1          = factor(mu1,        levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    all_results[[paste0("θ=", theta_s, "_μ1=", mu1_s)]] <- combined
  }
}

df_all_methods <- bind_rows(all_results, .id = "scenario_key")
saveRDS(df_all_methods, file = "data/resultados_sim.rds")
cat("✅ Simulación con KNN y LR (parametrizado) completada.\n")