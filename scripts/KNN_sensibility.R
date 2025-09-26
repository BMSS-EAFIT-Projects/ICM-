# --- 0) Librerías y paralelismo ---------------------------------------------
library(tidyr); library(dplyr); library(purrr); library(ggplot2)
library(future); library(future.apply); library(forcats)

workers_use <- max(1, parallel::detectCores(logical = FALSE) - 1)
plan(multisession, workers = workers_use)

# --- 1) Fuentes base ---------------------------------------------------------
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
source("R/montecarlo_helpers.R")
source("R/montecarlo_function.R")

# --- 2) Cargar KDE precomputado ---------------------------------------------
kde_bf_fixed <- readRDS("data/kde_bf_fixed.rds")

# --- 3) Escenarios y parámetros ---------------------------------------------
theta_vals <- c(100, 200)
mu1_vals   <- c(1, 1.5, 2)

# Umbrales a testear
h_vals     <- seq(1, 9, 0.5)

# Nº simulaciones por combinación
n_sim      <- 1000

# Longitud de entrenamiento y stream
m          <- 200
n_stream   <- 1000

# *** Barrido de K para KNN ***
K_grid     <- c(1, 7, 25, 50, 100, 150, 200)

# Betting functions (sin CBF)
bet_tbl <- tibble(
  bet_fun   = list(Constant_BF, Mixture_BF, kde_bf_fixed, histogram_betting_function),
  bet_lbl   = c("Constant BF", "Mixture BF", "Precomputed KDE BF", "Histogram BF"),
  params_bf = list(
    list(),                    # Constant
    list(),                    # Mixture
    list(),                    # KDE precomputada
    list(num_bins = 2)        # Histograma
  )
)

# --- 4) Simulación principal: NCM = KNN (con K_grid), MAD, IQR, LR ----------
t_total <- Sys.time()
all_results  <- list()
all_taus_tbl <- list()

for (theta_s in theta_vals) {
  for (mu1_s in mu1_vals) {
    t_scen <- Sys.time()
    message(sprintf("▶ Escenario θ=%s, μ1=%s ...", theta_s, mu1_s))
    
    # NCMs: KNN (con K), MAD, IQR, LR(mu1 dependiente)
    make_ncm_LR <- function(mu_shift) {
      function(xi, training_set, ...) {
        Non_conformity_LNR(xi, training_set, mu_r = mu_shift)
      }
    }
    
    ncm_tbl <- tibble(
      ncm_fun = list(
        Non_conformity_KNN,   # KNN
        Non_conformity_MAD,   # MAD
        Non_conformity_IQR,   # IQR
        make_ncm_LR(mu1_s)    # LR parametrizada con el salto actual
      ),
      ncm_lbl = c("KNN", "MAD", "IQR", "LR"),
      needs_k = c(TRUE, FALSE, FALSE, FALSE)
    )
    
    # Para cada NCM y cada BF, correr ICM clásico
    res_list <- expand_grid(ncm_tbl, bet_tbl) |>
      pmap(function(ncm_fun, ncm_lbl, needs_k, bet_fun, bet_lbl, params_bf) {
        
        # Si la NCM requiere K, barrer K_grid; si no, una sola corrida con k=NULL
        K_vec <- if (needs_k) K_grid else NA_integer_
        
        map(K_vec, function(K_try) {
          out_icm <- montecarlo_ICM(
            n_sim        = n_sim,
            h_vals       = h_vals,
            theta_stream = theta_s,
            mu1          = mu1_s,
            m            = m,
            ncm_fun      = ncm_fun,
            bet_fun      = bet_fun,
            k            = if (needs_k) K_try else NULL,
            params_bf    = params_bf,
            n_stream     = n_stream
          )
          
          sum_icm <- out_icm$summary |>
            mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM"),
                   ncm = ncm_lbl, bf = bet_lbl, K = if (needs_k) K_try else NA_integer_)
          tau_icm <- out_icm$taus |>
            mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM"),
                   ncm = ncm_lbl, bf = bet_lbl, K = if (needs_k) K_try else NA_integer_)
          
          list(summary = sum_icm, taus = tau_icm)
        }) |> 
          {\(.x) list(summary = bind_rows(map(.x, "summary")),
                      taus    = bind_rows(map(.x, "taus")))}()
      })
    
    scen_summary <- bind_rows(purrr::map(res_list, "summary")) |>
      mutate(
        theta_stream = factor(theta_s, levels = theta_vals),
        mu1          = factor(mu1_s,   levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    scen_taus <- bind_rows(purrr::map(res_list, "taus")) |>
      mutate(
        theta_stream = factor(theta_s, levels = theta_vals),
        mu1          = factor(mu1_s,   levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    key <- paste0("θ=", theta_s, "_μ1=", mu1_s)
    all_results[[key]]  <- scen_summary
    all_taus_tbl[[key]] <- scen_taus
    
    secs <- as.numeric(difftime(Sys.time(), t_scen, units = "secs"))
    message(sprintf("✓ Escenario θ=%s, μ1=%s listo en %.1f s (%.1f min)",
                    theta_s, mu1_s, secs, secs/60))
  }
}

df_ncm_compare_summary <- bind_rows(all_results,  .id = "scenario_key")
df_ncm_compare_taus    <- bind_rows(all_taus_tbl, .id = "scenario_key")

# Guardar resultados crudos
dir.create("data", showWarnings = FALSE, recursive = TRUE)
saveRDS(df_ncm_compare_summary, file = "data/ncm_compare_KNN_MAD_IQR_LR_summary.rds")
saveRDS(df_ncm_compare_taus,    file = "data/ncm_compare_KNN_MAD_IQR_LR_taus.rds")

secs_total <- as.numeric(difftime(Sys.time(), t_total, units = "secs"))
message(sprintf("⏱ Tiempo total: %.1f s (%.1f min)", secs_total, secs_total/60))

plan(sequential)
