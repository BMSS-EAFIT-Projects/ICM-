# --- 0) Librerías y plan de paralelismo --------------------------------------
library(tidyr); library(dplyr); library(purrr)
library(future); library(future.apply)

workers_use <- max(1, parallel::detectCores(logical = FALSE) - 1)
plan(multisession, workers = workers_use)

# --- 1) Cargar funciones base -------------------------------------------------
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
source("R/montecarlo_helpers.R")
source("R/montecarlo_function.R")

# --- 2) Cargar BF precomputada (si aplica) -----------------------------------
kde_bf_fixed <- readRDS("data/kde_bf_fixed.rds")

# Utilidad: NCM LR parametrizada
make_ncm_LR <- function(mu_shift) {
  function(xi, training_set, ...) {
    Non_conformity_LNR(xi, training_set, mu_r = mu_shift)
  }
}

# --- 3) Definir escenarios y parámetros --------------------------------------
# Escenarios multi: vector de cambios y niveles por tramo
scenarios_multi <- tribble(
  ~scenario_id,            ~theta_vec,          ~mu_levels,
  "S1: 2 changes up/down", c(500, 1300),         c(0, 1.5, 0),
  "S2: 3 changes up-up",   c(300, 1000, 1700),    c(0, 1.0, 2.0, 3.0)
)

# Umbrales y simulaciones
h_vals   <- seq(1, 6, 0.5)
n_sim    <- 2000
m        <- 200
n_stream <- 2000
k_par    <- 7

# Ventana de detección (elige una de dos modalidades)
window_mode <- "frac"    # "abs" = ventana en puntos; "frac" = fracción del segmento
window_abs  <- Inf
window_frac <- 1     

# Reentrenos (solo para ICM_multi_adaptive)
m_retrain  <- m
guard_band <- 20

# --- 5) Bucle principal: NCM × BF × Método × Escenario -----------------------
t_total <- Sys.time()
all_summary    <- list()
all_perchange  <- list()
all_alarms     <- list()

for (s in seq_len(nrow(scenarios_multi))) {
  sc <- scenarios_multi[s,]
  
  mu_levels_sc <- sc$mu_levels[[1]]
  jumps <- abs(diff(mu_levels_sc))
  delta_sc <- if (length(jumps)) median(jumps) else 1.0
  
  # --- NCMs: usa LR con Δ del escenario ------------------------------
  ncm_tbl <- tibble(
    ncm_fun  = list(
      Non_conformity_KNN,
      make_ncm_LR(delta_sc),
      Non_conformity_MAD,
      Non_conformity_IQR
    ),
    ncm_lbl  = c("KNN", "LR", "MAD", "IQR"),
    needs_k  = c(TRUE, FALSE, FALSE, FALSE)
  )
  
  # --- BF (puede quedarse igual) -------------------------------------
  bet_tbl <- tibble(
    bet_fun   = list(Constant_BF, Mixture_BF, kde_bf_fixed, histogram_betting_function),
    bet_lbl   = c("Constant BF", "Mixture BF", "Precomputed KDE BF", "Histogram BF"),
    params_bf = list(list(), list(), list(), list(num_bins = 2))
  )
  
  
  
  
  icm_multi_res <- expand_grid(ncm_tbl, bet_tbl) |>
    pmap(function(ncm_fun, ncm_lbl, needs_k, bet_fun, bet_lbl, params_bf) {
      params_bf$bf_name  <- bet_lbl
      params_bf$ncm_name <- ncm_lbl
      k_val <- if (needs_k) k_par else NULL
      
      # Método A: ICM_multi (sin reentrenos)
      out_m <- montecarlo_ICM_MULTI(
        n_sim      = n_sim,
        h_vals     = h_vals,
        m          = m,
        n_stream   = n_stream,
        theta_vec  = sc$theta_vec[[1]],
        mu_levels  = sc$mu_levels[[1]],
        ncm_fun    = ncm_fun,
        bet_fun    = bet_fun,
        params_bf  = params_bf,
        k          = k_val,
        window_mode = window_mode,
        window_abs  = window_abs,
        window_frac = window_frac
      )
      sum_m <- out_m$summary |>
        mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM-multi"),
               ncm = ncm_lbl, bf = bet_lbl,
               scenario_id = sc$scenario_id)
      per_m <- out_m$per_change |>
        mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM-multi"),
               ncm = ncm_lbl, bf = bet_lbl,
               scenario_id = sc$scenario_id)
      alm_m <- out_m$alarms |>
        mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM-multi"),
               ncm = ncm_lbl, bf = bet_lbl,
               scenario_id = sc$scenario_id)
      
      # Método B: ICM_multi_adaptive (con reentrenos)
      out_a <- montecarlo_ICM_MULTI_ADAPTIVE(
        n_sim      = n_sim,
        h_vals     = h_vals,
        m          = m,
        n_stream   = n_stream,
        theta_vec  = sc$theta_vec[[1]],
        mu_levels  = sc$mu_levels[[1]],
        ncm_fun    = ncm_fun,
        bet_fun    = bet_fun,
        params_bf  = params_bf,
        k          = k_val,
        m_retrain  = m_retrain,
        guard_band = guard_band,
        window_mode = window_mode,
        window_abs  = window_abs,
        window_frac = window_frac
      )
      sum_a <- out_a$summary |>
        mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM-multi-adaptive"),
               ncm = ncm_lbl, bf = bet_lbl,
               scenario_id = sc$scenario_id)
      per_a <- out_a$per_change |>
        mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM-multi-adaptive"),
               ncm = ncm_lbl, bf = bet_lbl,
               scenario_id = sc$scenario_id)
      alm_a <- out_a$alarms |>
        mutate(Method = paste(ncm_lbl, bet_lbl, "+ ICM-multi-adaptive"),
               ncm = ncm_lbl, bf = bet_lbl,
               scenario_id = sc$scenario_id)
      
      list(
        summary    = bind_rows(sum_m, sum_a),
        per_change = bind_rows(per_m, per_a),
        alarms     = bind_rows(alm_m, alm_a)
      )
    })
  
  all_summary[[sc$scenario_id]]   <- bind_rows(map(icm_multi_res, "summary"))
  all_perchange[[sc$scenario_id]] <- bind_rows(map(icm_multi_res, "per_change"))
  all_alarms[[sc$scenario_id]]    <- bind_rows(map(icm_multi_res, "alarms"))
  
  secs <- as.numeric(difftime(Sys.time(), t_total, units = "secs"))
  cat(sprintf("✅ %s listo. Acumulado: %.1f s (%.1f min)\n",
              sc$scenario_id, secs, secs/60))
}

df_multi_summary   <- bind_rows(all_summary,   .id = "scenario_group")
df_multi_perchange <- bind_rows(all_perchange, .id = "scenario_group")
df_multi_alarms    <- bind_rows(all_alarms,    .id = "scenario_group")

# --- 6) Guardar resultados ----------------------------------------------------
saveRDS(df_multi_summary,   file = "data/multi_summary.rds")
saveRDS(df_multi_perchange, file = "data/multi_perchange.rds")
saveRDS(df_multi_alarms,    file = "data/multi_alarms.rds.rds")

cat("Archivos guardados en data/: multi_summary.rds, multi_perchange.rds, multi_alarms.rds\n")

plan(sequential)