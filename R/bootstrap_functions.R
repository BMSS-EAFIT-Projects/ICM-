# 1) Cargar todas las funciones base
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")

## Ejemplo de training set in-control (Gauss N(0,1))
training_set <- rnorm(500)

## Generador in-control para el stream
r_incontrol <- function(n) {
  rnorm(n, mean = 0, sd = 1)
}

N_max <- 10000      # longitud máxima del stream (horizonte)

## ───── 1. Rutina de simulación in-control ────────────────────────
simulate_stream_incontrol <- function(h,
                                      training_set,
                                      N_max,
                                      ncm,
                                      bf,
                                      params_bf = list()) {
  
  stream <- r_incontrol(N_max)
  out <- ICM(training_set      = training_set,
             stream_data       = stream,
             non_conformity_measure = ncm,
             betting_function  = bf,
             th                = h,
             params_bf         = params_bf)
  
  as.integer(!is.na(out$change_point))  # 1 si hubo alarma
}

## ───── 2. Estimador de P_FA dado h ───────────────────────────────
estimate_P_FA <- function(h, M, seed_i,
                          training_set, N_max,
                          ncm, bf, params_bf = list()) {
  
  set.seed(seed_i)                       # semilla externa reproducible
  alarms <- replicate(
    M,
    simulate_stream_incontrol(h, training_set, N_max,
                              ncm, bf, params_bf)
  )
  mean(alarms)
}

## Calibrador de umbral

calibrate_h <- function(alpha, h_grid, M, seed_cal,
                        training_set, N_max,
                        ncm, bf, params_bf = list()) {
  
  set.seed(seed_cal)
  P_FA_vec <- sapply(
    h_grid,
    estimate_P_FA,
    M = M, seed_i = sample.int(1e9, 1),
    training_set = training_set,
    N_max = N_max,
    ncm = ncm,
    bf = bf,
    params_bf = params_bf
  )
  idx <- which(P_FA_vec <= alpha)
  h_hat <- if (length(idx)) h_grid[min(idx)] else max(h_grid)
  
  list(h_hat = h_hat,
       P_FA_est = P_FA_vec[match(h_hat, h_grid)],
       trace = data.frame(h = h_grid, P_FA = P_FA_vec))
}


## bootstrap 

bootstrap_calibration <- function(B = 100,
                                  alpha = 0.05,
                                  h_grid = seq(1, 6, 0.5),
                                  M = 5000,
                                  training_set,
                                  N_max,
                                  ncm, bf, params_bf = list(),
                                  seed_root = 20250710) {
  
  set.seed(seed_root)
  seeds <- sample.int(.Machine$integer.max, B)
  
  res <- lapply(
    seeds,
    calibrate_h,
    alpha = alpha,
    h_grid = h_grid,
    M = M,
    training_set = training_set,
    N_max = N_max,
    ncm = ncm,
    bf = bf,
    params_bf = params_bf
  )
  
  h_boot   <- vapply(res, `[[`, numeric(1), "h_hat")
  PFA_boot <- vapply(res, `[[`, numeric(1), "P_FA_est")
  
  list(h_boot = h_boot,
       PFA_boot = PFA_boot,
       alpha = alpha,
       M = M,
       h_grid = h_grid,
       seed_root = seed_root)
}