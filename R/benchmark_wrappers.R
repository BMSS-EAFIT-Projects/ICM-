#Wrapper for benchmark
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")

# --- Wrapper genérico para ICM -----------------------------------------------
# z               : vector numérico completo (incluye entrenamiento + stream)
# burn_in         : tamaño de training_set inicial
# threshold       : umbral de detección (th)
# non_conformity_measure : función que calcula el non-conformity (firma xi, training_set, ...)
# betting_function       : función g(p), toma p-value y devuelve apuesta
# ...             : parámetros adicionales para NCM o BF
run_icm <- function(z, burn_in, threshold, non_conformity_measure, betting_function, ...) {
  training_set <- z[1:burn_in]
  stream_data  <- z[(burn_in + 1):length(z)]
  
  res <- ICM(
    training_set = training_set,
    stream_data  = stream_data,
    non_conformity_measure = non_conformity_measure,
    betting_function       = betting_function,
    th = threshold,
    ...
  )
  
  # Devolver sólo el índice del cambio (o Inf si no ocurrió)
  return(res$change_point)
}

# --- Wrappers especializados ------------------------------------------------

run_icm <- function(z, burn_in, threshold, non_conformity_measure, betting_function, ...) {
  training_set <- z[1:burn_in]
  stream_data  <- z[(burn_in + 1):length(z)]
  
  res <- ICM(
    training_set = training_set,
    stream_data  = stream_data,
    non_conformity_measure = non_conformity_measure,
    betting_function       = betting_function,
    th = threshold,
    ...
  )
  
  # Devolver sólo el índice del cambio (o Inf si no ocurrió)
  return(res$change_point)
}

# --- Wrappers especializados ------------------------------------------------

# 1) KNN + Constant BF
run_icm_knn_constant <- function(z, burn_in = 100, threshold, k = 1) {
  run_icm(
    z, burn_in, threshold,
    non_conformity_measure = Non_conformity_KNN,
    betting_function       = Constant_BF,
    k = k
  )
}

# 2) KNN + Mixture BF
run_icm_knn_mixture <- function(z, burn_in = 100, threshold, k = 1) {
  run_icm(
    z, burn_in, threshold,
    non_conformity_measure = Non_conformity_KNN,
    betting_function       = Mixture_BF,
    k = k
  )
}

# 3) KNN + KDE BF (uso de función precomputada o default)
run_icm_knn_kde <- function(z, burn_in = 100, threshold, bet_fun = NULL) {
  # Si no se pasa bet_fun, se asume que el usuario cargó un objeto KDE precomputado via readRDS
  if (is.null(bet_fun)) {
    stop("Debes proporcionar 'bet_fun', la función KDE precomputada cargada desde tu .rds")
  }
  run_icm(
    z, burn_in, threshold,
    non_conformity_measure = Non_conformity_KNN,
    betting_function       = bet_fun
  )
}

# 4) LNR + Constant BF
run_icm_lnr_constant <- function(z, burn_in = 100, threshold, mu_r = 1) {
  run_icm(
    z, burn_in, threshold,
    non_conformity_measure = Non_conformity_LNR,
    betting_function       = Constant_BF,
    mu_r = mu_r
  )
}

# 5) LNR + Mixture BF
run_icm_lnr_mixture <- function(z, burn_in = 100, threshold, mu_r = 1) {
  run_icm(
    z, burn_in, threshold,
    non_conformity_measure = Non_conformity_LNR,
    betting_function       = Mixture_BF,
    mu_r = mu_r
  )
}

# 6) LNR + KDE BF (uso de función precomputada o default)
run_icm_lnr_kde <- function(z, burn_in = 100, threshold, bet_fun = NULL) {
  if (is.null(bet_fun)) {
    stop("Debes proporcionar 'bet_fun', la función KDE precomputada cargada desde tu .rds")
  }
  run_icm(
    z, burn_in, threshold,
    non_conformity_measure = Non_conformity_LNR,
    betting_function       = bet_fun
  )
}