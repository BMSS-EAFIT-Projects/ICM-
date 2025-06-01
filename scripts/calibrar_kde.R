# 1) Cargar funciones
source("R/non_conformity_measures.R")
source("R/betting_functions.R")

# 2) Definir parámetros de calibración
set.seed(2025)
m0          <- 200
n_calib     <- 1000
epsilon_cal <- 500
mu1_cal     <- 1
k0          <- ceiling(m0 / 2)

# 3) Generar training_set y calibration_data tal como en la Sección 2
train_for_kde <- rnorm(m0, 0, 1)
calib_stream  <- numeric(n_calib)
for (i in seq_len(n_calib)) {
  idx_global <- m0 + i
  if (idx_global < epsilon_cal) {
    calib_stream[i] <- rnorm(1, 0, 1)
  } else {
    calib_stream[i] <- rnorm(1, mu1_cal, 1)
  }
}

# 4) Ajustar la KDE y guardarla
kde_bf_fixed <- Precomputed_KDE_BF(
  training_set           = train_for_kde,
  calibration_data       = calib_stream,
  non_conformity_measure = Non_conformity_KNN,
  k                      = k0,
  n_grid                 = 512
)

# 5) Guardar la función en un archivo .rds para poder cargarla luego sin recalibrar
saveRDS(kde_bf_fixed, file = "data/kde_bf_fixed.rds")