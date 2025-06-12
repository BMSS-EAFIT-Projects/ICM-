#Bnechmarking

# 0) Cargar librerías y fuentes
library(microbenchmark)
library(dplyr)

# 0) Cargar librerías y fuentes
library(microbenchmark)
library(dplyr)

# Fuente de wrappers
source("R/benchmark_wrappers.R")

# 1) Cargar funciones KDE precomputadas (desde tus .rds)
kde_precomputed <- readRDS("data/kde_bf_fixed.rds")

# 2) Generar streams_list por simulación
# Parám. de simulación:
training_set_size <- 200  # número de puntos de entrenamiento (burn-in)
n_total          <- 1000 # longitud total de cada stream
# Escenarios de cambio (ajusta según tus necesidades):
thetas <- c(100, 200)
mu1s   <- c(1, 1.5, 2)

streams_list <- list()
for (theta in thetas) {
  for (mu1 in mu1s) {
    name <- paste0("theta", theta, "_mu1_", mu1)
    z <- c(
      rnorm(theta + training_set_size, mean = 0, sd = 1),
      rnorm(n_total - theta,         mean = mu1, sd = 1)
    )
    streams_list[[name]] <- z
  }
}

# 3) Parámetros de benchmarking
burn_in    <- training_set_size
thresholds <- seq(1, 6, 0.5)
times_each <- 20L  # repeticiones de microbenchmark

# 4) Contenedor para resultados crudos
timings <- list()

# 5) Iterar sobre escenarios y thresholds
for (scenario in names(streams_list)) {
  z <- streams_list[[scenario]]
  for (th in thresholds) {
    mb <- microbenchmark(
      KNN_const   = run_icm_knn_constant(z, burn_in, th),
      KNN_mixture = run_icm_knn_mixture(z, burn_in, th),
      KNN_kde     = run_icm_knn_kde(z, burn_in, th, bet_fun = kde_precomputed),
      LNR_const   = run_icm_lnr_constant(z, burn_in, th),
      LNR_mixture = run_icm_lnr_mixture(z, burn_in, th),
      LNR_kde     = run_icm_lnr_kde(z, burn_in, th, bet_fun = kde_precomputed),
      times       = times_each
    )
    df <- as.data.frame(mb)
    df$scenario  <- scenario
    df$threshold <- th
    timings[[paste0(scenario, "_th", th)]] <- df
  }
}

# 6) Combinar y guardar resultados crudos
all_timings <- bind_rows(timings)

output_dir <- "data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(all_timings, file.path(output_dir, "benchmarks_raw.rds"))
cat("Benchmark ejecutado y guardado")