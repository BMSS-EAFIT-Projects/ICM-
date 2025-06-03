set.seed(123)
source("R/non_conformity_measures.R")
source("R/icm_method.R")
source("R/betting_functions.R")

n <- 1000
theta <- 500

# Data frame base con tiempo
df_fig_2 <- data.frame(time = 1:n)

# Simulación con CP (LR, m = 1)
training_set_size <- 1
z <- c(rnorm(theta + training_set_size, mean = 0, sd = 1), rnorm(n - theta, mean = 1, sd = 1))

training_set <- z[1:training_set_size]
stream_data <- z[-(1:training_set_size)]

res_cp <- ICM(training_set, stream_data, Non_conformity_KNN, Constant_BF, k = 1)

df_fig_2$con_cp <- res_cp$Cn

# Simulaciones sin CP para m = 1:5
for (i in 1:5) {
  z_wcp <- rnorm(n + i, mean = 0, sd = 1)  # sin CP
  training_set_wcp <- z_wcp[1:i]
  stream_data_wcp  <- z_wcp[-(1:i)]
  
  res_m <- ICM(training_set_wcp, stream_data_wcp, Non_conformity_KNN, Constant_BF, k = i)
  
  df_fig_2[[paste0("m_", i)]] <- res_m$Cn
}

# Guardar para graficar en el notebook
saveRDS(df_fig_2, file = "data/figura2_icm.rds")