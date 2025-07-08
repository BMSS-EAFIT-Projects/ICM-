set.seed(123)
source("R/non_conformity_measures.R")
source("R/icm_method.R")
source("R/betting_functions.R")

n <- 1000
training_set_size <- 200
theta <- 500
W <- 50
z <- c(rnorm((theta + training_set_size), mean = 0, sd = 1),rnorm(n - theta, mean = 1, sd = 1))

training_set <- z[1:training_set_size]
stream_data <- z[-(1:training_set_size)]

res_knn_CBF <- ICM_CBF(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_KNN,
  betting_function = histogram_betting_function,
  W = W,
  epsilon = 0.01,
  params_bf = list(num_bins = 10),
  k = 7
)

res_lnr_CBF <- ICM_CBF(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_LNR,
  betting_function = histogram_betting_function,
  W = W,
  epsilon = 0.01,
  params_bf = list(num_bins = 10),
  mu_r = 1
)

res_knn_const_cbf <- ICM_CBF(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_KNN,
  betting_function = Constant_BF,
  W = W,
  epsilon = 0.01,
  k = 7
)

res_lnr_const_cbf <- ICM_CBF(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_LNR,
  betting_function = Constant_BF,
  W = W,
  epsilon = 0.01,
  mu_r = 1
)

res_knn_mix_cbf <- ICM_CBF(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_KNN,
  betting_function = Mixture_BF,
  W = W,
  epsilon = 0.01,
  k = 7
)

res_lnr_mix_cbf <- ICM_CBF(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_LNR,
  betting_function = Mixture_BF,
  W = W,
  epsilon = 0.01,
  mu_r = 1
)


plot((res_knn_CBF$S_vals), type = "l", col = "blue", lwd = 2,
     xlab = "Tiempo", ylab = expression(S[n]),
     main = "ICM_CBF - Comparación de funciones de apuesta")

lines((res_lnr_const_cbf$S_vals), col = "darkblue", lwd = 2, lty = 2)
lines((res_knn_mix_cbf$S_vals), col = "darkgreen", lwd = 2)
lines((res_lnr_mix_cbf$S_vals), col = "forestgreen", lwd = 2, lty = 2)
lines((res_knn_CBF$S_vals), col = "orange", lwd = 2)
lines((res_lnr_CBF$S_vals), col = "red", lwd = 2, lty = 2)
abline(v = 500, col = "gray", lty = 3)

#legend("bottomleft",
#       legend = c("Const-KNN", "Const-LNR", "Mix-KNN", "Mix-LNR", "Hist-KNN", "Hist-LNR"),
#       col = c("blue", "darkblue", "darkgreen", "forestgreen", "orange", "red"),
#       lwd = 2, lty = c(1,2,1,2,1,2))