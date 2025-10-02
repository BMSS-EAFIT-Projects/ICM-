source("R/non_conformity_measures.R")
source("R/icm_method.R")
source("R/betting_functions.R")

n <- 1000
training_set_size <- 200
theta <- 500

training_set <- rnorm(training_set_size, 0, 1)
z <- c(rnorm(theta, mean = 0, sd = 1),rnorm(n - theta, mean = 1, sd = 1))

stream_data <- z[-(1:training_set_size)]

res_knn <- ICM(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_KNN,
  betting_function = histogram_betting_function,
  th = 4.5,
  params_bf = list(num_bins = 10),
  k = 7
)

res_lnr <- ICM(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_LNR,
  betting_function = histogram_betting_function,
  th = 4.5,
  params_bf = list(num_bins = 10),
  mu_r = 1
)

res_MAD <- ICM(
  training_set = training_set,
  stream_data = stream_data,
  non_conformity_measure = Non_conformity_MAD,
  betting_function = histogram_betting_function,
  th = 4.5,
  params_bf = list(num_bins = 10)
)

plot(res_knn$Cn, type = "l", col = "blue", lwd = 2,
     ylim = range(c(res_knn$Cn, res_lnr$Cn)),
     xlab = "Tiempo", ylab = expression(C[n]),
     main = "Martingala ICM - Histograma")

lines(res_lnr$Cn, col = "darkgreen", lwd = 2)
abline(v = theta, lty = 2, col = "gray")

legend("topleft",
       legend = c("KNN", "MAD"),
       col = c("blue", "darkgreen"),
       lwd = 2)
