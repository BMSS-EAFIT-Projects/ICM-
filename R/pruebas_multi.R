set.seed(123)
source("R/non_conformity_measures.R")
source("R/icm_method.R")
source("R/betting_functions.R")

n <- 1000
training_set_size <- 100
theta <- 500

z <- c(rnorm((theta + training_set_size), mean = 0, sd = 1),rnorm(n - theta, mean = 2, sd = 1), rnorm(500, mean = 0, sd = 1))

training_set <- z[1:training_set_size]
stream_data <- z[-(1:training_set_size)]

res_adapt_multi <- ICM_multi_adaptive(stream_data,Non_conformity_KNN,Mixture_BF,th = 4, training_set, k = 7, m_retrain = 50)
res_multi <- ICM_multi(training_set, stream_data, Non_conformity_KNN, Mixture_BF, th = 4, k = 7)

res_adapt_multi$change_points_stream
res_multi$change_points

plot(res_adapt_multi$Cn_aligned, type = 'l')
plot(res_multi$Cn, type = 'l')