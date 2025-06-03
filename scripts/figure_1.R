set.seed(123)
source("R/non_conformity_measures.R")
source("R/icm_method.R")
source("R/betting_functions.R")

n <- 1000
training_set_size <- 1
theta <- 500

z_wcp <-c(rnorm(n + training_set_size, 0, 1))
z <- c(rnorm((theta + training_set_size), mean = 0, sd = 1),rnorm(n - theta, mean = 1, sd = 1))

training_set <- z[1:training_set_size]
stream_data <- z[-(1:training_set_size)]

training_set_wcp <-z_wcp[1:training_set_size]
stream_data_wcp <- z_wcp[-(1:training_set_size)]

res_1nn <- ICM(training_set, stream_data, Non_conformity_KNN, Mixture_BF, k=1)
res_ln <- ICM(training_set, stream_data, Non_conformity_LNR, Constant_BF)
res_wcp <- ICM(training_set_wcp, stream_data_wcp, Non_conformity_KNN, Mixture_BF)

df_fig_1 <- data.frame(
  time = 1:length(res_1nn$Cn),
  Cn_1NN = res_1nn$Cn,
  Cn_LR  = res_ln$Cn,
  Cn_wcp = res_wcp$Cn
)
saveRDS(df_fig_1, file = "data/figura1_icm.rds")