#All the functions here correspond to Non conformity measures

#NCM_KNN
#refer to this article "Inductive conforma martingales for change point detection"

Non_conformity_KNN <- function(xi, training_set, k = 1,...) {
  dist <- abs(xi - training_set)
  NN <- sort(dist)[1:k]
  return(mean(NN))
}

#Likelyhood ratio
#refer to this article "Inductive conforma martingales for change point detection".

Non_conformity_LNR <- function(xi, training_set, mu_r = 1, k=NULL,...) {
  # train_set: vector fijo con datos anteriores a z_1
  mu0_hat <- mean(training_set)
  sigma_f1 <- sqrt(1 + 1)
  
  numerator   <- dnorm(xi, mean = mu_r, sd = sigma_f1)
  denominator <- dnorm(xi, mean = mu0_hat, sd = sqrt(1))
  
  return(numerator / denominator)
}

#MAD
Non_conformity_MAD <- function(xi, training_set, K=NULL,...){
  median_data <- median(training_set)
  zscore <- abs(xi - median_data)
  return(zscore)
}


 