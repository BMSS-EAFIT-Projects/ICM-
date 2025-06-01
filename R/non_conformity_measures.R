#All the functions here correspond to Non conformity measures

#NCM_KNN
Non_conformity_KNN <- function(xi, training_set, k) {
  dist <- abs(xi - training_set)
  NN <- sort(dist)[1:k]
  return(mean(NN))
}

#Likelihood ratio
