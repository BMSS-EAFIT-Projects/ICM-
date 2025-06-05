# Oracles implementation

#Likelihood marginal without cp
log_marginal_likelihood <- function(z) {
  n <- length(z)
  z_bar <- mean(z)
  sum_z2 <- sum(z^2)
  return(
    -0.5 * n * log(2 * pi) +
      -0.5 * sum_z2 +
      0.5 * (n^2 * z_bar^2) / (n + 1) +
      -0.5 * log(n + 1)
  )
}

#likelihod marginal with cp
log_marginal_likelihood_cp <- function(z, theta) {
  n <- length(z)
  
  # Segmentos antes y después del cambio
  z1 <- z[1:(theta - 1)]
  z2 <- z[theta:n]
  
  # Sumas necesarias
  sum_z2_1 <- sum(z1^2)
  sum_z2_2 <- sum(z2^2)
  sum_z_1 <- sum(z1)
  sum_z_2 <- sum(z2)
  
  # Tamaños efectivos
  a <- theta - 1           # número de observaciones antes del cambio
  b <- n - theta + 1       # número de observaciones después del cambio
  
  # Log likelihood
  log_likelihood <- 
    -0.5 * n * log(2 * pi) +
    -0.5 * (sum_z2_1 + sum_z2_2) +
    0.5 * (sum_z_1^2) / (a + 1) +
    0.5 * (sum_z_2^2) / (b + 1) +
    -0.5 * log((a + 1) * (b + 1))
  
  return(log_likelihood)
}


#Cusum - oracle
cusum_oracle <- function(z, threshold) {
  n_total <- length(z)
  for (n in 2:n_total) {
    z_sub <- z[1:n]
    ll_total <- log_marginal_likelihood(z_sub)
    
    max_stat <- -Inf
    for (omega in 2:n) {  # omega = 1 no tiene parte anterior
      stat <- log_marginal_likelihood_cp(z_sub, omega) - ll_total
      if (stat > max_stat) max_stat <- stat
      if (stat >= threshold) break
    }
    
    if (max_stat >= threshold) {
      return(list(tau = n))
    }
  }
  return(list(tau = NA))
}

#Shiryaev-Roberts - oracle
sr_oracle <- function(z, threshold) {
  n_total <- length(z)
  for (n in 2:n_total) {
    z_sub <- z[1:n]
    ll_total <- log_marginal_likelihood(z_sub)
    
    sum_exp <- 0
    for (omega in 2:n) {
      stat <- log_marginal_likelihood_cp(z_sub, omega) - ll_total
      sum_exp <- sum_exp + exp(stat)
    }
    
    if (log(sum_exp) >= threshold) {
      return(list(tau = n))
    }
  }
  return(list(tau = NA))
}

# posterior probabilistic - oracle
posterior_oracle <- function(z, threshold, p = 1 / 100) {
  n_total <- length(z)
  
  for (n in 2:n_total) {
    z_sub <- z[1:n]
    ll_total <- log_marginal_likelihood(z_sub)
    
    weighted_sum <- 0
    for (omega in 2:n) {  # empieza en 2 porque omega=1 no tiene parte anterior
      prior <- p * (1 - p)^(omega - 1)
      
      # Likelihood con cambio en omega
      ll_cp <- log_marginal_likelihood_cp(z_sub, omega)
      
      weighted_sum <- weighted_sum + prior * exp(ll_cp - ll_total)
    }
    
    stat <- log(weighted_sum)
    if (stat >= threshold) {
      return(list(tau = n))
    }
  }
  
  return(list(tau = NA))
}
