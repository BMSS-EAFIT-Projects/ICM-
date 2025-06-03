# Oracles implementation

#Likelihood marginal
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

#Cusum - oracle
cusum_oracle <- function(z, threshold) {
  n_total <- length(z)
  for (n in 2:n_total) {
    z_sub <- z[1:n]
    ll_total <- log_marginal_likelihood(z_sub)
    
    max_stat <- -Inf
    for (omega in 1:n) {
      z_pre <- if (omega > 1) z_sub[1:(omega - 1)] else numeric(0)
      z_post <- z_sub[omega:n]
      
      ll_pre <- if (length(z_pre) == 0) 0 else log_marginal_likelihood(z_pre)
      ll_post <- log_marginal_likelihood(z_post)
      
      stat <- ll_pre + ll_post - ll_total
      
      if (stat > max_stat) {
        max_stat <- stat
      }
      
      # Early stopping interno opcional si supera threshold dentro del loop ω
      if (stat >= threshold) break
    }
    
    # Si el máximo ya supera h, guardamos y salimos
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
    for (omega in 1:n) {
      z_pre <- if (omega > 1) z_sub[1:(omega - 1)] else numeric(0)
      z_post <- z_sub[omega:n]
      
      ll_pre <- if (length(z_pre) == 0) 0 else log_marginal_likelihood(z_pre)
      ll_post <- log_marginal_likelihood(z_post)
      
      sum_exp <- sum_exp + exp(ll_pre + ll_post - ll_total)
    }
    
    stat <- log(sum_exp)
    if (stat >= threshold) {
      return(list(tau = n))
    }
  }
  
  return(list(tau = NA))
}

# posterior probabilistic - oracle
posterior_oracle <- function(z, threshold, p = 1/100) {
  n_total <- length(z)
  
  for (n in 2:n_total) {
    z_sub <- z[1:n]
    ll_total <- log_marginal_likelihood(z_sub)
    
    weighted_sum <- 0
    for (omega in 1:n) {
      prior <- p * (1 - p)^(omega - 1)
      
      z_pre <- if (omega > 1) z_sub[1:(omega - 1)] else numeric(0)
      z_post <- z_sub[omega:n]
      
      ll_pre <- if (length(z_pre) == 0) 0 else log_marginal_likelihood(z_pre)
      ll_post <- log_marginal_likelihood(z_post)
      
      weighted_sum <- weighted_sum + prior * exp(ll_pre + ll_post - ll_total)
    }
    
    stat <- log(weighted_sum)
    if (stat >= threshold) {
      return(list(tau = n))
    }
  }
  
  return(list(tau = NA))
}

#Method fot calling the oracle
evaluate_oracle <- function(sim_data, detector_fn, threshold, label, ...) {
  n_sims <- length(sim_data)
  res <- data.frame(
    sim = integer(n_sims),
    tau = integer(n_sims),
    omega_hat = integer(n_sims),
    delay = numeric(n_sims),
    false_alarm = logical(n_sims),
    method = label
  )
  
  for (i in seq_along(sim_data)) {
    z <- sim_data[[i]]$z
    omega_real <- sim_data[[i]]$omega_real
    
    out <- detector_fn(z, threshold = threshold, ...)
    
    res[i, "sim"] <- i
    res[i, "tau"] <- out$tau
    res[i, "omega_hat"] <- out$omega_hat
    res[i, "delay"] <- if (!is.na(out$tau) && out$tau > omega_real) out$tau - omega_real else NA
    res[i, "false_alarm"] <- ifelse(!is.na(out$tau) && out$tau < omega_real, TRUE, FALSE)
  }
  
  return(res)
}