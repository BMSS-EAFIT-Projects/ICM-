# Oracles implementation

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

#Likelihood marginal without cp
log_marginal_likelihood <- function(z) {
  n  <- length(z)
  s  <- sum(z)
  s2 <- sum(z^2)
  -0.5*n*log(2*pi) - 0.5*log(n+1) - 0.5*s2 + 0.5*(s^2)/(n+1)
}


#likelihod marginal with cp
.make_llcp_vec <- function(z, n) {
  cz  <- cumsum(z);       cz2 <- cumsum(z^2)
  s   <- cz[n];           s2  <- cz2[n]
  # omega = 1..n  => a = omega-1, b = n - a
  a   <- 0:(n-1);         b   <- n - a
  s1  <- c(0,  cz[1:(n-1)])     # sum_{1..a}
  ss1 <- c(0, cz2[1:(n-1)])
  s2s  <- s  - s1                 # sum_{omega..n}
  ss2s <- s2 - ss1
  (1 - n/2)*log(2*pi) - 0.5*log( (a+1)*(b+1) ) -
    0.5*(ss1 + ss2s) + 0.5*(s1^2)/(a+1) + 0.5*(s2s^2)/(b+1)
}


#Cusum - oracle
cusum_oracle <- function(z, threshold) {
  nT <- length(z)
  for (n in 2:nT) {
    ll_tot <- log_marginal_likelihood(z[1:n])
    ll_cp  <- .make_llcp_vec(z, n)
    stat   <- max(ll_cp - ll_tot)
    if (stat >= threshold) return(list(tau = n, stat = stat))
  }
  list(tau = NA)
}

#Shiryaev-Roberts - oracle
sr_oracle <- function(z, threshold) {
  nT <- length(z)
  for (n in 2:nT) {
    ll_tot <- log_marginal_likelihood(z[1:n])
    ll_cp  <- .make_llcp_vec(z, n)
    stat   <- logsumexp(ll_cp - ll_tot)
    if (stat >= threshold) return(list(tau = n, stat = stat))
  }
  list(tau = NA)
}

# posterior probabilistic - oracle
posterior_oracle <- function(z, threshold, p = 1/100) {
  nT <- length(z)
  for (n in 2:nT) {
    ll_tot <- log_marginal_likelihood(z[1:n])
    ll_cp  <- .make_llcp_vec(z, n)
    omega  <- 1:n
    logpri <- log(p) + (omega - 1) * log1p(-p)  # log[p (1-p)^{omega-1}]
    stat   <- logsumexp(logpri + (ll_cp - ll_tot)) - n*log1p(-p)
    if (stat >= threshold) return(list(tau = n, stat = stat))
  }
  list(tau = NA)
}