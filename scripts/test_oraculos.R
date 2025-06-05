set.seed(123)
n <- 1000
theta <- 500
mu1 <- 0
mu2 <- 1.5
z <- c(rnorm(theta - 1, mean = mu1), rnorm(n - theta + 1, mean = mu2))

for (th in 1:10){
  # Elegimos un threshold razonable (ajústalo según lo que observes)
  threshold <- 
  
  # Aplicar oráculos
  res_cusum <- cusum_oracle(z, th)
  res_sr <- sr_oracle(z, th)
  res_pp <- posterior_oracle(z, th)
  # Mostrar resultados
  cat("threshold", th, "\n")
  cat("Punto de cambio real: ", theta, "\n")
  cat("CUSUM Oracle detecta en: ", res_cusum$tau, "\n")
  cat("S-R Oracle detecta en: ", res_sr$tau, "\n")
  cat("Posterior Oracle detecta en: ", res_pp$tau, "\n")
}
