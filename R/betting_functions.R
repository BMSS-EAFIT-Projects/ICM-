#All the functions here correspond to Betting Functions

#Constant BF
#refer to this araticle "Inductive conforma martingales for change point detection".

Constant_BF <- function(p_values, new_p, i, ...) {
  if (new_p >= 0 && new_p < 0.5) {
    return(1.5)
  } else if (new_p >= 0.5 && new_p <= 1) {
    return(0.5)
  } 
}

#Mixture BF
#refer to this article "Inductive conforma martingales for change point detection".

Mixture_BF <- function(p_values, new_p, i, ...) {
  integrand <- function(epsilon) {
    epsilon * new_p^(epsilon - 1)
  }
  result <- integrate(integrand, lower = 0, upper = 1)
  return(result$value)
}

#KDE BF
KDE <- function(p_values, n_grid = 512) {
  p_values <- as.numeric(p_values)
  n        <- length(p_values)
  
  # Bandwidth: regla nrd0 + fallback
  h0 <- stats::bw.nrd0(p_values)
  if (!is.finite(h0) || h0 <= 0) {
    h0 <- max(diff(range(p_values)) * n^(-1/5), 1e-3)
  }
  
  # Reflect samples to mitigate bordes
  extended <- c(-p_values, p_values, 2 - p_values)
  
  # Ajustar KDE en [-1.5, 2.5] y luego recortar a [0,1]
  kde <- stats::density(
    x      = extended,
    kernel = "gaussian",
    bw     = h0,
    n      = n_grid,
    from   = -1.5,
    to     = 2.5
  )
  in01  <- (kde$x >= 0 & kde$x <= 1)
  x01   <- kde$x[in01]
  y01   <- kde$y[in01]
  
  # Normalizar para que ∫₀¹g(p)dp = 1
  delta <- x01[2] - x01[1]
  area  <- sum(y01) * delta
  if (area > 0) {
    y01 <- y01 / area
  } else {
    y01 <- rep(1, length(y01))
    y01 <- y01 / (sum(y01) * delta)
  }
  
  # Devuelve una función estática g(p)
  function(p) {
    p <- as.numeric(p)
    out <- rep(0, length(p))
    inside <- which(p >= 0 & p <= 1)
    if (length(inside)>0) {
      out[inside] <- approx(x01, y01, xout = p[inside], rule = 2)$y
    }
    out
  }
}

#PRECOMPUTED KDE BF
#refer to this article "Inductive conforma martingales for change point detection".

Precomputed_KDE_BF <- function(training_set,
                               calibration_data,
                               non_conformity_measure,
                               k = 1,
                               n_grid = 512) {
  m <- length(calibration_data)
  alphas  <- numeric(m)
  pvals   <- numeric(m)
  
  # 1) Calcular α₁…αₘ
  for (i in seq_len(m)) {
    alphas[i] <- non_conformity_measure(calibration_data[i], training_set, k)
  }
  
  # 2) Calcular p₁…pₘ de forma secuencial
  for (i in seq_len(m)) {
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    pvals[i] <- (greater + runif(1) * equal) / i
  }
  
  # 3) Ajustar KDE
  g_p <- KDE(pvals, n_grid = n_grid)
  
  # 4) Retornar función que se ajusta al framework esperado
  function(p_values, new_p, i, ...) {
    g_p(new_p)
  }
}

#HISTOGRAM BF
#refer to this article "A histogram based betting function for conformal martingales".
histogram_betting_function <- function(p_values, new_p, i, num_bins = 10L) {
  n <- length(p_values) + 1L
  if (n <= 1L) return(1)
  
  clamp01 <- function(x) pmin(1, pmax(0, x))
  new_p <- clamp01(new_p)
  p_prev <- clamp01(p_values)
  
  bin_of <- function(p) {
    b <- ceiling(p * k)
    pmax(1L, pmin(k, b))
  }
  
  counts <- tabulate(bin_of(p_prev), nbins = k)
  
  j <- bin_of(new_p)
  if (counts[j] == 0L) return(1)
  k * counts[j] / (n - 1L)
}
