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
  
  # 1) Calcular α₁…αₘ vs el mismo training_set
  for (i in seq_len(m)) {
    alphas[i] <- non_conformity_measure(calibration_data[i],
                                        training_set, k)
  }
  
  # 2) Calcular p₁…pₘ de forma secuencial
  for (i in seq_len(m)) {
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    pvals[i] <- (greater + runif(1)*equal) / i
  }
  
  # 3) Ajustar KDE una sola vez
  KDE(pvals, n_grid = n_grid)
}

#HISTOGRAM BF
#refer to this article "A histogram based betting function for conformal martingales".

histogram_betting_function <- function(p_values, new_p, i, num_bins = 10, ...) {
  if (length(p_values) == 0) return(1)  # No apostar si no hay historial
  
  breaks <- seq(0, 1, length.out = num_bins + 1)
  counts <- hist(p_values, breaks = breaks, plot = FALSE)$counts
  bin_index <- findInterval(new_p, breaks, rightmost.closed = TRUE)
  density <- if (sum(counts) == 0) 1 else (counts[bin_index] / sum(counts)) * num_bins
  return(max(density, 1e-10))
}

#CAUTIOUS BF
cautious_wrapper <- function(bf_function,             # función de apuesta base: g(p)
                             p_values,                # vector de p-valores previos
                             s_values,                # vector de martingala previa
                             new_p,                   # nuevo p-valor a transformar
                             W = 100,                 # tamaño de la ventana
                             epsilon = 100,           # umbral de cambio
                             ...) {                   # argumentos adicionales para bf_function
  
  n <- length(s_values)
  
  # Si no hay suficientes valores previos, no apostar aún
  if (n < W) return(1)
  
  # Calcular el ratio de crecimiento reciente
  s_window <- s_values[(n - W + 1):n]
  s_ratio <- s_values[n] / min(s_window)
  
  # Si no hay suficiente evidencia para apostar, no apuestes (retorna 1)
  if (s_ratio <= epsilon) {
    return(1)
  } else {
    # Si hay evidencia, usa la función base con los mismos argumentos
    return(bf_function(p_values = p_values, new_p = new_p, ...))
  }
}