# --- p-values estándar ICM a partir de alphas ---
pvals_from_alphas <- function(alphas) {
  n <- length(alphas); p <- numeric(n)
  for (i in seq_len(n)) {
    ai <- alphas[i]; prev <- alphas[seq_len(i)]
    gt <- sum(prev > ai); eq <- sum(prev == ai)
    p[i] <- (gt + runif(1) * eq) / i
  }
  p
}

# Calculate alphas for being re-used
alphas_from_ncm <- function(stream, training_set, ncm_fun, k = NULL, ...) {
  n <- length(stream); a <- numeric(n)
  # (opcional para replicar tu ICM): barajar el training_set
  training_set <- sample(training_set)
  for (i in seq_len(n)) {
    a[i] <- ncm_fun(xi = stream[i], training_set = training_set, k = k, ...)
  }
  a
}

# Calculate trayectories once per BF
icm_C_path_from_p <- function(p, betting_function, params_bf = list()) {
  n <- length(p); C <- numeric(n); C_prev <- 0
  for (i in seq_len(n)) {
    p_hist <- if (i == 1) numeric(0) else p[1:(i - 1)]
    g_i <- do.call(betting_function,
                   c(list(p_values = p_hist, new_p = p[i], i = i), params_bf))
    if (!is.finite(g_i) || g_i <= 0) g_i <- 1e-12
    C_prev <- max(0, C_prev + log(g_i))
    C[i] <- C_prev
  }
  C
}

icm_cbf_S_path_from_p <- function(p, bet_fun, W = 100, epsilon = 1.5, params_bf = list()) {
  n <- length(p); S <- numeric(n); S_prev <- 1.0
  step_bf <- function(prev_p, new_p, i) {
    fi <- do.call(bet_fun, c(list(p_values = prev_p, new_p = new_p, i = i), params_bf))
    if (!is.finite(fi) || fi <= 0) fi <- 1e-12
    fi
  }
  for (i in seq_len(n)) {
    if (i == 1) {
      ratio <- 0
    } else {
      w <- min(W, i - 1)
      min_window <- min(S[(i - w):(i - 1)])
      if (min_window <= 0) min_window <- 1e-12
      ratio <- (S_prev - min_window) / min_window
    }
    fi <- if (ratio <= epsilon) 1.0 else {
      prev_p <- if (i == 1) numeric(0) else p[seq_len(i - 1)]
      step_bf(prev_p, p[i], i)
    }
    S_prev <- S_prev * fi
    S[i] <- S_prev
  }
  S
}

# Determine the change point once per trayectories
first_cross_indices_linear <- function(path, h_vals) {
  ord <- order(h_vals)
  h_sorted <- h_vals[ord]
  idx_sorted <- rep(NA_integer_, length(h_sorted))
  j <- 1L
  for (i in seq_along(path)) {
    while (j <= length(h_sorted) && path[i] >= h_sorted[j]) {
      if (is.na(idx_sorted[j])) idx_sorted[j] <- i
      j <- j + 1L
    }
    if (j > length(h_sorted)) break
  }
  idx <- idx_sorted[order(ord)]
  idx
}