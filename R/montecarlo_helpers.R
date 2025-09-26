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
`%||%` <- function(a,b) if (is.null(a)) b else a

alphas_from_ncm <- function(stream, training_set, ncm_fun, k = NULL,
                            ..., shuffle_training = FALSE,
                            knn_stat = c("kth","mean")) {
  stopifnot(length(training_set) > 0)
  knn_stat <- match.arg(knn_stat)
  stream <- as.numeric(stream)
  training_set <- as.numeric(training_set)
  if (shuffle_training && length(training_set) > 1L) {
    training_set <- sample(training_set)
  }
  if (is.null(k)) k <- 7L
  
  # --- Caso 1: KNN (vectorizado con FNN) -------------------------------
  if (identical(ncm_fun, Non_conformity_KNN)) {
    if (!requireNamespace("FNN", quietly = TRUE)) {
      stop("Para KNN vectorizado instala 'FNN' (install.packages('FNN')).")
    }
    Xtr <- matrix(training_set, ncol = 1L)
    Xte <- matrix(stream,       ncol = 1L)
    nn  <- FNN::get.knnx(Xtr, Xte, k = k)$nn.dist  # n x k
    
    if (knn_stat == "kth") {
      # Si tu NCM_KNN usa el k-ésimo vecino:
      return(nn[, k])
    } else {
      # Si tu NCM_KNN usa promedio de k distancias:
      return(rowMeans(nn))
    }
  }
  
  # --- Caso 2: MAD (vectorizado) --------------------------------------
  if (identical(ncm_fun, Non_conformity_MAD)) {
    med  <- stats::median(training_set)
    s    <- stats::mad(training_set, constant = 1)  # usa tu constante si difiere
    s    <- if (s > 0) s else .Machine$double.eps
    return(abs(stream - med) / s)
  }
  
  # --- Caso 3: LR/LNR (vectorizado) -----------------------------------
  if (identical(ncm_fun, Non_conformity_LNR)) {
    dots <- list(...)
    mu_r <- dots$mu_r %||% 0
    mu0  <- mean(training_set)
    sd0  <- stats::sd(training_set); sd0 <- if (sd0 > 0) sd0 else .Machine$double.eps
    z    <- (stream - mu0 - mu_r) / sd0
    return(abs(z))  # ajusta si tu NCM difiere (e.g., z^2)
  }
    
  # --- Caso 4: IQR (NUEVO, vectorizado) --------------------------------
  if (identical(ncm_fun, Non_conformity_IQR)) {
    dots <- list(...)
    probs <- dots$probs %||% c(0.25, 0.50, 0.75)
    qtype <- dots$type  %||% 8L          # usas type = 8; lo respetamos
    c_iqr <- dots$c_iqr %||% 1.0         # si quisieras “normalizar” a sd, usa c_iqr = 1.349
    qs    <- as.numeric(stats::quantile(training_set, probs = probs, type = qtype, names = FALSE))
    med   <- qs[2]
    width <- max(qs[3] - qs[1], 1e-12)   # evita división por 0
    return(abs(stream - med) / (width * c_iqr))
  }
  
  # --- Fallback genérico (compatibilidad) ------------------------------
  vapply(stream, function(xi) {
    ncm_fun(xi = xi, training_set = training_set, k = k, ...)
  }, numeric(1))
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


## Helpers for multiple detection

make_stream_mean_shifts <- function(n_stream, theta_vec, mu_levels, sd = 1) {
  stopifnot(length(mu_levels) == length(theta_vec) + 1)
  if (length(theta_vec)) theta_vec <- sort(unique(as.integer(theta_vec)))
  seg_starts <- c(1L, theta_vec)
  seg_ends   <- c(if (length(theta_vec)) theta_vec - 1L else integer(0), n_stream)
  z <- numeric(0)
  for (j in seq_along(mu_levels)) {
    len_j <- seg_ends[j] - seg_starts[j] + 1L
    if (len_j > 0) z <- c(z, rnorm(len_j, mean = mu_levels[j], sd = sd))
  }
  list(stream = z, true_changes = as.integer(theta_vec))
}

# 2) Emparejamiento alarmas↔cambios con métricas multi -------------------------
# window_mode = "abs" usa ventana absoluta (puntos); "frac" usa fracción del segmento post-cambio
# Retorna: alarms_df (fila por alarma) y per_change (fila por cambio verdadero)
match_alarms_to_changes_multi <- function(
    alarms, true_changes, n_stream,
    window_mode = c("abs","frac"),
    window_abs = Inf, window_frac = 1.0
){
  window_mode <- match.arg(window_mode)
  alarms <- sort(unique(as.integer(alarms)))
  cj <- sort(as.integer(true_changes))
  J <- length(cj)
  
  # Segmentos [c_j, seg_end_j]
  seg_end <- if (J > 1) c(cj[-1] - 1L, n_stream) else n_stream
  seg_len <- if (J > 0) pmax(0L, seg_end - cj + 1L) else integer(0)
  
  # Ventanas (numéricas, sin convertir Inf a entero)
  if (window_mode == "abs") {
    win_raw <- rep(window_abs, J)
  } else {
    win_raw <- pmax(1, floor(as.numeric(seg_len) * window_frac))
  }
  # Ventana efectiva recortada por el segmento
  win_eff <- pmin(win_raw, as.numeric(seg_len))
  
  if (length(alarms) == 0L || J == 0L) {
    return(list(
      alarms_df = tibble(
        alarm = integer(0), is_tp = logical(0),
        change_id = integer(0), delay = numeric(0)
      ),
      per_change = tibble(
        change_id = seq_len(J), seg_len = seg_len,
        tp = 0L, delay = NA_real_, extras_after_tp = 0L, fa_in_segment = 0L
      ),
      tp_count = 0L, fa_count = length(alarms), miss_count = J,
      delays = numeric(0), extras_after_tp = 0L
    ))
  }
  
  used_alarm <- rep(FALSE, length(alarms))
  change_id_of_alarm <- rep(NA_integer_, length(alarms))
  delay_of_alarm <- rep(NA_real_, length(alarms))
  tp_per_change <- integer(J)
  delay_per_change <- rep(NA_real_, J)
  extras_after_tp <- integer(J)
  fa_in_segment <- integer(J)
  
  # Emparejar primer TP por cambio, dentro de [c_j, min(c_j + win_eff_j, seg_end_j)]
  for (j in seq_len(J)) {
    left_j  <- cj[j]
    right_j <- min(cj[j] + win_eff[j], seg_end[j]) # recorte por segmento
    idx <- which(!used_alarm & alarms >= left_j & alarms <= right_j)
    if (length(idx)) {
      a_idx <- idx[1]
      used_alarm[a_idx] <- TRUE
      change_id_of_alarm[a_idx] <- j
      delay_j <- alarms[a_idx] - cj[j]
      delay_of_alarm[a_idx] <- delay_j
      tp_per_change[j] <- 1L
      delay_per_change[j] <- delay_j
    }
  }
  
  # Contar extras y FAs por segmento
  for (j in seq_len(J)) {
    seg_l <- cj[j]
    seg_r <- seg_end[j]
    in_seg <- which(alarms >= seg_l & alarms <= seg_r)
    
    if (tp_per_change[j] == 1L) {
      # Posición del TP asignado a este cambio
      tp_pos <- alarms[which(used_alarm & change_id_of_alarm == j)][1]
      extras_after_tp[j] <- sum(alarms[in_seg] > tp_pos)
      fa_in_segment[j]   <- sum(!used_alarm[in_seg] & alarms[in_seg] < tp_pos)
    } else {
      extras_after_tp[j] <- 0L
      fa_in_segment[j]   <- length(in_seg) # todo es FA dentro del segmento
    }
  }
  
  is_tp_vec <- !is.na(change_id_of_alarm)
  
  # Aserciones defensivas (si algo viola, te avisa enseguida)
  bad_delay <- which(!is.na(delay_of_alarm) & delay_of_alarm < 0)
  if (length(bad_delay)) {
    stop(sprintf("Delay negativo en alarmas: %s", paste(alarms[bad_delay], collapse = ",")))
  }
  # delay no debe exceder ni la ventana efectiva ni el largo del segmento
  if (any(tp_per_change == 1L)) {
    ok <- TRUE
    for (j in which(tp_per_change == 1L)) {
      if (delay_per_change[j] > win_eff[j] + 1e-9 ||
          delay_per_change[j] > seg_len[j] + 1e-9) {
        ok <- FALSE
      }
    }
    if (!ok) stop("Se detectó un delay > win_eff o > seg_len (esto no debería pasar).")
  }
  
  alarms_df <- tibble(
    alarm     = alarms,
    is_tp     = is_tp_vec,
    change_id = change_id_of_alarm,
    delay     = delay_of_alarm
  )
  
  per_change <- tibble(
    change_id = seq_len(J),
    seg_len   = seg_len,
    tp        = tp_per_change,
    delay     = delay_per_change,
    extras_after_tp = extras_after_tp,
    fa_in_segment   = fa_in_segment
  )
  
  list(
    alarms_df       = alarms_df,
    per_change      = per_change,
    tp_count        = sum(tp_per_change),
    fa_count        = sum(!is_tp_vec),                # sin doble conteo
    miss_count      = J - sum(tp_per_change),
    delays          = alarms_df$delay[alarms_df$is_tp],
    extras_after_tp = sum(extras_after_tp)
  )
}

# --- DROP-IN: detector sin order() -----------------------------------
.detect_first_alarm_from_alphas <- function(alphas_seg, bet_fun, th, params_bf) {
  # 0) Chequeos de tipo (si falla, te dirá exactamente qué llegó mal)
  if (!is.numeric(alphas_seg)) {
    stop(sprintf(".detect_first_alarm_from_alphas: 'alphas_seg' no es numeric (clase: %s)",
                 paste(class(alphas_seg), collapse = "/")))
  }
  if (!is.function(bet_fun)) {
    stop(sprintf(".detect_first_alarm_from_alphas: 'bet_fun' no es function (clase: %s)",
                 paste(class(bet_fun), collapse = "/")))
  }
  if (!is.list(params_bf)) {
    stop(sprintf(".detect_first_alarm_from_alphas: 'params_bf' no es list (clase: %s)",
                 paste(class(params_bf), collapse = "/")))
  }
  if (!is.numeric(th) || length(th) != 1L || !is.finite(th)) {
    stop(".detect_first_alarm_from_alphas: 'th' debe ser numérico escalar finito")
  }
  
  n <- length(alphas_seg)
  if (n == 0L) return(NA_integer_)
  
  # 1) Evaluador de la BF sin 'order()' ni 'do.call' innecesario
  #    (si tus BFs SIEMPRE usan los mismos nombres, llama directo)
  bf_eval <- function(p_hist, new_p, i) {
    do.call(bet_fun, c(list(p_values = p_hist, new_p = new_p, i = i), params_bf))
  }
  
  # 2) Estado del tramo
  C_prev <- 0.0
  p_seg  <- numeric(n)
  # estructura ordenada incremental de alphas (sin llamar order() global)
  sorted_vals <- numeric(0L)
  
  for (i in seq_len  (n)) {
    ai  <- alphas_seg[i]
    pos <- findInterval(ai, sorted_vals)      # #anteriores <= ai
    gt_prev <- length(sorted_vals) - pos      # #anteriores > ai
    
    eq_prev <- 0L
    if (length(sorted_vals)) {
      L <- pos;       while (L >= 1L && sorted_vals[L] == ai) { eq_prev <- eq_prev + 1L; L <- L - 1L }
      R <- pos + 1L;  while (R <= length(sorted_vals) && sorted_vals[R] == ai) { eq_prev <- eq_prev + 1L; R <- R + 1L }
    }
    eq <- eq_prev + 1L   # <-- suma el actual
    
    u <- runif(1)
    if (u <= 0) u <- .Machine$double.eps
    if (u >= 1) u <- 1 - .Machine$double.eps
    p_i <- (gt_prev + u * eq) / i
    p_seg[i] <- p_i
    
    # insertar ai en sorted_vals (mantener ordenado sin order())
    if (pos == length(sorted_vals)) {
      sorted_vals <- c(sorted_vals, ai)
    } else if (pos == 0L) {
      sorted_vals <- c(ai, sorted_vals)
    } else {
      sorted_vals <- append(sorted_vals, ai, after = pos)
    }
    
    # apostar y acumular log-martingala truncada a 0
    g_i <- bf_eval(if (i == 1L) numeric(0) else p_seg[1:(i-1L)], p_i, i)
    if (!is.finite(g_i) || g_i <= 0) {
      # si la BF devuelve algo no válido, corta y reporta
      stop(sprintf("Betting function devolvió valor no positivo/finito en i=%d (g_i=%s)",
                   i, as.character(g_i)))
    }
    C_curr <- max(0, C_prev + log(g_i))
    if (C_curr > th) return(i)
    C_prev <- C_curr
  }
  
  NA_integer_
}

detect_alarms_from_pvals <- function(pvals, bet_fun, params_bf, h_vals) {
  n <- length(pvals)
  # g_t (factor de apuesta) en vector (sin suponer vectorización del BF)
  gvals <- vapply(pvals, function(p) {
    do.call(bet_fun, c(list(p = p), params_bf))
  }, numeric(1))
  
  # Escaneo por cada umbral: S se reinicia tras cruzar
  alarms_list <- vector("list", length(h_vals))
  for (j in seq_along(h_vals)) {
    h <- h_vals[j]
    S <- 1.0
    cps <- integer(0)
    for (t in seq_len(n)) {
      S <- S * gvals[t]
      if (S >= h) { cps <- c(cps, t); S <- 1.0 }
    }
    alarms_list[[j]] <- cps
  }
  list(gvals = gvals, alarms_by_h = alarms_list)
}


# === ICM_multi_fast: ENTRENAMIENTO FIJO (usa alphas precomputadas de TODO el stream) ===
ICM_multi_fast <- function(training_set,
                           stream_data,
                           non_conformity_measure,
                           betting_function,
                           th = 1,
                           params_bf = list(),
                           k = NULL,
                           alphas_full = NULL,
                           ...) {
  n <- length(stream_data)
  if (is.null(alphas_full)) {
    alphas_full <- alphas_from_ncm(stream_data, training_set,
                                   ncm_fun = non_conformity_measure, k = k, ...)
  }
  change_points <- integer(0); start_idx <- 1L
  while (start_idx <= n) {
    alphas_seg <- alphas_full[start_idx:n]
    tau_rel <- .detect_first_alarm_from_alphas(
      alphas_seg,
      bet_fun   = betting_function,
      th        = th,
      params_bf = params_bf
    )
    if (is.na(tau_rel)) break
    alarm_idx <- start_idx + tau_rel - 1L
    change_points <- c(change_points, alarm_idx)
    start_idx <- alarm_idx + 1L
  }
  list(change_points = change_points)
}

# === ICM_multi_adaptive_fast: REENTRENOS POR TRAMO (alphas por tramo) ===
ICM_multi_adaptive_fast <- function(stream_data,
                                    non_conformity_measure,
                                    betting_function,
                                    th = 1,
                                    training_set = NULL,
                                    training_size = NULL,
                                    m_retrain = NULL,
                                    guard_band = 0L,
                                    shuffle_training = TRUE,
                                    params_bf = list(),
                                    k = NULL,
                                    ...) {
  n <- length(stream_data)
  if (!is.numeric(th) || length(th) != 1L || !is.finite(th)) {
    stop("'th' debe ser numérico escalar finito")
  }
  if (is.null(k)) k <- 7L
  
  # --- entrenamiento inicial (igual lógica, con pequeñas defensas) ---
  if (!is.null(training_set)) {
    if (shuffle_training && length(training_set) > 1L) training_set <- sample(training_set)
    m0  <- length(training_set)
    pos <- 1L
  } else {
    if (is.null(training_size) || training_size <= 0L)
      stop("Especifica training_size > 0 si no pasas un training_set.")
    m0 <- as.integer(training_size)
    if (m0 >= n) stop("training_size no puede ser >= length(stream_data).")
    training_set <- stream_data[1:m0]
    if (shuffle_training && m0 > 1L) training_set <- sample(training_set)
    pos <- m0 + 1L
  }
  if (is.null(m_retrain)) m_retrain <- m0
  m_retrain <- as.integer(m_retrain)
  guard_band <- as.integer(guard_band)
  if (m_retrain <= 0L) stop("'m_retrain' debe ser > 0")
  
  # --- preasignación de alarmas y puntero ---
  cap   <- max(8L, min(n, ceiling(n / max(1L, m_retrain + guard_band + 1L))))
  cps   <- integer(cap)
  n_cps <- 0L
  
  # --- bucle por tramos ---
  while (pos <= n) {
    # 1) alphas del tramo actual (con tu NCM vectorizada si la tienes)
    #    Nota: si tu alphas_from_ncm ya es rápida, esto es el costo dominante.
    alphas_seg <- alphas_from_ncm(stream_data[pos:n], training_set,
                                  ncm_fun = non_conformity_measure, k = k, ...)
    
    # 2) detectar primera alarma relativa en el tramo
    tau_rel <- .detect_first_alarm_from_alphas(
      alphas_seg,
      bet_fun   = betting_function,
      th        = th,
      params_bf = params_bf
    )
    if (is.na(tau_rel)) break
    
    alarm_pos <- pos + tau_rel - 1L
    
    # 3) guardar alarma (sin concatenar)
    n_cps <- n_cps + 1L
    if (n_cps > length(cps)) cps <- c(cps, integer(length(cps)))  # crecimiento amortiguado
    cps[n_cps] <- alarm_pos
    
    # 4) ventana de reentrenamiento (con guard band)
    train_start <- alarm_pos + 1L + guard_band
    if (train_start > n) break
    train_end <- min(n, train_start + m_retrain - 1L)
    if (train_end < train_start) break
    
    new_train <- stream_data[train_start:train_end]
    if (shuffle_training && length(new_train) > 1L) new_train <- sample(new_train)
    training_set <- new_train
    
    # 5) siguiente tramo comienza después del reentrenamiento
    pos <- train_end + 1L
  }
  
  if (n_cps == 0L) {
    return(list(change_points_stream = integer(0)))
  } else {
    return(list(change_points_stream = cps[seq_len(n_cps)]))
  }
}

#------------------------ Helper outlier contamination-------------------------
rnorm_base <- function(n, mu, sd) stats::rnorm(n, mean = mu, sd = sd)

contaminate <- function(x, rate, model, params, base_mu = 0, base_sd = 1){
  stopifnot(rate >= 0, rate <= 1)
  n <- length(x)
  if (n == 0L || rate <= 0) return(x)
  
  idx <- which(stats::runif(n) < rate)
  if (length(idx) == 0L) return(x)
  
  delta  <- params$delta  %||% 6     # para "shift"
  lambda <- params$lambda %||% 5     # para "scale"
  
  if (model == "shift") {
    x[idx] <- stats::rnorm(length(idx), mean = base_mu + delta, sd = base_sd)
  } else if (model == "scale") {
    x[idx] <- stats::rnorm(length(idx), mean = base_mu, sd = base_sd * lambda)
  } else {
    stop("Modelo de contaminación no soportado: usa 'shift' o 'scale'.")
  }
  x
}

