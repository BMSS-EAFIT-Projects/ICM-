#Basic method of ICM

#icm method
#refer to this article "Inductive conforma martingales for change point detection".

ICM <- function(training_set,
                stream_data,
                non_conformity_measure,
                betting_function,
                th = 1,
                params_bf = list(),  # argumentos para la betting function
                ...) {               # argumentos para la NCM
  
  training_set <- sample(training_set)
  m <- length(training_set)
  n <- length(stream_data)
  
  alphas <- numeric(n)
  p_vals <- numeric(n)
  Cn     <- numeric(n)
  C_prev <- 0
  
  alarm<-NA
  for (i in seq_len(n)) {
    xi <- stream_data[i]
    
    alphas[i] <- do.call(non_conformity_measure,
                         c(list(xi = xi, training_set = training_set), list(...)))
    
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    p_vals[i] <- (greater + runif(1) * equal) / i
    
    p_hist <- if (i == 1) numeric(0) else p_vals[1:(i - 1)]
    
    g_i <- do.call(betting_function,
                   c(list(p_values = p_hist, new_p = p_vals[i], i = i), params_bf))
    
    val <- if (g_i > 0) (C_prev + log(g_i)) else -Inf
    Cn[i] <- max(0, val)
    C_prev <- Cn[i]
    
    if (Cn[i] > th) {
      alarm <- i
      break
    }
  }
  
  if (!is.na(alarm)) {
    alphas <- alphas[1:alarm]
    p_vals <- p_vals[1:alarm]
    Cn     <- Cn    [1:alarm]
  }
  
  list(
    Cn           = Cn,
    p_vals       = p_vals,
    change_point = alarm
  )
}


#ICM CAUTIOUS
#refer to this method in this article "A betting function for addresing concept drift with conformal martingales"

ICM_CBF <- function(training_set,
                    stream_data,
                    non_conformity_measure,
                    betting_function,
                    W        = 100,   # tamaño de la “ventana caliente”
                    epsilon  = 1.5,   # umbral de descongelado
                    th       = 20,    # ← umbral de alarma para S_n
                    params_bf = list(),   # args extra para BF
                    ...) {             # args extra para la NCM
  
  ## --- 1) Inicializar -------------------------------------------------------
  training_set <- sample(training_set)
  n  <- length(stream_data)
  
  alphas <- numeric(n)
  p_vals <- numeric(n)
  S_vals <- numeric(n)
  G_vals <- numeric(n)
  Ratio  <- rep(NA_real_, n)
  
  S_prev <- 1
  S_vals[1] <- S_prev
  alarm   <- NA             # posición de la alarma
  
  ## --- 2) Bucle principal ---------------------------------------------------
  for (i in seq_len(n)) {
    
    ## 2.1  NCM y p-valor
    xi <- stream_data[[i]]
    alphas[i] <- do.call(
      non_conformity_measure,
      c(list(xi = xi, training_set = training_set), list(...))
    )
    greater   <- sum(alphas[1:i] >  alphas[i])
    equal     <- sum(alphas[1:i] == alphas[i])
    p_vals[i] <- (greater + runif(1) * equal) / i
    
    ## 2.2  Betting function cautelosa
    if (i <= W) {
      g_i  <- 1
      ratio <- NA
    } else {
      s_hist <- S_vals[1:(i - 1)]
      wvals  <- tail(s_hist[s_hist > 0], W)
      min_w  <- min(wvals)
      ratio  <- S_prev / min_w
      
      if (ratio <= epsilon) {
        g_i <- 1            # congelado
      } else {
        g_i <- do.call(
          betting_function,
          c(list(
            p_values = p_vals[1:(i - 1)],
            new_p    = p_vals[i],
            i        = i
          ), params_bf)
        )
      }
    }
    
    ## 2.3  Actualizar martingala y registrar
    G_vals[i] <- g_i
    Ratio[i]  <- ratio
    S_vals[i] <- S_prev * g_i
    S_prev    <- S_vals[i]
    
    ## 2.4  ¿Se superó el umbral de alarma?
    if (S_vals[i] > th) {
      alarm <- i
      break            # corta ejecución en el primer disparo
    }
  }
  
  ## --- 3) Recorte si hubo alarma -------------------------------------------
  if (!is.na(alarm)) {
    alphas <- alphas[1:alarm]
    p_vals <- p_vals[1:alarm]
    S_vals <- S_vals[1:alarm]
    G_vals <- G_vals[1:alarm]
    Ratio  <- Ratio [1:alarm]
  }
  
  ## --- 4) Salida ------------------------------------------------------------
  list(
    S_vals       = S_vals,
    p_vals       = p_vals,
    G_vals       = G_vals,
    Ratio        = Ratio,
    change_point = alarm   # NA si nunca hubo disparo
  )
}

#ICM multiple
ICM_multi <- function(training_set,
                      stream_data,
                      non_conformity_measure,
                      betting_function,
                      th = 1,
                      params_bf = list(),  # argumentos para la betting function
                      ...) {               # argumentos para la NCM
  
  training_set <- sample(training_set)
  n <- length(stream_data)
  
  # Listas para almacenar los resultados completos
  all_Cn     <- numeric()
  all_p_vals <- numeric()
  change_points <- c()
  
  # Variables de trabajo internas
  start_idx <- 1
  C_prev <- 0
  
  while (start_idx <= n) {
    alphas <- numeric(n - start_idx + 1)
    p_vals <- numeric(n - start_idx + 1)
    Cn     <- numeric(n - start_idx + 1)
    
    for (i in seq_along(alphas)) {
      xi <- stream_data[start_idx + i - 1]
      
      alphas[i] <- do.call(non_conformity_measure,
                           c(list(xi = xi, training_set = training_set), list(...)))
      
      greater <- sum(alphas[1:i] > alphas[i])
      equal   <- sum(alphas[1:i] == alphas[i])
      p_vals[i] <- (greater + runif(1) * equal) / i
      
      p_hist <- if (i == 1) numeric(0) else p_vals[1:(i - 1)]
      
      g_i <- do.call(betting_function,
                     c(list(p_values = p_hist, new_p = p_vals[i], i = i), params_bf))
      
      val <- if (g_i > 0) (C_prev + log(g_i)) else -Inf
      Cn[i] <- max(0, val)
      C_prev <- Cn[i]
      
      if (Cn[i] > th) {
        alarm_idx <- start_idx + i - 1
        change_points <- c(change_points, alarm_idx)
        
        # Guardar resultados hasta la alarma
        all_Cn     <- c(all_Cn, Cn[1:i])
        all_p_vals <- c(all_p_vals, p_vals[1:i])
        
        # Reiniciar después de la alarma
        start_idx <- alarm_idx + 1
        C_prev <- 0
        break
      }
    }
    
    # Si no hubo alarma hasta el final
    if (start_idx <= n && length(Cn) + start_idx - 1 == n) {
      all_Cn     <- c(all_Cn, Cn)
      all_p_vals <- c(all_p_vals, p_vals)
      break
    }
  }
  
  list(
    Cn            = all_Cn,
    p_vals        = all_p_vals,
    change_points = change_points
  )
}


ICM_multi_adaptive <- function(stream_data,
                                        non_conformity_measure,
                                        betting_function,
                                        th = 1,
                                        training_set = NULL,
                                        training_size = NULL,
                                        m_retrain = NULL,         # tamaño del training post-alarma
                                        guard_band = 0,           # saltar G obs tras la alarma antes de re-entrenar
                                        shuffle_training = TRUE,
                                        params_bf = list(),
                                        ...) {
  n <- length(stream_data)
  
  # --- Resolver entrenamiento inicial y offset global (para mapear a 'z') ---
  if (!is.null(training_set)) {
    if (shuffle_training) training_set <- sample(training_set)
    m0 <- length(training_set)
    initial_offset <- m0           # global = stream + m0, ya que z = training(1:m0) + stream(m0+1:...)
    pos <- 1L                      # cursor absoluto en el stream
  } else {
    if (is.null(training_size) || training_size <= 0)
      stop("Especifica training_size > 0 si no pasas un training_set.")
    m0 <- training_size
    if (m0 >= n) stop("training_size no puede ser >= length(stream_data).")
    training_set   <- stream_data[1:m0]
    if (shuffle_training) training_set <- sample(training_set)
    initial_offset <- 0L
    pos <- m0 + 1L
  }
  mR <- if (is.null(m_retrain)) m0 else m_retrain
  
  # --- Buffers alineados a n (NA donde se re-entrena) ---
  all_Cn     <- rep(NA_real_, n)
  all_p_vals <- rep(NA_real_, n)
  
  # --- Salidas de puntos de cambio y eventos (para auditar) ---
  cp_stream <- integer(0)
  events <- list()
  
  # --- Estado del tramo actual (se reinicia tras cada alarma) ---
  C_prev <- 0
  i_seg <- 0L
  alphas_seg <- numeric(0)
  p_seg <- numeric(0)
  Cn_seg <- numeric(0)
  seg_start_pos <- pos
  
  while (pos <= n) {
    # Avance dentro del tramo
    i_seg <- i_seg + 1L
    xi <- stream_data[pos]
    
    # 1) NCM respecto al training_set vigente
    alpha_i <- do.call(non_conformity_measure,
                       c(list(xi = xi, training_set = training_set), list(...)))
    alphas_seg[i_seg] <- alpha_i
    
    # 2) p-valor ICM por ranking dentro del tramo
    greater <- sum(alphas_seg[1:i_seg] >  alpha_i)
    equal   <- sum(alphas_seg[1:i_seg] == alpha_i)
    p_i <- (greater + runif(1) * equal) / i_seg
    p_seg[i_seg] <- p_i
    
    # 3) Betting function
    p_hist <- if (i_seg == 1L) numeric(0) else p_seg[1:(i_seg - 1L)]
    g_i <- do.call(betting_function,
                   c(list(p_values = p_hist, new_p = p_i, i = i_seg), params_bf))
    
    # 4) Log-martingala con truncación a 0 (CUSUM-like en log)
    val <- if (is.finite(g_i) && g_i > 0) (C_prev + log(g_i)) else -Inf
    C_curr <- max(0, val)
    Cn_seg[i_seg] <- C_curr
    C_prev <- C_curr
    
    # 5) ¿Alarma?
    if (C_curr > th) {
      alarm_pos <- pos                          # ÍNDICE EN EL STREAM (1..n)
      cp_stream <- c(cp_stream, alarm_pos)
      
      # Escribir tramo detectado en los buffers alineados
      all_Cn[seg_start_pos:alarm_pos]     <- Cn_seg[1:i_seg]
      all_p_vals[seg_start_pos:alarm_pos] <- p_seg[1:i_seg]
      
      # 6) Re-entrenar con guard band y m_retrain
      train_start <- min(n, alarm_pos + 1L + as.integer(guard_band))
      if (train_start > n) break
      train_end <- min(n, train_start + mR - 1L)
      new_train <- stream_data[train_start:train_end]
      if (length(new_train) < 1L) break
      if (shuffle_training) new_train <- sample(new_train)
      training_set <- new_train
      
      # Marcar NA en la ventana de re-entrenamiento
      all_Cn[train_start:train_end]     <- NA_real_
      all_p_vals[train_start:train_end] <- NA_real_
      
      # Registrar evento para auditoría
      events[[length(events) + 1L]] <- data.frame(
        alarm_stream = alarm_pos,
        alarm_global = initial_offset + alarm_pos,
        seg_start    = seg_start_pos,
        train_start  = train_start,
        train_end    = train_end,
        next_start   = train_end + 1L
      )
      
      # 7) Reiniciar tramo y avanzar cursor a 'next_start'
      pos <- train_end + 1L
      seg_start_pos <- pos
      C_prev <- 0
      i_seg <- 0L
      alphas_seg <- numeric(0)
      p_seg <- numeric(0)
      Cn_seg <- numeric(0)
      next
    }
    
    # Si no hubo alarma en este paso, avanzar al siguiente punto del stream
    pos <- pos + 1L
  }
  
  # Volcar lo que haya quedado del último tramo sin alarma
  if (i_seg > 0L) {
    end_pos <- min(n, seg_start_pos + i_seg - 1L)
    all_Cn[seg_start_pos:end_pos]     <- Cn_seg[1:(end_pos - seg_start_pos + 1L)]
    all_p_vals[seg_start_pos:end_pos] <- p_seg[1:(end_pos - seg_start_pos + 1L)]
  }
  
  # Tabla de eventos y chequeo del offset (auditoría de indexación)
  ev <- if (length(events)) do.call(rbind, events) else
    data.frame(alarm_stream=integer(), alarm_global=integer(),
               seg_start=integer(), train_start=integer(),
               train_end=integer(), next_start=integer())
  if (nrow(ev) > 0L) {
    off <- unique(ev$alarm_global - ev$alarm_stream)
    # Si diste training_set inicial, el offset debe ser m0:
    if (!is.null(training_set) || initial_offset > 0L) {
      stopifnot(length(off) == 1L)
    }
  }
  
  list(
    Cn_aligned           = all_Cn,
    p_vals_aligned       = all_p_vals,
    change_points_stream = cp_stream,
    change_points_global = initial_offset + cp_stream,
    events               = ev,
    offset               = initial_offset
  )
}
