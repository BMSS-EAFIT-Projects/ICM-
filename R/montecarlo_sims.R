#montecarlo method for treshold simulations

montecarlo_ICM <- function(n_sim           = 200,
                           h_vals          = seq(1, 6, 0.5),
                           theta_stream    = 100,
                           mu1             = 1,
                           m               = 200,
                           ncm_fun         = Non_conformity_KNN,
                           bet_fun         = Constant_BF,
                           k               = NULL) {
  # Si k no se especifica, lo ponemos = ceiling(m/2)
  if (is.null(k)) k <- ceiling(m / 2)
  
  # Longitud del “stream” que ve ICM
  n_stream <- 1000
  
  # Calcular epsilon_global = posición global del cambio dentro de la serie completa
  # (sólo para referencia interna; en el detector realmente usamos pos_change = theta_stream).
  epsilon_global <- m + theta_stream
  
  # Preasignar matrices
  delays       <- matrix(NA, nrow = n_sim, ncol = length(h_vals))
  false_alarms <- matrix(0,  nrow = n_sim, ncol = length(h_vals))
  
  for (sim in seq_len(n_sim)) {
    # --------------------------------------------------
    # (1) Conjunto de entrenamiento de tamaño m (todos N(0,1))
    # --------------------------------------------------
    training_set <- rnorm(m, mean = 0, sd = 1)
    
    # --------------------------------------------------
    # (2) Generar “stream_data” de longitud n_stream con CP en theta_stream
    # --------------------------------------------------
    stream <- numeric(n_stream)
    for (i in seq_len(n_stream)) {
      # Índice global si lo contáramos desde el inicio de la serie completa sería m + i
      if (i < theta_stream) {
        # Antes del cambio (n < θ_stream) → N(0,1)
        stream[i] <- rnorm(1, mean = 0, sd = 1)
      } else {
        # i >= θ_stream → post-cambio → N(μ₁, 1)
        stream[i] <- rnorm(1, mean = mu1, sd = 1)
      }
    }
    
    # --------------------------------------------------
    # (3) Ejecutar ICM genérico con la función de apuesta “bet_fun”
    # --------------------------------------------------
    out_ICM <- ICM(
      training_set          = training_set,
      stream_data           = stream,
      non_conformity_measure = ncm_fun,
      betting_function      = bet_fun,
      k                     = k
    )
    Cn_vals <- out_ICM$Cn
    
    # --------------------------------------------------
    # (4) Para cada umbral h, hallamos τ y decidimos FA o delay
    # --------------------------------------------------
    pos_change <- theta_stream  # el cambio ocurre en el índice theta_stream dentro de “stream”
    
    for (j in seq_along(h_vals)) {
      h   <- h_vals[j]
      tau <- which(Cn_vals > h)[1]  # primer i tal que Cn_vals[i] > h
      
      if (is.na(tau)) {
        # nunca detectó → no hay FA y no hay delay posterior
        false_alarms[sim, j] <- 0
        delays[sim, j]       <- NA
      } else if (tau < pos_change) {
        # detectó antes del cambio → falsa alarma
        false_alarms[sim, j] <- 1
        delays[sim, j]       <- NA
      } else {
        # detectó en tau ≥ pos_change → retardo = tau - pos_change
        false_alarms[sim, j] <- 0
        delays[sim, j]       <- tau - pos_change
      }
    }
  }  # fin bucle simulaciones
  
  # --------------------------------------------------
  # (5) Calcular P_FA(h) y MeanDelay(h)
  # --------------------------------------------------
  p_FA       <- colMeans(false_alarms)                  # fracción de FA para cada h
  mean_delay <- colMeans(delays, na.rm = TRUE)          # retardo medio cond. no-FA
  
  # Empaquetar en data.frame
  df_res <- data.frame(
    theta_stream    = theta_stream,
    mu1             = mu1,
    threshold       = h_vals,
    p_false_alarm   = p_FA,
    mean_delay      = mean_delay
  ) %>%
    mutate(log_delay = log10(1 + mean_delay))
  
  return(df_res)
}