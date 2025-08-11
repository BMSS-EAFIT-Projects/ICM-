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
