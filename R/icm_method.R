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
  }
  
  change_point <- which(Cn > th)[1]
  
  list(
    Cn = Cn,
    p_vals = p_vals,
    change_point = change_point
  )
}


#ICM CAUTIOUS
#refer to this method in this article "A betting function for addresing concept drift with conformal martingales"

ICM_CBF <- function(training_set,
                    stream_data,
                    non_conformity_measure,
                    betting_function,
                    W        = 100,   # tamaño de la “ventana caliente”
                    epsilon  = 1.5,   # umbral de descongelado (suele usarse 1.2–2)
                    params_bf = list(),   # argumentos extra para la betting-function
                    ...) {             # argumentos extra para la NCM
  
  ## 1) Inicializar -----------------------------------------------------------
  training_set <- sample(training_set)     # baraja set de entrenamiento
  n  <- length(stream_data)                # longitud del stream
  
  alphas <- numeric(n)
  p_vals <- numeric(n)
  S_vals <- numeric(n)
  G_vals <- numeric(n)
  Ratio  <- numeric(n)
  change_point <- NA                       # (opcional) no se calcula aquí
  
  S_prev   <- 1                            # martingala arranca en 1
  S_vals[1] <- S_prev
  
  ## 2) Bucle principal -------------------------------------------------------
  for (i in seq_len(n)) {
    
    # 2.1) Medida de no-conformidad y p-valor
    xi <- stream_data[[i]]
    alphas[i] <- do.call(
      non_conformity_measure,
      c(list(xi = xi, training_set = training_set), list(...))
    )
    greater   <- sum(alphas[1:i] >  alphas[i])
    equal     <- sum(alphas[1:i] == alphas[i])
    p_vals[i] <- (greater + runif(1) * equal) / i
    
    # 2.2) Betting function cautelosa
    if (i <= W) {
      g_i  <- 1
      ratio <- NA
    } else {
      s_hist <- S_vals[1:(i-1)]
      wvals  <- tail(s_hist[s_hist > 0], W)   # últimos W valores positivos
      min_w  <- min(wvals)
      
      ratio <- S_prev / min_w                 # ← COCIENTE corregido
      
      if (ratio <= epsilon) {
        g_i <- 1                              # “congelado”
      } else {
        g_i <- do.call(
          betting_function,
          c(list(
            p_values = p_vals[1:(i-1)],
            new_p    = p_vals[i],
            i        = i
          ), params_bf)
        )
      }
    }
    
    # 2.3) Actualizar y almacenar
    G_vals[i] <- g_i
    Ratio[i]  <- ratio
    S_vals[i] <- S_prev * g_i
    S_prev    <- S_vals[i]
  }
  
  ## 3) Salida ----------------------------------------------------------------
  list(
    S_vals       = S_vals,
    p_vals       = p_vals,
    G_vals       = G_vals,
    Ratio        = Ratio,
    change_point = change_point
  )
}