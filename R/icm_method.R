#Basic method of ICM

#icm method
ICM <- function(training_set,
                stream_data,
                non_conformity_measure,
                betting_function,
                th = 1,
                ...) {
  # 1) Mezclamos el training_set para garantizar exchangeability
  training_set <- sample(training_set)
  
  m <- length(training_set)
  n <- length(stream_data)
  
  # Preasignamos vectores
  alphas    <- numeric(n)
  p_vals    <- numeric(n)
  Cn        <- numeric(n)
  
  C_prev <- 0  # C_0 = 0
  
  for (i in seq_len(n)) {
    # --- (a) Non-conformity score ω_i con respecto a 'training_set' fijo
    xi <- stream_data[i]
    alphas[i] <- do.call(non_conformity_measure, c(list(xi = xi, training_set = training_set), list(...)))
    
    # --- (b) p-value p_i según (1) en función de alphas[1:i]
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    p_vals[i] <- (greater + runif(1) * equal) / i
    
    # --- (c) betting value g_i = g(p_i)
    g_i <- betting_function(p_vals[i])
    # Si g_i = 0 genera ln(0) = -Inf; en tal caso definimos ln(g_i)= -Inf → C_i = 0 via max
    
    # --- (d) actualización de C_n: C_n = max{0, C_{n-1} + ln(g_i)}
    val <- if (g_i > 0) (C_prev + log(g_i)) else -Inf
    Cn[i] <- max(0, val)
    C_prev <- Cn[i]
  }
  
  change_point = which(Cn>th)[1]
  # Devolvemos el vector Cn y los p_vals (en caso de que quieras inspeccionarlos)
  list(
    Cn      = Cn,
    p_vals  = p_vals,
    change_point = change_point
  )
}
