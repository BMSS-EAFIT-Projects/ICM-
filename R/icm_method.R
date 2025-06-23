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
                    W = 100,
                    epsilon = 100,
                    params_bf = list(),  # ← aquí van los argumentos para la BF
                    ...) {               # ← aquí van solo los argumentos para la NCM
  
  training_set <- sample(training_set)
  m <- length(training_set)
  n <- length(stream_data)
  
  alphas <- numeric(n)
  p_vals <- numeric(n)
  S_vals <- numeric(n)
  S_prev <- 1
  
  for (i in seq_len(n)) {
    xi <- stream_data[i]
    
    alphas[i] <- do.call(non_conformity_measure, c(list(xi = xi, training_set = training_set), list(...)))
    
    greater <- sum(alphas[1:i] > alphas[i])
    equal   <- sum(alphas[1:i] == alphas[i])
    p_vals[i] <- (greater + runif(1) * equal) / i
    
    p_hist <- if (i == 1) numeric(0) else p_vals[1:(i - 1)]
    s_hist <- if (i == 1) numeric(0) else S_vals[1:(i - 1)]
    
    g_i <- if (i <= W || length(s_hist) < W || 
               (s_hist[length(s_hist)] / min(tail(s_hist, W))) <= epsilon) {
      1
    } else {
      do.call(betting_function, c(list(
        p_values = p_hist,
        new_p = p_vals[i],
        i = i
      ), params_bf))
    }
    
    S_vals[i] <- S_prev * g_i
    S_prev <- S_vals[i]
  }
  
  threshold <- 20
  change_point <- which(S_vals > threshold)[1]
  
  list(
    S_vals = S_vals,
    p_vals = p_vals,
    change_point = change_point
  )
}
