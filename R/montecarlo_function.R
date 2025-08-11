library(future.apply)
#montecarlo method for treshold simulations

montecarlo_ICM <- function(n_sim        = 200,
                              h_vals       = seq(1, 6, 0.5),
                              theta_stream = 100,
                              mu1          = 1,
                              m            = 200,
                              ncm_fun      = Non_conformity_KNN,
                              bet_fun      = Constant_BF,
                              k            = NULL,
                              params_bf    = list(),
                              n_stream     = 1000) {
  
  if (is.null(k)) k <- 7
  pos_change <- theta_stream
  
  simulate_once <- function(sim_id) {
    training_set <- rnorm(m, 0, 1)
    stream <- c(rnorm(theta_stream - 1, 0, 1),
                rnorm(n_stream - theta_stream + 1, mu1, 1))
    
    fa <- numeric(length(h_vals))
    dl <- rep(NA_real_, length(h_vals))
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      out <- ICM(training_set          = training_set,
                 stream_data           = stream,
                 non_conformity_measure = ncm_fun,
                 betting_function      = bet_fun,
                 params_bf             = params_bf,
                 k                     = k,
                 th                    = h)
      tau <- out$change_point
      if (is.na(tau)) {
        fa[j] <- 0
      } else if (tau < pos_change) {
        fa[j] <- 1
      } else {
        dl[j] <- tau - pos_change
      }
    }
    list(fa = fa, delay = dl)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(do.call(rbind, lapply(sims, `[[`, "fa"))),
    mean_delay    = colMeans(do.call(rbind, lapply(sims, `[[`, "delay")),
                             na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
}

#Monte carlo para oraculos
montecarlo_oraculo_par <- function(n_sim        = 200,
                                   h_vals       = seq(1, 6, 0.5),
                                   theta_stream = 100,
                                   mu1          = 1,
                                   m            = 200,      # se ignora, solo por simetría
                                   n_stream     = 1000,
                                   detector_fn,
                                   detector_label,
                                   ...) {
  
  pos_change <- theta_stream
  
  simulate_once <- function(sim_id) {
    stream <- c(rnorm(theta_stream - 1, 0, 1),
                rnorm(n_stream - theta_stream + 1, mu1, 1))
    
    fa <- numeric(length(h_vals))
    dl <- rep(NA_real_, length(h_vals))
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      res <- detector_fn(z = stream, threshold = h, ...)
      tau <- res$tau
      
      if (is.na(tau)) {
        fa[j] <- 0
      } else if (tau < pos_change) {
        fa[j] <- 1
      } else {
        dl[j] <- tau - pos_change
      }
    }
    list(fa = fa, delay = dl)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  
  data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(FALSE_ALARMS),
    mean_delay    = colMeans(DELAYS, na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1,
    Method        = detector_label
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
}

#Montecarlo para ICM_CBF
montecarlo_ICM_CBF_par <- function(n_sim        = 200,
                                   h_vals       = seq(1, 6, 0.5),
                                   theta_stream = 100,
                                   mu1          = 1,
                                   m            = 200,
                                   ncm_fun      = Non_conformity_KNN,
                                   bet_fun      = Constant_BF,
                                   k            = NULL,
                                   W            = 100,
                                   epsilon      = 0.01,
                                   params_bf    = list(),
                                   n_stream     = 1000) {
  
  if (is.null(k)) k <- 7
  pos_change <- theta_stream
  
  simulate_once <- function(sim_id) {
    ## Datos de la réplica -----------------------------------------------
    training_set <- rnorm(m, 0, 1)
    stream <- c(rnorm(theta_stream - 1, 0, 1),
                rnorm(n_stream - theta_stream + 1, mu1, 1))
    
    fa <- numeric(length(h_vals))
    dl <- rep(NA_real_, length(h_vals))
    
    ## Bucle de umbrales ---------------------------------------------------
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      
      out <- ICM_CBF(training_set          = training_set,
                     stream_data           = stream,
                     non_conformity_measure = ncm_fun,
                     betting_function      = bet_fun,
                     W                     = W,
                     epsilon               = epsilon,
                     params_bf             = params_bf,
                     k                     = k,
                     th                    = h)
      
      tau <- out$change_point   # o out$tau, según tu implementación
      
      if (is.na(tau)) {
        fa[j] <- 0
      } else if (tau < pos_change) {
        fa[j] <- 1
      } else {
        dl[j] <- tau - pos_change
      }
    }
    list(fa = fa, delay = dl)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  
  data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(FALSE_ALARMS),
    mean_delay    = colMeans(DELAYS, na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
}

