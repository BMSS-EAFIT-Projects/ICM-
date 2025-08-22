library(future.apply)
#montecarlo method for treshold simulations

source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/montecarlo_helpers.R")

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
    
    # (1) alphas y p una sola vez
    alphas <- alphas_from_ncm(stream, training_set, ncm_fun, k = k)
    p <- pvals_from_alphas(alphas)
    
    # (2) una sola trayectoria S_n para ESTA BF
    S <- icm_C_path_from_p(p, bet_fun, params_bf = params_bf)
    
    # (3) primeros cruces para TODOS los umbrales
    tau <- first_cross_indices_linear(S, h_vals)
    
    # (4) métricas por umbral
    fa <- as.integer(!is.na(tau) & tau < theta_stream)
    dl <- ifelse(is.na(tau) | tau < theta_stream, NA_real_, tau - theta_stream)
    
    list(fa = fa, delay = dl, tau = tau)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  TAUS <- do.call(rbind, lapply(sims, `[[`, "tau"))
  
  taus_wide <- as.data.frame(TAUS)
  colnames(taus_wide) <- paste0("h_", seq_along(h_vals))
  taus_wide$sim_id <- seq_len(nrow(taus_wide))
  
  # A formato largo y con el valor numérico de umbral
  taus_long <- tidyr::pivot_longer(
    taus_wide,
    cols = dplyr::starts_with("h_"),
    names_to = "h_idx",
    values_to = "tau"
  ) |>
    dplyr::mutate(
      threshold_idx = as.integer(sub("h_", "", h_idx)),
      threshold = h_vals[threshold_idx]
    ) |>
    dplyr::select(sim_id, threshold, tau)
  
  summary_df <- data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(FALSE_ALARMS),
    mean_delay    = colMeans(DELAYS, na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
  
  return(list(summary = summary_df, taus = taus_long))
}

#Monte carlo para oraculos
montecarlo_oraculo <- function(n_sim        = 200,
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
    tau_vec <- rep(NA_integer_, length(h_vals))
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      res <- detector_fn(z = stream, threshold = h, ...)
      tau <- res$tau
      tau_vec[j]<- tau
      
      if (is.na(tau)) {
        fa[j] <- 0
      } else if (tau < pos_change) {
        fa[j] <- 1
      } else {
        dl[j] <- tau - pos_change
      }
    }
    list(fa = fa, delay = dl, tau = tau_vec)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  TAUS <- do.call(rbind, lapply(sims, `[[`, "tau"))
  
  taus_wide <- as.data.frame(TAUS)
  colnames(taus_wide) <- paste0("h_", seq_along(h_vals))
  taus_wide$sim_id <- seq_len(nrow(taus_wide))
  
  # A formato largo y con el valor numérico de umbral
  taus_long <- tidyr::pivot_longer(
    taus_wide,
    cols = dplyr::starts_with("h_"),
    names_to = "h_idx",
    values_to = "tau"
  ) |>
    dplyr::mutate(
      threshold_idx = as.integer(sub("h_", "", h_idx)),
      threshold = h_vals[threshold_idx]
    ) |>
    dplyr::select(sim_id, threshold, tau)
  
  summary_df <- data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(FALSE_ALARMS),
    mean_delay    = colMeans(DELAYS, na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1,
    Method        = detector_label
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
  
  return(list(summary = summary_df, taus = taus_long))
}

#Montecarlo para ICM_CBF
montecarlo_ICM_CBF <- function(n_sim        = 200,
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
    training_set <- rnorm(m, 0, 1)
    stream <- c(rnorm(theta_stream - 1, 0, 1),
                rnorm(n_stream - theta_stream + 1, mu1, 1))
    
    # (1) alphas y p una sola vez
    alphas <- alphas_from_ncm(stream, training_set, ncm_fun, k = k)
    p <- pvals_from_alphas(alphas)
    
    # (2) una sola trayectoria S_n con CBF (sin log)
    S <- icm_cbf_S_path_from_p(p, bet_fun,
                               W = W, epsilon = epsilon,
                               params_bf = params_bf)
    
    # (3) primeros cruces para TODOS los umbrales
    tau <- first_cross_indices_linear(S, h_vals)
    
    # (4) métricas
    fa <- as.integer(!is.na(tau) & tau < theta_stream)
    dl <- ifelse(is.na(tau) | tau < theta_stream, NA_real_, tau - theta_stream)
    
    list(fa = fa, delay = dl, tau = tau)
  }
  
  sims <- future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  TAUS <- do.call(rbind, lapply(sims, `[[`, "tau"))
  
  taus_wide <- as.data.frame(TAUS)
  colnames(taus_wide) <- paste0("h_", seq_along(h_vals))
  taus_wide$sim_id <- seq_len(nrow(taus_wide))
  
  # A formato largo y con el valor numérico de umbral
  taus_long <- tidyr::pivot_longer(
    taus_wide,
    cols = dplyr::starts_with("h_"),
    names_to = "h_idx",
    values_to = "tau"
  ) |>
    dplyr::mutate(
      threshold_idx = as.integer(sub("h_", "", h_idx)),
      threshold = h_vals[threshold_idx]
    ) |>
    dplyr::select(sim_id, threshold, tau)
  
  summary_df <- data.frame(
    threshold     = h_vals,
    p_false_alarm = colMeans(FALSE_ALARMS),
    mean_delay    = colMeans(DELAYS, na.rm = TRUE),
    theta_stream  = theta_stream,
    mu1           = mu1
  ) |>
    dplyr::mutate(log_delay = log10(1 + mean_delay))
  
  return(list(summary = summary_df, taus = taus_long))
}

