library(future.apply)
library(data.table)
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



montecarlo_ICM_MULTI <- function(
    n_sim        = 200,
    h_vals       = seq(1, 6, 0.5),
    m            = 200,
    n_stream     = 1000,
    theta_vec    = c(300, 700),
    mu_levels    = c(0, 1.5, 0),
    ncm_fun      = Non_conformity_KNN,
    bet_fun      = Constant_BF,
    params_bf    = list(),
    k            = NULL,
    window_mode  = c("abs","frac"),
    window_abs   = Inf,
    window_frac  = 1.0
){
  if (is.null(k)) k <- 7
  window_mode <- match.arg(window_mode)
  J <- length(theta_vec)
  
  one_sim <- function(sim_id){
    # 1) stream + entrenamiento (una sola vez)
    training_set <- rnorm(m, 0, 1)
    gen <- make_stream_mean_shifts(n_stream, theta_vec, mu_levels, sd = 1)
    
    # 2) alphas UNA vez (para no recalcular NCM por cada h)
    alphas_full <- alphas_from_ncm(gen$stream, training_set,
                                   ncm_fun = ncm_fun, k = k)
    
    per_change_list <- vector("list", length(h_vals))
    alarms_list     <- vector("list", length(h_vals))
    stats_mat <- matrix(0, nrow = length(h_vals), ncol = 4) # TP, MISS, FA, OSEG
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      
      # *** Igual que antes: ICM_multi_fast con umbral ESCALAR (no NULL) ***
      out <- ICM_multi_fast(
        training_set           = training_set,
        stream_data            = gen$stream,
        non_conformity_measure = ncm_fun,
        betting_function       = bet_fun,
        th                     = h,                # <-- ya NO es NULL
        params_bf              = params_bf,
        k                      = k,
        alphas_full            = alphas_full       # <-- reutilizamos
      )
      
      eva <- match_alarms_to_changes_multi(
        alarms       = as.integer(out$change_points),
        true_changes = gen$true_changes,
        n_stream     = n_stream,
        window_mode  = window_mode,
        window_abs   = window_abs,
        window_frac  = window_frac
      )
      stats_mat[j,] <- c(eva$tp_count, eva$miss_count, eva$fa_count, eva$extras_after_tp)
      per_change_list[[j]] <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$per_change)
      alarms_list[[j]]     <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$alarms_df)
    }
    
    list(
      stats_mat  = stats_mat,
      per_change = dplyr::bind_rows(per_change_list),
      alarms     = dplyr::bind_rows(alarms_list)
    )
  }
  
  sims <- future.apply::future_lapply(seq_len(n_sim), one_sim, future.seed = TRUE)
  
  # apilar matrices TP/MISS/FA/OSEG
  stats_arr <- simplify2array(lapply(sims, `[[`, "stats_mat"))  # [length(h_vals) x 4 x n_sim]
  stats_sum <- apply(stats_arr, c(1,2), sum)
  colnames(stats_sum) <- c("TP","MISS","FA","OSEG")
  TP   <- stats_sum[,"TP"]
  MISS <- stats_sum[,"MISS"]
  FA   <- stats_sum[,"FA"]
  OSEG <- stats_sum[,"OSEG"]
  
  # Delays: sacarlos de per_change (promedios por h)
  per_change_all <- rbindlist(lapply(sims, `[[`, "per_change"))
  # resumen por h y cambio
  per_change_summ <- per_change_all[, .(
    seg_len          = mean(seg_len),
    recall_j         = mean(tp),
    miss_rate_j      = 1 - mean(tp),
    mean_delay_j     = mean(delay, na.rm = TRUE),
    median_delay_j   = median(delay, na.rm = TRUE),
    iqr_delay_j      = IQR(delay, na.rm = TRUE),
    norm_delay_j     = mean(delay / seg_len, na.rm = TRUE),
    extras_after_tp_j= mean(extras_after_tp),
    fa_in_segment_j  = mean(fa_in_segment)
  ),
  by = .(threshold, change_id)
  ]
  
  # mean_delay global por h (promediando delays válidos por simulación/cambio)
  delay_by_h <- per_change_all[, .(mean_delay = mean(delay, na.rm = TRUE)), by = .(threshold)]
  
  # summary por h
  dt_h <- data.table(threshold = h_vals)
  dt_h[, recall := TP / (J * n_sim)]
  dt_h[, miss_rate := 1 - recall]
  dt_h[, precision := ifelse(TP + FA > 0, TP / (TP + FA), NA_real_)]
  dt_h[, fa_per_1000 := FA / (n_stream * n_sim) * 1000]
  dt_h[, alarms_per_1000 := (TP + FA) / (n_stream * n_sim) * 1000]
  dt_h[, overseg_ratio := ifelse(TP > 0, OSEG / TP, NA_real_)]
  
  dt_h <- merge(dt_h, as.data.table(delay_by_h), by = "threshold", all.x = TRUE)
  dt_h[, log_delay := log10(1 + mean_delay)]
  
  list(
    summary    = as_tibble(dt_h),
    per_change = as_tibble(per_change_summ),
    alarms     = as_tibble(rbindlist(lapply(sims, `[[`, "alarms")))
  )
}


# 4) Wrapper: Monte Carlo para ICM_multi_adaptive (con reentrenos) -------------
montecarlo_ICM_MULTI_ADAPTIVE <- function(
    n_sim        = 200,
    h_vals       = seq(1, 6, 0.5),
    m            = 200,
    n_stream     = 1000,
    theta_vec    = c(300, 700),
    mu_levels    = c(0, 1.5, 0),
    ncm_fun      = Non_conformity_KNN,
    bet_fun      = Constant_BF,
    params_bf    = list(),
    k            = NULL,
    m_retrain    = NULL,
    guard_band   = 0L,
    window_mode  = c("abs","frac"),
    window_abs   = Inf,
    window_frac  = 1.0
){
  if (is.null(k)) k <- 7
  if (is.null(m_retrain)) m_retrain <- m
  window_mode <- match.arg(window_mode)
  J <- length(theta_vec)
  
  one_sim <- function(sim_id){
    gen <- make_stream_mean_shifts(n_stream, theta_vec, mu_levels, sd = 1)
    training_set0 <- rnorm(m, 0, 1)
    
    per_change_list <- vector("list", length(h_vals))
    alarms_list     <- vector("list", length(h_vals))
    stats_mat <- matrix(0, nrow = length(h_vals), ncol = 4)
    
    for (j in seq_along(h_vals)) {
      h <- h_vals[j]
      out <- ICM_multi_adaptive_fast(
        stream_data            = gen$stream,
        non_conformity_measure = ncm_fun,
        betting_function       = bet_fun,
        th                     = h,
        training_set           = training_set0,
        training_size          = NULL,
        m_retrain              = m_retrain,
        guard_band             = guard_band,
        shuffle_training       = TRUE,
        params_bf              = params_bf,
        k                      = k
      )
      eva <- match_alarms_to_changes_multi(
        alarms       = as.integer(out$change_points_stream),
        true_changes = gen$true_changes,
        n_stream     = n_stream,
        window_mode  = window_mode,
        window_abs   = window_abs,
        window_frac  = window_frac
      )
      stats_mat[j,] <- c(eva$tp_count, eva$miss_count, eva$fa_count, eva$extras_after_tp)
      per_change_list[[j]] <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$per_change)
      alarms_list[[j]]     <- tibble(sim_id = sim_id, threshold = h) |> dplyr::bind_cols(eva$alarms_df)
    }
    
    list(
      stats_mat = stats_mat,
      per_change = dplyr::bind_rows(per_change_list),
      alarms     = dplyr::bind_rows(alarms_list)
    )
  }
  
  sims <- future.apply::future_lapply(seq_len(n_sim), one_sim, future.seed = TRUE)
  
  # Agregación
  stats_arr <- simplify2array(lapply(sims, `[[`, "stats_mat"))
  stats_sum <- apply(stats_arr, c(1,2), sum)
  colnames(stats_sum) <- c("TP","MISS","FA","OSEG")
  TP   <- stats_sum[,"TP"]; MISS <- stats_sum[,"MISS"]; FA <- stats_sum[,"FA"]; OSEG <- stats_sum[,"OSEG"]
  
  per_change_all <- rbindlist(lapply(sims, `[[`, "per_change"))
  per_change_summ <- per_change_all[, .(
    seg_len          = mean(seg_len),
    recall_j         = mean(tp),
    miss_rate_j      = 1 - mean(tp),
    mean_delay_j     = mean(delay, na.rm = TRUE),
    median_delay_j   = median(delay, na.rm = TRUE),
    iqr_delay_j      = IQR(delay, na.rm = TRUE),
    norm_delay_j     = mean(delay / seg_len, na.rm = TRUE),
    extras_after_tp_j= mean(extras_after_tp),
    fa_in_segment_j  = mean(fa_in_segment)
  ),
  by = .(threshold, change_id)
  ]
  delay_by_h <- per_change_all[, .(mean_delay = mean(delay, na.rm = TRUE)), by = .(threshold)]
  
  dt_h <- data.table(threshold = h_vals)
  dt_h[, recall := TP / (J * n_sim)]
  dt_h[, miss_rate := 1 - recall]
  dt_h[, precision := ifelse(TP + FA > 0, TP / (TP + FA), NA_real_)]
  dt_h[, fa_per_1000 := FA / (n_stream * n_sim) * 1000]
  dt_h[, alarms_per_1000 := (TP + FA) / (n_stream * n_sim) * 1000]
  dt_h[, overseg_ratio := ifelse(TP > 0, OSEG / TP, NA_real_)]
  
  dt_h <- merge(dt_h, as.data.table(delay_by_h), by = "threshold", all.x = TRUE)
  dt_h[, log_delay := log10(1 + mean_delay)]
  
  list(
    summary    = as_tibble(dt_h),
    per_change = as_tibble(per_change_summ),
    alarms     = as_tibble(rbindlist(lapply(sims, `[[`, "alarms")))
  )
}

#----------------------montecarlo outliers--------------------------------------
montecarlo_ICM_contaminado <- function(
    n_sim        = 200,
    h_vals       = seq(1, 6, 0.5),
    theta_stream = 100,
    mu1          = 1,
    m            = 200,
    ncm_fun      = Non_conformity_KNN,
    bet_fun      = Constant_BF,
    k            = NULL,
    params_bf    = list(),
    n_stream     = 1000,
    # --- Contaminación entrenamiento (SOLO shift / scale) ---------------------
    train_contam_rate   = 0.05,
    train_contam_model  = c("shift","scale"),
    train_contam_params = list(delta = 6, lambda = 5),
    # --- Contaminación stream (pre / post θ) ---------------------------------
    pre_contam_rate     = 0.02,
    pre_contam_model    = c("shift","scale"),
    pre_contam_params   = list(delta = 6, lambda = 5),
    post_contam_rate    = 0.02,
    post_contam_model   = c("shift","scale"),
    post_contam_params  = list(delta = 6, lambda = 5)
){
  `%||%` <- function(a,b) if (is.null(a)) b else a
  if (is.null(k)) k <- 7L
  
  train_contam_model <- match.arg(train_contam_model)
  pre_contam_model   <- match.arg(pre_contam_model)
  post_contam_model  <- match.arg(post_contam_model)
  
  # -------- simulación de UNA réplica --------------------------------------
  simulate_once <- function(sim_id){
    # 1) Entrenamiento base + contaminación
    training_set <- rnorm_base(m, 0, 1)
    training_set <- contaminate(training_set, train_contam_rate, train_contam_model,
                                train_contam_params, base_mu = 0, base_sd = 1)
    
    # 2) Stream con cambio súbito en θ
    n1 <- theta_stream - 1L
    n2 <- n_stream - theta_stream + 1L
    
    pre  <- rnorm_base(n1, 0, 1)
    pre  <- contaminate(pre,  pre_contam_rate,  pre_contam_model,  pre_contam_params,
                        base_mu = 0,  base_sd = 1)
    
    post <- rnorm_base(n2, mu1, 1)
    post <- contaminate(post, post_contam_rate, post_contam_model, post_contam_params,
                        base_mu = mu1, base_sd = 1)
    
    stream <- c(pre, post)
    
    # 3) alphas y p
    alphas <- alphas_from_ncm(stream, training_set, ncm_fun, k = k)
    p <- pvals_from_alphas(alphas)
    
    # 4) Trayectoria martingala y cruces
    S <- icm_C_path_from_p(p, bet_fun, params_bf = params_bf)
    tau <- first_cross_indices_linear(S, h_vals)
    
    # 5) Métricas
    fa <- as.integer(!is.na(tau) & tau < theta_stream)
    dl <- ifelse(is.na(tau) | tau < theta_stream, NA_real_, tau - theta_stream)
    
    list(fa = fa, delay = dl, tau = tau)
  }
  
  sims <- future.apply::future_lapply(seq_len(n_sim), simulate_once, future.seed = TRUE)
  
  FALSE_ALARMS <- do.call(rbind, lapply(sims, `[[`, "fa"))
  DELAYS       <- do.call(rbind, lapply(sims, `[[`, "delay"))
  TAUS         <- do.call(rbind, lapply(sims, `[[`, "tau"))
  
  taus_wide <- as.data.frame(TAUS)
  colnames(taus_wide) <- paste0("h_", seq_along(h_vals))
  taus_wide$sim_id <- seq_len(nrow(taus_wide))
  
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
  
  list(summary = summary_df, taus = taus_long)
}

