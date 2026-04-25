source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(janitor)
library(ggplot2)
library(readxl)

wide_to_long <- function(df, estacion){
  names(df)[1] <- "anio"
  
  df %>%
    pivot_longer(
      cols = 2:13,
      names_to = "mes",
      values_to = "valor"
    ) %>%
    mutate(
      mes   = as.integer(gsub("^[xX]", "", mes)),
      valor = as.numeric(valor),
      fecha = as.Date(sprintf("%d-%02d-01", anio, mes)),
      estacion = estacion
    ) %>%
    relocate(estacion, anio, mes, fecha, valor)
}

path <- "datos_caudales/series_mensuales_por_estacion.xlsx"

sheets <- excel_sheets(path)

dfs_est <- setNames(lapply(sheets, function(sh){
  read_excel(path, sheet = sh)
}), sheets)

dfs_est_long <- imap(dfs_est, wide_to_long)


na_por_estacion <- imap_dfr(dfs_est_long, ~ tibble(
  estacion = .y,
  n_filas  = nrow(.x),
  n_na     = sum(is.na(.x$valor)),
  n_nan    = sum(is.nan(as.numeric(.x$valor))),
  n_na_nan = sum(is.na(.x$valor) | is.nan(as.numeric(.x$valor)))
)) %>%
  arrange(desc(n_na_nan))


clean_impute_station <- function(df){
  df <- df %>%
    arrange(fecha) %>%
    mutate(
      valor = as.numeric(valor),
      valor = ifelse(is.nan(valor), NA_real_, valor)
    )
  
  ok <- which(!is.na(df$valor))
  if (length(ok) == 0) return(df %>% slice(0))
  
  df <- df[min(ok):max(ok), , drop = FALSE]
  
  med_est <- median(df$valor, na.rm = TRUE)
  
  df %>%
    group_by(anio) %>%
    mutate(
      med_anual = median(valor, na.rm = TRUE),
      med_anual = ifelse(is.na(med_anual) | is.nan(med_anual), med_est, med_anual),
      valor     = ifelse(is.na(valor), med_anual, valor)
    ) %>%
    ungroup() %>%
    select(-med_anual)
}

dfs_est_filled <- map(dfs_est_long, clean_impute_station)

# Series de tiempo
for (i in sheets) {
  path_plot <- file.path("Series_de_tiempo_EH", paste0(i, ".", "png"))
  print(path_plot)
  df <- dfs_est_filled[[i]]
  p <- ggplot(df, aes(x=fecha, y=valor)) +
    geom_line() +
    theme_minimal()
  ggsave(path_plot,
         plot =p,
         bg = "white")
}


# Series de Tiempo con CP
cp_rows_wide <- list()

for (i in sheets) {
  path_plot <- file.path("Series_de_tiempo_CPEH", paste0(i, ".", "png"))
  df <- dfs_est_filled[[i]]
  
  n <- length(df$valor)
  ntrain <- round(n * 0.01, 0)
  
  train_Set   <- df$valor[1:ntrain]
  stream_flow <- df$valor[(ntrain + 1):n]
  
  res <- ICM_multi_adaptive(
    stream_flow, Non_conformity_MAD, Kernel_BF, th = 3, train_Set,
    m_retrain = ntrain/2, guard_band = 4,
    params_bf = list(L=12, n_grid=512, bw_floor=0.02, min_history=4)
  )
  
  cp_idx <- sort(unique(res$change_points + ntrain))
  cp_idx <- cp_idx[cp_idx >= 1 & cp_idx <= nrow(df)]
  cp_dates <- df$fecha[cp_idx]
  
  # --------- Guardar CPs en dataframes (wide y long) ---------
  # Wide: una fila por i, columnas CP1..CPk
  if (length(cp_dates) == 0) {
    wide_row <- tibble(i = i)  # sin CPs
  } else {
    wide_row <- tibble(i = i) %>%
      bind_cols(as_tibble_row(setNames(as.list(cp_dates), paste0("CP", seq_along(cp_dates)))))
  }
  cp_rows_wide[[as.character(i)]] <- wide_row
  
  starts <- c(1, cp_idx + 1)
  ends   <- c(cp_idx, nrow(df))
  
  reg <- tibble(i_start = starts, i_end = ends) %>%
    filter(i_start <= i_end) %>%
    mutate(
      xstart = df$fecha[i_start],
      xend   = df$fecha[i_end],
      mean_reg = map2_dbl(i_start, i_end, \(a,b) mean(df$valor[a:b], na.rm = TRUE))
    )
  
  p <- ggplot(df, aes(fecha, valor)) +
    geom_line(color = "grey70") +
    geom_segment(
      data = reg,
      aes(x = xstart, xend = xend, y = mean_reg, yend = mean_reg),
      color = "blue", linewidth = 1.1
    ) +
    theme_minimal()
  
  ggsave(path_plot, plot = p, bg = "white")
}

df_cp_wide <- bind_rows(cp_rows_wide)

write.csv(df_cp_wide,
          file = file.path("Series_de_tiempo_CPEH", "change_points_summary.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")