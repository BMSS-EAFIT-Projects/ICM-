library(readxl)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

# ── 1. Rutas ─────────────────────────────────────────────────
path_excel <- "datos_caudales/series_mensuales_por_estacion.xlsx"
path_cp    <- "Series_de_tiempo_CPEH/change_points_summary.csv"
out_dir    <- "Series_tendencia"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── 2. Leer y limpiar series (igual que antes) ───────────────
wide_to_long <- function(df, estacion) {
  names(df)[1] <- "anio"
  df %>%
    pivot_longer(cols = 2:13, names_to = "mes", values_to = "valor") %>%
    mutate(
      mes      = as.integer(gsub("^[xX]", "", mes)),
      valor    = as.numeric(valor),
      fecha    = as.Date(sprintf("%d-%02d-01", anio, mes)),
      estacion = estacion
    ) %>%
    relocate(estacion, anio, mes, fecha, valor)
}

clean_impute_station <- function(df) {
  df <- df %>%
    arrange(fecha) %>%
    mutate(valor = as.numeric(valor),
           valor = ifelse(is.nan(valor), NA_real_, valor))
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

sheets     <- excel_sheets(path_excel)
dfs_filled <- sheets %>%
  setNames(., .) %>%
  map(~ read_excel(path_excel, sheet = .x)) %>%
  imap(wide_to_long) %>%
  map(clean_impute_station)

# ── 3. Leer change points ────────────────────────────────────
cp_raw <- read.csv(path_cp, check.names = FALSE, stringsAsFactors = FALSE)
# La primera columna puede llamarse "i" o "Estacion" según tu CSV
cp_col <- names(cp_raw)[1]

get_cp_dates <- function(station_name) {
  row <- cp_raw[as.character(cp_raw[[cp_col]]) == as.character(station_name), ]
  if (nrow(row) == 0) return(as.Date(character(0)))
  vals <- row[1, grep("^CP", names(row))] %>% unlist() %>% na.omit()
  vals <- vals[nzchar(vals) & vals != "NA"]
  if (length(vals) == 0) return(as.Date(character(0)))
  as.Date(vals, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y"))
}

# ── 4. Construir tabla de regímenes a partir de CPs ──────────
build_regimes <- function(df, cp_dates) {
  n <- nrow(df)
  cp_idx <- integer(0)
  if (length(cp_dates) > 0) {
    cp_idx <- sort(unique(sapply(cp_dates, \(d) which.min(abs(df$fecha - d)))))
    cp_idx <- cp_idx[cp_idx >= 1 & cp_idx <= n]
  }
  starts <- c(1, cp_idx + 1)
  ends   <- c(cp_idx, n)
  tibble(i_start = starts, i_end = ends) %>%
    filter(i_start <= i_end) %>%
    mutate(
      fecha_inicio = df$fecha[i_start],
      fecha_fin    = df$fecha[i_end],
      mean_val     = map2_dbl(i_start, i_end, \(a, b) mean(df$valor[a:b], na.rm = TRUE)),
      regime_id    = row_number()
    )
}

# ── 5. Regresión lineal por régimen ──────────────────────────
lm_regime <- function(df, i_start, i_end, rid) {
  sub <- df[i_start:i_end, ]
  if (nrow(sub) < 3) {
    return(tibble(fecha = sub$fecha,
                  fitted = mean(sub$valor, na.rm = TRUE),
                  regime_id = rid))
  }
  sub <- mutate(sub, t = as.numeric(fecha))
  fit <- lm(valor ~ t, data = sub)
  tibble(fecha = sub$fecha, fitted = predict(fit), regime_id = rid)
}

# ── 6. Graficar y guardar ────────────────────────────────────
for (stn in sheets) {
  df       <- dfs_filled[[stn]]
  regimes  <- build_regimes(df, get_cp_dates(stn))
  
  lm_lines <- pmap_dfr(
    list(regimes$i_start, regimes$i_end, regimes$regime_id),
    \(a, b, rid) lm_regime(df, a, b, rid)
  )
  
  cp_dates_stn <- get_cp_dates(stn)
  cp_lines <- if (length(cp_dates_stn) > 0) {
    tibble(fecha = cp_dates_stn)
  } else {
    tibble(fecha = as.Date(character(0)))
  }
  
  p <- ggplot(df, aes(fecha, valor)) +
    # Serie original
    geom_line(color = "grey65", linewidth = 0.45) +
    # Líneas verticales en cada CP
    geom_vline(
      data = cp_lines,
      aes(xintercept = as.numeric(fecha)),
      color = "tomato", linetype = "dashed", linewidth = 0.55
    ) +
    # Media horizontal por régimen
    geom_segment(
      data = regimes,
      aes(x = fecha_inicio, xend = fecha_fin,
          y = mean_val,     yend = mean_val),
      color = "steelblue", linewidth = 1.0, lineend = "round"
    ) +
    # Tendencia lineal por régimen
    geom_line(
      data = lm_lines,
      aes(x = fecha, y = fitted, group = regime_id),
      color = "darkred", linewidth = 0.85
    ) +
    labs(
      title    = paste("Estación:", stn),
      subtitle = sprintf(
        "%d régimen(es)  |  azul = media  |  rojo = tendencia  |  línea punteada = CP",
        nrow(regimes)
      ),
      x = NULL,
      y = expression("Caudal (m"^3*"/s)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(color = "grey45", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    file.path(out_dir, paste0(stn, ".png")),
    plot = p, width = 11, height = 4, dpi = 150, bg = "white"
  )
  message("✓ ", stn)
}

message("\nListo. Imágenes guardadas en: ", out_dir)