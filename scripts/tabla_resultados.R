library(dplyr)
library(tidyr)
library(readr)
library(stringr)


df <- readRDS("data/prueba_modhist.rds")

# Columnas esperadas: Method, scenario_key, p_false_alarm, mean_delay
df <- df %>%
  mutate(p_false_alarm = as.numeric(p_false_alarm),
         mean_delay    = as.numeric(mean_delay))

# Columnas esperadas: Method, scenario_key, p_false_alarm, mean_delay
stopifnot(all(c("Method","scenario_key","p_false_alarm","mean_delay") %in% names(df)))

df <- df %>%
  mutate(p_false_alarm = as.numeric(p_false_alarm),
         mean_delay    = as.numeric(mean_delay))

# --- Normalización de strings para parsing robusto ---
norm <- function(x){
  x |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_replace_all("[-_]", " ") |>
    stringr::str_trim()
}

Method_clean <- norm(df$Method)
Method_up    <- toupper(Method_clean)

# --- Detectar NCM (incluye sinónimos) ---
detect_ncm <- function(UP) {
  dplyr::case_when(
    str_detect(UP, "\\bMAD\\b")                ~ "MAD",
    str_detect(UP, "\\bKNN\\b")                ~ "KNN",
    str_detect(UP, "\\bLNR\\b|LIKELIHOOD|\\bLR\\b") ~ "LR",
    str_detect(UP, "\\bIQR\\b") ~ "IQR",
    TRUE ~ NA_character_
  )
}

# --- Detectar BF (sinónimos y variantes) ---
detect_bf <- function(UP) {
  dplyr::case_when(
    str_detect(UP, "CONST(ANT)?\\s+BF")                          ~ "Constant BF",
    str_detect(UP, "MIX(TURE)?\\s+BF")                           ~ "Mixture BF",
    str_detect(UP, "HIST(OGRAM)?(\\s+BASED)?\\s+BF|\\bHIST\\s+BF\\b") ~ "Histogram BF",
    str_detect(UP, "KDE.*(PRE|PRE\\s*COMPUT|PRECOMPUT|PLUGIN|PLUG\\s*IN)") |
      str_detect(UP, "(PRE|PRE\\s*COMPUT|PRECOMPUT|PLUGIN|PLUG\\s*IN).*KDE") ~ "Precomputed KDE BF",
    TRUE ~ NA_character_
  )
}

df$NCM <- detect_ncm(Method_up)
df$BF  <- detect_bf(Method_up)

# --- Filtros solicitados ---
allowed_ncm <- c("MAD","LR","KNN", "IQR")
allowed_bf  <- c("Constant BF","Mixture BF","Histogram BF","Precomputed KDE BF")

df_f <- df %>%
  filter(NCM %in% allowed_ncm,
         BF  %in% allowed_bf,
         !str_detect(Method_up, "ORACLE"),
         !str_detect(Method_up, "\\bCBF\\b"),
         !str_detect(toupper(BF %||% ""), "\\bCBF\\b"))

# --- Retardo a PFA objetivo (interp con fallback al más cercano) ---
delay_at_interp <- function(p, y, target) {
  ok <- is.finite(p) & is.finite(y); p <- p[ok]; y <- y[ok]
  if (!length(p)) return(NA_real_)
  o <- order(p); p <- p[o]; y <- y[o]
  if (target >= min(p) && target <= max(p)) {
    as.numeric(approx(x = p, y = y, xout = target, ties = "ordered", rule = 1)$y)
  } else {
    y[which.min(abs(p - target))]
  }
}

bf_summary <- df_f %>%
  group_by(BF, scenario_key, NCM) %>%
  summarise(
    delay_5  = delay_at_interp(p_false_alarm, mean_delay, 0.05),
    delay_10 = delay_at_interp(p_false_alarm, mean_delay, 0.10),
    .groups = "drop"
  ) %>%
  mutate(delay_5 = round(delay_5, 2),
         delay_10 = round(delay_10, 2))

saveRDS(bf_summary, file = "data/tabla_delays.rds")