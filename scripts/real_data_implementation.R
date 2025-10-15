source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
library(tidychangepoint)
library(ggplot2)
library(dplyr)
library(readr)
library(trend)
library(tidyr)
library(lubridate)
library(zoo)
#------------------ Datos sinteticos--------------------------------------------

## Salto media 1
n_stream <- 1000
theta_stream <- 300
mu1 <- 1

make_ncm_LR <- function(mu_shift) {
  function(xi, training_set, ...) {
    Non_conformity_LNR(xi, training_set, mu_r = mu_shift)
  }
}
kde_bf_fixed <- readRDS("data/kde_bf_fixed.rds")

stream <- c(rnorm(theta_stream, 0,1),rnorm(n_stream - theta_stream, mu1 ,1))

train <- rnorm(200, 0, 1)

res_MAD <- ICM(train,stream, Non_conformity_MAD, kde_bf_fixed, th =2)

res_pelt <- segment(stream, method = 'pelt')

res_IQR <- ICM(train,stream, Non_conformity_IQR, kde_bf_fixed, th = 2)

res_KNN <- ICM(train, stream, Non_conformity_KNN, kde_bf_fixed, th = 2.5, k = 7)

cp_1 <- res_MAD$change_point
cp_2 <- res_pelt$model$tau 
cp_3 <- res_IQR$change_point
cp_4 <- res_KNN$change_point

df <- tibble::tibble(
  t = seq_along(stream),
  y = as.numeric(stream)
)

vref <- tibble::tibble(
  x   = as.integer(c(theta_stream, cp_1, cp_2[1], cp_3, cp_4)),
  ref = factor(
    c("θ verdadero", "τ estimado (MAD)", "τ estimado (PELT)", "τ estimado (IQR)", "τ estimado (KNN)"),
    levels = c("θ verdadero", "τ estimado (MAD)", "τ estimado (PELT)", "τ estimado (IQR)", "τ estimado (KNN)")
  )
)

obj_salto_mu_1 <- list(
  df   = df,
  vref = vref,
  meta = list(
    n_stream     = n_stream,
    theta_stream = theta_stream,
    mu1          = mu1,
    cp = list(
      MAD  = cp_1,
      PELT = cp_2,
      IQR  = cp_3,
      KNN  = cp_4
    )
  )
)

#esta linea esta comentada para asegurar los resultados
#saveRDS(obj_salto_mu_1, file = "data/figura_resultado_sintetico_salto_mu_1.rds")

ggplot(df, aes(t, y)) +
  geom_line(linewidth = 0.3, color = "grey7") +
  geom_vline(
    data = vref,
    mapping = aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"      = "#D55E00",
      "τ estimado (MAD)" = "#0072B2",
      "τ estimado (PELT)"  = "#009E73",
      "τ estimado (IQR)" = "#CC79A7",
      "τ estimado (KNN)" = "red"
    ),
    breaks = levels(vref$ref)
  ) +
  scale_linetype_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"      = "solid",
      "τ estimado (MAD)" = "solid",
      "τ estimado (PELT)"  = "dashed",
      "τ estimado (IQR)"  = "dashed",
      "τ estimado (KNN)"  = "dashed"
    ),
    breaks = levels(vref$ref),
    guide = "none"  # evita el warning de override duplicado
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
  labs(title = "Serie con posible punto de cambio",
       x = "Índice temporal", y = "Valor") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 20, 10, 10),
    legend.position  = "right"
  )

## Multiples puntos de cambio UP-UP
multi_stream <- c(rnorm(500,0,1), rnorm(500,1,1), rnorm(500,2,1))
train_multi <- rnorm(200,0,1)

res_multi_MAD <- ICM_multi(train_multi, multi_stream, Non_conformity_MAD, kde_bf_fixed, th = 2.5)

res_multi_IQR <- ICM_multi(train_multi, multi_stream, Non_conformity_IQR, kde_bf_fixed, th = 2.5)

res_multi_KNN <- ICM_multi(train_multi, multi_stream, Non_conformity_KNN, kde_bf_fixed, th = 2.5)
res_multi_KNN$change_points

res_multi_PELT <- segment(multi_stream, method = 'pelt')

cp_MAD <- res_multi_MAD$change_points
cp_IQR <- res_multi_IQR$change_points
cp_KNN <- res_multi_KNN$change_points
cp_PELT <- res_multi_PELT$model$tau

df_multi <- tibble::tibble(
  t = seq_along(multi_stream),
  y = as.numeric(multi_stream)
)


theta_true_multi <- c(500L, 1000L)

vref_multi <- dplyr::bind_rows(
  tibble::tibble(x = theta_true_multi, ref = "θ verdadero"),
  tibble::tibble(x = cp_MAD,           ref = "τ estimado (MAD)"),
  tibble::tibble(x = cp_PELT,          ref = "τ estimado (PELT)"),
  tibble::tibble(x = cp_IQR,           ref = "τ estimado (IQR)"),
  tibble::tibble(x = cp_KNN,           ref = "τ estimado (KNN)")
) %>%
  dplyr::arrange(x) %>%
  dplyr::mutate(
    ref = factor(
      ref,
      levels = c("θ verdadero", "τ estimado (MAD)", "τ estimado (PELT)", "τ estimado (IQR)", "τ estimado (KNN)")
    )
  )

obj_salto_mu_MULTI <- list(
  df   = df_multi,
  vref = vref_multi,
  meta = list(
    n_total        = length(multi_stream),
    theta_true     = theta_true_multi,  # c(500, 1000)
    cp = list(
      MAD  = cp_MAD,
      PELT = cp_PELT,
      IQR  = cp_IQR,
      KNN  = cp_KNN
    )
  )
)

#esta linea esta comentada para asegurar los resultados
#saveRDS(obj_salto_mu_MULTI, file = "data/salto_mu_MULTI_obj.rds")


ggplot(df_multi, aes(t, y)) +
  geom_line(linewidth = 0.3, color = "grey7") +
  geom_vline(
    data = vref_multi,
    mapping = aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"        = "#D55E00",
      "τ estimado (MAD)"   = "#0072B2",
      "τ estimado (PELT)"  = "#009E73",
      "τ estimado (IQR)"   = "#CC79A7",
      "τ estimado (KNN)"   = "red"
    ),
    breaks = levels(vref_multi$ref)
  ) +
  scale_linetype_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"        = "solid",
      "τ estimado (MAD)"   = "solid",
      "τ estimado (PELT)"  = "dashed",
      "τ estimado (IQR)"   = "dashed",
      "τ estimado (KNN)"   = "dashed"
    ),
    breaks = levels(vref_multi$ref),
    guide = "none"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
  labs(title = "Serie con múltiples puntos de cambio",
       x = "Índice temporal", y = "Valor") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 20, 10, 10),
    legend.position  = "right"
  )

## Multiples puntos de cambio UP-DN

multi_stream_updn<- c(rnorm(500,1,1), rnorm(500,0,1), rnorm(500,1,1))
train_multi_updn <- rnorm(200,1,1)


res_multiadapt_MAD <- ICM_multi_adaptive(multi_stream_updn, Non_conformity_MAD, kde_bf_fixed, th = 2.5, train_multi_updn, m_retrain = 200)
res_multiadapt_MAD$change_points

res_multiadapt_IQR <- ICM_multi_adaptive(multi_stream_updn, Non_conformity_IQR, kde_bf_fixed, th = 2.5, train_multi_updn, m_retrain = 200)
res_multiadapt_IQR$change_points

res_multiadapt_KNN <- ICM_multi_adaptive(multi_stream_updn, Non_conformity_KNN, kde_bf_fixed, th = 2.5, train_multi_updn, m_retrain = 200, K = 7)
res_multiadapt_KNN$change_points

res_multi_updn <-ICM_multi(train_multi_updn,multi_stream_updn, Non_conformity_MAD, kde_bf_fixed, th = 3)
res_multi_updn$change_points

res_multi_PELT <- segment(multi_stream_updn, method = 'pelt')


cp_MAD_adapt <- res_multiadapt_MAD$change_points
cp_MAD_UPDN <- res_multi_updn$change_points
cp_IQR_adapt <- res_multiadapt_IQR$change_points
cp_KNN_adapt <- res_multiadapt_KNN$change_points
cp_PELT_updn <- res_multi_PELT$model$tau

df_multiadapt <- tibble::tibble(
  t = seq_along(multi_stream_updn),
  y = as.numeric(multi_stream_updn)
)


theta_true_multi <- c(500L, 1000L)

vref_multiadapt <- dplyr::bind_rows(
  tibble::tibble(x = theta_true_multi, ref = "θ verdadero"),
  tibble::tibble(x = cp_MAD_adapt,           ref = "τ estimado (MAD)"),
  tibble::tibble(x = cp_PELT_updn,          ref = "τ estimado (PELT)"),
  tibble::tibble(x = cp_IQR_adapt,           ref = "τ estimado (IQR)"),
  tibble::tibble(x = cp_KNN_adapt,           ref = "τ estimado (KNN)"),
  tibble::tibble(x = cp_MAD_UPDN,           ref = "τ estimado (ICM multi)")
) %>%
  dplyr::arrange(x) %>%
  dplyr::mutate(
    ref = factor(
      ref,
      levels = c("θ verdadero", "τ estimado (MAD)", "τ estimado (PELT)", "τ estimado (IQR)", "τ estimado (KNN)", "τ estimado (ICM multi)")
    )
  )

obj_salto_mu_MULTIADAPT <- list(
  df   = df_multiadapt,
  vref = vref_multiadapt,
  meta = list(
    n_total        = length(multi_stream),
    theta_true     = theta_true_multi,  # c(500, 1000)
    cp = list(
      MAD  = cp_MAD,
      PELT = cp_PELT_updn,
      IQR  = cp_IQR,
      KNN  = cp_KNN,
      ICM_multi = cp_MAD_UPDN
    )
  )
)

#esta linea esta comentada para asegurar los resultados
saveRDS(obj_salto_mu_MULTIADAPT, file = "data/salto_mu_MULTIADAPT_obj.rds")


ggplot(df_multiadapt, aes(t, y)) +
  geom_line(linewidth = 0.3, color = "grey7") +
  geom_vline(
    data = vref_multiadapt,
    mapping = aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"        = "#D55E00",
      "τ estimado (MAD)"   = "#0072B2",
      "τ estimado (PELT)"  = "#009E73",
      "τ estimado (IQR)"   = "#CC79A7",
      "τ estimado (KNN)"   = "red",
      "τ estimado (ICM multi)" = "yellow"
    ),
    breaks = levels(vref_multiadapt$ref)
  ) +
  scale_linetype_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"        = "solid",
      "τ estimado (MAD)"   = "solid",
      "τ estimado (PELT)"  = "dashed",
      "τ estimado (IQR)"   = "dashed",
      "τ estimado (KNN)"   = "dashed",
      "τ estimado (ICM multi)" = "dashed"
    ),
    breaks = levels(vref_multiadapt$ref),
    guide = "none"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
  labs(title = "Serie con múltiples puntos de cambio",
       x = "Índice temporal", y = "Valor") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 20, 10, 10),
    legend.position  = "right"
  )



#------------------ Datos reales--------------------------------------------

## Datos de baseball

df_all <- tibble::tibble(
  year = as.integer(substr(mlb_diffs$yearID, 1, 4)),
  y    = as.numeric(mlb_diffs$hr_rate_diff)
) |>
  dplyr::filter(!is.na(year), !is.na(y))

test_size <- 5L
training_baseball <- df_all$y[1:(test_size - 1)]
stream_baseball   <- df_all$y[(test_size):length(df_all$y)]
year_stream       <- df_all$year[(test_size):length(df_all$year)]

paramsbf <- list(
  L = 20,
  n_grid = 256,
  bw_floor = 0.02,
  min_history = 2
)


res_baseball_MAD <- ICM(training_baseball, stream_baseball, Non_conformity_MAD, Kernel_BF, th=2, params_bf = paramsbf)
res_baseball_MAD$change_point

res_baseball_IQR <- ICM(training_baseball, stream_baseball, Non_conformity_IQR, Kernel_BF, th=2,params_bf = paramsbf)
res_baseball_IQR$change_point

res_baseball_KNN <- ICM(training_baseball, stream_baseball, Non_conformity_KNN, Kernel_BF, th=2,params_bf = paramsbf, k=3)
res_baseball_KNN$change_point

res_baseball_PELT <- segment(stream_baseball, method = "pelt")

cp_baseball_MAD <- res_baseball_MAD$change_point - 8
cp_baseball_IQR <- res_baseball_IQR$change_point - 8
cp_baseball_KNN <- res_baseball_KNN$change_point - 8
cp_baseball_PELT <- res_baseball_PELT$model$tau[2]
cp_real <- which(year_stream == 1973)[1]
 
year_cp_MAD  <- year_stream[cp_baseball_MAD]
year_cp_IQR  <- year_stream[cp_baseball_IQR]
year_cp_KNN  <- year_stream[cp_baseball_KNN]
year_cp_PELT <- year_stream[cp_baseball_PELT]

# --- Data frame con años en el eje X ---
df_baseball <- tibble::tibble(
  year = year_stream,
  y    = stream_baseball
)

# --- vref con años, incluyendo el año real del cambio (real_cp = 1973) ---
vref_baseball <- tibble::tibble(
  x   = as.integer(c(real_cp, year_cp_MAD, year_cp_PELT, year_cp_IQR, year_cp_KNN)),
  ref = factor(
    c("θ verdadero", "τ estimado (MAD)", "τ estimado (PELT)", "τ estimado (IQR)", "τ estimado (KNN)"),
    levels = c("θ verdadero", "τ estimado (MAD)", "τ estimado (PELT)", "τ estimado (IQR)", "τ estimado (KNN)")
  )
)

# --- (Opcional) guarda todo junto como .rds para reusar luego ---
obj_baseball <- list(
  df   = df_baseball,
  vref = vref_baseball,
  meta = list(
    test_size = test_size,
    real_cp   = real_cp,
    cp_idx    = list(MAD = cp_baseball_MAD, PELT = cp_baseball_PELT, IQR = cp_baseball_IQR, KNN = cp_baseball_KNN)
  )
)
#esta linea esta comentada para asegurar los resultados
#saveRDS(obj_baseball, "data/baseball_data_obj.rds")

# --- Gráfico con años en el eje X ---
ggplot(df_baseball, aes(year, y)) +
  geom_line(linewidth = 0.3, color = "grey7") +
  geom_vline(
    data = vref_baseball,
    mapping = aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 1
  ) +
  scale_color_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"        = "#D55E00",
      "τ estimado (MAD)"   = "#0072B2",
      "τ estimado (PELT)"  = "#009E73",
      "τ estimado (IQR)"   = "#CC79A7",
      "τ estimado (KNN)"   = "red"
    ),
    breaks = levels(vref_baseball$ref)
  ) +
  scale_linetype_manual(
    name   = "Referencia",
    values = c(
      "θ verdadero"        = "solid",
      "τ estimado (MAD)"   = "solid",
      "τ estimado (PELT)"  = "dashed",
      "τ estimado (IQR)"   = "dashed",
      "τ estimado (KNN)"   = "dashed"
    ),
    breaks = levels(vref_baseball$ref),
    guide  = "none"
  ) +
  scale_x_continuous(breaks = pretty) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
  labs(title = "HR rate diff: detección de cambio",
       x = "Año", y = "Diferencia en tasa de HR") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 20, 10, 10),
    legend.position  = "right"
  )

## Datos caudales

### Primer rio Anamichu

path <- file.path("datos_caudales", "CON_Q_MEDIA_D_22017010.csv")

MAX_GAP_INTERP <- 80

DO_WINSOR <- TRUE
WINSOR_PROBS <- c(0.005, 0.995)

rio_anamichu <- read_delim(
  file = path,
  delim = "|",
  skip = 1,
  col_names = TRUE,
  col_types = cols(
    Fecha = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
    Valor = col_double()
  ),
  locale = locale(encoding = "UTF-8"),
  trim_ws = TRUE,
  show_col_types = FALSE
)

#Regulariza frecuencia diaria
full_idx <- tibble(Fecha = seq(min(rio_anamichu$Fecha), max(rio_anamichu$Fecha), by = "1 day"))
tsd <- full_idx %>% left_join(rio_anamichu, by = "Fecha")

#Imputación
tsd <- tsd %>%
  mutate(Valor_fill = na.approx(Valor, x = Fecha,
                                maxgap = MAX_GAP_INTERP, na.rm = FALSE))

#Transformación log1p
tsd <- tsd %>% mutate(Y = log1p(Valor_fill))

#Climatología mensual robusta
clim_m <- tsd %>%
  mutate(mes = month(Fecha)) %>%
  group_by(mes) %>%
  summarise(
    med_m = median(Y, na.rm = TRUE),
    mad_m = mad(Y, constant = 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mad_m = ifelse(mad_m <= 1e-6 | is.na(mad_m), 1e-6, mad_m))

tsd <- tsd %>%
  mutate(mes = month(Fecha)) %>%
  left_join(clim_m, by = "mes") %>%
  mutate(
    Y_deseas = Y - med_m,
    Y_std    = Y_deseas / mad_m
  )

#Winsorización suave
winsor <- function(x, probs = c(0.005, 0.995)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
tsd <- tsd %>%
  mutate(
    stream = if (DO_WINSOR) winsor(Y_std, WINSOR_PROBS) else Y_std
  )


par(mfrow = c(3,1), mar = c(3,4,2,1))
plot(tsd$Fecha, tsd$Valor, type = "l", col = "gray30",
     main = "Caudal diario (raw)", xlab = "", ylab = "Valor")
plot(tsd$Fecha, tsd$Y, type = "l", col = "gray30",
     main = "log1p(Valor)", xlab = "", ylab = "log1p")
plot(tsd$Fecha, tsd$stream, type = "l", col = "gray30",
     main = "Anomalía robusta (deseasonalizada y estandarizada por mes)", xlab = "Fecha", ylab = "anomalía")

n  <- length(tsd$stream)
n_train <- 400

rio_anamichu_train <- tsd$stream[1:n_train]
rio_anamichu_stream<- tsd$stream[(n_train + 1):n]


res_rioanamichu <- ICM_multi_adaptive(rio_anamichu_stream, Non_conformity_MAD, Kernel_BF, th =5,
                                      rio_anamichu_train, m_retrain = 400,guard_band = 100, params_bf =list(
                                        L = 400,
                                        n_grid = 512,
                                        bw_floor = 0.02,
                                        min_history = 50
                                      ))
res_rioanamichu$change_points

res_pettitt_rioanimachu <- pettitt.test(tsd$stream)

#res_gacoen <- segment(tsd$stream, method = "ga-coen", maxiter = 50)

cp_icm_stream_idx <- as.integer(res_rioanamichu$change_points)
cp_icm_global_idx <- n_train + cp_icm_stream_idx
cp_pettitt <- as.integer(unname(res_pettitt_rioanimachu$estimate)[1])



df_rio_1 <- tibble(
  t = tsd$Fecha,
  y = as.numeric(tsd$Valor)
)

vref_rio_1 <- tibble(
  x   = c(tsd$Fecha[cp_icm_global_idx], tsd$Fecha[cp_pettitt]),
  ref = factor(
    c(rep("τ estimado (ICM)", length(cp_icm_global_idx)), "τ estimado (Pettitt)"),
    levels = c("τ estimado (ICM)", "τ estimado (Pettitt)")
  )
)

#esta linea esta comentada para asegurar los resultados 
#saveRDS(list(df = df_rio_1, vref = vref_rio_1), file = "data/rio_anamichu_df_vref.rds")


ggplot(df_rio_1, aes(t, y)) +
  geom_line(linewidth = 0.35, color = "grey7") +
  geom_vline(
    data = vref_rio_1,
    aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    "τ estimado (ICM)"     = "red",
    "τ estimado (Pettitt)" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "τ estimado (ICM)"     = "dashed",
    "τ estimado (Pettitt)" = "solid"
  )) +
  # Usa UNA de estas dos según la clase de t:
  scale_x_datetime(date_breaks = "5 years", date_labels = "%Y") +
  #scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(
    title = "Serie con puntos de cambio",
    x = "Fecha", y = "Valor"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank())

#--------------------------grafica segmenteada---------------------------------

df <- tibble(
  t = tsd$Fecha,
  y = tsd$Valor_fill
) %>% arrange(t)

cp_dates_icm <- sort(unique(tsd$Fecha[cp_icm_global_idx]))
cp_dates_icm <- cp_dates_icm[cp_dates_icm > min(df$t) & cp_dates_icm < max(df$t)]

breaks <- c(min(df$t), cp_dates_icm, max(df$t))
df <- df %>%
  mutate(regime_id = cut(t, breaks = breaks, labels = FALSE, include.lowest = TRUE))

reg_stats <- df %>%
  dplyr::group_by(regime_id) %>%
  dplyr::summarise(
    start   = min(t), end = max(t),
    n_days  = n(),
    mean_y  = mean(y, na.rm = TRUE),
    med_y   = median(y, na.rm = TRUE),
    .groups = "drop"
  )

mean_glob   <- mean(df$y, na.rm = TRUE)
median_glob <- median(df$y, na.rm = TRUE)

reg_stats <- reg_stats %>%
  dplyr::mutate(
    mean_glob   = mean_glob,
    median_glob = median_glob
  )
seg_mean <- reg_stats %>%
  dplyr::transmute(x = start, xend = end, y = mean_y, yend = mean_y)

seg_median <- reg_stats %>%
  dplyr::transmute(x = start, xend = end, y = med_y, yend = med_y)

#Serie con escalón de MEDIA por régimen
ggplot(df, aes(t, y)) +
  geom_line(linewidth = 0.35, color = "grey70") +
  geom_segment(data = seg_mean,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.8, color = "steelblue", inherit.aes = FALSE) +
  geom_hline(yintercept = unique(reg_stats$mean_glob),
             linetype = "dotted", linewidth = 0.7, color = "black") +
  labs(title = "Serie con media por régimen (ICM) + media global",
       x = "Fecha", y = tsd$Valor_fill) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())


vref_rioanamichu <- tibble::tibble(
  x   = as.POSIXct(cp_dates_icm),
  ref = factor(
    rep("τ estimado (ICM)", length(cp_dates_icm)),
    levels = c("τ estimado (ICM)")
  )
)

# Objeto compacto, sin funciones: solo datos para reusar en otro notebook
obj_rioanamichu <- list(
  df       = df %>% dplyr::select(t, y, regime_id),
  vref     = vref_rioanamichu,
  segments = list(
    mean   = seg_mean,    # columnas: x, xend, y, yend
    median = seg_median   # columnas: x, xend, y, yend
  ),
  stats    = reg_stats,   # start, end, n_days, mean_y, med_y, mean_glob, median_glob
  globals  = list(
    mean   = mean_glob,
    median = median_glob
  ),
  meta     = list(
    x_label = "Fecha",
    y_label = "Valor (rellenado)",
    title_mean   = "Serie con media por régimen (ICM) + media global",
    title_median = "Serie con mediana por régimen (ICM) + mediana global",
    cp_dates     = cp_dates_icm
  )
)

# Guarda como .rds (ajusta la ruta si quieres)
saveRDS(obj_rioanamichu, "data/icm_regimes_rioanamichu.rds")



### Segundo rio RIONEGRO

path <- file.path("datos_caudales", "CON_Q_MEDIA_D_23087830.csv")

MAX_GAP_INTERP <- 10000

DO_WINSOR <- TRUE
WINSOR_PROBS <- c(0.005, 0.995)

rio_rionegro <- read_delim(
  file = path,
  delim = "|",
  skip = 1,
  col_names = TRUE,
  col_types = cols(
    Fecha = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
    Valor = col_double()
  ),
  locale = locale(encoding = "UTF-8"),
  trim_ws = TRUE,
  show_col_types = FALSE
)

#Regulariza frecuencia diaria
full_idx <- tibble(Fecha = seq(min(rio_rionegro$Fecha), max(rio_rionegro$Fecha), by = "1 day"))
tsd <- full_idx %>% left_join(rio_rionegro, by = "Fecha")

#Imputación
tsd <- tsd %>%
  mutate(Valor_fill = na.approx(Valor, x = Fecha,
                                maxgap = MAX_GAP_INTERP, na.rm = FALSE))

#Transformación log1p
tsd <- tsd %>% mutate(Y = log1p(Valor_fill))

#Climatología mensual robusta
clim_m <- tsd %>%
  mutate(mes = month(Fecha)) %>%
  group_by(mes) %>%
  summarise(
    med_m = median(Y, na.rm = TRUE),
    mad_m = mad(Y, constant = 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mad_m = ifelse(mad_m <= 1e-6 | is.na(mad_m), 1e-6, mad_m))

tsd <- tsd %>%
  mutate(mes = month(Fecha)) %>%
  left_join(clim_m, by = "mes") %>%
  mutate(
    Y_deseas = Y - med_m,
    Y_std    = Y_deseas / mad_m
  )

#Winsorización suave
winsor <- function(x, probs = c(0.005, 0.995)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
tsd <- tsd %>%
  mutate(
    stream = if (DO_WINSOR) winsor(Y_std, WINSOR_PROBS) else Y_std
  )


par(mfrow = c(3,1), mar = c(3,4,2,1))
plot(tsd$Fecha, tsd$Valor_fill, type = "l", col = "gray30",
     main = "Caudal diario (raw)", xlab = "", ylab = "Valor")
plot(tsd$Fecha, tsd$Y, type = "l", col = "gray30",
     main = "log1p(Valor)", xlab = "", ylab = "log1p")
plot(tsd$Fecha, tsd$stream, type = "l", col = "gray30",
     main = "Anomalía robusta (deseasonalizada y estandarizada por mes)", xlab = "Fecha", ylab = "anomalía")


n  <- length(tsd$stream)
n_train <- 200

rio_rionegro_train <- tsd$stream[1:n_train]
rio_rionegro_stream<- tsd$stream[(n_train + 1):n]


res_rionegro <- ICM_multi_adaptive(rio_rionegro_stream, Non_conformity_MAD, Kernel_BF, th =5,
                                   rio_rionegro_train, m_retrain = 200,guard_band = 50, params_bf =list(
                                        L = 200,
                                        n_grid = 512,
                                        bw_floor = 0.02,
                                        min_history = 50
                                      ))
res_rionegro$change_points

res_pettitt_rionegro <- pettitt.test(tsd$stream)

#res_gacoen <- segment(tsd$stream, method = "ga-coen", maxiter = 50)

cp_icm_stream_idx <- as.integer(res_rionegro$change_points)
cp_icm_global_idx <- n_train + cp_icm_stream_idx
cp_pettitt <- as.integer(unname(res_pettitt_rionegro$estimate)[1])



df_rio_2 <- tibble(
  t = tsd$Fecha,
  y = as.numeric(tsd$Valor_fill)
)

vref_rio_2 <- tibble(
  x   = c(tsd$Fecha[cp_icm_global_idx], tsd$Fecha[cp_pettitt]),
  ref = factor(
    c(rep("τ estimado (ICM)", length(cp_icm_global_idx)), "τ estimado (Pettitt)"),
    levels = c("τ estimado (ICM)", "τ estimado (Pettitt)")
  )
)

#esta linea esta comentada para asegurar los resultados 
#saveRDS(list(df = df_rio_2, vref = vref_rio_2), file = "data/rio__df_vref.rds")


ggplot(df_rio_2, aes(t, y)) +
  geom_line(linewidth = 0.35, color = "grey7") +
  geom_vline(
    data = vref_rio_2,
    aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(
    "τ estimado (ICM)"     = "red",
    "τ estimado (Pettitt)" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "τ estimado (ICM)"     = "dashed",
    "τ estimado (Pettitt)" = "solid"
  )) +
  # Usa UNA de estas dos según la clase de t:
  scale_x_datetime(date_breaks = "5 years", date_labels = "%Y") +
  #scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(
    title = "Serie con puntos de cambio",
    x = "Fecha", y = "Valor"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank())

#---------------------segmentar rio cauca---------------------------------------

df <- tibble(
  t = tsd$Fecha,
  y = tsd$Valor_fill
) %>% arrange(t)

cp_dates_icm <- sort(unique(tsd$Fecha[cp_icm_global_idx]))
cp_dates_icm <- cp_dates_icm[cp_dates_icm > min(df$t) & cp_dates_icm < max(df$t)]

breaks <- c(min(df$t), cp_dates_icm, max(df$t))
df <- df %>%
  mutate(regime_id = cut(t, breaks = breaks, labels = FALSE, include.lowest = TRUE))

reg_stats <- df %>%
  dplyr::group_by(regime_id) %>%
  dplyr::summarise(
    start   = min(t), end = max(t),
    n_days  = n(),
    mean_y  = mean(y, na.rm = TRUE),
    med_y   = median(y, na.rm = TRUE),
    .groups = "drop"
  )

mean_glob   <- mean(df$y, na.rm = TRUE)
median_glob <- median(df$y, na.rm = TRUE)

reg_stats <- reg_stats %>%
  dplyr::mutate(
    mean_glob   = mean_glob,
    median_glob = median_glob
  )
seg_mean <- reg_stats %>%
  dplyr::transmute(x = start, xend = end, y = mean_y, yend = mean_y)

seg_median <- reg_stats %>%
  dplyr::transmute(x = start, xend = end, y = med_y, yend = med_y)

#Serie con escalón de MEDIA por régimen
ggplot(df, aes(t, y)) +
  geom_line(linewidth = 0.35, color = "grey70") +
  geom_segment(data = seg_mean,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.8, color = "steelblue", inherit.aes = FALSE) +
  geom_hline(yintercept = unique(reg_stats$mean_glob),
             linetype = "dotted", linewidth = 0.7, color = "black") +
  labs(title = "Serie con media por régimen (ICM) + media global",
       x = "Fecha", y = tsd$Valor_fill) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())


vref_riorionegro <- tibble::tibble(
  x   = as.POSIXct(cp_dates_icm),
  ref = factor(
    rep("τ estimado (ICM)", length(cp_dates_icm)),
    levels = c("τ estimado (ICM)")
  )
)


obj_riorionegro <- list(
  df       = df %>% dplyr::select(t, y, regime_id),
  vref     = vref_riorionegro,
  segments = list(
    mean   = seg_mean,    # columnas: x, xend, y, yend
    median = seg_median   # columnas: x, xend, y, yend
  ),
  stats    = reg_stats,   # start, end, n_days, mean_y, med_y, mean_glob, median_glob
  globals  = list(
    mean   = mean_glob,
    median = median_glob
  ),
  meta     = list(
    x_label = "Fecha",
    y_label = "Valor (rellenado)",
    title_mean   = "Serie con media por régimen (ICM) + media global",
    title_median = "Serie con mediana por régimen (ICM) + mediana global",
    cp_dates     = cp_dates_icm
  )
)


saveRDS(obj_riorionegro, "data/icm_regimes_riorionegro.rds")