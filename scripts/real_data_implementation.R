source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
library(tidychangepoint)
library(ggplot2)
library(dplyr)
library(readr)
library(trend)
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


stream <- c(rnorm(theta_stream, 0,1),rnorm(n_stream - theta_stream, mu1 ,1))

train <- rnorm(200, 0, 1)

res_MAD <- ICM(train,stream, Non_conformity_MAD, Mixture_BF, th =5)

res_pelt <- segment(stream, method = 'pelt')

res_IQR <- ICM(train,stream, Non_conformity_IQR, Mixture_BF, th = 5)

res_KNN <- ICM(train, stream, Non_conformity_KNN, Mixture_BF, th = 5, k = 7)

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

res_multi_MAD <- ICM_multi(train_multi, multi_stream, Non_conformity_MAD, Mixture_BF, th = 5)

res_multi_IQR <- ICM_multi(train_multi, multi_stream, Non_conformity_IQR, Mixture_BF, th = 5)

res_multi_KNN <- ICM_multi(train_multi, multi_stream, Non_conformity_KNN, Mixture_BF, th = 5)

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


res_multiadapt_MAD <- ICM_multi_adaptive(multi_stream_updn, Non_conformity_MAD, Mixture_BF, th = 5, train_multi_updn, m_retrain = 50)

res_multiadapt_IQR <- ICM_multi_adaptive(multi_stream_updn, Non_conformity_IQR, Mixture_BF, th = 5, train_multi_updn, m_retrain = 50)

res_multiadapt_KNN <- ICM_multi_adaptive(multi_stream_updn, Non_conformity_KNN, Mixture_BF, th = 5, train_multi_updn, m_retrain = 50, K = 7)

res_multi_updn <-ICM_multi(train_multi_updn,multi_stream_updn, Non_conformity_MAD, Mixture_BF, th = 5)

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
      PELT = cp_PELT,
      IQR  = cp_IQR,
      KNN  = cp_KNN,
      ICM_multi = cp_MAD_UPDN
    )
  )
)

#esta linea esta comentada para asegurar los resultados
#saveRDS(obj_salto_mu_MULTIADAPT, file = "data/salto_mu_MULTIADAPT_obj.rds")


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

baseball_data <- mlb_diffs$hr_rate_diff
real_cp <- 1973
test_size = 10

year_vector <- mlb_diffs$yearID[test_size:length(mlb_diffs$yearID)]
year_vector <- as.numeric(substr(year_vector, 1, 4))

training_baseball <- baseball_data[1:test_size-1]
stream_baseball <- baseball_data[test_size:length(baseball_data)]


res_baseball_MAD <- ICM(training_baseball, stream_baseball, Non_conformity_MAD, Mixture_BF, th=3)

res_baseball_IQR <- ICM(training_baseball, stream_baseball, Non_conformity_IQR, Mixture_BF, th=3)

res_baseball_KNN <- ICM(training_baseball, stream_baseball, Non_conformity_KNN, Mixture_BF, th=3, k=7)

res_baseball_PELT <- segment(stream_baseball, method = "pelt")

cp_baseball_MAD <- res_baseball_MAD$change_point
cp_baseball_IQR <- res_baseball_IQR$change_point
cp_baseball_KNN <- res_baseball_KNN$change_point
cp_baseball_PELT <- res_baseball_PELT$model$tau[2]

 
year_cp_MAD  <- year_vector[cp_baseball_MAD]
year_cp_IQR  <- year_vector[cp_baseball_IQR]
year_cp_KNN  <- year_vector[cp_baseball_KNN]
year_cp_PELT <- year_vector[cp_baseball_PELT]

# --- Data frame con años en el eje X ---
df_baseball <- tibble::tibble(
  year = year_vector,
  y    = as.numeric(stream_baseball)
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
  scale_x_continuous(breaks = pretty) +  # años "bonitos"
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

path <- file.path("datos_caudales", "CON_Q_MEDIA_D_22017010.csv")

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

ts <- rio_anamichu %>% arrange(Fecha) %>% filter(!is.na(Valor))
n  <- nrow(ts)
n_train <- floor(0.2 * n)

rio_anamichu_train <- ts[1:n_train, ]$Valor
rio_anamichu_stream<- ts[(n_train + 1):n, ]$Valor


res_rioanimachu <- ICM_multi_adaptive(rio_anamichu_stream, Non_conformity_MAD, histogram_betting_function, th =20,
                                      rio_anamichu_train, m_retrain = 400, params_bf =list(num_bins = 5))

res_pettitt_rioanimachu <- pettitt.test(rio_animachu_full)


cp_icm_stream_idx <- as.integer(res_rioanimachu$change_points)
cp_icm_global_idx <- n_train + cp_icm_stream_idx
cp_pettitt <- as.integer(unname(res_pettitt_rioanimachu$estimate)[1])



df_rio_1 <- tibble(
  t = seq_len(nrow(ts)),
  y = as.numeric(ts$Valor)
)

vref_rio_1 <- tibble(
  x   = c(cp_icm_global_idx, cp_pettitt),
  ref = factor(
    c(rep("τ estimado (ICM)", length(cp_icm_global_idx)), "τ estimado (Pettitt)"),
    levels = c("τ estimado (ICM)", "τ estimado (Pettitt)")
  )
)

#esta linea esta comentada para asegurar los resultados 
# saveRDS(list(df = df_rio_1, vref = vref_rio_1), file = "data/rio_anamichu_df_vref.rds")


ggplot(df_rio_1, aes(t, y)) +
  geom_line(linewidth = 0.35, color = "grey7") +
  geom_vline(
    data = vref_rio_1,
    aes(xintercept = x, color = ref, linetype = ref),
    linewidth = 0.7
  ) +
  scale_color_manual(values = c(
    "τ estimado (ICM)"    = "#0072B2",
    "τ estimado (Pettitt)" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "τ estimado (ICM)"    = "dashed",
    "τ estimado (Pettitt)" = "solid"
  )) +
  labs(
    title = "Serie con puntos de cambio",
    x = "Índice", y = "Valor"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank()
  )

idx <- c(2934, 3583 ,5901 ,6806 ,7472)

# construir data frame por índice
y  <- as.numeric(rio_anamichu_stream)
df <- tibble(t = seq_along(y), y = y)

# mantener solo índices válidos
idx <- idx[idx >= 1 & idx <= nrow(df)]

ggplot(df, aes(x = df, y = Valor)) +
  geom_line(linewidth = 0.35) +
  geom_vline(xintercept = idx, color = "firebrick", linetype = "dashed", linewidth = 0.5) +
  labs(
    title = "rio_animachu_stream",
    x = "Índice",
    y = "Valor"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())


