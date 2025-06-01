# 1) Cargar librerías y funciones si necesitas (por ej. ggplot2)
library(ggplot2)
library(dplyr)

# 2) Cargar los resultados pre‐calculados
df_all_methods <- readRDS("data/resultados_sim.rds")

# 3) Graficar para cada método por separado
para_metodo <- function(metodo, color, filename_pdf) {
  df_sub <- df_all_methods %>% filter(Method == metodo)
  
  p <- ggplot(df_sub, aes(x = p_false_alarm, y = log_delay)) +
    geom_line(color = color, size = 1) +
    geom_point(color = color, size = 1) +
    facet_grid(rows = vars(theta_stream), cols = vars(mu1),
               labeller = labeller(
                 theta_stream = function(x) paste0("θ = ", x),
                 mu1          = function(x) paste0("μ₁ = ", x)
               )) +
    theme_minimal(base_size = 13) +
    labs(
      x     = expression(P[0](tau <= theta)),
      y     = expression(log[10](1 + E[1][tau - theta ~|~ tau > theta])),
      title = paste("ICM con", metodo),
      subtitle = "Paneles: θ_stream = 100,200 y μ₁ = 1,1.5,2"
    ) +
    theme(
      strip.text.x     = element_text(size = 12, face = "bold"),
      strip.text.y     = element_text(size = 12, face = "bold"),
      panel.spacing    = unit(0.6, "lines"),
      axis.text        = element_text(size = 10),
      axis.title       = element_text(size = 12, face = "bold"),
      plot.title       = element_text(size = 14, face = "bold"),
      plot.subtitle    = element_text(size = 11),
      legend.position  = "none"
    )
  
  ggsave(filename = paste0("outputs/", filename_pdf),
         plot     = p,
         width    = 8,
         height   = 6)
}

# 4) Crear la carpeta outputs/ (si no existe) y guardar figuras
if (!dir.exists("outputs")) dir.create("outputs")

para_metodo("Constant BF",       "red",         "fig_constantBF.pdf")
para_metodo("Mixture BF",        "blue",        "fig_mixtureBF.pdf")
para_metodo("Precomputed KDE BF","darkgreen",   "fig_precomputedKDEBF.pdf")