# 1) Cargar todas las funciones base
source("R/non_conformity_measures.R")
source("R/betting_functions.R")
source("R/icm_method.R")
source("R/montecarlo_sims.R")

# 2) Cargar la betting function precomputada (Generada por calibrar_kde.R)
kde_bf_fixed <- readRDS("data/kde_bf_fixed.rds")

# 3) Definir escenarios y parámetros
theta_vals <- c(100, 200)
mu1_vals   <- c(1, 1.5, 2)
h_vals     <- seq(1, 6, 0.5)
n_sim      <- 200
m          <- 200
k_par      <- 7

# 4) Ejecutar Monte Carlo en cada escenario y cada método
all_results <- list()

for (theta_s in theta_vals) {
  for (mu1_s in mu1_vals) {
    df_CBF <- montecarlo_ICM(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      ncm_fun      = Non_conformity_KNN,
      bet_fun      = Constant_BF,
      k            = k_par
    ) %>% mutate(Method = "Constant BF")
    
    df_MBF <- montecarlo_ICM(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      ncm_fun      = Non_conformity_KNN,
      bet_fun      = Mixture_BF,
      k            = k_par
    ) %>% mutate(Method = "Mixture BF")
    
    df_KDEF <- montecarlo_ICM(
      n_sim        = n_sim,
      h_vals       = h_vals,
      theta_stream = theta_s,
      mu1          = mu1_s,
      m            = m,
      ncm_fun      = Non_conformity_KNN,
      bet_fun      = kde_bf_fixed,
      k            = k_par
    ) %>% mutate(Method = "Precomputed KDE BF")
    
    combined <- bind_rows(df_CBF, df_MBF, df_KDEF) %>%
      mutate(
        theta_stream = factor(theta_stream, levels = theta_vals),
        mu1          = factor(mu1,        levels = mu1_vals),
        scenario_id  = paste0("θ=", theta_s, "_μ1=", mu1_s)
      )
    
    all_results[[paste0("θ=", theta_s, "_μ1=", mu1_s)]] <- combined
  }
}

# 5) Concatenar y salvar
df_all_methods <- bind_rows(all_results, .id = "scenario_key")
saveRDS(df_all_methods, file = "data/resultados_sim.rds")