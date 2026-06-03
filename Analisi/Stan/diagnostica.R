rm(list=ls());gc();

library(ggplot2)
library(bayesplot)
library(dplyr)
library(ggridges) 
library(cmdstanr)

setwd("C:/Users/Utente/OneDrive/Universita/Magistrale/2025-2026/Iterazione/Progetto/Analisi")

load("dati_modificati.Rdata")
load("Stan/fit_mcmc_advanced.Rdata")
file_csv <- list.files(path="Stan/stan_output_advanced",
                       pattern="\\.csv",
                       full.names=T)
fit_ricaricato <- as_cmdstan_fit(file_csv)

# 1: CHECK NUMERICO RIGOROSO ----

# Estrazione dei diagnostici del campionatore (Divergenze, E-BFMI, Treedepth)
sampler_diagnostics <- fit$sampler_diagnostics(format = "df")

divergences <- sum(sampler_diagnostics$divergent__)
treedepth_hits <- sum(sampler_diagnostics$treedepth__ >= 12) # il tuo max_treedepth

cat("=== VERIFICA CRITICITÀ ===\n")
cat("Transizioni divergenti totali:", divergences, "\n")
cat("Iterazioni che hanno raggiunto il max_treedepth:", treedepth_hits, "\n")

# Controllo sistematico di Rhat ed ESS su TUTTI i parametri (inclusi gli effetti casuali)
summary_all <- fit_ricaricato$summary()

bad_rhat <- summary_all |> dplyr::filter(rhat > 1.05)
low_ess_bulk <- summary_all |> dplyr::filter(ess_bulk < 400)
low_ess_tail <- summary_all |> dplyr::filter(ess_tail < 400)

cat("\nParametri con Rhat > 1.05:", nrow(bad_rhat), "\n")
cat("Parametri con ESS Bulk scarso (< 400):", nrow(low_ess_bulk), "\n")
cat("Parametri con ESS Tail scarso (< 400):", nrow(low_ess_tail), "\n")

if(nrow(bad_rhat) > 0) print(bad_rhat |> select(variable, rhat, ess_bulk))

# 2: ISPEZIONE VISIVA AVANZATA ----

# Parametri strutturali principali da monitorare
target_params <- c("beta0_Asym", "beta0_b2", "beta0_b3", "sigma", "sigma_u_Asym", "sigma_u_b2", "sigma_u_b3", "rho")
#target_params <- c("beta0_Asym", "beta0_b2", "beta0_b3", "sigma", "sigma_u", "rho")

# 1. Rank Plots (Overlay di istogrammi dei ranghi per verificare l'uniformità)
mcmc_rank_overlay(fit_ricaricato$draws(target_params)) +
  theme_minimal() +
  labs(title = "Rank Plots: l'uniformità indica un buon mescolamento tra le catene")

# 2. Autocorrelazione (Fondamentale dato che hai inserito una struttura AR(1))
mcmc_acf(fit_ricaricato$draws(target_params), lags = 20) +
  theme_minimal() +
  labs(title = "Funzione di Autocorrelazione (ACF) per i parametri chiave")

# 3. Analisi delle Divergenze (Se presenti, capiamo dove si concentrano)
if (divergences > 0) {
  # Coppia critica tipica nei modelli misti: log(varianza effetto casuale) vs effetto casuale non centrato
  mcmc_pairs(fit$draws(), pars = c("sigma_u", "beta0_Asym"), 
             np = nuts_params(fit), max_treedepth = 12)
}

# 3: POSTERIOR PREDICTIVE CHECKS PER DATI LONGITUDINALI ----

Y_rep  <- fit_ricaricato$draws("Y_rep", format = "matrix")
dim(Y_rep)
# numero di righe = numero di catene
# numero di colonne = numero di u.s.

# Prepariamo un dataframe di supporto con la struttura del design sperimentale
# Assicurati che l'ordine di stan_data$time e stan_data$biomass corrisponda a Y_vec
df_check <- dati |> na.omit() |>
  dplyr::select(OD, tempo, id_biomassa) |> 
  dplyr::rename(Y_vec=OD,
         Settimana=tempo)
Y_vec <- dati |> na.omit() |> dplyr::pull(OD)
rm(dati_incrementi, dati_incrementi_rapporti, dati);gc();

# PPC RAGGRUPPATO PER SETTIMANA
# Questo ti mostra se il modello azzecca la distribuzione dei dati settimana per settimana.
# Vedrai visivamente se l'incertezza aumenta o diminuisce coerentemente con i dati reali.
ppc_violin_grouped(Y_vec, Y_rep, group = df_check$Settimana, probs = c(0.1, 0.5, 0.9)) +
  theme_minimal() +
  labs(title = "PPC: Distribuzione della Biomassa per ciascuna Settimana",
       subtitle = "I punti verdi/neri sono i dati reali, le bande sono gli intervalli predittivi del modello")

# VERIFICA SUI MASSIMI LOCALI
# Il modello riesce a intercettare il picco di crescita massimo delle biomasse?
ppc_stat(Y_vec, Y_rep, stat = "max") + 
  theme_minimal() +
  labs(title = "Verifica sul Valore Massimo (Picco di Biomassa)")
mean(apply(Y_rep, 1, max) > max(Y_vec))

# Confronto SE ----

ppc_stat(Y_vec, Y_rep, stat = "sd") + 
  theme_minimal() +
  labs(title = "Verifica sulla Deviazione Standard Globale")
mean(apply(Y_rep, 1, sd) < sd(Y_vec))

# 1. Recuperiamo i vettori dei tempi e dei dati reali
tempi <- stan_data$t_idx # Vettore lungo 369 con i tempi reali (1:7)
tempi_unici <- sort(unique(tempi)) # I 7 time point disponibili

# 2. Inizializziamo una matrice per le SD simulate [8000 iterazioni x 7 tempi]
sd_simulate <- matrix(NA, nrow = nrow(Y_rep), ncol = length(tempi_unici))
colnames(sd_simulate) <- paste0("Tempo_", tempi_unici)

# Inizializziamo un vettore per le SD reali
sd_reali <- numeric(length(tempi_unici))

# 3. Calcolo ciclico delle SD per ogni specifico tempo
for (i in seq_along(tempi_unici)) {
  t <- tempi_unici[i]
  colonne_del_tempo_t <- which(tempi == t) # Trova quali colonne di Y_rep sono al tempo t
  
  # SD dei dati reali al tempo t
  sd_reali[i] <- sd(stan_data$Y[colonne_del_tempo_t])
  
  # SD delle 8000 repliche del modello al tempo t
  sd_simulate[, i] <- apply(Y_rep[, colonne_del_tempo_t], 1, sd)
}

# 4. Riorganizzazione dei dati per ggplot
df_sim_long <- as.data.frame(sd_simulate) |> 
  tidyr::pivot_longer(cols = everything(), names_to = "Tempo", values_to = "SD_Simulata") |> 
  dplyr::mutate(Tempo = as.numeric(gsub("Tempo_", "", Tempo)))

df_real_long <- data.frame(
  Tempo = tempi_unici,
  SD_Reale = sd_reali
)

# 5. GRAFICO DI CONFRONTO
ggplot() +
  # Distribuzione delle SD simulate dal modello (Boxplot azzurri)
  geom_boxplot(data = df_sim_long, aes(x = factor(Tempo), y = SD_Simulata), 
               fill = "skyblue", alpha = 0.6, outlier.alpha = 0.1) +
  # Valore reale della SD nei dati (Rombi rossi)
  geom_point(data = df_real_long, aes(x = factor(Tempo), y = SD_Reale), 
             color = "red", size = 4, shape = 18) +
  theme_minimal() +
  labs(
    title = "Verifica della Deviazione Standard (SD) lungo il Tempo",
    subtitle = "I boxplot indicano l'incertezza del modello (Y_rep). I rombi rossi indicano i dati reali (Y_real).",
    x = "Tempo (Rilevazione da 1 a 7)",
    y = "Deviazione Standard (Variabilità tra biomasse)"
  )

# FASE 4: CRITICITÀ DEI DATI (LOOCV) ----

# Calcolo LOO sfruttando l'oggetto di CmdStanR per maggiore stabilità
loo_obj <- fit_ricaricato$loo(cores = 6)
print(loo_obj)

# Plot dei valori di Pareto k
plot(loo_obj, label_points = TRUE) +
  theme_minimal() +
  labs(title = "Diagnostica dei valori di Pareto k (LOO-CV)",
       subtitle = "Punti con k > 0.7 indicano osservazioni altamente influenti o outlier")

# Identificazione geometrica delle osservazioni problematiche
k_table <- data.frame(obs_index = 1:stan_data$S, k = loo_obj$diagnostics$pareto_k)
problematic_obs <- k_table |> filter(k > 0.7)

if(nrow(problematic_obs) > 0) {
  cat("\nAttenzione: ci sono", nrow(problematic_obs), "osservazioni con Pareto k > 0.7\n")
  print(problematic_obs)
}

stan_data$Y[stan_data$id_biomassa==45]

# 5. Densità a Posteriori dei Parametri (con intervalli di credibilità) ----

plot_c_density <- function(coef){
  mcmc_areas(
    fit_ricaricato$draws(coef),
    prob = 0.5,       # Area interna (più scura): intervallo di credibilità al 50%
    prob_outer = 0.9, # Area esterna (più chiara): intervallo di credibilità al 90%
    point_est = "median" # Linea verticale sul valore mediano
  ) +
    geom_vline(aes(xintercept=0, col="red")) +
    theme(legend.position = "none") +
    theme_minimal() +
    labs(
      title = "Distribuzioni a Posteriori dei Parametri Chiave",
      subtitle = "I punti indicano le mediane, le aree colorate gli intervalli di credibilità al 50% e 90%"
    )
}

fit_ricaricato$metadata()$model_params |> head(20)

#target_params <- c("sigma", "sigma_u", "rho")
#plot_c_density(target_params)

target_params <- c("sigma", "sigma_u_Asym", "sigma_u_b2", "sigma_u_b3", "rho")
plot_c_density(target_params)

target_params <- grep(fit_ricaricato$metadata()$model_params, pattern="beta[0-9]_Asym", value = T)
plot_c_density(target_params)

target_params <- grep(fit_ricaricato$metadata()$model_params, pattern="beta[0-9]_b2", value = T)
plot_c_density(target_params)

target_params <- grep(fit_ricaricato$metadata()$model_params, pattern="beta[0-9]_b3", value = T)
plot_c_density(target_params)

# Test: tempo 6 vs tempo 7 ----

mu_new <- fit_ricaricato$draws("mu_new", format = "matrix")
dim(mu_new)

new_data_grid <- expand.grid(
  I_cov = c(-1, 0, 1),
  D_cov = c(-1, 0, 1),
  P_cov = c(-1, 0, 1)
) %>% mutate(id_combinazione = row_number())

matrice_differenze <- matrix(NA, nrow = nrow(mu_new), ncol = 27)
colnames(matrice_differenze) <- paste0("Scenario_", 1:27)

for (n in 1:27) {
  col_t6 <- paste0("mu_new[", n, ",6]")
  col_t7 <- paste0("mu_new[", n, ",7]")
  
  # Calcoliamo la differenza mantenendo il vettore di 8000 draw
  matrice_differenze[, n] <- mu_new[, col_t7] - mu_new[, col_t6]
}

df_distribuzioni <- as.data.frame(matrice_differenze) %>%
  mutate(id_draw_mcmc = row_number()) %>% 
  tidyr::pivot_longer(cols = starts_with("Scenario_"), 
                      names_to = "Scenario", 
                      values_to = "Differenza_T7_T6") %>%
  mutate(id_combinazione = as.numeric(gsub("Scenario_", "", Scenario))) %>%
  inner_join(new_data_grid, by = "id_combinazione") %>%
  select(id_draw_mcmc, I_cov, D_cov, P_cov, Differenza_T7_T6)
dim(df_distribuzioni)
# 216.000 righe (8000 draw * 27 scenari)
head(df_distribuzioni)

df_plot <- df_distribuzioni %>%
  mutate(Scenario_Label = paste0("I:", I_cov, " | D:", D_cov, " | P:", P_cov))

ggplot(df_plot, aes(x = Differenza_T7_T6, y = reorder(Scenario_Label, Differenza_T7_T6))) +
  # Disegna le densità colorandole in base al fatto che siano maggiori o minori di zero
  stat_density_ridges(aes(fill = stat(x) > 0), geom = "density_ridges_gradient", 
                      calc_ecdf = TRUE, alpha = 0.7, scale = 1.5) +
  # Linea verticale sullo zero: separa visivamente crescita da decrescita
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  scale_fill_manual(values = c("#d95f02", "#1b9e77"), name = "Incertezza MCMC",
                    labels = c("Decrescita (T7 < T6)", "Crescita (T7 > T6)")) +
  theme_minimal() +
  labs(
    title = "Distribuzione delle Differenze di Crescita (Tempo 7 - Tempo 6)",
    subtitle = "Analisi della forma della distribuzione su 8000 draw MCMC",
    x = "Valore della Differenza (mu_new[T7] - mu_new[T6])",
    y = "Scenari (Combinazioni Covariate)"
  ) +
  theme(panel.grid.minor = element_blank(), text = element_text(size = 12))

head(df_distribuzioni)
res <- df_distribuzioni |>  
  group_by(I_cov, D_cov, P_cov) |>  
  summarise(
    Mediana_Diff = median(Differenza_T7_T6),
    Prob_Decrescita = (mean(Differenza_T7_T6 < 0)*100) |> round(2),
    Prob_Crescita = (mean(Differenza_T7_T6 > 0)*100) |> round(2),
    .groups = "drop"
  ) |> 
  arrange(Mediana_Diff)
print(res, n=Inf)

# Confronto al tempo 6 ----

# matrice che contiene le previsioni al tempo 6
mu_t6 <- matrix(NA, nrow = nrow(mu_new), ncol = 27)
colnames(mu_t6) <- paste0("Scenario_", 1:27)
for (n in 1:27) {
  mu_t6[, n] <- mu_new[, paste0("mu_new[", n, ",6]")]
}

matrice_dominanza_t6 <- matrix(NA, nrow = 27, ncol = 27)
for (i in 1:27) {
  for (j in 1:27) {
    if (i == j) {
      matrice_dominanza_t6[i, j] <- 0.5 # Contro se stesso è un pareggio teorico
    } else {
      # Quante volte l'MCMC dello scenario i supera lo scenario j al tempo 7?
      matrice_dominanza_t6[i, j] <- mean(mu_t6[, i] > mu_t6[, j])
    }
  }
}

scenari_battuti <- rowSums(matrice_dominanza_t6 > 0.50)

res_confronto <- new_data_grid %>%
  mutate(
    Scenari_Battuti = scenari_battuti,
    Percentuale_Dominio = round((Scenari_Battuti / 26) * 100, 1)
  ) %>%
  arrange(desc(Scenari_Battuti))

print(res_confronto)

# Prepariamo le etichette per rendere il grafico leggibile
etichette <- new_data_grid %>%
  mutate(Label = paste0("I:", I_cov, " | D:", D_cov, " | P:", P_cov)) %>%
  select(id_combinazione, Label)

# Trasformiamo la matrice 27x27 in formato lungo per ggplot
df_heatmap <- as.data.frame(matrice_dominanza_t6) %>%
  mutate(Scenario_A = row_number()) %>%
  tidyr::pivot_longer(cols = starts_with("V"), names_to = "Scenario_B", values_to = "Prob_Vittoria") %>%
  mutate(Scenario_B = as.numeric(gsub("V", "", Scenario_B))) %>%
  # Attacchiamo le etichette descrittive
  inner_join(etichette, by = c("Scenario_A" = "id_combinazione")) %>%
  rename(Label_A = Label) %>%
  inner_join(etichette, by = c("Scenario_B" = "id_combinazione")) %>%
  rename(Label_B = Label)

# Uniamo il punteggio di vittorie per ordinare l'asse Y in modo logico
df_heatmap <- df_heatmap %>%
  inner_join(select(res_confronto, id_combinazione, Scenari_Battuti), by = c("Scenario_A" = "id_combinazione"))

# GRAFICO
ggplot(df_heatmap, aes(x = Label_B, y = reorder(Label_A, Scenari_Battuti), fill = Prob_Vittoria)) +
  geom_tile(color = "white", size = 0.1) +
  # Scala divergente: viola/blu se perde (<0.5), bianca al pareggio (0.5), verde/oro se vince (>0.5)
  scale_fill_gradient2(low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c", 
                       midpoint = 0.5, name = "Probabilità\nSotto > Destra") +
  theme_minimal() +
  labs(
    title = "Matrice di Dominanza Bayesiana al Tempo 7",
    subtitle = "Legenda: Uno scenario sull'asse Y domina quello sull'asse X se la cella tende al rosso (> 0.5)",
    x = "Scenario Avversario (Confronto)",
    y = "Scenario in Esame (Ordinati dal migliore in alto al peggiore in basso)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

# Confronto al tempo 7 ----

# matrice che contiene le previsioni al tempo 7
mu_t7 <- matrix(NA, nrow = nrow(mu_new), ncol = 27)
colnames(mu_t7) <- paste0("Scenario_", 1:27)
for (n in 1:27) {
  mu_t7[, n] <- mu_new[, paste0("mu_new[", n, ",7]")]
}

matrice_dominanza_t7 <- matrix(NA, nrow = 27, ncol = 27)
for (i in 1:27) {
  for (j in 1:27) {
    if (i == j) {
      matrice_dominanza_t7[i, j] <- 0.5 # Contro se stesso è un pareggio teorico
    } else {
      # Quante volte l'MCMC dello scenario i supera lo scenario j al tempo 7?
      matrice_dominanza_t7[i, j] <- mean(mu_t7[, i] > mu_t7[, j])
    }
  }
}

scenari_battuti <- rowSums(matrice_dominanza_t7 > 0.50)

res_confronto <- new_data_grid %>%
  mutate(
    Scenari_Battuti = scenari_battuti,
    Percentuale_Dominio = round((Scenari_Battuti / 26) * 100, 1)
  ) %>%
  arrange(desc(Scenari_Battuti))

print(res_confronto)

# Prepariamo le etichette per rendere il grafico leggibile
etichette <- new_data_grid %>%
  mutate(Label = paste0("I:", I_cov, " | D:", D_cov, " | P:", P_cov)) %>%
  select(id_combinazione, Label)

# Trasformiamo la matrice 27x27 in formato lungo per ggplot
df_heatmap <- as.data.frame(matrice_dominanza_t7) %>%
  mutate(Scenario_A = row_number()) %>%
  tidyr::pivot_longer(cols = starts_with("V"), names_to = "Scenario_B", values_to = "Prob_Vittoria") %>%
  mutate(Scenario_B = as.numeric(gsub("V", "", Scenario_B))) %>%
  # Attacchiamo le etichette descrittive
  inner_join(etichette, by = c("Scenario_A" = "id_combinazione")) %>%
  rename(Label_A = Label) %>%
  inner_join(etichette, by = c("Scenario_B" = "id_combinazione")) %>%
  rename(Label_B = Label)

# Uniamo il punteggio di vittorie per ordinare l'asse Y in modo logico
df_heatmap <- df_heatmap %>%
  inner_join(select(res_confronto, id_combinazione, Scenari_Battuti), by = c("Scenario_A" = "id_combinazione"))

# GRAFICO
ggplot(df_heatmap, aes(x = Label_B, y = reorder(Label_A, Scenari_Battuti), fill = Prob_Vittoria)) +
  geom_tile(color = "white", size = 0.1) +
  # Scala divergente: viola/blu se perde (<0.5), bianca al pareggio (0.5), verde/oro se vince (>0.5)
  scale_fill_gradient2(low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c", 
                       midpoint = 0.5, name = "Probabilità\nSotto > Destra") +
  theme_minimal() +
  labs(
    title = "Matrice di Dominanza Bayesiana al Tempo 7",
    subtitle = "Legenda: Uno scenario sull'asse Y domina quello sull'asse X se la cella tende al rosso (> 0.5)",
    x = "Scenario Avversario (Confronto)",
    y = "Scenario in Esame (Ordinati dal migliore in alto al peggiore in basso)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )
