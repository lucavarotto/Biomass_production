# ==============================================================================
# ottimizzazione.R
# Ottimizzazione della superficie di risposta
#
# Prerequisito: modello.R deve essere stato eseguito e il fit salvato in
#   Modello finale/fits/re_b3_het_pars.rds
#
# Per ogni orizzonte temporale t* si cerca la combinazione (xI, xD, xP) in
# [-1,1]^3 che massimizza l'OD predetto:
#   OD(t*, x) = exp(eta_A) * exp(-exp(eta_B) * plogis(eta_C)^t*)
# con v = 0 (previsione per nuova unità sperimentale, RE marginato).
# L'ottimizzazione viene fatta su griglia densa per ogni campione posteriore.
# ==============================================================================

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

SCRIPT_DIR  <- dirname(rstudioapi::getActiveDocumentContext()$path)
GRAFICI_DIR <- file.path(SCRIPT_DIR, "grafici")
dir.create(GRAFICI_DIR, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# Caricamento fit
# ==============================================================================

fit <- readRDS(file.path(SCRIPT_DIR, "fits", "modello.rds"))
cat("Fit caricato da fits/modello.rds\n")

# ==============================================================================
# Sezione 1 — Estrai campioni posteriori (subsample a 2000 per velocità)
# ==============================================================================

cat("\n=== OTTIMIZZAZIONE ===\n")

draws_opt <- fit$draws(format = "matrix")

beta_A_all <- draws_opt[, grep("^beta_A\\[", colnames(draws_opt))]
beta_B_all <- draws_opt[, grep("^beta_B\\[", colnames(draws_opt))]
beta_C_all <- draws_opt[, grep("^beta_C\\[", colnames(draws_opt))]

set.seed(42)
idx_sub <- sample(nrow(beta_A_all), 2000)
beta_A  <- beta_A_all[idx_sub, ]   # 2000 × 10
beta_B  <- beta_B_all[idx_sub, ]
beta_C  <- beta_C_all[idx_sub, ]

# ==============================================================================
# Sezione 2 — Griglia densa su [-1, 1]^3 (30 punti per lato = 27.000 punti)
# ==============================================================================

n_grid   <- 30
grid_seq <- seq(-1, 1, length.out = n_grid)
grid_df  <- expand.grid(xI = grid_seq, xD = grid_seq, xP = grid_seq)

X_grid <- model.matrix(
  ~ xI + xD + xP + xI:xD + xI:xP + xD:xP + I(xI^2) + I(xD^2) + I(xP^2),
  data = grid_df
)  # 27000 × 10

cat("Griglia:", nrow(grid_df), "punti |",
    round(2000 * nrow(grid_df) * 8 / 1e6), "MB per matrice eta\n")

# ==============================================================================
# Sezione 3 — Predittori lineari (campioni × punti griglia)
# ==============================================================================

eta_A <- beta_A %*% t(X_grid)   # 2000 × 27000
eta_B <- beta_B %*% t(X_grid)
eta_C <- beta_C %*% t(X_grid)
# RE su b3: v = 0 per nuova unità (previsione marginata)

# ==============================================================================
# Sezione 4 — Ottimo per ogni t*
# ==============================================================================

t_star_vec <- c(5, 6, 7, 8, 9)   # scala 1-7 osservata; 8-9 estrapolazione

results_list <- lapply(t_star_vec, function(t_star) {
  OD_mat   <- exp(eta_A) * exp(-exp(eta_B) * plogis(eta_C)^t_star)
  best_idx <- max.col(OD_mat)
  opt_df   <- grid_df[best_idx, ]
  opt_df$OD_opt <- OD_mat[cbind(seq_len(nrow(OD_mat)), best_idx)]
  opt_df$t_star <- t_star
  opt_df
})

results <- bind_rows(results_list)

# Back-transform in unità originali
results <- results |>
  mutate(
    I_opt = xI * 190 + 280,
    D_opt = xD * 6   + 18,
    P_opt = xP * 24.5 + 25.5
  )

# ==============================================================================
# Sezione 5 — Tabella riassuntiva: mediana + IC 90% per ogni t*
# ==============================================================================

summary_opt <- results |>
  group_by(t_star) |>
  summarise(
    xI_med = median(xI), xI_q05 = quantile(xI, 0.05), xI_q95 = quantile(xI, 0.95),
    xD_med = median(xD), xD_q05 = quantile(xD, 0.05), xD_q95 = quantile(xD, 0.95),
    xP_med = median(xP), xP_q05 = quantile(xP, 0.05), xP_q95 = quantile(xP, 0.95),
    I_med  = median(I_opt),
    D_med  = median(D_opt),
    P_med  = median(P_opt),
    OD_med = median(OD_opt),
    OD_q05 = quantile(OD_opt, 0.05),
    OD_q95 = quantile(OD_opt, 0.95),
    .groups = "drop"
  )

cat("\n--- Tabella ottimo per t* ---\n")
print(summary_opt, width = 120)

# ==============================================================================
# Sezione 6 — Plot A: distribuzione posteriore dell'ottimo per variabile e t*
# ==============================================================================

results_long <- results |>
  select(t_star, xI, xD, xP) |>
  pivot_longer(cols = c(xI, xD, xP), names_to = "variabile", values_to = "valore")

p_violin <- ggplot(results_long, aes(x = factor(t_star), y = valore)) +
  geom_violin(fill = "steelblue", alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ variabile, labeller = labeller(variabile = c(
    xI = "Intensita' luminosa (xI)",
    xD = "Durata esposizione (xD)",
    xP = "Fosforo (xP)"
  ))) +
  scale_x_discrete(name = "Orizzonte temporale t* (settimane)") +
  scale_y_continuous(name = "Valore ottimale (scala RSM [-1, 1])") +
  labs(title = "Distribuzione posteriore dell'ottimo per ogni t*") +
  theme_minimal()

print(p_violin)
ggsave(file.path(GRAFICI_DIR, "ottimo_distribuzione.png"), plot = p_violin,
       width = 10, height = 5, dpi = 300)

# ==============================================================================
# Sezione 7 — Plot B: scatter triplette ottimali per ogni t*
# ==============================================================================

plot_scatter_ottimo <- function(data, t_val) {
  d  <- filter(data, t_star == t_val)
  p1 <- ggplot(d, aes(x = xI, y = xD)) +
    geom_point(alpha = 0.1, color = "steelblue", size = 0.8) +
    geom_density_2d(color = "navy", linewidth = 0.4) +
    labs(x = "xI*", y = "xD*") + theme_minimal()
  p2 <- ggplot(d, aes(x = xI, y = xP)) +
    geom_point(alpha = 0.1, color = "steelblue", size = 0.8) +
    geom_density_2d(color = "navy", linewidth = 0.4) +
    labs(x = "xI*", y = "xP*") + theme_minimal()
  p3 <- ggplot(d, aes(x = xD, y = xP)) +
    geom_point(alpha = 0.1, color = "steelblue", size = 0.8) +
    geom_density_2d(color = "navy", linewidth = 0.4) +
    labs(x = "xD*", y = "xP*") + theme_minimal()
  (p1 | p2 | p3) +
    plot_annotation(title = paste0("Scatter ottimo posteriore — t* = ", t_val, " settimane"))
}

for (tv in t_star_vec) {
  p_sc <- plot_scatter_ottimo(results, tv)
  print(p_sc)
  ggsave(file.path(GRAFICI_DIR, paste0("ottimo_scatter_t", tv, ".png")),
         plot = p_sc, width = 10, height = 4, dpi = 300)
}

# ==============================================================================
# Sezione 8 — Plot C: OD ottimale in funzione di t*
# ==============================================================================

p_od_tstar <- ggplot(summary_opt, aes(x = t_star, y = OD_med)) +
  geom_ribbon(aes(ymin = OD_q05, ymax = OD_q95), fill = "steelblue", alpha = 0.3) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  labs(
    x        = "Orizzonte temporale t* (settimane)",
    y        = "OD ottimale (mediana posteriore)",
    title    = "OD massimo raggiungibile in funzione del tempo",
    subtitle = "Banda: IC 90%"
  ) +
  theme_minimal()

print(p_od_tstar)
ggsave(file.path(GRAFICI_DIR, "od_ottimale_per_t.png"), plot = p_od_tstar,
       width = 7, height = 5, dpi = 300)

cat("\n=== OTTIMIZZAZIONE COMPLETATA ===\n")
