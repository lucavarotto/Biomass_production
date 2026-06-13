rm(list=ls());

SCRIPT_DIR  <- dirname(rstudioapi::getActiveDocumentContext()$path)
GRAFICI_DIR <- file.path(SCRIPT_DIR, "grafici")
dir.create(GRAFICI_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# modello.R
#
# Sezioni:
#   1 — Preparazione dati
#   2 — Stima
#   3 — Diagnostiche numeriche
#   4 — Rank plot e traceplot
#   5 — Posterior predictive check
#   6 — Posteriori coefficienti RSM
#   7 — Posteriori parametri strutturali (sigma, delta, sigma_C)
#   7b— Residui standardizzati vs valori predetti (verifica eteroschedasticita')
# ==============================================================================

library(cmdstanr)
library(posterior)
library(bayesplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

color_scheme_set("blue")

# ==============================================================================
# Sezione 1 — Preparazione dati
# ==============================================================================

load(file.path(SCRIPT_DIR, "dati_modificati.Rdata"))

dati <- dati |>
  mutate(
    xI      = (I - 280) / 190,
    xD      = (D - 18)  / 6,
    xP      = (P - 25.5) / 24.5,
    tempo_t = tempo + 1
  ) |>
  filter(!is.na(OD)) |>
  arrange(id_biomassa, tempo_t)

X_full <- model.matrix(
  ~ xI + xD + xP + xI:xD + xI:xP + xD:xP + I(xI^2) + I(xD^2) + I(xP^2),
  data = dati
)

curve_ids   <- unique(dati$id_biomassa)
K           <- length(curve_ids)
X_per_curve <- X_full[match(curve_ids, dati$id_biomassa), ]

start_idx <- integer(K)
end_idx   <- integer(K)
for (i in seq_len(K)) {
  idx <- which(dati$id_biomassa == curve_ids[i])
  start_idx[i] <- min(idx)
  end_idx[i]   <- max(idx)
}

stan_data <- list(
  N = nrow(dati), K = K, P = ncol(X_per_curve),
  OD = dati$OD, tempo = dati$tempo_t,
  X = X_per_curve, start_idx = start_idx, end_idx = end_idx
)

cat("Dati: N =", stan_data$N, "| K =", stan_data$K, "| P =", stan_data$P, "\n")

# ==============================================================================
# Sezione 2 — Stima
# ==============================================================================

mod <- cmdstan_model(file.path(SCRIPT_DIR, "modello.stan"))

fit <- mod$sample(
  data            = stan_data,
  chains          = 4,
  iter_warmup     = 1000,
  iter_sampling   = 2000,
  max_treedepth   = 15,
  seed            = 42,
  parallel_chains = 4
)

dir.create(file.path(SCRIPT_DIR, "fits"), showWarnings = FALSE)
fit$save_object(file.path(SCRIPT_DIR, "fits", "modello.rds"))
cat("Fit salvato in fits/modello.rds\n")

# ==============================================================================
# Sezione 3 — Diagnostiche numeriche
# ==============================================================================

cat("\n=== DIAGNOSTICHE ===\n")
diag <- fit$diagnostic_summary()
cat("Divergenze:       ", sum(diag$num_divergent), "\n")
cat("Max treedepth:    ", sum(diag$num_max_treedepth), "\n")
cat("E-BFMI:           ", round(diag$ebfmi, 3), "\n")

summ <- fit$summary(variables = c("beta_A", "beta_B", "beta_C",
                                   "sigma_C", "sigma", "delta"))
cat("R-hat max:        ", round(max(summ$rhat,      na.rm = TRUE), 4), "\n")
cat("ESS bulk min:     ", round(min(summ$ess_bulk,  na.rm = TRUE), 0), "\n")
cat("ESS tail min:     ", round(min(summ$ess_tail,  na.rm = TRUE), 0), "\n")
print(summ |> dplyr::select(variable, mean, median, sd, q5, q95, rhat, ess_bulk), n = 35)

# ==============================================================================
# Sezione 4 — Rank plot e traceplot
# ==============================================================================

beta_labels <- c("intercetta", "xI", "xD", "xP",
                 "xI:xD", "xI:xP", "xD:xP", "xI^2", "xD^2", "xP^2")

p_rank <- mcmc_rank_overlay(fit$draws(
    c("beta_A[1]", "beta_B[1]", "beta_C[1]", "sigma_C", "sigma", "delta"),
    format = "array")) +
    ggtitle("Rank plot") + theme_minimal()
print(p_rank)
ggsave(file.path(GRAFICI_DIR, "rank_plot.png"), plot = p_rank,
       width = 8, height = 5, dpi = 300)

draws_A <- fit$draws(paste0("beta_A[", 1:10, "]"), format = "array")
dimnames(draws_A)$variable <- paste0("A: ", beta_labels)
p_trace_A <- mcmc_trace(draws_A) + ggtitle("Traceplot beta_A") + theme_minimal(base_size = 8)
print(p_trace_A)
ggsave(file.path(GRAFICI_DIR, "traceplot_beta_A.png"), plot = p_trace_A,
       width = 10, height = 7, dpi = 300)

draws_B <- fit$draws(paste0("beta_B[", 1:10, "]"), format = "array")
dimnames(draws_B)$variable <- paste0("B: ", beta_labels)
p_trace_B <- mcmc_trace(draws_B) + ggtitle("Traceplot beta_B") + theme_minimal(base_size = 8)
print(p_trace_B)
ggsave(file.path(GRAFICI_DIR, "traceplot_beta_B.png"), plot = p_trace_B,
       width = 10, height = 7, dpi = 300)

draws_C <- fit$draws(paste0("beta_C[", 1:10, "]"), format = "array")
dimnames(draws_C)$variable <- paste0("C: ", beta_labels)
p_trace_C <- mcmc_trace(draws_C) + ggtitle("Traceplot beta_C") + theme_minimal(base_size = 8)
print(p_trace_C)
ggsave(file.path(GRAFICI_DIR, "traceplot_beta_C.png"), plot = p_trace_C,
       width = 10, height = 7, dpi = 300)

p_trace_str <- mcmc_trace(fit$draws(c("sigma_C", "sigma", "delta"), format = "array")) +
    ggtitle("Traceplot parametri strutturali") + theme_minimal()
print(p_trace_str)
ggsave(file.path(GRAFICI_DIR, "traceplot_strutturali.png"), plot = p_trace_str,
       width = 8, height = 5, dpi = 300)

# ==============================================================================
# Sezione 5 — Posterior predictive check
# ==============================================================================

Y_vec     <- dati$OD
Y_rep     <- fit$draws("OD_rep", format = "matrix")
Y_rep_sub <- Y_rep[sample(nrow(Y_rep), 200, replace = FALSE), ]

p_ppc_dens <- ppc_dens_overlay(Y_vec, Y_rep_sub) +
  labs(title = "PPC — distribuzione marginale") + theme_minimal()
print(p_ppc_dens)
ggsave(file.path(GRAFICI_DIR, "ppc_densita.png"), plot = p_ppc_dens,
       width = 7, height = 5, dpi = 300)

p_ppc_sd <- ppc_stat_grouped(Y_vec, Y_rep_sub, group = dati$tempo_t, stat = "sd") +
  labs(title = "PPC — SD per settimana") + theme_minimal()
print(p_ppc_sd)
ggsave(file.path(GRAFICI_DIR, "ppc_sd_settimana.png"), plot = p_ppc_sd,
       width = 8, height = 5, dpi = 300)

p_ppc_med <- ppc_stat_grouped(Y_vec, Y_rep_sub, group = dati$tempo_t, stat = "median") +
  labs(title = "PPC — Mediana per settimana") + theme_minimal()
print(p_ppc_med)
ggsave(file.path(GRAFICI_DIR, "ppc_mediana_settimana.png"), plot = p_ppc_med,
       width = 8, height = 5, dpi = 300)

mu_hat <- colMeans(Y_rep)
acf_df <- dati |>
  mutate(residuo = OD - mu_hat) |>
  group_by(id_biomassa) |>
  arrange(tempo_t) |>
  summarise(acf1 = if (n() > 2) cor(residuo[-n()], residuo[-1],
                                     use = "complete") else NA_real_,
            .groups = "drop")

p_acf <- ggplot(acf_df, aes(x = acf1)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick") +
  labs(title = "ACF lag-1 residui per curva",
       x = "Autocorrelazione lag-1", y = "Curve") + theme_minimal()
print(p_acf)
ggsave(file.path(GRAFICI_DIR, "acf_lag1_residui.png"), plot = p_acf,
       width = 6, height = 5, dpi = 300)
cat("Mediana ACF lag-1:", round(median(acf_df$acf1, na.rm = TRUE), 3), "\n")

# ==============================================================================
# Sezione 6 — Posteriori coefficienti RSM
# ==============================================================================

posterior_summary <- function(draws_mat, labels, component) {
  as.data.frame(t(apply(draws_mat, 2,
                        quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))) |>
    setNames(c("q05", "q25", "med", "q75", "q95")) |>
    mutate(parametro = labels, componente = component)
}

d_A <- fit$draws(paste0("beta_A[", 1:10, "]"), format = "matrix")
d_B <- fit$draws(paste0("beta_B[", 1:10, "]"), format = "matrix")
d_C <- fit$draws(paste0("beta_C[", 1:10, "]"), format = "matrix")

lbl_main <- c("xI", "xD", "xP")
lbl_int  <- c("xI:xD", "xI:xP", "xD:xP")
lbl_quad <- c("xI^2", "xD^2", "xP^2")

for (grp in list(list(2:4, lbl_main, "Effetti principali",   "coeff_principali"),
                 list(5:7, lbl_int,  "Interazioni",          "coeff_interazioni"),
                 list(8:10, lbl_quad, "Quadratici",           "coeff_quadratici"))) {
  idx <- grp[[1]]; lbl <- grp[[2]]; titolo <- grp[[3]]; fname <- grp[[4]]
  df <- bind_rows(
    posterior_summary(d_A[, idx], lbl, "log(Asym)"),
    posterior_summary(d_B[, idx], lbl, "log(b2)"),
    posterior_summary(d_C[, idx], lbl, "logit(b3)")
  )
  p_coeff <- ggplot(df, aes(y = parametro, x = med, color = componente)) +
    geom_pointrange(aes(xmin = q05, xmax = q95),
                    position = position_dodge(0.5), linewidth = 0.7) +
    geom_pointrange(aes(xmin = q25, xmax = q75),
                    position = position_dodge(0.5), linewidth = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("log(Asym)" = "#fd8d3c",
                                  "log(b2)"   = "#31a354",
                                  "logit(b3)" = "#756bb1")) +
    labs(title = titolo,
         subtitle = "Spesso: IC 50% | Sottile: IC 90%",
         x = NULL, y = NULL, color = NULL) +
    theme_minimal() + theme(legend.position = "bottom")
  print(p_coeff)
  ggsave(file.path(GRAFICI_DIR, paste0(fname, ".png")), plot = p_coeff,
         width = 7, height = 5, dpi = 300)
}

# ==============================================================================
# Sezione 7 — Posteriori parametri strutturali
# ==============================================================================

d_str <- fit$draws(c("sigma_C", "sigma", "delta"), format = "df")

p_sC <- ggplot(d_str, aes(x = sigma_C)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = median(d_str$sigma_C), linetype = "dashed") +
  annotate("text", x = median(d_str$sigma_C), y = Inf,
           label = paste0("med=", round(median(d_str$sigma_C), 3)),
           vjust = 1.5, hjust = -0.1, size = 3.5) +
  labs(title = "sigma_C (SD random effect su b3)") + theme_minimal()

p_si <- ggplot(d_str, aes(x = sigma)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = median(d_str$sigma), linetype = "dashed") +
  annotate("text", x = median(d_str$sigma), y = Inf,
           label = paste0("med=", round(median(d_str$sigma), 3)),
           vjust = 1.5, hjust = -0.1, size = 3.5) +
  labs(title = "sigma (scala SD residua)") + theme_minimal()

p_de <- ggplot(d_str, aes(x = delta)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = median(d_str$delta), linetype = "dashed") +
  geom_vline(xintercept = 2, color = "firebrick", linetype = "dotted") +
  annotate("text", x = 2, y = Inf, label = "delta=2\n(CV cost.)",
           vjust = 1.5, hjust = 1.1, size = 3, color = "firebrick") +
  labs(title = "delta (esponente power-law)") + theme_minimal()

p_strutturali <- p_sC / p_si / p_de
print(p_strutturali)
ggsave(file.path(GRAFICI_DIR, "parametri_strutturali.png"), plot = p_strutturali,
       width = 6, height = 9, dpi = 300)

cat("\n=== MODELLO B COMPLETATO ===\n")

# ==============================================================================
# Sezione 7b — Residui standardizzati vs valori predetti
#
# Verifica che il modello eteroschedastico abbia catturato correttamente la
# struttura della varianza. Usiamo le medie posteriori di tutti i parametri
# per calcolare mu_ij, poi standardizziamo per sigma * mu_ij^(delta/2).
# Se il modello e' ben specificato, i residui standardizzati devono avere
# dispersione costante al variare di mu (nessun imbuto residuo).
# ==============================================================================

# Estraiamo le medie posteriori dei parametri
draws_df <- fit$draws(format = "df")

get_post_mean <- function(pattern) {
  cols <- grep(pattern, names(draws_df), value = TRUE)
  colMeans(draws_df[, cols, drop = FALSE])
}

beta_A_pm <- get_post_mean("^beta_A\\[")
beta_B_pm <- get_post_mean("^beta_B\\[")
beta_C_pm <- get_post_mean("^beta_C\\[")
v_pm      <- get_post_mean("^v\\[")
sigma_pm  <- mean(draws_df$sigma)
delta_pm  <- mean(draws_df$delta)

# Ricalcoliamo mu_ij per ogni osservazione usando i valori posteriori medi
df_resid_std <- dati |>
  mutate(
    row_k   = match(id_biomassa, curve_ids),
    eta_A   = as.numeric(X_full[match(row_number(), seq_len(nrow(dati))), ] %*% beta_A_pm),
    eta_B   = as.numeric(X_full[match(row_number(), seq_len(nrow(dati))), ] %*% beta_B_pm),
    eta_C   = as.numeric(X_full[match(row_number(), seq_len(nrow(dati))), ] %*% beta_C_pm),
    Asym_i  = exp(eta_A),
    b2_i    = exp(eta_B),
    b3_i    = plogis(eta_C + v_pm[row_k]),
    mu_ij   = Asym_i * exp(-b2_i * b3_i^tempo_t),
    sd_ij   = sigma_pm * mu_ij^(delta_pm / 2),
    resid_std = (OD - mu_ij) / sd_ij
  )

p_resid_std <- ggplot(df_resid_std, aes(x = mu_ij, y = resid_std)) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_hline(yintercept = 0,  linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept =  2, linetype = "dotted", color = "gray50",    linewidth = 0.6) +
  geom_hline(yintercept = -2, linetype = "dotted", color = "gray50",    linewidth = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick",
              linewidth = 1, span = 1) +
  labs(
    title    = "Residui standardizzati vs valori predetti",
    subtitle = "Media posteriore dei parametri | Linee punteggiate: +/-2 SD | Se piatto -> varianza catturata",
    x = expression(hat(mu)[ij] ~ "(media posteriore)"),
    y = expression(frac(OD[ij] - hat(mu)[ij], hat(sigma) %.% hat(mu)[ij]^{hat(delta)/2}))
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

print(p_resid_std)
ggsave(file.path(GRAFICI_DIR, "residui_std_vs_fitted.png"), plot = p_resid_std,
       width = 7, height = 5, dpi = 300)

# ==============================================================================
# Sezione 8 — Confronto tra combinazioni osservate
#
# Per ogni combinazione sperimentale osservata (I, D, P) si calcola la
# distribuzione posteriore di OD a t=6 e t=7 con v=0 (curva media).
# Per ogni coppia di combinazioni (i, j) si calcola:
#   - differenza media posteriore:  mean(OD_i - OD_j)
#   - probabilita' che i batta j:   P(OD_i > OD_j) = mean(OD_i > OD_j)
# Risultato: due heatmap 17x17 per ogni t* (prob e differenza media).
# ==============================================================================

cat("\n=== CONFRONTO COMBINAZIONI OSSERVATE ===\n")

# --- 8.1 Estrai campioni posteriori (tutti gli 8000) ---

draws_conf  <- fit$draws(format = "matrix")
bA_conf     <- draws_conf[, grep("^beta_A\\[", colnames(draws_conf))]
bB_conf     <- draws_conf[, grep("^beta_B\\[", colnames(draws_conf))]
bC_conf     <- draws_conf[, grep("^beta_C\\[", colnames(draws_conf))]

# --- 8.2 Identifica le combinazioni uniche osservate ---
# Ogni combinazione unica di (I, D, P) ha un vettore x identico.
# Estrae le righe uniche da X_per_curve e la loro etichetta leggibile.

comb_df <- dati |>
  distinct(I, D, P, xI, xD, xP) |>
  arrange(I, D, P) |>
  mutate(
    label = paste0("I=", I, "\nD=", D, "\nP=", P),
    label_short = paste0(I, "/", D, "/", P)
  )

n_comb <- nrow(comb_df)
cat("Combinazioni uniche osservate:", n_comb, "\n")

# Costruisce X per le combinazioni uniche
X_comb <- model.matrix(
  ~ xI + xD + xP + xI:xD + xI:xP + xD:xP + I(xI^2) + I(xD^2) + I(xP^2),
  data = comb_df
)  # n_comb × 10

# --- 8.3 Calcola distribuzione posteriore di OD per ogni combinazione e t* ---
# eta_X: 8000 campioni × n_comb combinazioni

eta_A_c <- bA_conf %*% t(X_comb)   # 8000 × n_comb
eta_B_c <- bB_conf %*% t(X_comb)
eta_C_c <- bC_conf %*% t(X_comb)   # v = 0

plot_heatmap_confronto <- function(t_star) {

  # OD_mat: 8000 × n_comb
  OD_mat <- exp(eta_A_c) * exp(-exp(eta_B_c) * plogis(eta_C_c)^t_star)

  # --- 8.4 Calcola prob e diff media per ogni coppia (i, j) ---
  # prob[i,j]  = P(OD_i > OD_j)
  # mdiff[i,j] = mean(OD_i - OD_j)

  prob_mat  <- matrix(NA, n_comb, n_comb)
  mdiff_mat <- matrix(NA, n_comb, n_comb)

  for (i in seq_len(n_comb)) {
    for (j in seq_len(n_comb)) {
      if (i == j) next
      diff_ij          <- OD_mat[, i] - OD_mat[, j]
      prob_mat[i, j]   <- mean(diff_ij > 0)
      mdiff_mat[i, j]  <- mean(diff_ij)
    }
  }

  lbl <- comb_df$label_short

  # --- 8.5 Heatmap probabilita' ---
  prob_df <- as.data.frame(prob_mat) |>
    setNames(lbl) |>
    mutate(riga = lbl) |>
    pivot_longer(-riga, names_to = "colonna", values_to = "prob") |>
    mutate(
      riga    = factor(riga,    levels = lbl),
      colonna = factor(colonna, levels = lbl)
    )

  p_prob <- ggplot(prob_df, aes(x = colonna, y = riga, fill = prob)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(prob), "",
                                 sprintf("%.2f", prob))),
              size = 2.2, color = "white") +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#d6604d",
      midpoint = 0.5,
      na.value = "grey85",
      limits   = c(0, 1),
      name     = "P(riga > col)"
    ) +
    labs(
      title    = paste0("P(OD_riga > OD_colonna) — t* = ", t_star),
      subtitle = "Rosso: riga batte colonna | Blu: colonna batte riga",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())

  # --- 8.6 Heatmap differenza media ---
  mdiff_df <- as.data.frame(mdiff_mat) |>
    setNames(lbl) |>
    mutate(riga = lbl) |>
    pivot_longer(-riga, names_to = "colonna", values_to = "mdiff") |>
    mutate(
      riga    = factor(riga,    levels = lbl),
      colonna = factor(colonna, levels = lbl)
    )

  lim_d <- max(abs(mdiff_df$mdiff), na.rm = TRUE)

  p_diff <- ggplot(mdiff_df, aes(x = colonna, y = riga, fill = mdiff)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(mdiff), "",
                                 sprintf("%.1f", mdiff))),
              size = 2.2, color = "white") +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#d6604d",
      midpoint = 0,
      na.value = "grey85",
      limits   = c(-lim_d, lim_d),
      name     = "OD_riga - OD_col"
    ) +
    labs(
      title    = paste0("Differenza media OD (riga - colonna) — t* = ", t_star),
      subtitle = "Rosso: riga produce piu' OD | Blu: colonna produce piu' OD",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())

  p_heatmap <- p_prob / p_diff
  print(p_heatmap)
  ggsave(file.path(GRAFICI_DIR, paste0("heatmap_confronto_t", t_star, ".png")),
         plot = p_heatmap, width = 10, height = 14, dpi = 300)

  # --- 8.7 Tabella ranking: OD medio per combinazione ---
  ranking <- data.frame(
    combinazione = lbl,
    OD_med       = colMeans(OD_mat),
    OD_q05       = apply(OD_mat, 2, quantile, 0.05),
    OD_q95       = apply(OD_mat, 2, quantile, 0.95)
  ) |>
    arrange(desc(OD_med))

  cat("\n--- Ranking combinazioni a t* =", t_star, "---\n")
  print(ranking, digits = 3)

  # --- 8.8 Curve medie stimate, un colore per combinazione ---
  # Ordine legenda: dalla combinazione con OD piu' alto a t* (in cima)
  # alla piu' bassa (in fondo). Colori discreti, uno per combinazione.

  eta_A_med <- apply(eta_A_c, 2, median)
  eta_B_med <- apply(eta_B_c, 2, median)
  eta_C_med <- apply(eta_C_c, 2, median)

  t_seq <- seq(1, 7, by = 0.1)

  # Ordine delle combinazioni per OD decrescente a t*
  ord_lbl <- ranking$combinazione   # gia' ordinato desc per OD_med

  curve_df <- lapply(seq_len(n_comb), function(k) {
    mu_t <- exp(eta_A_med[k]) *
            exp(-exp(eta_B_med[k]) * plogis(eta_C_med[k])^t_seq)
    data.frame(
      t            = t_seq,
      OD           = mu_t,
      combinazione = comb_df$label_short[k]
    )
  }) |>
    bind_rows() |>
    mutate(combinazione = factor(combinazione, levels = ord_lbl))

  # Palette discreta blu → rosso (coerente con heatmap):
  # rosso = OD più alto (in cima alla legenda), blu = OD più basso
  pal <- colorRampPalette(c("#2166ac", "#d6604d"))(n_comb)
  pal <- rev(pal)   # primo = rosso (OD max), ultimo = blu (OD min)
  names(pal) <- ord_lbl

  p_curve <- ggplot(curve_df, aes(x = t, y = OD,
                                   group = combinazione,
                                   color = combinazione)) +
    geom_line(linewidth = 0.9, alpha = 0.85) +
    scale_color_manual(
      values = pal,
      name   = "Combinazione\n(I/D/P)",
      guide  = guide_legend(reverse = FALSE,
                            keyheight = unit(0.5, "cm"),
                            ncol = 1)
    ) +
    scale_x_continuous(breaks = 1:7) +
    labs(
      title    = "Curve medie stimate per combinazione",
      subtitle = "Legenda ordinata da OD piu' alto (in cima) a piu' basso (in fondo)",
      x = "Tempo (settimane)", y = "OD predetto (mediana posteriore)"
    ) +
    theme_minimal() +
    theme(legend.position  = "right",
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 8))

  print(p_curve)
  ggsave(file.path(GRAFICI_DIR, paste0("curve_medie_t", t_star, ".png")),
         plot = p_curve, width = 9, height = 6, dpi = 300)

  invisible(list(prob = prob_mat, mdiff = mdiff_mat, ranking = ranking))
}

res_t6 <- plot_heatmap_confronto(t_star = 6)
res_t7 <- plot_heatmap_confronto(t_star = 7)

cat("\n=== CONFRONTO COMBINAZIONI COMPLETATO ===\n")
