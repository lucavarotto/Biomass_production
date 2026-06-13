#rm(list = ls()); gc()

library(readxl)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(minpack.lm)   # nlsLM: più stabile di nls per Gompertz

# Percorsi automatici basati sulla posizione dello script
SCRIPT_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Cartella output grafici
GRAFICI_DIR <- file.path(SCRIPT_DIR, "grafici")
dir.create(GRAFICI_DIR, showWarnings = FALSE, recursive = TRUE)

# Costanti grafiche condivise
PALETTE <- "plasma"
DPI     <- 300

# Palette per covariata (chiaro → scuro = livello basso → alto)
COL_I <- c("#fdd0a2", "#fd8d3c", "#d94801")  # arancione: I = 90, 280, 470 lux
COL_D <- c("#a1d99b", "#31a354", "#006d2c")  # verde:     D = 12, 18, 24 ore
COL_P <- c("#bcbddc", "#756bb1", "#54278f")  # viola:     P = 1, 25.5, 50 mg/L
COL_REF <- "firebrick"   # colore unico per linee di riferimento
COL_PT  <- "steelblue"   # colore unico per punti singoli
COL_GREY <- "gray50"     # grigio per linee accessorie

# =============================================================================
# 1. CARICAMENTO E PREPARAZIONE DEI DATI
# =============================================================================

dati <- readxl::read_xlsx(file.path(SCRIPT_DIR, "nostoc.xlsx"), skip = 2)
colnames(dati) <- c("row_id", "I", "D", "P", "tempo", "OD")

dati$OD <- as.numeric(dati$OD)

# id univoco per le 53 biomasse
dati$id_biomassa <- rep(1:53, each = 7)

# Etichetta combinazione fattori
dati$condizione_sperimentale <- paste(dati$I, dati$D, dati$P, sep = "_")

# Consumo energetico totale (I x D)
dati$consumo <- as.factor(dati$I * dati$D)

# Traslazione tempo da [0-6] a [1-7] per evitare b3^0 = 1 (non identificabile)
if (min(dati$tempo, na.rm = TRUE) == 0) dati$tempo <- dati$tempo + 1

# Valore iniziale N0 (OD a t=1 dopo traslazione)
dati$N0 <- rep(dati[dati$tempo == 1, ]$OD, each = 7)

# Incrementi assoluti e rapporti settimanali
dati_incrementi <- dati |>
  group_by(id_biomassa) |>
  arrange(tempo) |>
  mutate(
    incremento_OD = OD - lag(OD),
    intervallo    = paste0(lag(tempo), "-", tempo)
  ) |>
  filter(!is.na(incremento_OD)) |>
  ungroup()

dati_incrementi_rapporti <- dati |>
  group_by(id_biomassa) |>
  arrange(tempo) |>
  mutate(
    incremento_OD = OD / lag(OD),
    intervallo    = paste0(lag(tempo), "/", tempo)
  ) |>
  filter(!is.na(incremento_OD)) |>
  ungroup()

save(dati, dati_incrementi, dati_incrementi_rapporti,
     file = file.path(SCRIPT_DIR, "dati_modificati.Rdata"))

cat("Distribuzioni marginali:\n")
print(table(dati$condizione_sperimentale) / 7)
print(table(dati$I) / 7)
print(table(dati$D) / 7)
print(table(dati$P) / 7)

# =============================================================================
# 2. PANORAMICA DELLE CURVE
# =============================================================================

# Tutte le curve (senza colore, pulito)
p_tutte <- ggplot(dati, aes(x = tempo, y = OD, group = id_biomassa)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE, alpha = 1, size = 1) +
  labs(x = "Tempo (settimane)", y = "Resa di Biomassa (OD)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(GRAFICI_DIR, "tutte_curve.png"),
       plot = p_tutte, width = 6, height = 5, dpi = DPI)

# Tutte le curve + media globale in grassetto
p_media_globale <- ggplot(dati, aes(x = tempo, y = OD)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  labs(
    title    = "Traiettorie di crescita — tutte le curve",
    subtitle = "Linea spessa: media globale",
    x = "Tempo (settimane)",
    y = "Resa di Biomassa (OD)"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(GRAFICI_DIR, "media_globale.png"),
       plot = p_media_globale, width = 6, height = 5, dpi = DPI)

# Tutte le curve colorate per condizione sperimentale
p_curve_condizione <- ggplot(dati, aes(x = tempo, y = OD,
                                       color = as.factor(condizione_sperimentale),
                                       group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  scale_color_viridis_d(option = PALETTE, end = 0.9) +
  labs(
    x     = "Tempo (settimane)",
    y     = "Optical Density (OD)",
    color = "Condizione Sperimentale"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(GRAFICI_DIR, "tutte_curve_condizione.png"),
       plot = p_curve_condizione, width = 6, height = 5, dpi = DPI)

# Curve medie per condizione (sfondo grigio + media colorata + ribbon ±SE)
p_curve_medie <- ggplot(dati, aes(x = tempo, y = OD,
                                   color = as.factor(condizione_sperimentale),
                                   fill  = as.factor(condizione_sperimentale))) +
  geom_line(aes(group = id_biomassa), color = "gray85", alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(group = id_biomassa), color = "gray85", alpha = 0.5, size = 1, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line",  linewidth = 0.7, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 1, na.rm = TRUE) +
  scale_color_viridis_d(option = PALETTE, end = 0.9) +
  scale_fill_viridis_d(option = PALETTE,  end = 0.9) +
  labs(
    x     = "Tempo (settimane)",
    y     = "Optical Density (OD)",
    color = "Condizione Sperimentale",
    fill  = "Condizione Sperimentale"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(GRAFICI_DIR, "Curve_medie_condizione.png"),
       plot = p_curve_medie, width = 6, height = 5, dpi = DPI)

# =============================================================================
# 3. GRAFICI CONDIZIONATI
# =============================================================================

grafico_condizionate <- function(condizionante, f1, f2 = NULL, colori) {
  if (!is.null(f2)) {
    dati$gruppo_colore <- paste(dati[[f1]], dati[[f2]], sep = "_")
    label_colore <- paste(f1, f2, sep = "_")
  } else {
    dati$gruppo_colore <- as.factor(dati[[f1]])
    label_colore <- f1
  }

  ggplot(dati, aes(x = tempo, y = OD,
                   color = gruppo_colore,
                   fill  = gruppo_colore)) +
    geom_line(aes(group = id_biomassa), alpha = 0.2, na.rm = TRUE) +
    stat_summary(fun.data = mean_se, geom = "ribbon",
                 alpha = 0.15, color = NA, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "line",  linewidth = 1.2, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point", size = 2,        na.rm = TRUE) +
    facet_wrap(vars(.data[[condizionante]])) +
    scale_color_manual(values = colori) +
    scale_fill_manual(values  = colori) +
    labs(
      title = paste("Dinamica di crescita condizionata a", condizionante),
      x     = "Tempo (settimane)",
      y     = "Optical Density (OD)",
      color = label_colore,
      fill  = label_colore
    ) +
    theme_minimal()
}

p_cond_DI <- grafico_condizionate("I", "D", colori = COL_D)
ggsave(file.path(GRAFICI_DIR, "condizionata_D_I.png"),
       plot = p_cond_DI, width = 9, height = 5, dpi = DPI)

p_cond_DP <- grafico_condizionate("P", "D", colori = COL_D)
ggsave(file.path(GRAFICI_DIR, "condizionata_D_P.png"),
       plot = p_cond_DP, width = 9, height = 5, dpi = DPI)

p_cond_PI <- grafico_condizionate("I", "P", colori = COL_P)
ggsave(file.path(GRAFICI_DIR, "condizionata_P_I.png"),
       plot = p_cond_PI, width = 9, height = 5, dpi = DPI)

# =============================================================================
# 4. GRAFICI MARGINALI (I, D, P)
# =============================================================================

marginal_plot <- function(fattore, label, colori) {
  ggplot(dati, aes(x = tempo, y = OD,
                   color = as.factor(.data[[fattore]]),
                   fill  = as.factor(.data[[fattore]]))) +
    geom_line(aes(group = id_biomassa), color = "gray75", alpha = 0.6, na.rm = TRUE) +
    geom_point(aes(group = id_biomassa), color = "gray75", alpha = 0.6, size = 0.8, na.rm = TRUE) +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.25, color = NA, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "line",  linewidth = 1.2, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point", size = 2,        na.rm = TRUE) +
    scale_color_manual(values = colori) +
    scale_fill_manual(values  = colori) +
    labs(x = "Tempo (settimane)", y = "Optical Density (OD)",
         color = label, fill = label) +
    theme_minimal()
}

marginale_I <- marginal_plot("I", "I (lux)",  COL_I)
marginale_D <- marginal_plot("D", "D (ore)",  COL_D)
marginale_P <- marginal_plot("P", "P (mg/L)", COL_P)

ggsave(
  file.path(GRAFICI_DIR, "marginali.png"),
  plot  = arrangeGrob(marginale_I, marginale_D, marginale_P, ncol = 3),
  width = 18, height = 5, dpi = DPI
)

# =============================================================================
# 5. DATI MANCANTI
# =============================================================================

na_ids <- dati$id_biomassa[is.na(dati$OD)]

p_na_zoom <- ggplot(
  subset(dati, id_biomassa %in% na_ids),
  aes(x = tempo, y = OD, group = id_biomassa)
) +
  geom_line(color = COL_REF, linewidth = 1.0, alpha = 0.9, na.rm = TRUE) +
  geom_point(color = COL_REF, alpha = 0.7, na.rm = TRUE) +
  labs(
    title = "Unità con dati mancanti",
    x     = "Tempo (settimane)",
    y     = "Optical Density (OD)"
  ) +
  ylim(0, max(dati$OD, na.rm = TRUE)) +
  theme_minimal()

p_na_contesto <- ggplot(dati, aes(x = tempo, y = OD, group = id_biomassa)) +
  geom_line(data = subset(dati, !(id_biomassa %in% na_ids)),
            color = "gray70", linewidth = 0.4, alpha = 0.6, na.rm = TRUE) +
  geom_line(data = subset(dati, id_biomassa %in% na_ids),
            color = COL_REF, linewidth = 1.0, alpha = 0.9, na.rm = TRUE) +
  geom_point(data = subset(dati, id_biomassa %in% na_ids),
             color = COL_REF, alpha = 0.7, na.rm = TRUE) +
  labs(
    title = "Posizione nel dataset completo",
    x     = "Tempo (settimane)",
    y     = "Optical Density (OD)",
    color = "Unità"
  ) +
  ylim(0, max(dati$OD, na.rm = TRUE)) +
  theme_minimal()

ggsave(
  file.path(GRAFICI_DIR, "dati_mancanti.png"),
  plot  = arrangeGrob(p_na_zoom, p_na_contesto, ncol = 2),
  width = 12, height = 5, dpi = DPI
)

# =============================================================================
# 6. ETEROSCHEDASTICITÀ E RELAZIONE MEDIA-VARIANZA
# =============================================================================

# Impatto del consumo energetico (I x D) sulla crescita
p_consumo <- ggplot(dati, aes(x = tempo, y = OD, color = consumo)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  scale_color_viridis_d(option = PALETTE, end = 0.9) +
  labs(
    title    = "Impatto del consumo energetico (I × D) sulla crescita",
    subtitle = "Linee spesse: media per livello di consumo",
    x        = "Tempo (settimane)",
    y        = "Resa di Biomassa (OD)",
    color    = "Consumo (I×D)"
  ) +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave(file.path(GRAFICI_DIR, "consumo_energetico.png"),
       plot = p_consumo, width = 6, height = 5, dpi = DPI)

# Effetto del valore iniziale N0 sulla traiettoria
p_n0 <- ggplot(dati, aes(x = tempo, y = OD, colour = N0)) +
  geom_line(aes(group = id_biomassa), na.rm = TRUE, linewidth = 1.5) +
  facet_wrap(~I) +
  scale_color_viridis_c(option = PALETTE, end = 0.9) +
  labs(
    title    = "Effetto del valore iniziale (N0) sulla traiettoria di crescita",
    subtitle = "Pannelli per intensità luminosa I | Colore: OD a t=1",
    x        = "Tempo (settimane)",
    y        = "Resa di Biomassa (OD)",
    color    = "N0 (OD iniziale)"
  ) +
  theme_minimal() +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave(file.path(GRAFICI_DIR, "effetto_n0.png"),
       plot = p_n0, width = 8, height = 5, dpi = DPI)

# Relazione media-varianza (scala originale)
df_summary <- dati |>
  group_by(tempo) |>
  summarise(media_OD = mean(OD, na.rm = TRUE), varianza_OD = var(OD, na.rm = TRUE))

t_var <- ggplot(df_summary, aes(x = tempo, y = varianza_OD)) +
  geom_point(size = 2.5, color = COL_PT, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = COL_GREY, linewidth = 0.8) +
  geom_smooth(method = "loess", color = COL_REF, linewidth = 1.2,
              se = FALSE, span = 1) +
  labs(x = "Tempo", y = "Varianza di OD") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

m_var <- ggplot(df_summary, aes(x = media_OD, y = varianza_OD)) +
  geom_point(size = 2.5, color = COL_PT, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = COL_GREY, linewidth = 0.8) +
  geom_smooth(method = "loess", color = COL_REF, linewidth = 1.2,
              se = FALSE, span = 1) +
  labs(x = "Media di OD", y = "Varianza di OD") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(
  file.path(GRAFICI_DIR, "relazione_m_var.png"),
  plot  = arrangeGrob(t_var, m_var, ncol = 2),
  width = 12, height = 5, dpi = DPI
)

# Stima empirica di delta in scala log-log
# Se Var = sigma^2 * mu^delta, allora log(Var) = cost + delta * log(mu).
df_loglog <- df_summary |>
  filter(media_OD > 0, varianza_OD > 0) |>
  mutate(log_media = log(media_OD), log_var = log(varianza_OD))

fit_loglog <- lm(log_var ~ log_media, data = df_loglog)
delta_emp  <- round(coef(fit_loglog)[2], 2)
cat("Stima empirica di delta (pendenza log-log):", delta_emp, "\n")

p_loglog <- ggplot(df_loglog, aes(x = log_media, y = log_var)) +
  geom_point(size = 2.5, color = COL_PT, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = COL_REF,
              linewidth = 1.2, fill = COL_REF, alpha = 0.15) +
  annotate("text", x = min(df_loglog$log_media), y = max(df_loglog$log_var),
           label = paste0("delta stimato = ", delta_emp),
           hjust = 0, vjust = 1, size = 4, color = COL_REF) +
  geom_abline(intercept = coef(fit_loglog)[1], slope = 1,
              linetype = "dashed", color = COL_GREY, linewidth = 0.8) +
  geom_abline(intercept = coef(fit_loglog)[1], slope = 2,
              linetype = "dotted", color = COL_GREY, linewidth = 0.8) +
  annotate("text", x = max(df_loglog$log_media) * 0.85,
           y = coef(fit_loglog)[1] + max(df_loglog$log_media) * 0.85 + 0.1,
           label = "delta=1", size = 3, color = "gray50") +
  annotate("text", x = max(df_loglog$log_media) * 0.7,
           y = coef(fit_loglog)[1] + 2 * max(df_loglog$log_media) * 0.7 + 0.1,
           label = "delta=2", size = 3, color = "gray50") +
  labs(
    title    = "Relazione media-varianza in scala log-log",
    subtitle = paste0("Pendenza stimata (delta) = ", delta_emp,
                      " | Tratteggio: delta=1 | Puntinato: delta=2"),
    x = "log(Media OD)",
    y = "log(Varianza OD)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(GRAFICI_DIR, "relazione_m_var_loglog.png"),
       plot = p_loglog, width = 7, height = 5, dpi = DPI)

# Residui grezzi vs fitted (Gompertz per curva, stima indipendente)
df_resid_raw <- do.call(rbind, lapply(unique(dati$id_biomassa), function(id) {
  sub <- dati[dati$id_biomassa == id & !is.na(dati$OD), ]
  if (nrow(sub) < 4) return(NULL)
  tryCatch({
    fit_tmp <- nlsLM(
      OD ~ Asym * exp(-b2 * b3^tempo),
      data  = sub,
      start = list(Asym = max(sub$OD, na.rm = TRUE), b2 = 5, b3 = 0.7),
      lower = c(0, 0, 0),
      upper = c(Inf, Inf, 1),
      control = nls.lm.control(maxiter = 200)
    )
    data.frame(id_biomassa = id,
               fitted      = fitted(fit_tmp),
               residuo     = residuals(fit_tmp))
  }, error = function(e) NULL)
}))

p_resid_raw <- ggplot(df_resid_raw, aes(x = fitted, y = residuo)) +
  geom_point(alpha = 0.4, size = 1.5, color = COL_PT) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF, linewidth = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = COL_REF,
              linewidth = 1, span = 1) +
  labs(
    title = "Residui grezzi vs valori predetti",
    x = "Valori predetti (fitted)",
    y = "Residuo (OD − fitted)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(GRAFICI_DIR, "residui_grezzi_vs_fitted.png"),
       plot = p_resid_raw, width = 7, height = 5, dpi = DPI)

# =============================================================================
# 7. INCREMENTI SETTIMANALI
# =============================================================================

# Differenze assolute per condizione sperimentale
diff <- ggplot(dati_incrementi,
               aes(x = tempo, y = incremento_OD,
                   color = condizione_sperimentale,
                   fill  = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF, linewidth = 0.8) +
  scale_color_viridis_d(option = PALETTE) +
  scale_fill_viridis_d(option = PALETTE) +
  labs(
    x     = "Tempo (settimane)",
    y     = bquote(Delta ~ "OD"),
    color = "Condizione sperimentale",
    fill  = "Condizione sperimentale"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Rapporti per condizione sperimentale
rapp <- ggplot(dati_incrementi_rapporti,
               aes(x = tempo, y = incremento_OD,
                   color = condizione_sperimentale,
                   fill  = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = COL_REF, linewidth = 0.8) +
  scale_color_viridis_d(option = PALETTE) +
  scale_fill_viridis_d(option = PALETTE) +
  labs(
    x     = "Tempo (settimane)",
    y     = bquote(frac(OD[t], OD[t - 1])),
    color = "Condizione sperimentale",
    fill  = "Condizione sperimentale"
  ) +
  theme_minimal()

ggsave(
  file.path(GRAFICI_DIR, "incrementi.png"),
  plot  = arrangeGrob(diff, rapp, ncol = 2),
  width = 12, height = 5, dpi = DPI
)

# Differenze per consumo energetico (I x D)
p_incrementi_consumo <- ggplot(dati_incrementi,
                                aes(x = tempo, y = incremento_OD, color = consumo)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = COL_REF, linewidth = 1) +
  scale_color_viridis_d(option = PALETTE) +
  labs(
    title    = "Incrementi settimanali di OD per consumo energetico (I × D)",
    subtitle = "Linea spessa: media per livello di consumo | Tratteggio rosso: incremento zero",
    x        = "Tempo (settimane)",
    y        = expression(Delta ~ "OD"),
    color    = "Consumo (I × D)"
  ) +
  theme_minimal()

ggsave(file.path(GRAFICI_DIR, "incrementi_consumo.png"),
       plot = p_incrementi_consumo, width = 6, height = 5, dpi = DPI)

# =============================================================================
# 8. ISPEZIONE OUTLIER
# =============================================================================

p_outlier <- ggplot(subset(dati, condizione_sperimentale == "470_18_25.5"),
                    aes(x = tempo, y = OD, group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  labs(
    title    = "Ispezione outlier: condizione I=470, D=18, P=25.5",
    subtitle = "Curve replicate della stessa condizione sperimentale",
    x = "Tempo (settimane)",
    y = "Optical Density (OD)"
  ) +
  theme_minimal()

ggsave(file.path(GRAFICI_DIR, "outlier_470_18_255.png"),
       plot = p_outlier, width = 6, height = 5, dpi = DPI)

cat("\nGrafici esplorativi completati. Output in:", GRAFICI_DIR, "\n")
