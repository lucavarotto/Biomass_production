################################################################################
#                                                                              #
#   ANALISI_ESPLORATIVA.R                                                      #
#                                                                              #
#   Unione dei tre script di analisi del progetto, che prima erano separati:   #
#                                                                              #
#     Parte 1  <- Esplorative.R   preparazione dei dati e grafici esplorativi  #
#     Parte 2  <- Marco.R         approfondimenti, incrementi, normalita'      #
#     Parte 3  <- Riccardo.R      incrementi, GAM e superfici di risposta      #
#                                                                              #
#   (Giorgia.R era vuoto, quindi non ha contribuito nulla.)                    #
#                                                                              #
#   Le parti vanno eseguite nell'ordine: la Parte 1 costruisce gli oggetti     #
#   dati, dati_incrementi e dati_incrementi_rapporti che le altre due usano,   #
#   e li salva in ../Dati non serve, restano in memoria. La modellazione in    #
#   Stan vive a parte, in Stan/modello.R.                                      #
#                                                                              #
#   Librerie, percorsi e costanti grafiche sono dichiarati una volta sola qui  #
#   sotto: nel corpo delle tre parti non e' stata tolta nessuna altra riga.    #
#                                                                              #
################################################################################

rm(list = ls()); gc()

# --- Librerie -----------------------------------------------------------------
# dalla Parte 1
library(readxl)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(minpack.lm)   # nlsLM: piu' stabile di nls per Gompertz
# aggiunte dalla Parte 2
library(tidyverse)
library(MuMIn)
library(nlme)
# aggiunte dalla Parte 3
library(mgcv)
library(plotly)
# le librerie usate solo in singoli blocchi (brms, boot, ggforce,
# marginaleffects, tseries, viridis) restano caricate dove servono

library(readxl)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(minpack.lm)   # nlsLM: più stabile di nls per Gompertz

# Percorsi automatici basati sulla posizione dello script
SCRIPT_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Cartella output grafici
# i grafici esplorativi stanno nella cartella Figures/ del progetto
GRAFICI_DIR <- file.path(dirname(SCRIPT_DIR), "Figures")
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


################################################################################
#  PARTE 1 - PREPARAZIONE DEI DATI E ANALISI ESPLORATIVA
#  (dal file Esplorative.R)
################################################################################

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


################################################################################
#  PARTE 2 - APPROFONDIMENTI, INCREMENTI E VERIFICHE
#  (dal file Marco.R)
################################################################################

# Nota: qui c'era un load() di dati_modificati.Rdata, non piu' necessario
# perche' la Parte 1 lascia gli stessi oggetti in memoria.


# # Caricamento e pulizia
# dati <- readxl::read_xlsx("Codice/nostoc.xlsx", skip = 2)
# # Rinominiamo le colonne (basandoci sulla struttura del file)
# colnames(dati) <- c("row_id", "I", "D", "P", "tempo", "OD")
# 
# # Trasformiamo OD in numerico (gestisce eventuali "-" o caratteri non validi)
# dati$OD <- as.numeric(dati$OD)
# 
# # Creiamo l'id_biomassa univoco per le 53 biomasse
# dati$id_biomassa <- rep(1:53, each = 7)
# 
# # Creiamo l'indicatore di gruppo (combinazione fattori)
# dati$condizione_sperimentale <- paste(dati$I,
#                                       dati$D,
#                                       dati$P,
#                                       sep = "_")
# 
# table(dati$condizione_sperimentale)/7
# 
# # Analisi esplorativa ----
# 
## Media globale ----

ggplot(dati, aes(x = tempo, y = OD)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  scale_color_viridis_d(option = "plasma", end = 0.9) +
  labs(
    title = "Impatto del Consumo Energetico (I x D) sulla Crescita",
    subtitle = "Le linee spesse indicano la media per ogni livello di intensità luminosa totale",
    x = "Tempo (giorni)",
    y = "Resa di Biomassa (OD)",
    color = "Consumo (I*D)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

## Marginale ----

library(ggplot2)

ggplot(dati, aes(x = tempo, y = OD, color = as.factor(I), group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per il fattore I",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattore I"
  ) + theme_minimal()
table(dati$I)/7

ggplot(dati, aes(x = tempo, y = OD, color = as.factor(D), group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per il fattore D",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattore D"
  ) + theme_minimal()
table(dati$D)/7

ggplot(dati, aes(x = tempo, y = OD, color = as.factor(P), group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per il fattore I",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattore P"
  ) + theme_minimal()
table(dati$P)/7

## Congiunta ----

ggplot(dati, aes(x = tempo, y = OD, color = condizione_sperimentale, group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per combinazione di fattori",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  theme_minimal()

ggplot(dati, aes(x = tempo, y = OD, color = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.2, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 2, na.rm = TRUE) +
  labs(
    title = "Traiettorie medie per condizione sperimentale",
    subtitle = "Le linee spesse rappresentano la media di gruppo; le linee chiare le singole biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  theme_minimal()

## Condizionate ----

dati$condizione_sperimentale2 <- paste(dati$I, dati$D, sep = "_")

ggplot(dati, aes(x = tempo, y = OD, color = condizione_sperimentale2)) +
  geom_line(aes(group = id_biomassa), alpha = 0.2, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 2, na.rm = TRUE) +
  labs(
    title = "Traiettorie medie per condizione sperimentale",
    subtitle = "Le linee spesse rappresentano la media di gruppo; le linee chiare le singole biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) + theme_minimal()

ggplot(dati, aes(x = tempo, y = OD, color = condizione_sperimentale2)) +
  geom_line(aes(group = id_biomassa), alpha = 0.2, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 2, na.rm = TRUE) +
  facet_wrap(~I)+
  labs(
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  theme_minimal()

ggplot(dati, aes(x = tempo, y = OD, color = condizione_sperimentale2)) +
  geom_line(aes(group = id_biomassa), alpha = 0.2, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 2, na.rm = TRUE) +
  facet_wrap(~P)+
  labs(
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  theme_minimal()

## Esplorazione outliers ----

ggplot(subset(dati, condizione_sperimentale=="470_18_25.5"), aes(x = tempo, y = OD, color = condizione_sperimentale, group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per combinazione di fattori",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  theme_minimal()

## Esplorazione NA ----

ggplot(subset(dati, id_biomassa %in% dati$id_biomassa[is.na(dati$OD)]), aes(x = tempo, y = OD, color = condizione_sperimentale, group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per combinazione di fattori",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  ylim(0,max(dati$OD, na.rm=T))+
  theme_minimal()

ggplot(dati, aes(x = tempo, y = OD, color = id_biomassa %in% dati$id_biomassa[is.na(dati$OD)], group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per combinazione di fattori",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  ylim(0,max(dati$OD, na.rm=T))+
  theme_minimal()

ggplot(subset(dati, I==90 & D==24 & (P == 1 | P==25.5) ), aes(x = tempo, y = OD, color = condizione_sperimentale, group = id_biomassa)) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) + # Opzionale: aggiunge i punti per chiarezza
  labs(
    title = "Traiettorie di crescita per combinazione di fattori",
    subtitle = "Ogni linea rappresenta una delle 53 biomasse",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattori (I_D_P)"
  ) +
  theme_minimal()






library(viridis)
dati$consumo <- as.factor(dati$I * dati$D)

ggplot(dati, aes(x = tempo, y = OD, color = consumo)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  scale_color_viridis_d(option = "plasma", end = 0.9) +
  labs(
    title = "Impatto del Consumo Energetico (I x D) sulla Crescita",
    subtitle = "Le linee spesse indicano la media per ogni livello di intensità luminosa totale",
    x = "Tempo (giorni)",
    y = "Resa di Biomassa (OD)",
    color = "Consumo (I*D)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Dati incrementali: differenze ----

library(dplyr)

dati_incrementi <- dati %>%
  group_by(id_biomassa) %>%
  arrange(tempo) %>%
  mutate(
    incremento_OD = OD - lag(OD),
    intervallo = paste0(lag(tempo), "-", tempo)
  ) |>
  filter(!is.na(incremento_OD)) |> # Rimuoviamo il tempo 0 (che non ha incremento)
  ungroup()

# Visualizzazione della struttura del nuovo dataset (6 righe per biomassa)
# str(dati_incrementi)

# 2. Grafico degli incrementi con palette Plasma
ggplot(dati_incrementi, aes(x = tempo, y = incremento_OD, color = consumo)) +
  # Linee sottili per le singole biomasse
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  
  # Linea spessa della media per livello di consumo
  stat_summary(fun = mean, geom = "line", size = 1.2, na.rm = TRUE) +
  
  # Palette Plasma (da viridis)
  scale_color_viridis_d(option = "plasma") +
  
  labs(
    title = "Velocità di crescita: Incrementi giornalieri di OD",
    subtitle = "L'incremento indica quanto la biomassa aumenta tra un giorno e il successivo",
    x = "Giorno (t)",
    y = "Delta OD (OD[t] - OD[t-1])",
    color = "Consumo (I*D)"
  ) +
  theme_minimal()


# Dati incrementali: rapporti ----

dati_incrementi_rapporti <- dati %>%
  group_by(id_biomassa) %>%
  arrange(tempo) %>%
  mutate(
    incremento_OD = OD / lag(OD),
    intervallo = paste0(lag(tempo), "/", tempo)
  ) |>
  filter(!is.na(incremento_OD)) |> # Rimuoviamo il tempo 0 (che non ha incremento)
  ungroup()

# Visualizzazione della struttura del nuovo dataset (6 righe per biomassa)
# str(dati_incrementi)

# 2. Grafico degli incrementi con palette Plasma
ggplot(dati_incrementi_rapporti, aes(x = tempo, y = incremento_OD, color = consumo)) +
  # Linee sottili per le singole biomasse
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  
  # Linea spessa della media per livello di consumo
  stat_summary(fun = mean, geom = "line", size = 1.2, na.rm = TRUE) +
  
  # Palette Plasma (da viridis)
  scale_color_viridis_d(option = "plasma") +
  
  labs(
    title = "Velocità di crescita: Incrementi giornalieri di OD",
    subtitle = "L'incremento indica quanto la biomassa aumenta tra un giorno e il successivo",
    x = "Giorno (t)",
    y = "Delta OD (OD[t] - OD[t-1])",
    color = "Consumo (I*D)"
  ) +
  theme_minimal()

# aggiungo alcune linee di commenti per testare il successo del commit + push

dati_m_gompertz <- na.omit(dati)
dati_m_gompertz <- dati_m_gompertz |> mutate(
  I = I - 280,
  D = D - 18,
  P = P - 25.5
)

modello_gompertz_fisso <- gnls(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  params = list(
    Asym ~ (I + D + P)^2,  # Interazione per l'asintoto
    b2 ~ (I + D + P)^2,    # Interazione per la latenza (richiesta)
    b3 ~ I+D+P     # Interazione per la rapidità
  ),
  # Dobbiamo aggiornare gli start: ora abbiamo 8 parametri per OGNI componente
  # (Intercetta, I, D, P, I:D, I:P, D:P, I:D:P) x 3 parametri della Gompertz = 24 coefficienti
  start = c(
    Asym = c(3, rep(0, 6)),  
    b2   = c(1.5, rep(0, 6)),               
    b3   = c(0.8, rep(0, 3))
  )
)

set.seed(123) # Per riproducibilità
ids_da_mostrare <- 1:53

dati_m_gompertz$predicted_fisso <- predict(modello_gompertz_fisso)

dati_plot <- dati_m_gompertz %>%
  filter(id_biomassa %in% ids_da_mostrare)



library(ggforce)
ggplot(dati_plot, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = OD), alpha = 0.6) +
  geom_line(aes(y = OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted_fisso, group = id_biomassa), linewidth = 1.2)+
  # SOSTITUIAMO facet_wrap con facet_wrap_paginate
  facet_wrap_paginate(~ id_biomassa, scales = "free_y", ncol = 4, nrow = 3, page = 2) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Confronto Dati Reali vs Stima Gompertz",
    subtitle = "I punti indicano le osservazioni reali; le linee continue rappresentano il modello NLME",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)"
  ) +
  theme_minimal()

modello_gompertz <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+P,
    b2 ~ (I+D)^2+P,           
    b3 ~ (I+D)^2+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  method = 'REML',
  # Nota: i valori di start devono ora includere gli zeri per i nuovi termini di interazione
  # In un modello I*D*P ci sono 8 coefficienti (intercetta + 3 principali + 3 doppie + 1 tripla)
  start = c(
    Asym = c(5, rep(0, 4)), 
    b2 = c(2, rep(0, 4)),     
    b3 = c(0.8, rep(0, 4)) )
)

modello_gompertz_ml <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+P,
    b2 ~ (I+D)^2+P,           
    b3 ~ (I+D)^2+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  method = 'ML',
  # Nota: i valori di start devono ora includere gli zeri per i nuovi termini di interazione
  # In un modello I*D*P ci sono 8 coefficienti (intercetta + 3 principali + 3 doppie + 1 tripla)
  start = c(
    Asym = c(5, rep(0, 4)), 
    b2 = c(2, rep(0, 4)),     
    b3 = c(0.8, rep(0, 4)) )
)
?nlme::nlme
modello_gompertz %>% AIC
modello_gompertz_ml %>% AIC
summary(modello_gompertz)
summary(modello_gompertz_ml)

modello_gompertz_ml

# modello_gompertz_fisso %>% AIC
# dati_m_gompertz$predicted_fisso <- predict(modello_gompertz)
dati_m_gompertz$predicted <- predict(modello_gompertz)

res_norm <- residuals(modello_gompertz, type = "normalized")
res_raw  <- residuals(modello_gompertz, type = "response")
par(mfrow = c(1, 2))
# Q-Q plot
qqnorm(res_norm, main = "Q-Q Plot residui normalizzati")
qqline(res_norm, col = "red", lwd = 2)

# Istogramma con curva normale
hist(res_norm, freq = FALSE, breaks = 20,
     main = "Distribuzione residui", xlab = "Residui")
curve(dnorm(x, mean = mean(res_norm), sd = sd(res_norm)),
      add = TRUE, col = "red", lwd = 2)

# Shapiro-Wilk (ottimo per n < 5000)
shapiro.test(res_norm)

# Kolmogorov-Smirnov
ks.test(res_norm, "pnorm", mean = mean(res_norm), sd = sd(res_norm))

# Jarque-Bera (sensibile a skewness e curtosi)
# install.packages("tseries")
library(tseries)
jarque.bera.test(res_norm)

modello_gompertz_log <- nlme(
  model = log(OD) ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+P,
    b2 ~ (I+D)^2+P,           
    b3 ~ (I+D)^2+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  # Nota: i valori di start devono ora includere gli zeri per i nuovi termini di interazione
  # In un modello I*D*P ci sono 8 coefficienti (intercetta + 3 principali + 3 doppie + 1 tripla)
  start = c(
    Asym = c(5, rep(0, 4)), 
    b2 = c(2, rep(0, 4)),     
    b3 = c(0.8, rep(0, 4)) )
)

modello_gompertz_log %>% AIC
# modello_gompertz_fisso %>% AIC
# dati_m_gompertz$predicted_fisso <- predict(modello_gompertz)
dati_m_gompertz$predicted_log <- predict(modello_gompertz_log)

res_norm_log <- residuals(modello_gompertz_log, type = "normalized")
res_raw_log  <- residuals(modello_gompertz_log, type = "response")
par(mfrow = c(1, 2))
# Q-Q plot
qqnorm(res_norm_log,main = "Q-Q Plot residui normalizzati")
qqline(res_norm_log, col = "red", lwd = 2)

# Istogramma con curva normale
hist(res_norm_log, freq = FALSE, breaks = 20,
     main = "Distribuzione residui", xlab = "Residui")
curve(dnorm(x, mean = mean(res_norm), sd = sd(res_norm)),
      add = TRUE, col = "red", lwd = 2)

# Shapiro-Wilk (ottimo per n < 5000)
shapiro.test(res_norm_log)

# Kolmogorov-Smirnov
ks.test(res_norm_log, "pnorm", mean = mean(res_norm_log), sd = sd(res_norm_log))

# Jarque-Bera (sensibile a skewness e curtosi)
# install.packages("tseries")
library(tseries)
jarque.bera.test(res_norm_log)

modello_gompertz_sqrt <- nlme(
  model = sqrt(OD) ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+P,
    b2 ~ (I+D)^2+P,           
    b3 ~ (I+D)^2+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  # Nota: i valori di start devono ora includere gli zeri per i nuovi termini di interazione
  # In un modello I*D*P ci sono 8 coefficienti (intercetta + 3 principali + 3 doppie + 1 tripla)
  start = c(
    Asym = c(5, rep(0, 4)), 
    b2 = c(2, rep(0, 4)),     
    b3 = c(0.8, rep(0, 4)) )
)

modello_gompertz_sqrt %>% AIC
# modello_gompertz_fisso %>% AIC
# dati_m_gompertz$predicted_fisso <- predict(modello_gompertz)
dati_m_gompertz$predicted_sqrt <- predict(modello_gompertz_sqrt)

res_norm_sqrt <- residuals(modello_gompertz_sqrt, type = "normalized")
res_raw_sqrt  <- residuals(modello_gompertz_sqrt, type = "response")
par(mfrow = c(1, 2))
# Q-Q plot
qqnorm(res_norm_sqrt,main = "Q-Q Plot residui normalizzati")
qqline(res_norm_sqrt, col = "red", lwd = 2)

# Istogramma con curva normale
hist(res_norm_sqrt, freq = FALSE, breaks = 20,
     main = "Distribuzione residui", xlab = "Residui")
curve(dnorm(x, mean = mean(res_norm_sqrt), sd = sd(res_norm_sqrt)),
      add = TRUE, col = "red", lwd = 2)

# Shapiro-Wilk (ottimo per n < 5000)
shapiro.test(res_norm_sqrt)

# Kolmogorov-Smirnov
ks.test(res_norm_sqrt, "pnorm", mean = mean(res_norm_sqrt), sd = sd(res_norm_sqrt))

# Jarque-Bera (sensibile a skewness e curtosi)
# install.packages("tseries")
library(tseries)
jarque.bera.test(res_norm_sqrt)


set.seed(123) # Per riproducibilità
ids_da_mostrare <- 1:53

dati_plot <- dati_m_gompertz %>%
  filter(id_biomassa %in% ids_da_mostrare)

library(ggforce)
x11()
ggplot(dati_plot, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = OD), alpha = 0.6) +
  geom_line(aes(y = OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted, group = id_biomassa), linewidth = 1.2)+
  facet_wrap_paginate(~ id_biomassa, scales = "free_y", ncol = 5, nrow = 5, page = 1) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Confronto Dati Reali vs Stima Gompertz",
    subtitle = "I punti indicano le osservazioni reali; le linee continue rappresentano il modello NLME",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)"
  ) +
  theme_minimal()

summary(modello_gompertz)
# Controlla i residui per unità
plot(modello_gompertz)

# Vedi quali unità hanno residui grandi
residuals(modello_gompertz, type = "normalized") |> 
  abs() |> 
  sort(decreasing = TRUE) |> 
  head(10)

# Controlla i BLUP stimati
random.effects(modello_gompertz)
# Visualizza le curve delle unità problematiche
dati_m_gompertz |> 
  filter(id_biomassa %in% c(45, 36, 37, 52)) |> 
  ggplot(aes(x = tempo, y = OD, color = as.factor(id_biomassa))) +
  geom_point() +
  geom_line(aes(y = predict(modello_gompertz)[
    dati_m_gompertz$id_biomassa %in% c(45, 36, 37, 52)
  ])) +
  facet_wrap(~ id_biomassa, scales = "free_y")

# previsione delle nuove curve

# Valori osservati di I, D, P nel disegno sperimentale
library(marginaleffects)
pred_curve <- predictions(
  modello_gompertz,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  group_split()

# previsione delle nuove curve, vedendo se all'aggiunta di nuove settimane si raggiunge plateau
pred_curve_10s <- predictions(
  modello_gompertz,
  newdata = datagrid(
    tempo = 0:9,  # i 7 punti osservati
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_10s |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  group_split()

pred_curve_tutto_sc <- predictions(
  modello_gompertz,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = seq(-190,190,by = 5),
    D = seq(-6,6,by = .5),
    P = seq(-24.5,24.5,by = 1),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_tutto_sc |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  slice_tail(n = 1) |> 
  ungroup() |>
  arrange(desc(estimate))



# guardiamo comparisons
library(marginaleffects)

# Effetto di I tenendo D e P ai loro valori unici
comp_I <- comparisons(
  modello_gompertz,
  newdata = datagrid(
    tempo = Inf,
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  variables = "I",
  re.form = NA
)

comp_I |> as_tibble() |> select(I, D, P, contrast, estimate, conf.low, conf.high, p.value) |> print(n = 27)

library(marginaleffects)

# Effetto di D tenendo I e P ai loro valori unici
comp_D <- comparisons(
  modello_gompertz,
  newdata = datagrid(
    tempo = Inf,
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  variables = "D",
  re.form = NA
)

comp_D |> as_tibble() |> select(I, D, P, contrast, estimate, conf.low, conf.high, p.value)|> print(n = 27)

library(marginaleffects)

# Effetto di P tenendo D e I ai loro valori unici
comp_P <- comparisons(
  modello_gompertz,
  newdata = datagrid(
    tempo = Inf,
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  variables = "P",
  re.form = NA
)

comp_P |> as_tibble() |> select(I, D, P, contrast, estimate, conf.low, conf.high, p.value) |> print(n = 27)

# avg comparisons
avg_comparisons(
  modello_gompertz,
  newdata = datagrid(
    tempo = Inf,
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  variables = c("I", "D", "P"),
  re.form = NA
)

View(dati)

# previsione su griglia di valori out of sample
pred_curve_nuove <- predictions(
  modello_gompertz,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = seq(570,670,by = 10) - 280,
    D = seq(18,27,by=1) - 18,
    P = seq(60,75,by=3) - 25.5,
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_nuove |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  slice_tail(n = 1) |> 
  ungroup() |>
  arrange(desc(estimate)) |> 
  print(n = 50)

pred_curve |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  filter(is.finite(tempo)) |>  # esclude tempo = Inf
  ggplot(aes(x = tempo, y = estimate, color = as.factor(I))) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(I)),
              alpha = 0.15, color = NA) +
  facet_grid(D ~ P, labeller = label_both) +
  labs(
    x = "Tempo",
    y = "OD predetto",
    color = "I",
    fill = "I",
    title = "Curve predette per condizione sperimentale"
  ) +
  theme_minimal()

pred_marginale_I <- pred_curve_10s |>
  as_tibble() |>
  filter(is.finite(tempo)) |>
  mutate(I = I + 280) |>  # o qualunque sia il valore di centratura
  group_by(I, tempo) |>
  summarise(estimate = mean(estimate), .groups = "drop")

ggplot(dati, aes(x = tempo, y = OD, color = as.factor(I), group = id_biomassa)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.3) +
  geom_line(data = pred_marginale_I,
            aes(x = tempo, y = estimate, color = as.factor(I), group = as.factor(I)),
            linewidth = 1.5, inherit.aes = FALSE) +
  labs(
    title = "Traiettorie di crescita per il fattore I",
    subtitle = "Linee spesse = curva predetta marginale; linee sottili = dati osservati",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)",
    color = "Fattore I"
  ) +
  theme_minimal()

unique(dati$I)
unique(pred_marginale_I$I)

# modello gompertz con effetto quadratico dei predittori sui parametri della gompertz
modello_gompertz_quad <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+I(I^2) + I(D^2) + P + I(P^2),
    b2 ~ (I+D)^2+I(I^2) + I(D^2) + P + I(P^2),           
    b3 ~ (I+D)^2+I(I^2) + I(D^2) + P + I(P^2)     
  ),
  random = Asym ~ 1 | id_biomassa,
  # Nota: i valori di start devono ora includere gli zeri per i nuovi termini di interazione
  # In un modello I*D*P ci sono 8 coefficienti (intercetta + 3 principali + 3 doppie + 1 tripla)
  start = c(
    Asym = c(5, rep(0, 7)), 
    b2 = c(2, rep(0, 7)),     
    b3 = c(0.8, rep(0, 7)) )
)
AIC(modello_gompertz,modello_gompertz_quad)

pred_curve_quad <- predictions(
  modello_gompertz_quad,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = unique(dati_m_gompertz$I),
    D = unique(dati_m_gompertz$D),
    P = unique(dati_m_gompertz$P),
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_quad |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  group_split()

pred_curve_nuove_quad <- predictions(
  modello_gompertz_quad,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = seq(570,670,by = 10) - 280,
    D = seq(18,27,by=1) - 18,
    P = seq(60,75,by=3) - 25.5,
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_nuove_quad |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  slice_tail(n = 1) |> 
  ungroup() |>
  arrange(desc(estimate)) |> 
  print(n = 20)

predicted_quad <- predict(modello_gompertz_quad)

ggplot(dati_m_gompertz, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = OD), alpha = 0.6) +
  geom_line(aes(y = OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted_quad, group = id_biomassa), linewidth = 1.2)+
  facet_wrap_paginate(~ id_biomassa, scales = "free_y", ncol = 5, nrow = 5, page = 3) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Confronto Dati Reali vs Stima Gompertz con termini quadratici",
    subtitle = "I punti indicano le osservazioni reali; le linee continue rappresentano il modello NLME",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)"
  ) +
  theme_minimal()

# modelli standardizzati
dati_m_gompertz_st <- na.omit(dati)
dati_m_gompertz_st <- dati_m_gompertz_st |> mutate(
  I = (I - 280) / 190,
  D = (D - 18) / 6,
  P = (P - 25.5)/24.5
)

modello_gompertz_st <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz_st,
  fixed = list(
    Asym ~ (I+D)^2+P,
    b2 ~ (I+D)^2+P,           
    b3 ~ (I+D)^2+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  # Nota: i valori di start devono ora includere gli zeri per i nuovi termini di interazione
  # In un modello I*D*P ci sono 8 coefficienti (intercetta + 3 principali + 3 doppie + 1 tripla)
  start = c(
    Asym = c(5, rep(0, 4)), 
    b2 = c(2, rep(0, 4)),     
    b3 = c(0.8, rep(0, 4)) )
)

modello_gompertz_st %>% AIC
# modello_gompertz_fisso %>% AIC
# dati_m_gompertz$predicted_fisso <- predict(modello_gompertz)
dati_m_gompertz_st$predicted <- predict(modello_gompertz_st)

library(marginaleffects)
pred_curve_st <- predictions(
  modello_gompertz_st,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = unique(dati_m_gompertz_st$I),
    D = unique(dati_m_gompertz_st$D),
    P = unique(dati_m_gompertz_st$P),
    id_biomassa = dati_m_gompertz_st$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_st |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  group_split()

pred_curve_nuove_st <- predictions(
  modello_gompertz_st,
  newdata = datagrid(
    tempo = 0:6,  # i 7 punti osservati
    I = (seq(570,670,by = 10) - 280)/190,
    D = (seq(18,27,by=1) - 18)/6,
    P = (seq(60,75,by=3) - 25.5)/24.5,
    id_biomassa = dati_m_gompertz$id_biomassa[1]
  ),
  re.form = NA
)

pred_curve_nuove_st |> 
  as_tibble() |>
  select(I, D, P, tempo, estimate, conf.low, conf.high) |>
  group_by(I, D, P) |>
  slice_tail(n = 1) |> 
  ungroup() |>
  arrange(desc(estimate)) |> 
  print(n = 5)

# approccio bayesiano
# install.packages('brms')
# library(brms)

# prior <- c(
#   # Intercette
#   prior(normal(14, 5),    class = "b", coef = "Intercept", nlpar = "Asym"),
#   prior(normal(5, 2),     class = "b", coef = "Intercept", nlpar = "b2"),
#   prior(normal(0.7, 0.3), class = "b", coef = "Intercept", nlpar = "b3"),
#   
#   # Coefficienti delle covariate (tutti gli altri b)
#   prior(normal(0, 2), class = "b", nlpar = "Asym"),
#   prior(normal(0, 2), class = "b", nlpar = "b2"),
#   prior(normal(0, 2), class = "b", nlpar = "b3"),
#   
#   # Varianza random
#   prior(cauchy(0, 2), class = "sd", nlpar = "Asym"),
#   
#   # Nu della Student-t (già ha un default ma puoi restringerlo)
#   prior(gamma(2, 0.1), class = "nu")
# )
# brm(bf(OD ~ Asym * exp(-exp(-(tempo - b2) / b3)),
#          Asym ~ (I + D)^2 + P + (1 | id_biomassa),
#          b2 ~ (I + D)^2 + P,
#          b3 ~ (I + D)^2 + P,
#          nl = TRUE,family = student()), 
#     data = dati_m_gompertz,prior = prior, chains = 0)  # chains = 0 controlla senza campionare
# 
# fit_brms <- brm(
#   bf(OD ~ Asym * exp(-exp(-(tempo - b2) / b3)),
#      Asym ~ (I + D)^2 + P + (1 | id_biomassa),
#      b2 ~ (I + D)^2 + P,
#      b3 ~ (I + D)^2 + P,
#      nl = TRUE,family = student()),
#   data = dati_m_gompertz,
#   prior = prior,
#   chains = 4,
#   iter = 2000,
#   cores = 4
# )
# 
# default_prior(
#   bf(OD ~ Asym * exp(-exp(-(tempo - b2) / b3)),
#      Asym ~ (I + D)^2 + P + (1 | id_biomassa),
#      b2 ~ (I + D)^2 + P,
#      b3 ~ (I + D)^2 + P,
#      nl = TRUE),
#   family = student(),
#   data = dati_m_gompertz
# )

#### Baranyi - Roberts pipeline

# bootstrap
coef(modello_gompertz_ml)
summary(modello_gompertz_ml)

dati_m_gompertz <- na.omit(dati)
dati_m_gompertz <- dati_m_gompertz |> mutate(
  I = I - 280,
  D = D - 18,
  P = P - 25.5
)

Sys.time()
modello_gompertz_ml <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+P, # cambio
    b2 ~ (I+D)^2+P, # cambio          
    b3 ~ (I+D)^2+P  # cambio   
  ),
  random = Asym ~ 1 | id_biomassa,
  method = 'ML',
  start = c(
    Asym = c(5, rep(0, 4)), 
    b2 = c(2, rep(0, 4)),     
    b3 = c(0.8, rep(0, 4)) )
)
Sys.time()
fixef(modello_gompertz_ml)
start_vals <- c(
  Asym = fixef(modello_gompertz_ml)[1:5],
  b2   = fixef(modello_gompertz_ml)[6:10],
  b3   = fixef(modello_gompertz_ml)[11:15]
)

library(boot)
?boot

bm_ids <- unique(dati_m_gompertz$id_biomassa)
df_idbiomass <- data.frame(id_biomassa = bm_ids) 

fun.stat <- function(data,i){
  
  sampled_biomass <- data[i]
  
  # Genera il dataset di bootstrap combinando i dati campionati in un unico data frame
  dati_boot <- map_dfr(seq_along(sampled_biomass), ~ {
    # rappresenta l'indice corrente del ciclo (es. 1, 2, 3...)
    dati_m_gompertz %>%
      #  Isola i dati della biomassa estratta in questa specifica iterazione
      filter(id_biomassa == sampled_biomass[.x]) %>%
      # Riassegna un ID univoco basato sulla posizione del clone ("cl_1", "cl_2", ecc.).
      # Questo è fondamentale se lo stesso ID originale è stato estratto più volte.
      mutate(id_biomassa = paste0("cl_", .x))
  })
  
  mod <- tryCatch(
    nlme(
      model = OD ~ SSgompertz(tempo, Asym, b2, b3),
      data = dati_boot,
      fixed = list(
        Asym ~ (I+D)^2+P, # cambio
        b2 ~ (I+D)^2+P, # cambio          
        b3 ~ (I+D)^2+P  # cambio   
      ),
      random = Asym ~ 1 | id_biomassa,
      method = 'ML',
      start = start_vals,
      control = nlmeControl(
        maxIter = 500,
        msMaxIter = 500,
        tolerance = 1e-3
      )
    ),
    error = function(e) NULL
  )
  if (is.null(mod)) return(rep(NA, 15))
  return(fixef(mod))
}

Sys.time()
res.boot <- boot(bm_ids,statistic = fun.stat,R = 100)
Sys.time()
cat("Fallimenti:", sum(is.na(res.boot$t[, 1])), "/", 100, "\n")
res.boot

nomi <- names(fixef(modello_gompertz_ml))

ic_list <- lapply(1:15, function(k) {
  ci <- boot.ci(res.boot, index = k, type = c("perc", "bca"))
  data.frame(
    parametro = nomi[k],
    lower_perc = ci$percent[4],
    upper_perc = ci$percent[5],
    lower_bca  = ci$bca[4],
    upper_bca  = ci$bca[5]
  )
})

ic_df <- do.call(rbind, ic_list)
print(ic_df)

# bootstrap sui residui
fitted_vals <- fitted(modello_gompertz_ml, level = 0)  # solo effetti fissi
residui <- residuals(modello_gompertz_ml, level = 0)   # residui rispetto agli effetti fissi

fun.stat.resid <- function(data, i) {
  
  sampled_biomass <- data[i]  # data è bm_ids, un vettore
  
  dati_boot <- map_dfr(seq_along(sampled_biomass), ~ {
    idx <- which(dati_m_gompertz$id_biomassa == sampled_biomass[.x])  # riferimento esplicito
    sub <- dati_m_gompertz[idx, ]
    sub$OD <- fitted_vals[idx] + residui[idx]
    sub$id_biomassa <- paste0("cl_", .x)
    sub
  })
  
  mod <- tryCatch(
    nlme(
      model = OD ~ SSgompertz(tempo, Asym, b2, b3),
      data = dati_boot,
      fixed = list(
        Asym ~ (I+D)^2+P,
        b2 ~ (I+D)^2+P,           
        b3 ~ (I+D)^2+P     
      ),
      random = Asym ~ 1 | id_biomassa,
      method = 'ML',
      start = start_vals,
      control = nlmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-3)
    ),
    error = function(e) NULL
  )
  
  if (is.null(mod)) return(rep(NA, 15))
  return(fixef(mod))
}

res.boot.resid <- boot(bm_ids, statistic = fun.stat.resid, R = 100)

cat("Fallimenti:", sum(is.na(res.boot.resid$t[, 1])), "/", 100, "\n")

ic_list2 <- lapply(1:15, function(k) {
  ci <- boot.ci(res.boot.resid, index = k, type = c("perc", "bca"))
  data.frame(
    parametro = nomi[k],
    lower_perc = ci$percent[4],
    upper_perc = ci$percent[5],
    lower_bca  = ci$bca[4],
    upper_bca  = ci$bca[5]
  )
})

ic_df2 <- do.call(rbind, ic_list)
print(ic_df2)


################################################################################
#  PARTE 3 - INCREMENTI, GAM E SUPERFICI DI RISPOSTA
#  (dal file Riccardo.R)
################################################################################

# Grafico incrementi ------------------------------------------------------

# Differenze per consumo
ggplot(dati_incrementi, aes(x = tempo, y = incremento_OD, color = consumo)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", size = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Velocità di crescita: Incrementi giornalieri di OD",
    subtitle = "L'incremento indica quanto la biomassa aumenta tra un giorno e il successivo",
    x = "Giorno (t)",
    y = "Delta OD (OD[t] - OD[t-1])",
    color = "Consumo (I*D)"
  ) +
  theme_minimal()

# Differenze per combinazione di condizioni sperimentali
ggplot(dati_incrementi, aes(x = tempo, y = incremento_OD, color = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", size = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Velocità di crescita: Incrementi giornalieri di OD",
    subtitle = "L'incremento indica quanto la biomassa aumenta tra un giorno e il successivo",
    x = "Settimana (t)",
    y = "Delta OD (OD[t] - OD[t-1])",
    color = "Condizioni Sperimentali"
  ) +
  theme_minimal()

# Rapporti per comb. condiz. sperimentali
ggplot(dati_incrementi_rapporti, aes(x = tempo, y = incremento_OD, color = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", size = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Velocità di crescita: Incrementi giornalieri di OD",
    subtitle = "L'incremento indica quanto la biomassa aumenta tra un giorno e il successivo",
    x = "Settimana (t)",
    y = "Delta OD (OD[t] - OD[t-1])",
    color = "Condizioni Sperimentali"
  ) +
  theme_minimal()


# Caricamento librerie necessarie
library(plotly)
library(dplyr)

# In questo caso fisso il tempo (per esempio a 7, poi andrebbe trovato il t ottimale in cui fermarsi)
# e fisso anche il fosforo, quindi faccio variare solo I e D

# Creiamo una griglia densa per I, D e fissiamo P (o creiamo più slice)
i_seq <- seq(90, 470, length.out = 50)
d_seq <- seq(12, 24, length.out = 50)
p_valore_fisso <- 25.5  # Esempio: visualizziamo la sezione centrale di P

# 1. Creiamo la griglia con i valori CENTRATI (come il modello si aspetta)
grid <- expand.grid(
  I = seq(90 - 280, 470 - 280, length.out = 50), # Range: -190 a +190
  D = seq(12 - 18, 24 - 18, length.out = 50),    # Range: -6 a +6
  P = 0,                               # Fissato a 0 (che è il vecchio 25.5)
  tempo = max(dati_m_gompertz$tempo)
)

# 2. Predizione
grid$Risposta_Predetta <- predict(modello_gompertz, newdata = grid, level = 0)

# 3. Trasformazione per il grafico (Matrice Z)
z_matrix <- matrix(grid$Risposta_Predetta, nrow = 50, ncol = 50)

# 4. GRAFICO: Usiamo i valori REALI per gli assi X e Y per leggibilità
plot_ly(x = ~seq(90, 470, length.out = 50), # Valori originali per l'asse X
        y = ~seq(12, 24, length.out = 50),  # Valori originali per l'asse Y
        z = ~z_matrix) %>%
  add_surface() %>%
  layout(
    scene = list(
      xaxis = list(title = "Intensità Luminosa (I) REALE"),
      yaxis = list(title = "Durata (D) REALE"),
      zaxis = list(title = "Risposta (OD)")
    ),
    title = paste0("Superficie di Risposta P = ", p_valore_fisso)
  )


# -------------------------------------------------------------------------
# Al variare di P

# Parametri fissi
tempo_target <- max(dati_m_gompertz$tempo)
p_levels <- c(1, 25.5, 50)
i_seq_real <- seq(90, 470, length.out = 40) # Diminuito leggermente per fluidità
d_seq_real <- seq(12, 24, length.out = 40)

# 1. Funzione per generare i dati di ogni superficie
get_surface_data <- function(val_p) {
  grid_p <- expand.grid(
    I = i_seq_real - 280,
    D = d_seq_real - 18,
    P = val_p - 25.5,
    tempo = tempo_target
  )
  grid_p$Risposta_Predetta <- predict(modello_gompertz, newdata = grid_p, level = 0)
  return(matrix(grid_p$Risposta_Predetta, nrow = 40, ncol = 40))
}

# 2. Creazione del grafico con scene multiple
plot_ly() %>%
  # Superficie 1 (P = 1)
  add_surface(z = ~get_surface_data(1), x = ~i_seq_real, y = ~d_seq_real, 
              scene = "scene1", name = "P=1", colorscale = "Viridis", showscale = FALSE) %>%
  # Superficie 2 (P = 25.5)
  add_surface(z = ~get_surface_data(25.5), x = ~i_seq_real, y = ~d_seq_real, 
              scene = "scene2", name = "P=25.5", colorscale = "Viridis", showscale = FALSE) %>%
  # Superficie 3 (P = 50)
  add_surface(z = ~get_surface_data(50), x = ~i_seq_real, y = ~d_seq_real, 
              scene = "scene3", name = "P=50", colorscale = "Viridis", showscale = TRUE) %>%
  layout(
    title = "Confronto Superfici di Risposta al variare di P",
    # Configurazione spaziale delle tre scene
    grid = list(rows = 1, columns = 3, pattern = 'independent'),
    scene = list(domain = list(column = 0), xaxis = list(title = "I"), yaxis = list(title = "D"), zaxis = list(title = "OD"), title = "P=1"),
    scene2 = list(domain = list(column = 1), xaxis = list(title = "I"), yaxis = list(title = "D"), zaxis = list(title = "OD"), title = "P=25.5"),
    scene3 = list(domain = list(column = 2), xaxis = list(title = "I"), yaxis = list(title = "D"), zaxis = list(title = "OD"), title = "P=50")
  )


# GAM :incrementi ---------------------------------------------------------

library(mgcv)

dati_gam <- dati_incrementi %>%
  mutate(
    I_c = I - 280,
    D_c = D - 18,
    P_c = P - 25.5,
    id_biomassa = as.factor(id_biomassa) # Necessario per l'effetto random
  )

modello_gam <- gam(
  incremento_OD ~ s(tempo, k = 4) +           # Andamento temporale comune
    s(I_c, D_c, k = 4) +                # Superficie di risposta (interazione I e D)
    sP_c +                        # Effetto lineare di P
    s(id_biomassa, bs = "re"),   # Effetto random sulle biomasse (id_biomassa)
  data = dati_gam,
  method = "REML"                              # Metodo più robusto per la stima
)

summary(modello_gam)

plot(modello_gam, terms = "s(I_c, D_c, k = 4)", se = T)

# Creiamo la griglia
grid_gam <- expand.grid(
  tempo = 6,                   # Scegliamo il giorno del picco osservato nei grafici
  I_c = seq(-190, 190, length.out = 50),
  D_c = seq(-6, 6, length.out = 50),
  P_c = 0,                     # Punto centrale (P=25.5)
  id_biomassa = dati_gam$id_biomassa[1] # Un ID arbitrario (verrà ignorato con exclude)
)

# Predizione escludendo l'effetto random (per avere la media)
grid_gam$Pred_Incremento <- predict(modello_gam, newdata = grid_gam, exclude = "s(id_biomassa)")

# Trovare il massimo
ottimo_gam <- grid_gam[which.max(grid_gam$Pred_Incremento), ]

# Creazione di un nuovo dataframe con i valori reali
ottimo_reale <- ottimo_gam %>%
  mutate(
    I_reale = I_c + 280,
    D_reale = D_c + 18,
    P_reale = P_c + 25.5,
    Tempo_reale = tempo
  ) %>%
  select(Tempo_reale, I_reale, D_reale, P_reale, Pred_Incremento)

print(ottimo_reale)

par(mfrow=c(2,2))
plot(modello_gam, all.terms = T)

# -------------------------------------------------------------------------
gam_biv <- gam(incremento_OD~lo(I, D, P, span = 0.1) + s(tempo, k = 4), data = dati_incrementi)

I_grid <- seq(min(dati_incrementi$I), max(dati_incrementi$I), l = 50)
D_grid <- seq(min(dati_incrementi$D), max(dati_incrementi$D), l = 50)
Xg <- expand.grid(I = I_grid, D = D_grid, tempo = 6, P = 0)
y_hat <- predict(gam_biv, Xg)
z_matrix <- matrix(y_hat, nrow = length(I_grid), ncol = length(D_grid))

persp(I_grid, D_grid, z_matrix)
# -------------------------------------------------------------------------

modello_gam2 <- gam(
  incremento_OD ~ s(tempo, k = 4) +           # Andamento temporale comune
    s(I_c, D_c, k = 4) +                # Superficie di risposta (interazione I e D)
    s(P_c, k = 3) +                        # Effetto lineare di P
    s(id_biomassa, bs = "re"),   # Effetto random sulle biomasse (id_biomassa)
  data = dati_gam,
  method = "REML"                              # Metodo più robusto per la stima
)

summary(modello_gam)

par(mfrow=c(2,2))
plot(modello_gam2, se = T)

# Creiamo la griglia
grid_gam <- expand.grid(
  tempo = 6,                   # Scegliamo il giorno del picco osservato nei grafici
  I_c = seq(-190, 190, length.out = 50),
  D_c = seq(-6, 6, length.out = 50),
  P_c = seq(-25.5, 25.5, length.out = 50),                     
  id_biomassa = dati_gam$id_biomassa[1] # Un ID arbitrario (verrà ignorato con exclude)
)

# Predizione escludendo l'effetto random (per avere la media)
grid_gam$Pred_Incremento <- predict(modello_gam2, newdata = grid_gam, exclude = "s(id_biomassa)")

# Trovare il massimo
ottimo_gam <- grid_gam[which.max(grid_gam$Pred_Incremento), ]

# Creazione di un nuovo dataframe con i valori reali
ottimo_reale <- ottimo_gam %>%
  mutate(
    I_reale = I_c + 280,
    D_reale = D_c + 18,
    P_reale = P_c + 25.5,
    Tempo_reale = tempo
  ) %>%
  select(Tempo_reale, I_reale, D_reale, P_reale, Pred_Incremento)

print(ottimo_reale)