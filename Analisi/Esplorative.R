rm(list=ls());gc();

library(ggplot2)
library(gridExtra)

load("Analisi/dati_modificati.Rdata")

# Palette unica e risoluzione per tutti i grafici
PALETTE <- "plasma"
DPI     <- 300

# Cartella di destinazione dei grafici (percorso assoluto)
#OUT_DIR <- "C:/Users/Antonio/Desktop/UNI/MAGISTRALE/SECONDO ANNO/2_SEMESTRE/ITERAZIONE/progetto/BiomassProduction/Grafici/"
 OUT_DIR <- "C:/Users/Utente/OneDrive/Universita/Magistrale/2025-2026/Iterazione/Progetto/Grafici"
dati$tempo <- dati$tempo + 1

# Tutte le curve --------------------------------------------------------------

p_tutte <- ggplot(dati, aes(x = tempo, y = OD, group = id_biomassa)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE, alpha = 1, size = 1) +
  labs(
    x = "Tempo (settimane)",
    y = "Resa di Biomassa (OD)"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "tutte_curve.png"), plot = p_tutte,
       width = 6, height = 5, dpi = DPI)


# Curve per condizione sperimentale (singole) ---------------------------------

p1 <- ggplot(dati, aes(x = tempo, y = OD,
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
  theme(legend.position = "none")

ggsave(
  file.path(OUT_DIR, "tutte_curve_condizione.png"),
  plot  = p1,
  width = 8, height = 5, dpi = DPI
)

# Curve medie per condizione sperimentale (sfondo grigio + media colorata) ----

p2 <- ggplot(dati, aes(x = tempo, y = OD,
                       color = as.factor(condizione_sperimentale),
                       fill  = as.factor(condizione_sperimentale))) +
  geom_line(aes(group = id_biomassa), color = "gray85", alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(group = id_biomassa), color = "gray85", alpha = 0.5, size = 1, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line",  linewidth = 0.7, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 1,        na.rm = TRUE) +
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

ggsave(
  file.path(OUT_DIR, "Curve_medie_condizione.png"),
  plot  = p2,
  width = 7, height = 5, dpi = DPI
)


# Grafici condizionati --------------------------------------------------------

grafico_condizionate <- function(condizionante, f1, f2 = NULL) {
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
    # .data[[]] è la forma robusta per indicizzare colonne dinamiche in vars()
    facet_wrap(vars(.data[[condizionante]])) +
    scale_color_viridis_d(option = PALETTE, end = 0.9) +
    scale_fill_viridis_d(option = PALETTE,  end = 0.9) +
    labs(
      x     = "Tempo (settimane)",
      y     = "OD",
      color = label_colore,
      fill  = label_colore
    ) +
    theme_minimal()
}

p1 <- grafico_condizionate("I", "D")
p2 <- grafico_condizionate("P", "D")
p3 <- grafico_condizionate("I", "P")

ggsave(file.path(OUT_DIR, "condizionata_D_I.png"),
       plot = p1,
       width = 9, height = 4, dpi = DPI)

ggsave(file.path(OUT_DIR, "condizionata_D_P.png"),
       plot = p2,
       width = 9, height = 4, dpi = DPI)

ggsave(file.path(OUT_DIR, "condizionata_P_I.png"),
       plot = p3,
       width = 9, height = 4, dpi = DPI)

library(patchwork)
# Unisci i grafici in un'unica colonna
plot_combinato <- p1 / p2 / p3 +
  plot_annotation(
  ) &
  theme(plot.margin = margin(10, 10, 10, 10))

# Esporta impostando le dimensioni esatte della pagina (es. A4, tenendo conto dei margini)
ggsave("Grafici/condizionate_combinate.pdf",
       plot = plot_combinato,
       width = 9,
       height = 12, # Dimensioni A4, modificale in base a quanto spazio vuoi occupare
       units = "in")


# Grafici marginali -----------------------------------------------------------
# Singole curve sullo sfondo (alpha basso) + media per livello in grassetto.

# marginale_I <- ggplot(dati, aes(x = tempo, y = OD,
#                                 color = as.factor(I), fill = as.factor(I))) +
#   geom_line(aes(group = id_biomassa), alpha = 0.3, na.rm = TRUE) +
#   geom_point(aes(group = id_biomassa), alpha = 0.3, size = 0.8, na.rm = TRUE) +
#   stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
#   stat_summary(fun = mean, geom = "line",  linewidth = 1.2, na.rm = TRUE) +
#   stat_summary(fun = mean, geom = "point", size = 2,        na.rm = TRUE) +
#   scale_color_viridis_d(option = PALETTE, end = 0.9) +
#   scale_fill_viridis_d(option = PALETTE,  end = 0.9) +
#   labs(x = "Tempo (settimane)", y = "Optical Density (OD)",
#        color = "Fattore I", fill = "Fattore I") +
#   theme_minimal()
# 
# marginale_D <- ggplot(dati, aes(x = tempo, y = OD,
#                                 color = as.factor(D), fill = as.factor(D))) +
#   geom_line(aes(group = id_biomassa), alpha = 0.3, na.rm = TRUE) +
#   geom_point(aes(group = id_biomassa), alpha = 0.3, size = 0.8, na.rm = TRUE) +
#   stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
#   stat_summary(fun = mean, geom = "line",  linewidth = 1.2, na.rm = TRUE) +
#   stat_summary(fun = mean, geom = "point", size = 2,        na.rm = TRUE) +
#   scale_color_viridis_d(option = PALETTE, end = 0.9) +
#   scale_fill_viridis_d(option = PALETTE,  end = 0.9) +
#   labs(x = "Tempo (settimane)", y = "Optical Density (OD)",
#        color = "Fattore D", fill = "Fattore D") +
#   theme_minimal()
# 
# marginale_P <- ggplot(dati, aes(x = tempo, y = OD,
#                                 color = as.factor(P), fill = as.factor(P))) +
#   geom_line(aes(group = id_biomassa), alpha = 0.3, na.rm = TRUE) +
#   geom_point(aes(group = id_biomassa), alpha = 0.3, size = 0.8, na.rm = TRUE) +
#   stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
#   stat_summary(fun = mean, geom = "line",  linewidth = 1.2, na.rm = TRUE) +
#   stat_summary(fun = mean, geom = "point", size = 2,        na.rm = TRUE) +
#   scale_color_viridis_d(option = PALETTE, end = 0.9) +
#   scale_fill_viridis_d(option = PALETTE,  end = 0.9) +
#   labs(x = "Tempo (settimane)", y = "Optical Density (OD)",
#        color = "Fattore P", fill = "Fattore P") +
#   theme_minimal()


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

marginale_I <- marginal_plot("I", "I",  COL_I)
marginale_D <- marginal_plot("D", "D",  COL_D)
marginale_P <- marginal_plot("P", "P", COL_P)


ggsave(
  file.path(OUT_DIR, "Grafici/marginali.png"),
  plot  = arrangeGrob(marginale_I, marginale_D, marginale_P, ncol = 3),
  width = 18, height = 5, dpi = DPI
)


# Dati mancanti ---------------------------------------------------------------
# Le due unità con OD mancante all'ultima settimana appartengono alla
# condizione I=90, D=24, con P diverso (1 e 25.5): il pattern non è casuale
# rispetto alla condizione sperimentale (non MCAR).

na_ids <- dati$id_biomassa[is.na(dati$OD)]

p_na_zoom <- ggplot(
  subset(dati, id_biomassa %in% na_ids),
  aes(x = tempo, y = OD,
      color = as.factor(condizione_sperimentale),
      group = id_biomassa)
) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +

  scale_color_viridis_d(option = PALETTE, end = 0.9) +
  labs(
    x     = "Tempo (settimane)",
    y     = "Optical Density (OD)",
    color = "Condizione sperimentale"
  ) +
  ylim(0, max(dati$OD, na.rm = TRUE)) +
  theme_minimal()

p_na_contesto <- ggplot(
  dati,
  aes(x = tempo, y = OD,
      color = id_biomassa %in% na_ids,
      group = id_biomassa)
) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  scale_color_manual(
    values = c("FALSE" = "gray70", "TRUE" = "#D64E12"),
    labels = c("Completa", "Con NA")
  ) +
  labs(
    x     = "Tempo (settimane)",
    y     = "Optical Density (OD)",
    color = "Unità"
  ) +
  ylim(0, max(dati$OD, na.rm = TRUE)) +
  theme_minimal()

ggsave(
  file.path(OUT_DIR, "dati_mancanti.png"),
  plot  = arrangeGrob(p_na_zoom, p_na_contesto, ncol = 2),
  width = 14, height = 5, dpi = DPI
)


# Grafici incrementi ----------------------------------------------------------
# "consumo" è il prodotto I x D (energia luminosa totale ricevuta).
# La variabile deve essere presente in dati_incrementi e dati_incrementi_rapporti
# prima di eseguire questo blocco (cfr. script di preparazione dati).

diff <- ggplot(dati_incrementi,
               aes(x = tempo, y = incremento_OD, color = condizione_sperimentale, fill = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
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

rapp <- ggplot(dati_incrementi_rapporti,
               aes(x = tempo, y = incremento_OD, color = condizione_sperimentale, fill = condizione_sperimentale)) +
  geom_line(aes(group = id_biomassa), alpha = 0.15, na.rm = TRUE) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA, na.rm = TRUE) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +

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
  file.path(OUT_DIR, "incrementi.png"),
  plot  = arrangeGrob(diff, rapp, ncol = 2),
  width = 12, height = 5, dpi = DPI
)


# Relazione media-varianza ----

library(dplyr)
# 1. Calcolo di media e varianza di OD per ogni punto temporale
df_summary <- dati |>
  group_by(tempo) |>
  summarise(
    media_OD = mean(OD, na.rm = TRUE),
    varianza_OD = var(OD, na.rm = TRUE)
  )

t_var <- ggplot(df_summary, aes(x = tempo, y = varianza_OD)) +
  # Punti con leggera trasparenza per gestire l'overplotting
  geom_point(size = 2.5, color = "#2c3e50", alpha = 0.6) +
  # Linea di riferimento di Poisson (Varianza = Media) per vedere l'iperdispersione
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_smooth(method = "loess", color = "#e74c3c", linewidth = 1.2, se = FALSE, span=1, fill = "#e74c3c")+
  labs(
    x = "Tempo",
    y = "Varianza di OD"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

m_var <- ggplot(df_summary, aes(x = media_OD, y = varianza_OD)) +
  # Punti con leggera trasparenza per gestire l'overplotting
  geom_point(size = 2.5, color = "#2c3e50", alpha = 0.6) +
  # Linea di riferimento di Poisson (Varianza = Media) per vedere l'iperdispersione
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_smooth(method = "loess", color = "#e74c3c", linewidth = 1.2, se = FALSE, span=1, fill = "#e74c3c")+
  labs(
    x = "Media di OD",
    y = "Varianza di OD"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave(
  file.path(OUT_DIR, "relazione_m_var.png"),
  plot  = arrangeGrob(t_var, m_var, ncol = 2),
  width = 12, height = 5, dpi = DPI
)

fit <- lm(varianza_OD ~ media_OD + I(media_OD^2) - 1,
          data=df_summary)
fit
cbind(df_summary[,2]*coef(fit)[1] + (df_summary[,2]^2)*coef(fit)[2],
      df_summary[,3])

ggplot(df_summary, aes(x = media_OD, y = varianza_OD)) +
  geom_point(size = 3, color = "#0d0887", alpha = 0.7) + # Blu scuro plasma

  # 1. Spline Cubica (Giallo plasma)
  geom_smooth(
    method = "lm",
    formula = y ~ splines::bs(x, df = 3),
    color = "#f0f921",
    se = FALSE,
    linewidth = 1.2
  ) +

  # 2. Modello Mu + Beta*Mu^2 (Arancione plasma)
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x^2) -1,
    color = "#f89540",
    se = FALSE,
    linewidth = 1.2,
    linetype = "dashed" # Distingue visivamente i due modelli
  ) +
  labs(
    title = "Relazione Media-Varianza",
    subtitle = "Giallo: Spline | Arancione tratteggiato: Modello mu + beta*mu^2",
    x = "Media di OD",
    y = "Varianza di OD"
  ) +
  theme_minimal(base_size = 14)
