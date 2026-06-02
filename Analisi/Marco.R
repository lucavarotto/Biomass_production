rm(list=ls());gc();

load("Analisi/dati_modificati.Rdata")
library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)
library(MuMIn)
library(nlme)

# # Caricamento e pulizia
# dati <- readxl::read_xlsx("Analisi/nostoc.xlsx", skip = 2)
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
