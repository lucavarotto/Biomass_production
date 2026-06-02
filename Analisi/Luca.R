rm(list=ls());gc();

load("Analisi/dati_modificati.Rdata")

# bilanciamento disegno ----

library(ggplot2)
library(dplyr)

conteggi <- dati %>%
  group_by(I, D, P) %>%
  summarise(n_biomasse = n_distinct(id_biomassa), .groups = 'drop')

ggplot(conteggi, aes(x = as.factor(I), y = as.factor(D))) +
  geom_tile(aes(fill = n_biomasse), color = "white") +
  # Aggiungiamo il numero scritto dentro ogni cella
  geom_text(aes(label = n_biomasse), color = "white", fontface = "bold") +
  facet_wrap(~ P, labeller = label_both) +
  scale_fill_viridis_c(option = "magma", end=0.8) +
  labs(
    x = "Intensità Luminosa (I)",
    y = "Durata Esposizione (D)",
    fill = "N. Repliche"
  ) +
  theme_minimal()
ggsave("Grafici/disegno_sperimentale.pdf")
print("Immagine salvata")

# Modelli biologici ----

# install.packages("MicrobialGrowth")

#colnames(dati)
#library(MicrobialGrowth)
#g <- MicrobialGrowth(x, y, model = "gompertz")
#g


# Regressione gompertz con effetti casuali ----

# Asym rappresenta la capacità portante o il valore massimo teorico

# b2 (Parametro di posizione / Latenza
# Controlla lo spostamento della curva lungo l'asse X (il tempo).
# Un valore di b2 più alto "sposta" la curva verso destra,
# indicando che la crescita inizia più tardi.
# b2 è legato al valore della variabile y quando il tempo t = 0

# b3 rappresenta la velocità o il tasso di crescita intrinseco.
# Se b3 è alto, la curva sarà molto ripida (crescita esplosiva); se b3 è basso, la curva salirà molto dolcemente.

dati_m_gompertz <- na.omit(dati)
dati_m_gompertz <- dati_m_gompertz |> mutate(
  I = (I - 280) / ((max(dati$I) - min(dati$I))/2),
  D = (D - 18 / ((max(dati$D) - min(dati$D))/2)),
  P = (P - 25.5) / ((max(dati$P) - min(dati$P))/2)
)

library(nlme)

modello_gompertz <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)^2+P,
    b2 ~ (I+D)^2+P,           
    b3 ~ (I+D)^2+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  start = c(
    Asym = c(14, rep(0, 4)), 
    b2 = c(5.5, rep(0, 4)),     
    b3 = c(0.75, rep(0, 4))),
  control = nlmeControl(maxIter = 350, # Il numero massimo di iterazioni durante la
                        # fase di ottimizzazione della stima di massima verosimiglianza
                        pnlsMaxIter = 25, # il modello tiene fissi
                        # gli effetti casuali e cerca di aggiornare
                        # i parametri fissi per minimizzare l'errore.
                        msMaxIter = 200, # si riferisce all'ottimizzatore numerico
                        msMaxEval = 500, # Aumentiamo il numero di valutazioni della funzione obiettivo
                        niterEM = 50, # importante per trovare buoni punti di partenza per le varianze
                        opt = "nlminb")
)
modello_gompertz2 <- update(modello_gompertz,
                            random = b2 ~ 1 | id_biomassa)
modello_gompertz3 <- update(modello_gompertz,
                            random = b3 ~ 1 | id_biomassa)
modello_gompertz4 <- update(modello_gompertz,
                            random = Asym + b3 ~ 1 | id_biomassa)
#modello_gompertz6 <- update(modello_gompertz,
#                            random = list(id_biomassa = pdDiag(Asym + b2 ~ 1)))
AIC(modello_gompertz, modello_gompertz2, modello_gompertz3,
    modello_gompertz4)

?MuMIn::dredge()

summary(modello_gompertz)

dati_m_gompertz$predicted <- predict(modello_gompertz)
library(ggforce)
x11()
ggplot(dati_m_gompertz, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = OD), alpha = 0.6) +
  geom_line(aes(y = OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted, group = id_biomassa), linewidth = 1.2)+
  facet_wrap_paginate(~ id_biomassa, scales = "free_y",
                      ncol = 3, nrow = 3, page = 1) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Confronto Dati Reali vs Stima Gompertz",
    subtitle = "I punti indicano le osservazioni reali; le linee continue rappresentano il modello NLME",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)"
  ) +
  theme_minimal()

X_m <- model.matrix(~(I+D)^2+P, data=dati_m_gompertz[1,])
dim(X_m)

f <- matrix(modello_gompertz$coefficients$fixed, nrow=3, byrow=T)
dim(f)

re <- modello_gompertz$coefficients$random$id_biomassa

eval_gomp <- function(m, t, re=0){
  Asym <- m[1]
  b2 <- m[2]
  b3 <- m[3]
  (Asym  + re)*exp(-b2*(b3^t))
}

eval_gomp(f %*% t(X_m), 5, re[1])

## Regressione gompertz no interazione ----

gompertz_no_int <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)+P,
    b2 ~ (I+D)+P,           
    b3 ~ (I+D)+P     
  ),
  random = Asym ~ 1 | id_biomassa,
  start = c(
    Asym = c(14, rep(0, 3)), 
    b2 = c(5.5, rep(0, 3)),     
    b3 = c(0.75, rep(0, 3)))
)

gompertz_no_int_AR <- nlme(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  fixed = list(
    Asym ~ (I+D)+P,
    b2 ~ (I+D)+P,           
    b3 ~ (I+D)+P     
  ),
  correlation = corAR1(form = ~ tempo | id_biomassa),
  random = Asym ~ 1 | id_biomassa,
  start = c(
    Asym = c(14, rep(0, 3)), 
    b2 = c(5.5, rep(0, 3)),     
    b3 = c(0.75, rep(0, 3)))
)

anova(gompertz_no_int, gompertz_no_int_AR)

2*(-gompertz_no_int$logLik + gompertz_no_int_AR$logLik)

1-pchisq(2*(-gompertz_no_int$logLik + gompertz_no_int_AR$logLik), df = 1)

AIC(modello_gompertz, gompertz_no_int, gompertz_no_int_AR)
BIC(modello_gompertz, gompertz_no_int, gompertz_no_int_AR)

# Regressione gompertz solo fisso ----

modello_gompertz_fisso <- gnls(
  model = OD ~ SSgompertz(tempo, Asym, b2, b3),
  data = dati_m_gompertz,
  params = list(
    Asym ~ (I + D)^2+P, # Interazione per l'asintoto
    b2 ~ (I + D)^2+P, # Interazione per la latenza
    b3 ~ (I + D)^2+P # Interazione per la rapidità
  ),
  start = c(
    Asym = c(3, rep(0, 4)),  
    b2   = c(1.5, rep(0, 4)),               
    b3   = c(0.8, rep(0, 4))
  )
)

AIC(modello_gompertz, modello_gompertz_fisso)

anova(modello_gompertz_fisso)

summary(modello_gompertz_fisso)
coef(modello_gompertz_fisso) |> round(4)

dati_m_gompertz$predicted_fisso <- predict(modello_gompertz)

library(ggforce)
ggplot(dati_m_gompertz, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = OD), alpha = 0.6) +
  geom_line(aes(y = OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted_fisso, group = id_biomassa), linewidth = 1.2)+
  facet_wrap_paginate(~ id_biomassa, scales = "free_y", ncol = 3, nrow = 3, page = 2) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Confronto Dati Reali vs Stima Gompertz",
    subtitle = "I punti indicano le osservazioni reali; le linee continue rappresentano il modello NLME",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)"
  ) + theme_minimal()

matrix(coef(modello_gompertz_fisso), nrow=3, byrow=T) %*%
  t(as.matrix(model.matrix(~(I + D)^2+P, data=dati[1,])))


# Modello parametrico sugli incrementi ----

dati_lme <- na.omit(dati_incrementi)
dati_lme <- dati_incrementi |> mutate(
  I = (I - 280) / ((max(dati$I) - min(dati$I))/2),
  D = (D - 18 / ((max(dati$D) - min(dati$D))/2)),
  P = (P - 25.5) / ((max(dati$P) - min(dati$P))/2)
)
dati_lme$tempo <- dati_lme$tempo + 1
dati_lme$tempo2 <- dati_lme$tempo^2

lm0 <- lm(OD ~ tempo + tempo2 +
            I*D + P + I(I^2):I(D^2) + I(P^2),
          data=dati_lme)
summary(lm0)
car::vif(lm0, type="predictor")

library(lme4)
colnames(dati_lme)
m0 <- lme(incremento_OD ~ tempo + (I + D)^2+P + tempo:I + tempo:D + tempo:P,
          random = ~ 1 + tempo | as.factor(id_biomassa), data=dati_lme,
          control = lmeControl(maxIter = 1000, # Il numero massimo di iterazioni durante la
                               # fase di ottimizzazione della stima di massima verosimiglianza
                               msMaxIter = 100, # si riferisce all'ottimizzatore numerico
                               niterEM = 50,
                               msMaxEval = 400,
                               opt = "optim")
)
m1 <- lme(incremento_OD ~ tempo + I((log(tempo))^2) + (I + D)^2+P +
            I((log(tempo))^2):I + I((log(tempo))^2):D + I((log(tempo))^2):P,
          random = ~ 1 + tempo | as.factor(id_biomassa), data=dati_lme,
          control = lmeControl(maxIter = 1000, # Il numero massimo di iterazioni durante la
                               # fase di ottimizzazione della stima di massima verosimiglianza
                               msMaxIter = 100, # si riferisce all'ottimizzatore numerico
                               niterEM = 50,
                               msMaxEval = 400,
                               opt = "optim")
)
m2 <- lme(incremento_OD ~ tempo + I(tempo^2) + I(tempo^3) +
            (I + D)^2+P +
            I(tempo):I + I(tempo):D + I(tempo):P +
            I(tempo^2):I + I(tempo^2):D + I(tempo^2):P +
            I(tempo^3):I + I(tempo^3):D + I(tempo^3):P,
          random = ~ 1 + tempo | as.factor(id_biomassa), data=dati_lme,
          control = lmeControl(maxIter = 1000, # Il numero massimo di iterazioni durante la
                               # fase di ottimizzazione della stima di massima verosimiglianza
                               msMaxIter = 100, # si riferisce all'ottimizzatore numerico
                               niterEM = 50,
                               msMaxEval = 400,
                               opt = "optim")
)
m3 <- lme(incremento_OD ~ tempo + I(tempo^2) + I(tempo^3) +
            (I + D)^2+P,
          random = ~ 1 + tempo | as.factor(id_biomassa), data=dati_lme,
          control = lmeControl(maxIter = 1000, # Il numero massimo di iterazioni durante la
                               # fase di ottimizzazione della stima di massima verosimiglianza
                               msMaxIter = 100, # si riferisce all'ottimizzatore numerico
                               niterEM = 50,
                               msMaxEval = 400,
                               opt = "optim")
)
summary(m3)
AIC(m0, m1, m2, m3)



library(ggforce)
library(ggplot2)

dati_lme$predicted <- fitted(m3)
ggplot(dati_lme, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = incremento_OD), alpha = 0.6) +
  geom_line(aes(y = incremento_OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted, group = id_biomassa), linewidth = 1.2) +
  facet_wrap_paginate(~ id_biomassa, scales = "free_y", ncol = 3, nrow = 3, page = 3) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Confronto incrementi veri vs stimati",
    subtitle = "Punti: Osservazioni reali | Linea continua: Predizione specifica per biomassa",
    x = "Tempo (giorni)",
    y = "Optical Density (OD)"
  ) + 
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))


N0 <- 0.01
dati_lme <- dati_lme %>%
  group_by(id_biomassa) %>%
  arrange(tempo) %>%
  mutate(
    # Somma cumulata degli incrementi predetti + il valore iniziale reale N0
    predicted_OD_cum = N0 + cumsum(predicted)
  ) %>%
  ungroup()

ggplot(dati_lme, aes(x = tempo, color = as.factor(condizione_sperimentale))) +
  geom_point(aes(y = OD), alpha = 0.6) +
  geom_line(aes(y = OD, group = id_biomassa), linetype = "dashed", alpha = 0.4) +
  geom_line(aes(y = predicted_OD_cum, group = id_biomassa), linewidth = 1.2) +
  facet_wrap_paginate(~ id_biomassa, scales = "free_y", ncol = 3, nrow = 3, page = 5) +
  scale_color_viridis_d(option = "plasma", guide = "none") +
  labs(
    title = "Ricostruzione della Curva di Crescita (da modello LME su incrementi)",
    subtitle = "Punti: OD osservata | Linea continua: Somma cumulata degli incrementi stimati",
    x = "Tempo (giorni)",
    y = "Optical Density (OD Totale)"
  ) + 
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))




# modello normale ----

# Modello gerarchico non lineare
# Y_i ~ N_{n_i}(mu_i(theta), V_i), i = 1,...,53
# mu_ij = Asym_i * exp(-b2_i * b3_i^tempo)
# Asym_i = eta_i(I,D,P) + u_i,  u_i ~ N(0,1)
# b2_i   = eta_i(I,D,P)
# b3_i   = eta_i(I,D,P)
# eta_i  = beta0 + beta1*I_i + beta2*D_i + beta3*P_i + beta4*I_i*D_i

dati <- na.omit(dati)
dati$n_i <- 7
dati$n_i[dati$id_biomassa == 12 | dati$id_biomassa == 15] <- 6
if (min(dati$tempo) == 0){ dati$tempo <- dati$tempo + 1}

library(nlme)

mean_fun <- function(tempo, Asym, b2, b3) {
  Asym * exp(-b2 * b3^tempo)
}

eta <- function(beta, I, D, P) {
  beta[1] + beta[2]*I + beta[3]*D + beta[4]*P + beta[5]*I*D
}

# Log-verosimiglianza marginale
# (effetti fissi + effetto casuale su Asym)
# u_i ~ N(0,1) => integrazione numerica su u_i
log_lik <- function(par, data) {
  beta_Asym <- par[1:5]   # coefficienti per Asym
  beta_b2   <- par[6:10]  # coefficienti per b2
  beta_b3   <- par[11:15] # coefficienti per b3
  log_sigma <- par[16]    # log(sigma) residuo (per gruppo)
  sigma     <- exp(log_sigma)
  
  groups <- unique(data$id_biomassa)
  ll <- 0
  
  for (i in groups) {
    
    sub  <- data[data$id_biomassa == i, ]
    
    I_i  <- sub$I[1]; D_i <- sub$D[1]; P_i <- sub$P[1]
    
    eta_Asym <- eta(beta_Asym, I_i, D_i, P_i)
    b2_i     <- eta(beta_b2,   I_i, D_i, P_i)
    b3_i     <- eta(beta_b3,   I_i, D_i, P_i)
    
    # Integrazione su u_i ~ N(0,1) con quadratura di Gauss-Hermite
    gh        <- statmod::gauss.quad.prob(20, dist = "normal")
    nodes     <- gh$nodes          # quantili di N(0,1)
    weights   <- gh$weights
    
    ll_i <- 0
    for (k in seq_along(nodes)) {
      Asym_ik <- eta_Asym + nodes[k]   # Asym_i = eta + u_i
      
      mu_ik <- mean_fun(sub$tempo, Asym_ik, b2_i, b3_i)
      resid <- sub$OD - mu_ik
      # densità normale del vettore Y_i | u_i (V_i = sigma^2 * I)
      log_f <- sum(dnorm(resid, mean = 0, sd = sigma, log = TRUE))
      ll_i  <- ll_i + weights[k] * exp(log_f)
    }
    ll <- ll + log(ll_i + .Machine$double.eps)
  }
  return(-ll)   # restituisce la neg-log-lik per optim()
}

# 5.  Stima con optim()
par_init <- rep(0, 16)
par_init[16] <- log(1)   # sigma iniziale = 1

fit <- optim(
  par     = par_init,
  fn      = log_lik,
  data    = dati,
  method  = "BFGS",
  control = list(maxit = 2000, trace = 1),
  hessian = TRUE
)

parametri <- matrix(fit$par[1:15], nrow=3, byrow=T)
dim(parametri)

X_m <- model.matrix(~I*D+P, data=dati)
X_m <- X_m[!duplicated(X_m),]

apply(X_m, 1, function(x) c(x %*% t(parametri))) |> View()

# errori standard da hessiano
se <- sqrt(diag(solve(fit$hessian)))

cat("Stime:", fit$par, "\n")
cat("SE:   ", se,      "\n")

# 6.  Alternativa con nlme (se V_i = sigma^2 * I e random solo su Asym)
fit_nlme <- nlme(
  model  = OD ~ Asym * exp(-b2 * b3^tempo),
  data   = df,
  fixed  = Asym + b2 + b3 ~ I + D + P + I:D,
  random = Asym ~ 1 | id_biomassa,     # u_i ~ N(0, psi^2) stimato
  start  = c(
    Asym = c(0,0,0,0,0),
    b2   = c(0,0,0,0,0),
    b3   = c(0,0,0,0,0)
  )
)
summary(fit_nlme)







# modello normale senza RE ----

# Modello non lineare con errori AR(1) — senza effetti casuali
# Y_i ~ N_{n_i}(mu_i(theta), V_i),   V_i = sigma^2 * R_i(rho)
# mu_ij = Asym_i * exp(-b2_i * b3_i^tempo)
# Asym_i = eta_i(I,D,P)
# b2_i   = eta_i(I,D,P)
# b3_i   = eta_i(I,D,P)
# eta_i  = beta0 + beta1*I_i + beta2*D_i + beta3*P_i

library(statmod)

dati <- na.omit(dati)
dati$n_i <- 7
dati$n_i[dati$id_biomassa == 12 | dati$id_biomassa == 15] <- 6
if (min(dati$tempo) == 0){ dati$tempo <- dati$tempo + 1}

# Funzioni di supporto

mean_fun <- function(tempo, Asym, b2, b3) {
  Asym * exp(-b2 * b3^tempo)
}

eta <- function(beta, I, D, P) {
  beta[1] + beta[2]*I + beta[3]*D + beta[4]*P
}

# Matrice di correlazione AR(1) di dimensione n
ar1_cor <- function(n, rho) {
  exp_mat <- abs(outer(1:n, 1:n, "-"))
  rho^exp_mat
}

# Log-verosimiglianza
log_lik <- function(par, data) {
  
  beta_Asym <- par[1:4]
  beta_b2   <- par[5:8]
  beta_b3   <- par[9:12]
  sigma     <- exp(par[13])          # sigma > 0
  rho       <- tanh(par[14])         # rho in (-1, 1)
  
  groups <- unique(data$id_biomassa)
  ll     <- 0
  
  for (i in groups) {
    
    sub <- data[data$id_biomassa == i, ]
    n_i <- nrow(sub)
    
    I_i <- sub$I[1];  D_i <- sub$D[1];  P_i <- sub$P[1]
    
    Asym_i <- eta(beta_Asym, I_i, D_i, P_i)   # nessun u_i
    b2_i   <- eta(beta_b2,   I_i, D_i, P_i)
    b3_i   <- eta(beta_b3,   I_i, D_i, P_i)
    
    mu_i  <- mean_fun(sub$tempo, Asym_i, b2_i, b3_i)
    resid <- sub$OD - mu_i
    
    # V_i = sigma^2 * R_i(rho)
    R_i   <- ar1_cor(n_i, rho)
    V_i   <- sigma^2 * R_i
    
    # log-densità N_{n_i}(0, V_i) valutata nei residui
    ll_i  <- mvtnorm::dmvnorm(resid, mean = rep(0, n_i), sigma = V_i, log = TRUE)
    ll    <- ll + ll_i
  }
  
  return(-ll)
}

# Stima con optim()

par_init        <- rep(0, 14)
par_init[13]    <- log(0.5)     # sigma = 1
par_init[14]    <- atanh(0.3) # rho   = 0.3

fit <- optim(
  par     = par_init,
  fn      = log_lik,
  data    = dati,
  method  = "BFGS",
  control = list(maxit = 2000, trace = 1),
  hessian = TRUE
)

parametri <- matrix(fit$par[1:12], nrow=3, byrow=T)
dim(parametri)

X_m <- model.matrix(~I+D+P, data=dati)
X_m <- X_m[!duplicated(X_m),]

apply(X_m, 1, function(x) c(x %*% t(parametri))) |> View()

# Errori standard da hessiano (scala originale via delta method se serve)
se <- sqrt(diag(solve(fit$hessian)))

# Recupera parametri in scala leggibile
stime        <- fit$par
stime[13]    <- exp(fit$par[13])          # sigma
stime[14]    <- tanh(fit$par[14])         # rho

nomi <- c(
  paste0("Asym_beta", 0:3),
  paste0("b2_beta",   0:3),
  paste0("b3_beta",   0:3),
  "sigma", "rho"
)

result <- data.frame(
  parametro = nomi,
  stima     = round(stime, 4),
  se_raw    = round(se, 4)
)
print(result)




fit_nlme <- nlme(
  model       = OD ~ Asym * exp(-b2 * b3^tempo),
  data        = dati,
  fixed       = Asym + b2 + b3 ~ I + D + P,
  correlation = corAR1(form = ~ tempo | id_biomassa),
  random    = list(id_biomassa = pdDiag(Asym ~ 0)),
  method      = "ML", # <--- essenziale per LRT e AIC confrontabili
  start       = c(
    Asym = c(0,0,0,0),
    b2   = c(0,0,0,0),
    b3   = c(0,0,0,0)
  )
)

