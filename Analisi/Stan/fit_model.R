# =============================================================
#  Script R per preparare i dati e stimare il modello Stan
# =============================================================
library(cmdstanr)   # oppure: library(rstan)
#install_cmdstan(overwrite = TRUE) # se è la prima volta
library(posterior)
library(dplyr)

# -------------------------------------------------------------
# 1. Prepara i dati
#    Adatta questa sezione al tuo data.frame reale.
#    Si assume un data.frame `df` con colonne:
#      subj  : indice soggetto (1..53)
#      j     : indice temporale (1..7, oppure 1..6 per i=12,15)
#      Y     : risposta
#      I_cov, D_cov, P_cov : covariate (costanti per soggetto)
# -------------------------------------------------------------

rm(list=ls());
setwd("C:/Users/Utente/OneDrive/Universita/Magistrale/2025-2026/Iterazione/Progetto/Analisi")
load("dati_modificati.Rdata")
rm(dati_incrementi, dati_incrementi_rapporti);gc();

str(dati)

S <- length(unique(dati$id_biomassa))
short <- unique(dati$id_biomassa[is.na(dati$OD)])
n_vec <- ifelse(seq_len(S) %in% short, 6L, 7L)

dati <- na.omit(dati)

subj_vec <- dati$id_biomassa
tidx_vec <- dati$tempo + 1

N <- sum(n_vec)
N

Y_vec <- dati |> pull(OD)

cov_df <- dati |> select(id_biomassa, I, D, P) |> mutate(
  I = (I - 280) / ((max(dati$I) - min(dati$I))/2),
  D = (D - 18) / ((max(dati$D) - min(dati$D))/2),
  P = (P - 25.5) / ((max(dati$P) - min(dati$P))/2)
) |> 
  rename(I_cov = I,
         D_cov = D,
         P_cov = P)
cov_df <- cov_df[(cumsum(c(1, n_vec))[-(S+1)]),]

# Creiamo la griglia con le 27 combinazioni (3 * 3 * 3)
new_data_grid <- expand.grid(
  I_cov = c(-1, 0, 1),
  D_cov = c(-1, 0, 1),
  P_cov = c(-1, 0, 1)
)
N_new <- nrow(new_data_grid) # Numero di nuove combinazioni (27)

rm(dati); gc();

stan_data <- list(
  N      = N,
  S      = S,
  n_i    = n_vec,
  id_biomassa = subj_vec,
  t_idx  = tidx_vec,
  Y      = Y_vec,
  I_cov  = cov_df$I_cov,
  D_cov  = cov_df$D_cov,
  P_cov  = cov_df$P_cov,
  I_new  = new_data_grid$I_cov,
  D_new  = new_data_grid$D_cov,
  P_new  = new_data_grid$P_cov,
  N_new  = N_new
)

set_cmdstan_path()

modello_semplice <- function(){
  mod <- cmdstan_model("Stan/model_base.stan")
  
  dir <- "Stan/stan_output_base"
  dir.create(dir, showWarnings = FALSE)
  fit <- mod$sample(
    # Passa la lista contenente i dati strutturati richiesti dal file .stan
    data            = stan_data,
    seed            = 123,
    chains          = 6, # Numero di catene indipendenti da far girare
    parallel_chains = 6, # Numero di core della CPU da allocare
    output_dir      = dir, # I file si salvano qui dentro in tempo reale
    iter_warmup     = 1000, # Iterazioni iniziali di "warm-up", che vengono scartate
    iter_sampling   = 2000, # Iterazioni di campionamento effettive per catena
    adapt_delta     = 0.95, # Target di accettazione
    max_treedepth   = 12, # Profondità massima dell'albero NUTS
    refresh         = 200 # Ogni quanto mostrare i log di avanzamento nella console di R
  )
  
  save(fit, stan_data, file="Stan/fit_mcmc_base.Rdata")
}

modello_avanzato <- function(){
  mod <- cmdstan_model("Stan/model_advanced.stan")
  
  dir <- "Stan/stan_output_advanced"
  dir.create(dir, showWarnings = FALSE)
  fit <- mod$sample(
    # Passa la lista contenente i dati strutturati richiesti dal file .stan
    data            = stan_data,
    seed            = 123,
    chains          = 6, # Numero di catene indipendenti da far girare
    parallel_chains = 6, # Numero di core della CPU da allocare
    output_dir      = dir, # I file si salvano qui dentro in tempo reale
    iter_warmup     = 1000, # Iterazioni iniziali di "warm-up", che vengono scartate
    iter_sampling   = 2000, # Iterazioni di campionamento effettive per catena
    adapt_delta     = 0.95, # Target di accettazione
    max_treedepth   = 12, # Profondità massima dell'albero NUTS
    refresh         = 200 # Ogni quanto mostrare i log di avanzamento nella console di R
  )
  save(fit, stan_data, file="Stan/fit_mcmc_advanced.Rdata")
}

modello_semplice(); modello_avanzato();
