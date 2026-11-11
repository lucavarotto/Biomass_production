// =============================================================
//  Modello nonlineare misto con varianza eteroschedastica a due parametri
//  Y_i ~ prodotto di normali indipendenti con Var dipendente da mu
//
//  Struttura della media:
//    mu_ij = Asym_i * exp(-b2_i * b3_i^j)
//    Asym_i = exp(eta_Asym_i + u_i),   u_i ~ N(0, sigma2_u)
//    b2_i   = exp(eta_b2_i)
//    b3_i   = inv_logit(eta_b3_i)
//
//  Predittori lineari ESTESI con termini quadratici:
//    eta_param_i = beta0 + beta1*I + beta2*D + beta3*P + beta4*I*D
//                  + beta5*I^2 + beta6*D^2 + beta7*P^2
//    (i termini quadratici hanno prior N(0, 0.1): regularizzazione)
//
//  Struttura di varianza ETEROSCHEDASTICA con parametri stimati:
//    Var(Y_ij) = phi1 * mu_ij + phi2 * mu_ij^2  +  1e-5  (nugget)
//
//    dove phi1 >= 0 e phi2 >= 0 sono parametri liberi stimati dall'MCMC.
//
//    Interpretazione dei casi limite:
//      phi1 > 0, phi2 = 0  => varianza proporzionale alla media   (Poisson-like)
//      phi1 = 0, phi2 > 0  => varianza proporzionale a mu^2       (CV costante)
//      phi1 > 0, phi2 > 0  => forma quadratica generale           (NB generalizzata)
//
//    Il modello precedente (phi1 = phi2 = 1) è un caso speciale annidato.
//
//    Il nugget 1e-5 sotto la radice impedisce sd = 0 se mu_ij -> 0,
//    garantendo stabilità numerica della Cholesky interna di normal_lpdf.
//
//  Effetto casuale: solo su Asym (parametrizzazione non centrata)
// =============================================================

data {
  int<lower=1> N;                     // numero totale osservazioni (= 53*7 - 2 = 369)
  int<lower=1> S;                     // numero soggetti (= 53)
  array[S] int<lower=1> n_i;          // n_i: 6 per i=12,15; 7 per gli altri
  array[N] int<lower=1> id_biomassa;  // indice soggetto per ogni osservazione
  array[N] int<lower=1> t_idx;        // indice temporale j (1..n_i) dentro ogni soggetto
  vector[N] Y;                        // variabile risposta (in ordine soggetto, poi tempo)

  // Covariate (una per soggetto)
  vector[S] I_cov;                    // covariata I
  vector[S] D_cov;                    // covariata D
  vector[S] P_cov;                    // covariata P

  // --- DATI PER LA PREVISIONE ---
  int<lower=0> N_new;                 // numero nuove unità (sarà 27)
  vector[N_new] I_new;
  vector[N_new] D_new;
  vector[N_new] P_new;
}

// Viene eseguito una sola volta prima dell'inizio del campionamento MCMC.
// Calcola gli indici di inizio e fine di ogni soggetto nel vettore globale Y.
transformed data {
  array[S] int start;
  array[S] int end;

  {
    int pos = 1;
    for (s in 1:S) {
      start[s] = pos;
      end[s]   = pos + n_i[s] - 1;
      pos      = pos + n_i[s];
    }
  }
}

// Parametri MCMC da stimare.
//
// Rispetto al modello precedente (Var = mu + mu^2 fissa):
//  - AGGIUNTI: phi1, phi2 >= 0, i coefficienti della funzione di varianza
//    Var = phi1*mu + phi2*mu^2. Renderli parametri liberi (invece di fissarli
//    a 1) consente al modello di stimare la forma della dispersione
//    direttamente dai dati, senza imporre a priori né Poisson né CV costante.
//
// Tutti gli altri parametri (beta*, z, sigma_u) sono invariati.
parameters {
  // --- EFFETTI FISSI: predittore di Asym (scala log) ---
  real beta0_Asym;   // Intercetta
  real beta1_Asym;   // Effetto lineare covariata I
  real beta2_Asym;   // Effetto lineare covariata D
  real beta3_Asym;   // Effetto lineare covariata P
  real beta4_Asym;   // Effetto interazione I * D
  real beta5_Asym;   // Effetto quadratico I^2  (prior fortemente regularizzante)
  real beta6_Asym;   // Effetto quadratico D^2  (prior fortemente regularizzante)
  real beta7_Asym;   // Effetto quadratico P^2  (prior fortemente regularizzante)

  // --- EFFETTI FISSI: predittore di b2 (scala log) ---
  real beta0_b2;     // Intercetta
  real beta1_b2;     // Effetto lineare covariata I
  real beta2_b2;     // Effetto lineare covariata D
  real beta3_b2;     // Effetto lineare covariata P
  real beta4_b2;     // Effetto interazione I * D
  real beta5_b2;     // Effetto quadratico I^2  (prior fortemente regularizzante)
  real beta6_b2;     // Effetto quadratico D^2  (prior fortemente regularizzante)
  real beta7_b2;     // Effetto quadratico P^2  (prior fortemente regularizzante)

  // --- EFFETTI FISSI: predittore di b3 (scala logit, b3 in (0,1)) ---
  real beta0_b3;     // Intercetta
  real beta1_b3;     // Effetto lineare covariata I
  real beta2_b3;     // Effetto lineare covariata D
  real beta3_b3;     // Effetto lineare covariata P
  real beta4_b3;     // Effetto interazione I * D
  real beta5_b3;     // Effetto quadratico I^2  (prior fortemente regularizzante)
  real beta6_b3;     // Effetto quadratico D^2  (prior fortemente regularizzante)
  real beta7_b3;     // Effetto quadratico P^2  (prior fortemente regularizzante)

  // --- EFFETTO CASUALE su Asym (parametrizzazione non centrata) ---
  // z[s] ~ N(0,1); l'effetto random effettivo è u[s] = sigma_u * z[s]
  vector[S] z;

  // --- PARAMETRO DI VARIANZA DELL'EFFETTO CASUALE ---
  // sigma_u: deviazione standard dell'effetto random su Asym
  real<lower=0> sigma_u;

  // --- PARAMETRI DELLA FUNZIONE DI VARIANZA ETEROSCHEDASTICA ---
  // Var(Y_ij) = phi1 * mu_ij + phi2 * mu_ij^2
  //
  // Entrambi vincolati a essere non negativi: una varianza negativa
  // non ha senso fisico e causerebbe sd < 0, invalidando normal_lpdf.
  //
  // phi1 ~ Gamma(2, 1): prior debolmente informativa, media = 2, moda = 1.
  //   Favorisce valori dell'ordine di 1 (come nel modello precedente),
  //   ma lascia libertà di allontanarsene se i dati lo supportano.
  // phi2 ~ Gamma(2, 1): stessa logica per il termine quadratico.
  //
  // Nota: usare Gamma (supporto > 0) invece di half-Normal è più naturale
  // qui perché il gradiente di Gamma a zero è finito (con shape > 1),
  // riducendo il rischio di catene bloccate vicino al boundary.
  real<lower=0> phi1;   // coefficiente del termine lineare nella Var (effetto Poisson-like)
  real<lower=0> phi2;   // coefficiente del termine quadratico nella Var (effetto CV costante)
}

// Calcola i predittori lineari e i parametri individuali.
// La struttura è identica al modello base, con l'aggiunta dei termini
// quadratici I^2, D^2, P^2 nei predittori lineari.
transformed parameters {

  // Predittori lineari della popolazione
  vector[S] eta_Asym;
  vector[S] eta_b2;
  vector[S] eta_b3;

  // Parametri individuali finali (con vincoli di segno/ampiezza)
  vector<lower=0>[S] Asym_i;
  vector<lower=0>[S] b2_i;
  vector<lower=0, upper=1>[S] b3_i;

  // Effetto random su scala originale
  vector[S] u;   // u[s] = sigma_u * z[s]

  for (s in 1:S) {
    real id  = I_cov[s] * D_cov[s];      // termine interazione I*D
    real I2  = square(I_cov[s]);          // termine quadratico I^2
    real D2  = square(D_cov[s]);          // termine quadratico D^2
    real P2  = square(P_cov[s]);          // termine quadratico P^2

    // 1. Predittori lineari estesi (effetti fissi di popolazione)
    eta_Asym[s] = beta0_Asym
                + beta1_Asym * I_cov[s]  + beta2_Asym * D_cov[s]
                + beta3_Asym * P_cov[s]  + beta4_Asym * id
                + beta5_Asym * I2        + beta6_Asym * D2
                + beta7_Asym * P2;

    eta_b2[s]   = beta0_b2
                + beta1_b2 * I_cov[s]    + beta2_b2 * D_cov[s]
                + beta3_b2 * P_cov[s]    + beta4_b2 * id
                + beta5_b2 * I2          + beta6_b2 * D2
                + beta7_b2 * P2;

    eta_b3[s]   = beta0_b3
                + beta1_b3 * I_cov[s]    + beta2_b3 * D_cov[s]
                + beta3_b3 * P_cov[s]    + beta4_b3 * id
                + beta5_b3 * I2          + beta6_b3 * D2
                + beta7_b3 * P2;

    // 2. Effetto random (non centrato -> centrato)
    u[s] = sigma_u * z[s];

    // 3. Parametri individuali
    Asym_i[s] = exp(eta_Asym[s] + u[s]);   // log-link + effetto casuale
    b2_i[s]   = exp(eta_b2[s]);             // log-link
    b3_i[s]   = inv_logit(eta_b3[s]);       // logit-link
  }
}

model {


  // Asym (scala log)
  target += normal_lpdf(beta0_Asym | 2,   1.0);
  target += normal_lpdf(beta1_Asym | 0,   0.1);
  target += normal_lpdf(beta2_Asym | 0,   0.1);
  target += normal_lpdf(beta3_Asym | 0,   0.1);
  target += normal_lpdf(beta4_Asym | 0,   0.1);
  target += normal_lpdf(beta5_Asym | 0,   0.1);
  target += normal_lpdf(beta6_Asym | 0,   0.1);
  target += normal_lpdf(beta7_Asym | 0,   0.1);

  // b2 (scala log)
  target += normal_lpdf(beta0_b2 | 1.61, 1.0);
  target += normal_lpdf(beta1_b2 | 0,    0.1);
  target += normal_lpdf(beta2_b2 | 0,    0.1);
  target += normal_lpdf(beta3_b2 | 0,    0.1);
  target += normal_lpdf(beta4_b2 | 0,    0.1);
  target += normal_lpdf(beta5_b2 | 0,    0.1);
  target += normal_lpdf(beta6_b2 | 0,    0.1);
  target += normal_lpdf(beta7_b2 | 0,    0.1);

  // b3 (scala logit)
  target += normal_lpdf(beta0_b3 | 2,    0.5);
  target += normal_lpdf(beta1_b3 | 0,    0.1);
  target += normal_lpdf(beta2_b3 | 0,    0.1);
  target += normal_lpdf(beta3_b3 | 0,    0.1);
  target += normal_lpdf(beta4_b3 | 0,    0.1);
  target += normal_lpdf(beta5_b3 | 0,    0.1);
  target += normal_lpdf(beta6_b3 | 0,    0.1);
  target += normal_lpdf(beta7_b3 | 0,    0.1);

  // -------------------------------------------------------
  // PRIOR — varianza dell'effetto casuale
  // Half-Normal(0, 0.3): mette la massa principale in [0, 0.6],
  // debolmente informativo ma riduce sovra-stima della varianza.
  // -------------------------------------------------------
  target += normal_lpdf(sigma_u | 0, 0.3) - normal_lccdf(0 | 0, 0.3);

  // -------------------------------------------------------
  // PRIOR — effetti random standardizzati (parametrizzazione non centrata)
  // z[s] ~ N(0,1); implicitamente u[s] ~ N(0, sigma_u^2)
  // -------------------------------------------------------
  target += std_normal_lpdf(z);

  // -------------------------------------------------------
  // PRIOR — parametri della funzione di varianza (phi1, phi2)
  //
  // Gamma(2, 1): media = shape/rate = 2, varianza = shape/rate^2 = 2.
  // La scelta shape = 2 (> 1) garantisce che la densità sia zero in zero
  // e abbia un massimo (moda = 1) lontano dal boundary, migliorando
  // la geometria della posterior per HMC rispetto a Exponential o
  // half-Normal che hanno la moda esattamente in zero.
  // -------------------------------------------------------
  target += gamma_lpdf(phi1 | 2, 1);
  target += gamma_lpdf(phi2 | 2, 1);

  // -------------------------------------------------------
  // VEROSIMIGLIANZA: normale eteroschedastica indipendente
  //
  // Ogni osservazione Y_ij ~ N(mu_ij, sd_ij) con:
  //   Var_ij = phi1 * mu_ij + phi2 * mu_ij^2 + 1e-5
  //   sd_ij  = sqrt(Var_ij)
  //
  // Il nugget 1e-5 è una costante additiva sotto la radice che:
  //   (a) impedisce sd = 0 quando mu_ij -> 0 (stabilità numerica);
  //   (b) ha effetto trascurabile per mu_ij >> sqrt(1e-5) ≈ 0.003,
  //       che è il caso tipico per biomasse.
  //
  // La matrice V_i è DIAGONALE: le osservazioni sono
  // condizionatamente indipendenti dato mu_i (nessuna struttura AR).
  // -------------------------------------------------------
  for (s in 1:S) {
    int ni = n_i[s];

    for (j in 1:ni) {
      // Media del soggetto s al time-point j
      real mu_ij = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));

      // Varianza eteroschedastica a due parametri: Var = phi1*mu + phi2*mu^2 + 1e-5
      // Il nugget 1e-5 garantisce sd > 0 anche se mu_ij -> 0.
      real sd_ij = sqrt(phi1 * mu_ij + phi2 * square(mu_ij) + 1e-5);

      // Contributo alla log-verosimiglianza
      target += normal_lpdf(Y[start[s] + j - 1] | mu_ij, sd_ij);
    }
  }
}

// -----------------------------------------------------------------------------
// BLOCCO GENERATED QUANTITIES
//
// Calcola:
//  - log_lik[s]: log-verosimiglianza per soggetto (per LOO-CV via loo())
//  - Y_rep[n]:   repliche predittive posteriori (per PPC)
//  - mu_new, Y_new: previsioni sulle 27 nuove combinazioni covariate
//
// La struttura eteroschedastica (Var = phi1*mu + phi2*mu^2 + 1e-5) è
// replicata fedelmente sia per log_lik che per il campionamento predittivo:
// i valori di phi1 e phi2 usati sono quelli campionati dall'MCMC in ogni
// iterazione, propagando così la loro incertezza nella predittiva posteriore.
// -----------------------------------------------------------------------------
generated quantities {

  vector[S] log_lik;  // Log-verosimiglianza aggregata per soggetto (usata da loo())
  vector[N] Y_rep;    // Repliche PPC sui punti osservati storici

  { // Blocco locale: isola le variabili temporanee dall'output finale
    for (s in 1:S) {
      int ni = n_i[s];

      // Accumula la log-verosimiglianza del soggetto s
      // come somma dei contributi indipendenti di ogni time-point
      real ll_s = 0;

      for (j in 1:ni) {
        real mu_ij = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));

        // Stessa formula di varianza del blocco model: coerenza fondamentale
        // per LOO e PPC. Usare una formula diversa qui invaliderebbe log_lik.
        real sd_ij = sqrt(phi1 * mu_ij + phi2 * square(mu_ij) + 1e-5);

        // Log-verosimiglianza dell'osservazione (j, s)
        ll_s += normal_lpdf(Y[start[s] + j - 1] | mu_ij, sd_ij);

        // Campione predittivo posteriore per l'osservazione (j, s)
        Y_rep[start[s] + j - 1] = normal_rng(mu_ij, sd_ij);
      }

      log_lik[s] = ll_s;
    }
  }

  // ---- PREVISIONE SULLE NUOVE 27 COMBINAZIONI (TRAIETTORIE TEMPORALI) ----
  // Per ogni nuova unità si campiona un effetto random da N(0, sigma_u^2)
  // per propagare l'incertezza inter-individuale nelle previsioni.
  matrix[N_new, 7] mu_new;  // Media predittiva per le 7 settimane
  matrix[N_new, 7] Y_new;   // Campione predittivo (mu + rumore eteroschedastico)

  {
    int n_t = 7; // numero di time-point previsti

    for (n in 1:N_new) {
      real id_new = I_new[n] * D_new[n];
      real I2_new = square(I_new[n]);
      real D2_new = square(D_new[n]);
      real P2_new = square(P_new[n]);

      // Predittori lineari di popolazione per la nuova unità (con termini quadratici)
      real eta_Asym_new = beta0_Asym
                        + beta1_Asym * I_new[n]  + beta2_Asym * D_new[n]
                        + beta3_Asym * P_new[n]  + beta4_Asym * id_new
                        + beta5_Asym * I2_new    + beta6_Asym * D2_new
                        + beta7_Asym * P2_new;

      real eta_b2_new   = beta0_b2
                        + beta1_b2 * I_new[n]    + beta2_b2 * D_new[n]
                        + beta3_b2 * P_new[n]    + beta4_b2 * id_new
                        + beta5_b2 * I2_new      + beta6_b2 * D2_new
                        + beta7_b2 * P2_new;

      real eta_b3_new   = beta0_b3
                        + beta1_b3 * I_new[n]    + beta2_b3 * D_new[n]
                        + beta3_b3 * P_new[n]    + beta4_b3 * id_new
                        + beta5_b3 * I2_new      + beta6_b3 * D2_new
                        + beta7_b3 * P2_new;

      // Campionamento effetto random per la nuova unità (variabilità inter-individuale)
      real z_new       = normal_rng(0, 1);
      real Asym_new    = exp(eta_Asym_new + sigma_u * z_new);
      real b2_new      = exp(eta_b2_new);
      real b3_new      = inv_logit(eta_b3_new);

      for (t in 1:n_t) {
        // Media al time-point t per la nuova unità
        real mu_nt = Asym_new * exp(-b2_new * pow(b3_new, t));

        // Deviazione standard eteroschedastica con parametri phi1, phi2 stimati
        real sd_nt = sqrt(phi1 * mu_nt + phi2 * square(mu_nt) + 1e-5);

        // Salva media e campione predittivo
        mu_new[n, t] = mu_nt;
        Y_new[n, t]  = normal_rng(mu_nt, sd_nt);
      }
    }
  }
}
