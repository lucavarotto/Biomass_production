// =============================================================
//  Modello nonlineare misto con covarianza AR(1)
//  Y_i ~ N_{n_i}(mu_i(theta), V_i)
//  mu_ij = Asym_i * exp(-b2_i * b3_i^j)
//  Asym_i = exp(eta_Asym_i + u_Asym_i),  u_Asym ~ N(0, sigma_u_Asym^2)
//  b2_i   = exp(eta_b2_i   + u_b2_i),    u_b2   ~ N(0, sigma_u_b2^2)   [NUOVO]
//  b3_i   = inv_logit(eta_b3_i + u_b3_i), u_b3  ~ N(0, sigma_u_b3^2)   [NUOVO]
//  eta_*_i = beta0_* + beta1_*·I_i + beta2_*·D_i + beta3_*·P_i + beta4_*·I_i·D_i
//  (V_i)_{jk} = sigma^2 * rho^|j-k|
// =============================================================

data {
  int<lower=1> N;                     // numero totale osservazioni (= 53*7 - 2 = 369)
  int<lower=1> S;                     // numero soggetti (= 53)
  array[S] int<lower=1> n_i;          // n_i: 6 per i=12,15; 7 per gli altri
  array[N] int<lower=1> id_biomassa;  // indice soggetto per ogni osservazione
  array[N] int<lower=1> t_idx;        // indice temporale j (1..n_i) dentro ogni soggetto
  vector[N] Y;                        // variabile risposta (in ordine soggetto, poi tempo)

  // Covariate (una per soggetto)
  vector[S] I_cov;
  vector[S] D_cov;
  vector[S] P_cov;

  // --- DATI PER LA PREVISIONE ---
  int<lower=0> N_new;       // Sarà 27
  vector[N_new] I_new;
  vector[N_new] D_new;
  vector[N_new] P_new;
}

transformed data {
  array[S] int start;
  array[S] int end;
  {
    int pos = 1;
    for (s in 1:S){
      start[s] = pos;
      end[s]   = pos + n_i[s] - 1;
      pos      = pos + n_i[s];
    }
  }
}

parameters {
  // --- EFFETTI FISSI ---
  real beta0_Asym; real beta1_Asym; real beta2_Asym; real beta3_Asym; real beta4_Asym;
  real beta0_b2;   real beta1_b2;   real beta2_b2;   real beta3_b2;   real beta4_b2;
  real beta0_b3;   real beta1_b3;   real beta2_b3;   real beta3_b3;   real beta4_b3;

  // --- EFFETTI CASUALI (parametrizzazione non centrata) ---
  vector[S] z_Asym;   // effetti raw per Asym
  vector[S] z_b2;     // effetti raw per b2  [NUOVO]
  vector[S] z_b3;     // effetti raw per b3  [NUOVO]

  // --- PARAMETRI DI VARIANZA E CORRELAZIONE ---
  real<lower=0> sigma;         // SD errore residuo
  real<lower=0> sigma_u_Asym;  // SD effetto random su Asym
  real<lower=0> sigma_u_b2;    // SD effetto random su b2  [NUOVO]
  real<lower=0> sigma_u_b3;    // SD effetto random su b3  [NUOVO]
  real<lower=-1, upper=1> rho; // correlazione AR(1)
}

transformed parameters {
  vector[S] eta_Asym;
  vector[S] eta_b2;
  vector[S] eta_b3;

  vector<lower=0>[S] Asym_i;
  vector<lower=0>[S] b2_i;
  vector<lower=0, upper=1>[S] b3_i;

  // Effetti random su scala originale (per interpretazione e diagnostica)
  vector[S] u_Asym;
  vector[S] u_b2;    // [NUOVO]
  vector[S] u_b3;    // [NUOVO]

  for (s in 1:S) {
    real id = I_cov[s] * D_cov[s];

    eta_Asym[s] = beta0_Asym + beta1_Asym * I_cov[s] + beta2_Asym * D_cov[s]
                + beta3_Asym * P_cov[s] + beta4_Asym * id;
    eta_b2[s]   = beta0_b2   + beta1_b2   * I_cov[s] + beta2_b2   * D_cov[s]
                + beta3_b2   * P_cov[s] + beta4_b2   * id;
    eta_b3[s]   = beta0_b3   + beta1_b3   * I_cov[s] + beta2_b3   * D_cov[s]
                + beta3_b3   * P_cov[s] + beta4_b3   * id;

    // Ricostruzione effetti random (da non centrata a centrata)
    u_Asym[s] = sigma_u_Asym * z_Asym[s];
    u_b2[s]   = sigma_u_b2   * z_b2[s];   // [NUOVO]
    u_b3[s]   = sigma_u_b3   * z_b3[s];   // [NUOVO]

    // Parametri individuali con effetto random su tutti e tre
    Asym_i[s] = exp(eta_Asym[s] + u_Asym[s]);
    b2_i[s]   = exp(eta_b2[s]   + u_b2[s]);    // [NUOVO]
    b3_i[s]   = inv_logit(eta_b3[s] + u_b3[s]); // [NUOVO]
  }
}

model {
  real sigma2 = square(sigma);

  // -------------------------------------------------------
  // Prior — effetti fissi
  // -------------------------------------------------------
  // Asym (scala log)
  target += normal_lpdf(beta0_Asym | 2,   0.8);
  target += normal_lpdf(beta1_Asym | 0.5, 0.4);
  target += normal_lpdf(beta2_Asym | 0,   0.4);
  target += normal_lpdf(beta3_Asym | 0,   0.4);
  target += normal_lpdf(beta4_Asym | 0,   0.4);

  // b2 (scala log)
  target += normal_lpdf(beta0_b2 | 1.61, 0.8);
  target += normal_lpdf(beta1_b2 | 0.5,  0.4);
  target += normal_lpdf(beta2_b2 | 0,    0.4);
  target += normal_lpdf(beta3_b2 | 0,    0.4);
  target += normal_lpdf(beta4_b2 | 0,    0.4);

  // b3 (scala logit)
  target += normal_lpdf(beta0_b3 | 2,   0.4);
  target += normal_lpdf(beta1_b3 | 0.5, 0.2);
  target += normal_lpdf(beta2_b3 | 0,   0.2);
  target += normal_lpdf(beta3_b3 | 0,   0.2);
  target += normal_lpdf(beta4_b3 | 0,   0.2);

  // -------------------------------------------------------
  // Prior — parametri di varianza
  // Sostituiamo exponential(1) con half-Normal(0, 0.3):
  // mette massa principalmente in [0, 0.6], riducendo la
  // sovra-stima della varianza pur restando debolmente informativa.
  // -------------------------------------------------------
  target += normal_lpdf(sigma       | 0, 0.3) - normal_lccdf(0 | 0, 0.3);
  target += normal_lpdf(sigma_u_Asym| 0, 0.3) - normal_lccdf(0 | 0, 0.3);
  target += normal_lpdf(sigma_u_b2  | 0, 0.3) - normal_lccdf(0 | 0, 0.3); // [NUOVO]
  target += normal_lpdf(sigma_u_b3  | 0, 0.3) - normal_lccdf(0 | 0, 0.3); // [NUOVO]

  // -------------------------------------------------------
  // Prior — correlazione rho: Beta(2,2) ricentrata su (-1,1)
  // Penalizza |rho| vicino a 1 per evitare quasi-singolarità di V_i.
  // -------------------------------------------------------
  target += beta_lpdf((rho + 1) / 2 | 2, 2);

  // -------------------------------------------------------
  // Prior — effetti random (parametrizzazione non centrata)
  // -------------------------------------------------------
  target += std_normal_lpdf(z_Asym);
  target += std_normal_lpdf(z_b2);   // [NUOVO]
  target += std_normal_lpdf(z_b3);   // [NUOVO]

  // -------------------------------------------------------
  // Verosimiglianza: MVN per soggetto con covarianza AR(1)
  // -------------------------------------------------------
  for (s in 1:S) {
    int ni = n_i[s];
    vector[ni] mu_s;
    matrix[ni, ni] Sigma_s;

    for (j in 1:ni)
      mu_s[j] = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));

    for (j in 1:ni) {
      Sigma_s[j, j] = sigma2 + 1e-9;
      for (k in (j+1):ni) {
        Sigma_s[j, k] = sigma2 * pow(rho, k - j);
        Sigma_s[k, j] = Sigma_s[j, k];
      }
    }

    matrix[ni, ni] L_s = cholesky_decompose(Sigma_s);
    vector[ni] Y_s = Y[start[s]:end[s]];
    target += multi_normal_cholesky_lpdf(Y_s | mu_s, L_s);
  }
}

// -----------------------------------------------------------------------------
// BLOCCO GENERATED QUANTITIES
// Eseguito rigidamente una volta per iterazione di sampling (escluso il warm-up).
// Calcola le metriche di validazione (LOO-CV), le repliche sui dati storici (PPC)
// e proietta le traiettorie predittive sulle 27 combinazioni di interesse.
// -----------------------------------------------------------------------------
generated quantities {
  
// OUTPUT SALVATI IN R: Strutturati come [Iterazioni MCMC x Dimensione Variabile]
  vector[S] log_lik;   // Log-verosimiglianza aggregata per SOGGETTO (vettore multivariato)
  vector[N] Y_rep;     // Repliche predittive posteriori sui singoli punti-osservazione storici

  { // Blocco locale {}: isola le matrici temporanee per non intasare l'output di R
    real sigma2 = square(sigma);

    // LOOP 1: Diagnostica e Repliche sul dataset storico (53 Biomasse)
    for (s in 1:S) {
      int ni = n_i[s];
      vector[ni] mu_s;
      matrix[ni, ni] Sigma_s;

      // Ricalcolo di mu_s e Sigma_s (identico al blocco model)
      for (j in 1:ni) {
        mu_s[j] = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));
      }

    // Costruzione della matrice di covarianza teorica AR(1) per il soggetto s
    for (j in 1:ni) {
      // Aggiunta di 1e-9 (nugget) per stabilità numerica
      Sigma_s[j, j] = sigma2 + 1e-9; 
      for (k in (j+1):ni) {
        Sigma_s[j, k] = sigma2 * pow(rho, k - j);
        Sigma_s[k, j] = Sigma_s[j, k]; // Forza la simmetria perfetta
      }
    }

    // Usa nuovamente Cholesky per coerenza con il blocco model
    matrix[ni, ni] L_s = cholesky_decompose(Sigma_s);

    vector[ni] Y_s = Y[start[s]:end[s]];

    // Log-verosimiglianza soggetto-specifica
    log_lik[s] = multi_normal_cholesky_lpdf(Y_s | mu_s, L_s);

    // Campione dalla predittiva a posteriori per ogni soggetto
    // multi_normal_cholesky_rng genera un campione da N(mu_s, L_s L_s^T)
    vector[ni] Y_rep_s = to_vector(multi_normal_cholesky_rng(mu_s, L_s));
    Y_rep[start[s]:end[s]] = Y_rep_s;
    }
  }

  // ---- PREVISIONE SULLE NUOVE N_new COMBINAZIONI ----
  //
  // Per una nuova unità non osservata occorre:
  //   1. Campionare l'effetto random z_new ~ N(0,1) per ognuno dei
  //      tre parametri (Asym, b2, b3), così da propagare la
  //      variabilità inter-individuale nella predizione.
  //   2. Costruire mu_new con i parametri individuali campionati.
  //   3. Generare Y_new dalla MVN con covarianza AR(1) completa,
  //      non da normali indipendenti (corregge il bug originale).
  //
  matrix[N_new, 7] mu_new;
  matrix[N_new, 7] Y_new;

  {
    real sigma2 = square(sigma);
    int n_t = 7; // numero di time-point previsti

    for (n in 1:N_new) {
      real id_new = I_new[n] * D_new[n];

      // Predittori lineari di popolazione per la nuova unità
      real eta_Asym_new = beta0_Asym + beta1_Asym * I_new[n] + beta2_Asym * D_new[n]
                        + beta3_Asym * P_new[n] + beta4_Asym * id_new;
      real eta_b2_new   = beta0_b2   + beta1_b2   * I_new[n] + beta2_b2   * D_new[n]
                        + beta3_b2   * P_new[n] + beta4_b2   * id_new;
      real eta_b3_new   = beta0_b3   + beta1_b3   * I_new[n] + beta2_b3   * D_new[n]
                        + beta3_b3   * P_new[n] + beta4_b3   * id_new;

      // Campionamento effetto random per la nuova unità (variabilità inter-individuale)
      real z_Asym_new = normal_rng(0, 1);
      real z_b2_new   = normal_rng(0, 1);   // [NUOVO]
      real z_b3_new   = normal_rng(0, 1);   // [NUOVO]

      // Parametri individuali della nuova unità
      real Asym_new_val = exp(eta_Asym_new + sigma_u_Asym * z_Asym_new);
      real b2_new_val   = exp(eta_b2_new   + sigma_u_b2   * z_b2_new);   // [NUOVO]
      real b3_new_val   = inv_logit(eta_b3_new + sigma_u_b3 * z_b3_new); // [NUOVO]

      // Vettore di media non lineare per i 7 time-point
      vector[n_t] mu_s;
      for (t in 1:n_t)
        mu_s[t] = Asym_new_val * exp(-b2_new_val * pow(b3_new_val, t));

      // Matrice di covarianza AR(1) per i 7 time-point
      matrix[n_t, n_t] Sigma_s;
      for (j in 1:n_t) {
        Sigma_s[j, j] = sigma2 + 1e-9;
        for (k in (j+1):n_t) {
          Sigma_s[j, k] = sigma2 * pow(rho, k - j);
          Sigma_s[k, j] = Sigma_s[j, k];
        }
      }

      // Cholesky e campionamento MVN (corregge la normale indipendente originale)
      matrix[n_t, n_t] L_s = cholesky_decompose(Sigma_s);
      vector[n_t] Y_s = to_vector(multi_normal_cholesky_rng(mu_s, L_s));

      // Salva media e predizione per questa nuova unità
      mu_new[n] = mu_s';
      Y_new[n]  = Y_s';
    }
  }
}
