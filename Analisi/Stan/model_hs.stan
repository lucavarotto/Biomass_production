// =============================================================================
//  Modello nonlineare misto con covarianza AR(1) — Versione Horseshoe
//
//  STRUTTURA DELLA MEDIA:
//    Y_ij ~ MVN(mu_i, V_i)
//    mu_ij = Asym_i * exp(-b2_i * b3_i^j)
//
//  PARAMETRI INDIVIDUALI:
//    Asym_i = exp(eta_Asym_i + u_Asym_i),  u_Asym ~ N(0, sigma_u_Asym^2)
//    b2_i   = exp(eta_b2_i)                 [SOLO EFFETTO FISSO — niente random]
//    b3_i   = inv_logit(eta_b3_i + u_b3_i), u_b3  ~ N(0, sigma_u_b3^2)
//
//  PREDITTORI LINEARI (ora con termini quadratici per I, D, P):
//    eta_*_i = beta0_* + beta1_* I + beta2_* D + beta3_* P
//            + beta4_* I*D                  (interazione originale)
//            + beta5_* I^2 + beta6_* D^2 + beta7_* P^2
//
//  PRIOR — VARIABLE SELECTION:
//    Horseshoe regolarizzata (RHS) di Piironen & Vehtari (2017) per tutti
//    gli effetti fissi (beta1..beta7, intercette a parte) e per le deviazioni
//    standard degli effetti casuali.
//    Parametrizzazione non centrata per z_Asym e z_b3.
//
//  COVARIANZA ERRORE:
//    (V_i)_{jk} = sigma^2 * rho^|j-k|    [AR(1) omoschedastica]
// =============================================================================

// =============================================================================
// BLOCCO DATA
// =============================================================================
data {

  // Dimensioni del dataset
  int<lower=1> N;                     // n. totale osservazioni (= 53*7 - 2 = 369)
  int<lower=1> S;                     // n. soggetti (= 53)
  array[S] int<lower=1> n_i;          // lunghezze serie individuali (6 o 7)
  array[N] int<lower=1> id_biomassa;  // indice soggetto per ogni osservazione
  array[N] int<lower=1> t_idx;        // indice temporale j (1..n_i[s])
  vector[N] Y;                        // risposta (ordinata per soggetto, poi tempo)

  // Covariate soggetto-specifiche (standardizzate prima di passarle al modello)
  vector[S] I_cov;   // Irrigazione
  vector[S] D_cov;   // Densità
  vector[S] P_cov;   // Fertilizzazione (Phosphorus o simile)

  // --- DATI PER LA PREVISIONE FUORI-CAMPIONE ---
  int<lower=0> N_new;       // n. nuove combinazioni (= 27)
  vector[N_new] I_new;
  vector[N_new] D_new;
  vector[N_new] P_new;

  // --- IPERPARAMETRI HORSESHOE (passati come dati per flessibilità) ---
  // slab_scale e slab_df controllano la coda regolarizzata (RHS).
  // Valori tipici: slab_scale = 2, slab_df = 4.
  real<lower=0> slab_scale;   // scala dello "slab" gaussiano
  real<lower=0> slab_df;      // gradi di libertà della distribuzione t dello slab

  // tau0: scala a priori per il parametro di shrinkage globale tau.
  // Scelta: tau0 = (p0/p) / sqrt(N), dove p0 è il n. di predittori
  // "rilevanti" atteso a priori e p è il totale dei predittori.
  // Passa-lo come dato così puoi ricalibrarla senza ricompilare.
  real<lower=0> tau0;
}

// =============================================================================
// BLOCCO TRANSFORMED DATA
// Calcola gli offset start/end per ogni soggetto e i termini quadratici.
// =============================================================================
transformed data {

  // Offset per indicizzare Y per soggetto
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

  // Termini quadratici pre-calcolati per i soggetti osservati
  // (evita ricalcoli ridondanti nel loop del blocco model)
  vector[S] I2_cov = square(I_cov);
  vector[S] D2_cov = square(D_cov);
  vector[S] P2_cov = square(P_cov);

  // Termini quadratici per le nuove unità
  vector[N_new] I2_new = square(I_new);
  vector[N_new] D2_new = square(D_new);
  vector[N_new] P2_new = square(P_new);

  // Numero totale di predittori "covariate" per parametro
  // (escluse le intercette, che hanno prior separate)
  // beta1..beta7 => 7 predittori per ciascuno dei 3 predittori non-lineari
  int p_fx = 7;   // per Asym, b2, b3 ciascuno
  int p_re = 2;   // sigma_u_Asym, sigma_u_b3: 2 parametri di varianza random
}

// =============================================================================
// BLOCCO PARAMETERS
// Parametrizzazione non centrata dove possibile per migliorare il mixing MCMC.
// Parametri horseshoe: tau (shrinkage globale) e lambda (shrinkage locale).
// =============================================================================
parameters {

  // ---------------------------------------------------------------------------
  // INTERCETTE degli effetti fissi (prior separata, informativa debole)
  // Restano fuori dalla horseshoe: l'intercetta è quasi sempre non-zero e
  // includerla nello shrinkage distocerebbe l'intera scala del predittore.
  // ---------------------------------------------------------------------------
  real beta0_Asym;   // intercetta log(Asym)
  real beta0_b2;     // intercetta log(b2)
  real beta0_b3;     // intercetta logit(b3)

  // ---------------------------------------------------------------------------
  // EFFETTI FISSI (slopes) — soggetti a horseshoe
  // Organizzati come vettori per semplificare l'applicazione della prior RHS.
  // Ordine per Asym, b2, b3:
  //   [1] I, [2] D, [3] P, [4] I*D, [5] I^2, [6] D^2, [7] P^2
  // ---------------------------------------------------------------------------
  vector[p_fx] beta_Asym_raw;   // slopes non scalate per Asym
  vector[p_fx] beta_b2_raw;     // slopes non scalate per b2
  vector[p_fx] beta_b3_raw;     // slopes non scalate per b3

  // ---------------------------------------------------------------------------
  // PARAMETRI HORSESHOE — effetti fissi
  // Ogni predittore ha il proprio shrinkage locale lambda_ij^2 > 0.
  // tau è lo shrinkage globale condiviso tra tutti i predittori dello stesso
  // parametro non lineare; qui usiamo un tau comune a tutti e tre (Asym/b2/b3)
  // per semplicità, ma potrebbe essere separato.
  // ---------------------------------------------------------------------------
  real<lower=0> tau_fx;                    // shrinkage globale effetti fissi
  vector<lower=0>[p_fx] lambda_Asym;       // shrinkage locale per beta_Asym
  vector<lower=0>[p_fx] lambda_b2;         // shrinkage locale per beta_b2
  vector<lower=0>[p_fx] lambda_b3;         // shrinkage locale per beta_b3

  // Parametri ausiliari per la regolarizzazione "slab" (distribuzione t):
  // nella RHS, ogni lambda_j^2 viene moltiplicato per c^2 = slab_scale^2 * v/(v+...)
  // Lo implementiamo tramite variabili chi-squared ausiliarie (parametrizzazione
  // di Piironen & Vehtari 2017, eq. 2-3).
  vector<lower=0>[p_fx] caux_Asym;   // variabili ausiliarie chi^2 per lo slab
  vector<lower=0>[p_fx] caux_b2;
  vector<lower=0>[p_fx] caux_b3;

  // ---------------------------------------------------------------------------
  // EFFETTI CASUALI — parametrizzazione non centrata z ~ N(0,1)
  // b2 NON ha più effetto random (rimosso su richiesta).
  // ---------------------------------------------------------------------------
  vector[S] z_Asym;   // componenti raw per u_Asym
  vector[S] z_b3;     // componenti raw per u_b3

  // ---------------------------------------------------------------------------
  // PARAMETRI DI VARIANZA DEGLI EFFETTI CASUALI — soggetti a horseshoe
  // sigma_u_Asym e sigma_u_b3 ricevono una prior half-horseshoe:
  //   sigma_u ~ half-N(0, tau_re * lambda_re)
  // dove tau_re è lo shrinkage globale per le varianze random.
  // ---------------------------------------------------------------------------
  real<lower=0> tau_re;                   // shrinkage globale effetti casuali
  vector<lower=0>[p_re] lambda_re;        // shrinkage locale per le varianze random
  vector<lower=0>[p_re] caux_re;          // ausiliari chi^2 per lo slab random

  // ---------------------------------------------------------------------------
  // PARAMETRI DI ERRORE E CORRELAZIONE
  // ---------------------------------------------------------------------------
  real<lower=0> sigma;          // SD errore residuo
  real<lower=-1, upper=1> rho;  // correlazione AR(1)
}

// =============================================================================
// BLOCCO TRANSFORMED PARAMETERS
// Ricostruisce i parametri interpretabili dalla parametrizzazione non centrata
// e applica la regolarizzazione horseshoe (slab) agli slopes.
// =============================================================================
transformed parameters {

  // ---------------------------------------------------------------------------
  // VARIANZE LOCALI SCALATE (RHS: slab che smorza le code)
  // c_j^2 = slab_scale^2 * (slab_df/2) / caux_j   [da Piironen & Vehtari 2017]
  // lambda_tilde_j^2 = lambda_j^2 * c_j^2 / (c_j^2 + tau^2 * lambda_j^2)
  // Questo "mozzica" i lambda troppo grandi, rendendo la prior propriamente regolare.
  // ---------------------------------------------------------------------------

  // -- Effetti fissi Asym --
  vector<lower=0>[p_fx] c2_Asym = square(slab_scale) * (slab_df / 2.0) ./ caux_Asym;
  vector<lower=0>[p_fx] lam2_tilde_Asym =
      square(lambda_Asym) .* c2_Asym
      ./ (c2_Asym + square(tau_fx) * square(lambda_Asym));
  vector[p_fx] beta_Asym = sqrt(lam2_tilde_Asym) .* beta_Asym_raw * tau_fx;

  // -- Effetti fissi b2 --
  vector<lower=0>[p_fx] c2_b2 = square(slab_scale) * (slab_df / 2.0) ./ caux_b2;
  vector<lower=0>[p_fx] lam2_tilde_b2 =
      square(lambda_b2) .* c2_b2
      ./ (c2_b2 + square(tau_fx) * square(lambda_b2));
  vector[p_fx] beta_b2 = sqrt(lam2_tilde_b2) .* beta_b2_raw * tau_fx;

  // -- Effetti fissi b3 --
  vector<lower=0>[p_fx] c2_b3 = square(slab_scale) * (slab_df / 2.0) ./ caux_b3;
  vector<lower=0>[p_fx] lam2_tilde_b3 =
      square(lambda_b3) .* c2_b3
      ./ (c2_b3 + square(tau_fx) * square(lambda_b3));
  vector[p_fx] beta_b3 = sqrt(lam2_tilde_b3) .* beta_b3_raw * tau_fx;

  // ---------------------------------------------------------------------------
  // VARIANZE DEGLI EFFETTI CASUALI con horseshoe half-normal
  // sigma_u_k = |tau_re * lambda_tilde_re_k|  (half-normal risultante)
  // ---------------------------------------------------------------------------
  vector<lower=0>[p_re] c2_re = square(slab_scale) * (slab_df / 2.0) ./ caux_re;
  vector<lower=0>[p_re] lam2_tilde_re =
      square(lambda_re) .* c2_re
      ./ (c2_re + square(tau_re) * square(lambda_re));
  // sigma_u_Asym = indice 1, sigma_u_b3 = indice 2
  real<lower=0> sigma_u_Asym = tau_re * sqrt(lam2_tilde_re[1]);
  real<lower=0> sigma_u_b3   = tau_re * sqrt(lam2_tilde_re[2]);

  // ---------------------------------------------------------------------------
  // PREDITTORI LINEARI E PARAMETRI INDIVIDUALI
  // eta_*[s] = intercetta + parte lineare + quadratica + interazione
  // ---------------------------------------------------------------------------
  vector[S] eta_Asym;
  vector[S] eta_b2;
  vector[S] eta_b3;

  vector<lower=0>[S]           Asym_i;   // asintoto individuale (> 0)
  vector<lower=0>[S]           b2_i;     // tasso individuale    (> 0, solo fisso)
  vector<lower=0, upper=1>[S]  b3_i;     // persistenza individuale (0,1)

  // Effetti random su scala originale (utili per diagnostica e post-proc.)
  vector[S] u_Asym;
  vector[S] u_b3;

  for (s in 1:S) {
    real id = I_cov[s] * D_cov[s];   // termine di interazione I*D

    // Predittore lineare per Asym (log-scala)
    //   [1]*I + [2]*D + [3]*P + [4]*I*D + [5]*I^2 + [6]*D^2 + [7]*P^2
    eta_Asym[s] = beta0_Asym
                + beta_Asym[1] * I_cov[s]  + beta_Asym[2] * D_cov[s]
                + beta_Asym[3] * P_cov[s]  + beta_Asym[4] * id
                + beta_Asym[5] * I2_cov[s] + beta_Asym[6] * D2_cov[s]
                + beta_Asym[7] * P2_cov[s];

    // Predittore lineare per b2 (log-scala) — SOLO effetto fisso, niente random
    eta_b2[s] = beta0_b2
              + beta_b2[1] * I_cov[s]  + beta_b2[2] * D_cov[s]
              + beta_b2[3] * P_cov[s]  + beta_b2[4] * id
              + beta_b2[5] * I2_cov[s] + beta_b2[6] * D2_cov[s]
              + beta_b2[7] * P2_cov[s];

    // Predittore lineare per b3 (logit-scala)
    eta_b3[s] = beta0_b3
              + beta_b3[1] * I_cov[s]  + beta_b3[2] * D_cov[s]
              + beta_b3[3] * P_cov[s]  + beta_b3[4] * id
              + beta_b3[5] * I2_cov[s] + beta_b3[6] * D2_cov[s]
              + beta_b3[7] * P2_cov[s];

    // Ricostruzione effetti random (non centrata -> centrata)
    u_Asym[s] = sigma_u_Asym * z_Asym[s];
    u_b3[s]   = sigma_u_b3   * z_b3[s];

    // Parametri individuali
    Asym_i[s] = exp(eta_Asym[s] + u_Asym[s]);
    b2_i[s]   = exp(eta_b2[s]);                       // nessun random su b2
    b3_i[s]   = inv_logit(eta_b3[s] + u_b3[s]);
  }
}

// =============================================================================
// BLOCCO MODEL
// Prior + verosimiglianza.
// =============================================================================
model {
  real sigma2 = square(sigma);

  // ---------------------------------------------------------------------------
  // PRIOR — INTERCETTE (weakly informative, fuori horseshoe)
  // Centrate su valori ragionevoli a priori sulla scala di eta.
  // ---------------------------------------------------------------------------
  target += normal_lpdf(beta0_Asym | 2,    0.8);   // log(Asym) ≈ 7.4 kg/m²
  target += normal_lpdf(beta0_b2   | 1.61, 0.8);   // log(b2)   ≈ 5
  target += normal_lpdf(beta0_b3   | 2,    0.4);   // logit(b3) ≈ 0.88

  // ---------------------------------------------------------------------------
  // PRIOR HORSESHOE REGOLARIZZATA — effetti fissi slopes
  //
  // Struttura RHS (Piironen & Vehtari 2017):
  //   beta_j | lambda_j, tau ~ N(0, lambda_tilde_j^2 * tau^2)
  //   lambda_j ~ half-Cauchy(0, 1)
  //   tau      ~ half-Cauchy(0, tau0)   [tau0 passato come dato]
  //   caux_j   ~ Inv-Gamma(slab_df/2, slab_df/2)  =>  chi^2 equivalente
  //
  // Nota: beta_*_raw ~ N(0,1) e beta_* = sqrt(lam2_tilde) .* beta_*_raw * tau
  //       è la parametrizzazione non centrata dell'horseshoe. Questo schema
  //       migliora drammaticamente il mixing rispetto alla centrata quando
  //       molti beta sono vicini a zero (il caso tipico di variable selection).
  // ---------------------------------------------------------------------------

  // Shrinkage globale
  target += cauchy_lpdf(tau_fx | 0, tau0) - cauchy_lccdf(0 | 0, tau0);

  // Shrinkage locale e ausiliari per Asym
  target += cauchy_lpdf(lambda_Asym | 0, 1) - cauchy_lccdf(0 | 0, 1);
  target += inv_gamma_lpdf(caux_Asym | slab_df / 2.0, slab_df / 2.0);
  target += std_normal_lpdf(beta_Asym_raw);   // prior non centrata sui raw

  // Shrinkage locale e ausiliari per b2
  target += cauchy_lpdf(lambda_b2 | 0, 1) - cauchy_lccdf(0 | 0, 1);
  target += inv_gamma_lpdf(caux_b2 | slab_df / 2.0, slab_df / 2.0);
  target += std_normal_lpdf(beta_b2_raw);

  // Shrinkage locale e ausiliari per b3
  target += cauchy_lpdf(lambda_b3 | 0, 1) - cauchy_lccdf(0 | 0, 1);
  target += inv_gamma_lpdf(caux_b3 | slab_df / 2.0, slab_df / 2.0);
  target += std_normal_lpdf(beta_b3_raw);

  // ---------------------------------------------------------------------------
  // PRIOR HORSESHOE REGOLARIZZATA — varianze effetti casuali
  //
  // Stessa struttura RHS applicata a sigma_u_Asym e sigma_u_b3.
  // Questo consente che la varianza random collassi a zero se i dati
  // non supportano eterogeneità inter-soggetto, ottenendo selezione
  // automatica della struttura random.
  // ---------------------------------------------------------------------------
  target += cauchy_lpdf(tau_re | 0, tau0) - cauchy_lccdf(0 | 0, tau0);
  target += cauchy_lpdf(lambda_re | 0, 1) - cauchy_lccdf(0 | 0, 1);
  target += inv_gamma_lpdf(caux_re | slab_df / 2.0, slab_df / 2.0);

  // ---------------------------------------------------------------------------
  // PRIOR — effetti random (parametrizzazione non centrata)
  // z_Asym e z_b3 ~ N(0,1); la scala viene incorporata tramite sigma_u_*.
  // ---------------------------------------------------------------------------
  target += std_normal_lpdf(z_Asym);
  target += std_normal_lpdf(z_b3);

  // ---------------------------------------------------------------------------
  // PRIOR — errore residuo e correlazione AR(1)
  //
  // sigma: half-Cauchy(0, 0.5) — più conservativa di half-N(0, 0.3),
  //        ma con coda più pesante per gestire deviazioni inattese.
  // rho:   Beta(2,2) ricentrata su (-1,1) come nel modello originale.
  //        Penalizza |rho| vicino a 1 (quasi-singolarità della V_i).
  // ---------------------------------------------------------------------------
  target += cauchy_lpdf(sigma | 0, 0.5) - cauchy_lccdf(0 | 0, 0.5);
  target += beta_lpdf((rho + 1) / 2 | 2, 2);

  // ---------------------------------------------------------------------------
  // VEROSIMIGLIANZA: MVN per soggetto con covarianza AR(1)
  // L_s è la fattorizzazione di Cholesky di Sigma_s (più stabile di inv()).
  // ---------------------------------------------------------------------------
  for (s in 1:S) {
    int ni = n_i[s];
    vector[ni] mu_s;
    matrix[ni, ni] Sigma_s;

    // Media non lineare: modello di Chapman-Richards / Gompertz discreto
    for (j in 1:ni)
      mu_s[j] = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));

    // Matrice di covarianza AR(1) con nugget 1e-9 per stabilità Cholesky
    for (j in 1:ni) {
      Sigma_s[j, j] = sigma2 + 1e-9;
      for (k in (j+1):ni) {
        Sigma_s[j, k] = sigma2 * pow(rho, k - j);
        Sigma_s[k, j] = Sigma_s[j, k];   // simmetria esplicita
      }
    }

    matrix[ni, ni] L_s = cholesky_decompose(Sigma_s);
    vector[ni] Y_s = Y[start[s]:end[s]];
    target += multi_normal_cholesky_lpdf(Y_s | mu_s, L_s);
  }
}

// =============================================================================
// BLOCCO GENERATED QUANTITIES
// Calcola: log-verosimiglianza per LOO-CV, repliche PPC, previsioni fuori campione.
// =============================================================================
generated quantities {

  // Log-verosimiglianza per soggetto (usata da loo() in R)
  vector[S] log_lik;

  // Repliche predittive posteriori sui punti osservati (per PPC grafica)
  vector[N] Y_rep;

  { // Blocco locale: le matrici temporanee non vengono salvate nell'output MCMC
    real sigma2 = square(sigma);

    for (s in 1:S) {
      int ni = n_i[s];
      vector[ni] mu_s;
      matrix[ni, ni] Sigma_s;

      // Ricalcolo di mu_s per ogni soggetto
      for (j in 1:ni)
        mu_s[j] = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));

      // Ricalcolo di Sigma_s (identico al blocco model — necessario perché
      // le matrici definite in model non sono visibili qui)
      for (j in 1:ni) {
        Sigma_s[j, j] = sigma2 + 1e-9;
        for (k in (j+1):ni) {
          Sigma_s[j, k] = sigma2 * pow(rho, k - j);
          Sigma_s[k, j] = Sigma_s[j, k];
        }
      }

      matrix[ni, ni] L_s = cholesky_decompose(Sigma_s);
      vector[ni] Y_s = Y[start[s]:end[s]];

      // Log-verosimiglianza soggetto-specifica
      log_lik[s] = multi_normal_cholesky_lpdf(Y_s | mu_s, L_s);

      // Campionamento dalla predittiva posteriore
      Y_rep[start[s]:end[s]] = to_vector(multi_normal_cholesky_rng(mu_s, L_s));
    }
  }

  // ---------------------------------------------------------------------------
  // PREVISIONE FUORI CAMPIONE: N_new nuove combinazioni di covariate
  //
  // Per ogni nuova unità:
  //   1. Calcola eta_* usando gli effetti fissi campionati.
  //   2. Campiona z_new ~ N(0,1) per Asym e b3 per propagare la
  //      variabilità inter-individuale. b2 non ha random: usa solo la media.
  //   3. Genera un intero vettore di 7 time-point dalla MVN AR(1)
  //      (non da normali indipendenti: cattura la dipendenza temporale).
  // ---------------------------------------------------------------------------
  matrix[N_new, 7] mu_new;    // medie predittive per le N_new unità
  matrix[N_new, 7] Y_new;     // campioni predittivi

  {
    real sigma2 = square(sigma);
    int n_t = 7;   // numero fisso di time-point previsti

    for (n in 1:N_new) {
      real id_new = I_new[n] * D_new[n];

      // Predittori lineari di popolazione per la nuova unità
      // (stessa formula del blocco transformed parameters)
      real eta_Asym_new = beta0_Asym
                        + beta_Asym[1] * I_new[n]  + beta_Asym[2] * D_new[n]
                        + beta_Asym[3] * P_new[n]  + beta_Asym[4] * id_new
                        + beta_Asym[5] * I2_new[n] + beta_Asym[6] * D2_new[n]
                        + beta_Asym[7] * P2_new[n];

      real eta_b2_new   = beta0_b2
                        + beta_b2[1] * I_new[n]  + beta_b2[2] * D_new[n]
                        + beta_b2[3] * P_new[n]  + beta_b2[4] * id_new
                        + beta_b2[5] * I2_new[n] + beta_b2[6] * D2_new[n]
                        + beta_b2[7] * P2_new[n];

      real eta_b3_new   = beta0_b3
                        + beta_b3[1] * I_new[n]  + beta_b3[2] * D_new[n]
                        + beta_b3[3] * P_new[n]  + beta_b3[4] * id_new
                        + beta_b3[5] * I2_new[n] + beta_b3[6] * D2_new[n]
                        + beta_b3[7] * P2_new[n];

      // Effetti random campionati per la nuova unità
      // (b2 non ha random: si usa direttamente eta_b2_new)
      real z_Asym_new = normal_rng(0, 1);
      real z_b3_new   = normal_rng(0, 1);

      // Parametri individuali della nuova unità
      real Asym_new_val = exp(eta_Asym_new + sigma_u_Asym * z_Asym_new);
      real b2_new_val   = exp(eta_b2_new);                              // solo fisso
      real b3_new_val   = inv_logit(eta_b3_new + sigma_u_b3 * z_b3_new);

      // Media non lineare per i 7 time-point
      vector[n_t] mu_s;
      for (t in 1:n_t)
        mu_s[t] = Asym_new_val * exp(-b2_new_val * pow(b3_new_val, t));

      // Covarianza AR(1) per i 7 time-point
      matrix[n_t, n_t] Sigma_s;
      for (j in 1:n_t) {
        Sigma_s[j, j] = sigma2 + 1e-9;
        for (k in (j+1):n_t) {
          Sigma_s[j, k] = sigma2 * pow(rho, k - j);
          Sigma_s[k, j] = Sigma_s[j, k];
        }
      }

      // Campionamento MVN (propaga correlazione AR)
      matrix[n_t, n_t] L_s = cholesky_decompose(Sigma_s);
      vector[n_t] Y_s = to_vector(multi_normal_cholesky_rng(mu_s, L_s));

      mu_new[n] = mu_s';
      Y_new[n]  = Y_s';
    }
  }
}
