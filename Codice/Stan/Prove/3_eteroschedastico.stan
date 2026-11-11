// =============================================================
//  Modello nonlineare misto con covarianza ETEROSCHEDASTICA
//  (varianza crescente nel tempo, struttura AR rimossa)
//
//  Modello statistico complessivo:
//    Y_i ~ N_{n_i}(mu_i(theta), Sigma_i)
//
//  Media non lineare (curva di crescita asintotica):
//    mu_ij = Asym_i * exp(-b2_i * b3_i^j)
//
//  Parametri individuali (effetti fissi + random dove previsto):
//    Asym_i = exp(eta_Asym_i + u_Asym_i),  u_Asym ~ N(0, sigma_u_Asym^2)
//    b2_i   = exp(eta_b2_i)                 [NESSUN effetto random su b2]
//    b3_i   = inv_logit(eta_b3_i + u_b3_i), u_b3  ~ N(0, sigma_u_b3^2)
//
//  Predittori lineari (con termini quadratici aggiunti):
//    eta_*_i = beta0_* + beta1_*·I_i + beta2_*·D_i + beta3_*·P_i
//            + beta4_*·I_i·D_i
//            + beta5_*·I_i^2 + beta6_*·D_i^2 + beta7_*·P_i^2
//
//  Struttura di covarianza — ETEROSCHEDASTICA ESPONENZIALE:
//    Sigma_i = D_t * R * D_t
//    dove:
//      R   = I_{n_i}  (matrice identità: residui non correlati)
//      D_t = diag(sigma_1, sigma_2, ..., sigma_{n_i})
//      sigma_t = exp(alpha_v + beta_v * t)   [modello esponenziale]
//
//  SCELTA DEL MODELLO DI VARIANZA — GIUSTIFICAZIONE:
//  Si è scelto il modello ESPONENZIALE sigma_t = exp(alpha_v + beta_v * t)
//  rispetto al modello LINEARE sigma_t = alpha_v + beta_v * t per due
//  ragioni fondamentali:
//    1. GEOMETRIA HMC: i parametri alpha_v e beta_v sono non vincolati
//       (appartengono a R), così il log-posterior è definito su tutto R^2
//       e i gradienti rimangono ben comportati. Il modello lineare richiede
//       vincoli <lower=0> che creano barriere rigide incompatibili con NUTS.
//    2. POSITIVITÀ GARANTITA: exp(.) è intrinsecamente positivo per qualsiasi
//       valore dei parametri, anche durante il warmup quando il campionatore
//       esplora regioni lontane dalla moda posteriore.
//
//  --- COME PASSARE AL MODELLO LINEARE (sigma_t = alpha_v + beta_v * t) ---
//  Se si vuole comunque usare il modello lineare, apportare le seguenti
//  modifiche al codice:
//
//  1. Blocco parameters: sostituire
//       real alpha_v;
//       real beta_v;
//     con:
//       real<lower=0> alpha_v;   // intercetta >= 0 per garantire sigma > 0
//       real<lower=0> beta_v;    // pendenza >= 0 per varianza non decrescente
//     NOTA: il vincolo su beta_v è opzionale se si accetta varianza
//     decrescente nel tempo; quello su alpha_v è invece obbligatorio.
//
//  2. Nel blocco transformed parameters (o nel blocco model), sostituire
//     la riga:
//       sigma_t[t] = exp(alpha_v + beta_v * t);
//     con:
//       sigma_t[t] = alpha_v + beta_v * t;
//
//  3. Prior: i prior half-Normal devono essere rimossi (i parametri sono già
//     vincolati). Sostituire con prior esponenziali o normali troncati, es.:
//       target += exponential_lpdf(alpha_v | 5);   // media a priori ≈ 0.2
//       target += exponential_lpdf(beta_v  | 10);  // piccola crescita
// =============================================================

data {
  int<lower=1> N;                     // Numero totale di osservazioni (= 369)
  int<lower=1> S;                     // Numero di soggetti (= 53)
  array[S] int<lower=1> n_i;          // Lunghezza serie per soggetto (6 o 7)
  array[N] int<lower=1> id_biomassa;  // Indice soggetto per ogni osservazione
  array[N] int<lower=1> t_idx;        // Indice temporale j (1..n_i) per soggetto
  vector[N] Y;                        // Variabile risposta (biomassa osservata)

  // Covariate (una per soggetto): Irraggiamento, Densità, Percentuale
  vector[S] I_cov;
  vector[S] D_cov;
  vector[S] P_cov;

  // --- DATI PER LA PREVISIONE SU NUOVE UNITÀ ---
  int<lower=0> N_new;       // Numero di nuove combinazioni (= 27)
  vector[N_new] I_new;
  vector[N_new] D_new;
  vector[N_new] P_new;
}

// -----------------------------------------------------------------------------
// TRANSFORMED DATA
// Precalcola gli indici di inizio e fine di ogni soggetto nel vettore Y.
// Precalcola anche i termini quadratici e l'interazione per evitare
// ricalcoli ridondanti nei blocchi successivi.
// -----------------------------------------------------------------------------
transformed data {
  // Indici di slicing nel vettore Y per ogni soggetto
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

  // Termini quadratici e di interazione per i soggetti osservati
  // Calcolati una volta sola qui per efficienza computazionale
  vector[S] I2_cov;   // I^2
  vector[S] D2_cov;   // D^2
  vector[S] P2_cov;   // P^2
  vector[S] ID_cov;   // I*D (interazione)
  for (s in 1:S) {
    I2_cov[s] = square(I_cov[s]);
    D2_cov[s] = square(D_cov[s]);
    P2_cov[s] = square(P_cov[s]);
    ID_cov[s] = I_cov[s] * D_cov[s];
  }

  // Termini quadratici e di interazione per le nuove unità (previsione)
  vector[N_new] I2_new;
  vector[N_new] D2_new;
  vector[N_new] P2_new;
  vector[N_new] ID_new;
  for (n in 1:N_new) {
    I2_new[n] = square(I_new[n]);
    D2_new[n] = square(D_new[n]);
    P2_new[n] = square(P_new[n]);
    ID_new[n] = I_new[n] * D_new[n];
  }

  // Numero massimo di time-point (usato nei blocchi generated quantities)
  int T_max = max(n_i);
}

// -----------------------------------------------------------------------------
// PARAMETERS
// Tutti i parametri sono campionati su scala non vincolata dove possibile.
// Le quantità derivate (Asym_i, b2_i, b3_i) sono costruite in
// transformed parameters tramite trasformazioni analitiche (exp, inv_logit).
// -----------------------------------------------------------------------------
parameters {

  // =========================================================
  // EFFETTI FISSI — Asym (asintoto, scala log)
  // Modello lineare esteso con termini quadratici:
  //   eta_Asym = beta0 + beta1*I + beta2*D + beta3*P + beta4*I*D
  //            + beta5*I^2 + beta6*D^2 + beta7*P^2
  // =========================================================
  real beta0_Asym;   // Intercetta (media di log(Asym) al riferimento)
  real beta1_Asym;   // Effetto lineare di I su log(Asym)
  real beta2_Asym;   // Effetto lineare di D su log(Asym)
  real beta3_Asym;   // Effetto lineare di P su log(Asym)
  real beta4_Asym;   // Interazione I*D su log(Asym)
  real beta5_Asym;   // Effetto quadratico di I su log(Asym) [NUOVO]
  real beta6_Asym;   // Effetto quadratico di D su log(Asym) [NUOVO]
  real beta7_Asym;   // Effetto quadratico di P su log(Asym) [NUOVO]

  // =========================================================
  // EFFETTI FISSI — b2 (tasso di decadimento, scala log)
  // NESSUN effetto random: b2 è ora puramente a effetti fissi.
  // La variabilità inter-individuale su b2 viene assorbita dal
  // modello di varianza eteroschedastica e dall'effetto random su Asym.
  // =========================================================
  real beta0_b2;     // Intercetta per log(b2)
  real beta1_b2;     // Effetto lineare di I su log(b2)
  real beta2_b2;     // Effetto lineare di D su log(b2)
  real beta3_b2;     // Effetto lineare di P su log(b2)
  real beta4_b2;     // Interazione I*D su log(b2)
  real beta5_b2;     // Effetto quadratico di I su log(b2) [NUOVO]
  real beta6_b2;     // Effetto quadratico di D su log(b2) [NUOVO]
  real beta7_b2;     // Effetto quadratico di P su log(b2) [NUOVO]

  // =========================================================
  // EFFETTI FISSI — b3 (fattore di scala temporale, scala logit)
  // =========================================================
  real beta0_b3;     // Intercetta per logit(b3)
  real beta1_b3;     // Effetto lineare di I su logit(b3)
  real beta2_b3;     // Effetto lineare di D su logit(b3)
  real beta3_b3;     // Effetto lineare di P su logit(b3)
  real beta4_b3;     // Interazione I*D su logit(b3)
  real beta5_b3;     // Effetto quadratico di I su logit(b3) [NUOVO]
  real beta6_b3;     // Effetto quadratico di D su logit(b3) [NUOVO]
  real beta7_b3;     // Effetto quadratico di P su logit(b3) [NUOVO]

  // =========================================================
  // EFFETTI CASUALI — parametrizzazione non centrata (NCP)
  // z ~ N(0,1), u = sigma_u * z  =>  u ~ N(0, sigma_u^2)
  // La NCP migliora il mixing dell'HMC disaccoppiando la scala
  // dell'effetto random dal suo valore standardizzato.
  // =========================================================
  vector[S] z_Asym;     // Effetti raw non centrati per Asym
  // z_b2 RIMOSSO: b2 non ha più effetto random (vedi nota sopra)
  vector[S] z_b3;       // Effetti raw non centrati per b3

  // =========================================================
  // VARIANZE DI GRUPPO
  // =========================================================
  real<lower=0> sigma_u_Asym;  // SD effetto random su Asym
  // sigma_u_b2 RIMOSSA: b2 ora ha solo effetti fissi
  real<lower=0> sigma_u_b3;    // SD effetto random su b3

  // =========================================================
  // PARAMETRI DELLA VARIANZA ETEROSCHEDASTICA (modello esponenziale)
  //
  //   sigma_t = exp(alpha_v + beta_v * t)
  //
  // alpha_v: intercetta su scala log della SD al tempo t=0
  //          (in pratica, log della SD basale)
  // beta_v:  tasso di crescita logaritmica della SD nel tempo
  //          beta_v > 0 => varianza crescente (atteso per biomasse)
  //          beta_v < 0 => varianza decrescente (non escluso a priori)
  //          beta_v = 0 => omoschedasticità (caso limite)
  //
  // Entrambi i parametri sono NON VINCOLATI (appartengono a R):
  // la positività di sigma_t è garantita dall'esponenziale.
  // =========================================================
  real alpha_v;   // Intercetta log-SD (non vincolata)
  real beta_v;    // Pendenza log-SD nel tempo (non vincolata)
}

// -----------------------------------------------------------------------------
// TRANSFORMED PARAMETERS
// Costruisce le quantità derivate usate nel blocco model e generated quantities.
// Separare le trasformazioni qui permette a Stan di calcolare automaticamente
// il Jacobian aggiustamento dove necessario, e riduce il codice duplicato.
// -----------------------------------------------------------------------------
transformed parameters {

  // Predittori lineari (su scala di trasformazione) per i parametri individuali
  vector[S] eta_Asym;
  vector[S] eta_b2;
  vector[S] eta_b3;

  // Parametri individuali su scala originale
  vector<lower=0>[S] Asym_i;          // Asintoto (> 0)
  vector<lower=0>[S] b2_i;            // Tasso di decadimento (> 0)
  vector<lower=0, upper=1>[S] b3_i;   // Fattore di scala temporale (in (0,1))

  // Effetti random su scala originale (utili per diagnosi e output)
  vector[S] u_Asym;
  // u_b2 RIMOSSO
  vector[S] u_b3;

  // Deviazioni standard per ogni time-point (condivise tra soggetti)
  // Dimensione T_max: la varianza dipende solo dal tempo, non dal soggetto.
  // Per soggetti con n_i < T_max si useranno solo i primi n_i elementi.
  vector[T_max] sigma_t;

  // ------------------------------------------------------------------
  // Calcolo sigma_t per ogni time-point tramite modello esponenziale
  // ------------------------------------------------------------------
  for (t in 1:T_max) {
    // sigma_t = exp(alpha_v + beta_v * t)
    // La varianza al tempo t è quindi sigma_t^2 = exp(2*(alpha_v + beta_v*t))
    sigma_t[t] = exp(alpha_v + beta_v * t);
    // [ALTERNATIVA LINEARE]: sigma_t[t] = alpha_v + beta_v * t;
    //  (richiede alpha_v, beta_v con <lower=0> nel blocco parameters)
  }

  // ------------------------------------------------------------------
  // Calcolo dei parametri individuali per ogni soggetto
  // ------------------------------------------------------------------
  for (s in 1:S) {

    // --- Predittori lineari estesi (con termini quadratici) ---
    eta_Asym[s] = beta0_Asym
                + beta1_Asym * I_cov[s]  + beta2_Asym * D_cov[s]
                + beta3_Asym * P_cov[s]  + beta4_Asym * ID_cov[s]
                + beta5_Asym * I2_cov[s]                           // quadratico I [NUOVO]
                + beta6_Asym * D2_cov[s]                           // quadratico D [NUOVO]
                + beta7_Asym * P2_cov[s];                          // quadratico P [NUOVO]

    eta_b2[s]   = beta0_b2
                + beta1_b2 * I_cov[s]  + beta2_b2 * D_cov[s]
                + beta3_b2 * P_cov[s]  + beta4_b2 * ID_cov[s]
                + beta5_b2 * I2_cov[s]                             // quadratico I [NUOVO]
                + beta6_b2 * D2_cov[s]                             // quadratico D [NUOVO]
                + beta7_b2 * P2_cov[s];                            // quadratico P [NUOVO]

    eta_b3[s]   = beta0_b3
                + beta1_b3 * I_cov[s]  + beta2_b3 * D_cov[s]
                + beta3_b3 * P_cov[s]  + beta4_b3 * ID_cov[s]
                + beta5_b3 * I2_cov[s]                             // quadratico I [NUOVO]
                + beta6_b3 * D2_cov[s]                             // quadratico D [NUOVO]
                + beta7_b3 * P2_cov[s];                            // quadratico P [NUOVO]

    // --- Ricostruzione effetti random (da NCP a scala originale) ---
    u_Asym[s] = sigma_u_Asym * z_Asym[s];
    // u_b2 RIMOSSO
    u_b3[s]   = sigma_u_b3   * z_b3[s];

    // --- Parametri individuali su scala originale ---
    // b2 NON ha più effetto random: dipende solo dagli effetti fissi eta_b2
    Asym_i[s] = exp(eta_Asym[s] + u_Asym[s]);          // random su Asym
    b2_i[s]   = exp(eta_b2[s]);                         // solo fissi su b2
    b3_i[s]   = inv_logit(eta_b3[s] + u_b3[s]);        // random su b3
  }
}

// -----------------------------------------------------------------------------
// MODEL
// Specifica la log-densità posteriore (prior + verosimiglianza).
// Stan massimizza target += log p(...) durante il sampling.
// -----------------------------------------------------------------------------
model {

  // ===========================================================
  // PRIOR — effetti fissi di Asym (scala log)
  // Centrati attorno a valori biologicamente ragionevoli.
  // I termini quadratici ricevono prior più stretti (0.2) per
  // favore la parsimonia e prevenire overfitting sui dati limitati.
  // ===========================================================
  target += normal_lpdf(beta0_Asym | 2,   0.8);   // log(Asym) base ≈ 7.4 g
  target += normal_lpdf(beta1_Asym | 0.5, 0.4);   // effetto lineare I
  target += normal_lpdf(beta2_Asym | 0,   0.4);   // effetto lineare D
  target += normal_lpdf(beta3_Asym | 0,   0.4);   // effetto lineare P
  target += normal_lpdf(beta4_Asym | 0,   0.4);   // interazione I*D
  target += normal_lpdf(beta5_Asym | 0,   0.2);   // quadratico I [NUOVO]
  target += normal_lpdf(beta6_Asym | 0,   0.2);   // quadratico D [NUOVO]
  target += normal_lpdf(beta7_Asym | 0,   0.2);   // quadratico P [NUOVO]

  // ===========================================================
  // PRIOR — effetti fissi di b2 (scala log)
  // ===========================================================
  target += normal_lpdf(beta0_b2 | 1.61, 0.8);   // log(b2) base ≈ 5 (Gompertz)
  target += normal_lpdf(beta1_b2 | 0.5,  0.4);
  target += normal_lpdf(beta2_b2 | 0,    0.4);
  target += normal_lpdf(beta3_b2 | 0,    0.4);
  target += normal_lpdf(beta4_b2 | 0,    0.4);
  target += normal_lpdf(beta5_b2 | 0,    0.2);   // quadratico I [NUOVO]
  target += normal_lpdf(beta6_b2 | 0,    0.2);   // quadratico D [NUOVO]
  target += normal_lpdf(beta7_b2 | 0,    0.2);   // quadratico P [NUOVO]

  // ===========================================================
  // PRIOR — effetti fissi di b3 (scala logit)
  // b3 ∈ (0,1): logit(b3) = 2 => b3 ≈ 0.88, sensato per crescita
  // ===========================================================
  target += normal_lpdf(beta0_b3 | 2,   0.4);
  target += normal_lpdf(beta1_b3 | 0.5, 0.2);
  target += normal_lpdf(beta2_b3 | 0,   0.2);
  target += normal_lpdf(beta3_b3 | 0,   0.2);
  target += normal_lpdf(beta4_b3 | 0,   0.2);
  target += normal_lpdf(beta5_b3 | 0,   0.1);   // quadratico I [NUOVO]
  target += normal_lpdf(beta6_b3 | 0,   0.1);   // quadratico D [NUOVO]
  target += normal_lpdf(beta7_b3 | 0,   0.1);   // quadratico P [NUOVO]

  // ===========================================================
  // PRIOR — varianze di gruppo (half-Normal)
  // half-Normal(0, 0.3) mette la maggior parte della massa in
  // [0, 0.6], scoraggiando iper-varianza ma restando debolmente
  // informativo rispetto ai dati.
  // ===========================================================
  target += normal_lpdf(sigma_u_Asym | 0, 0.3) - normal_lccdf(0 | 0, 0.3);
  // sigma_u_b2 RIMOSSA
  target += normal_lpdf(sigma_u_b3   | 0, 0.3) - normal_lccdf(0 | 0, 0.3);

  // ===========================================================
  // PRIOR — parametri della varianza eteroschedastica esponenziale
  //
  //   sigma_t = exp(alpha_v + beta_v * t)
  //
  // alpha_v ~ Normal(-1, 1):
  //   exp(-1) ≈ 0.37 come SD basale (al t=0), debolmente informativo.
  //   Esclude a priori SD basali > exp(2) ≈ 7.4 (biologicamente implausibili).
  //
  // beta_v ~ Normal(0, 0.3):
  //   Centrato su zero: a priori non si assume crescita né decrescita.
  //   SD = 0.3: con 7 time-point, beta_v = 0.3 implica che sigma al t=7
  //   è exp(0.3*7) ≈ 8x la SD basale — prior già generoso.
  //   Valori |beta_v| > 1 sono fortemente penalizzati (crescita esplosiva).
  // ===========================================================
  target += normal_lpdf(alpha_v | -1, 1);    // Intercetta log-SD
  target += normal_lpdf(beta_v  |  0, 0.3);  // Pendenza log-SD nel tempo

  // ===========================================================
  // PRIOR — effetti random standardizzati (NCP)
  // z ~ N(0,1) per costruzione della parametrizzazione non centrata
  // ===========================================================
  target += std_normal_lpdf(z_Asym);
  // z_b2 RIMOSSO
  target += std_normal_lpdf(z_b3);

  // ===========================================================
  // VEROSIMIGLIANZA
  //
  // Per ogni soggetto s, Y_s ~ MVN(mu_s, Sigma_s) dove:
  //   Sigma_s = D_s * I * D_s = diag(sigma_t^2)
  //
  // Poiché R = I (matrice identità), la MVN si fattorizza in
  // n_i normali indipendenti eteroschedastiche:
  //   Y_{sj} ~ Normal(mu_{sj}, sigma_{t_j})  per j = 1..n_i
  //
  // Questo è EQUIVALENTE alla forma MVN con covarianza diagonale,
  // ma computazionalmente molto più efficiente: O(n_i) invece di O(n_i^3).
  //
  // La matrice Sigma_s = D_s * R * D_s con R = I si riduce esattamente
  // a diag(sigma_1^2, sigma_2^2, ..., sigma_{n_i}^2).
  // ===========================================================
  for (s in 1:S) {
    int ni = n_i[s];

    for (j in 1:ni) {
      // Media non lineare al time-point j per il soggetto s
      real mu_sj = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));

      // Deviazione standard al time-point j (eteroschedastica, condivisa tra soggetti)
      // sigma_t[j] è già calcolato in transformed parameters
      // Il jitter 1e-9 è aggiunto qui per sicurezza numerica (previene sd=0
      // in casi estremi del campionamento durante il warmup)
      real sd_j = sigma_t[j] + 1e-9;

      // Log-verosimiglianza contributo scalare (efficiente: evita Cholesky)
      target += normal_lpdf(Y[start[s] + j - 1] | mu_sj, sd_j);
    }
  }
}

// -----------------------------------------------------------------------------
// GENERATED QUANTITIES
// Calcola diagnostiche (log_lik per LOO-CV), repliche predittive (PPC),
// e previsioni su nuove unità.
// -----------------------------------------------------------------------------
generated quantities {

  // Log-verosimiglianza aggregata per soggetto (vettore usato da loo())
  vector[S] log_lik;

  // Repliche predittive posteriori sui dati osservati (per PPC)
  vector[N] Y_rep;

  // ------------------------------------------------------------------
  // LOOP 1: Diagnostica e repliche sul dataset storico (53 biomasse)
  // ------------------------------------------------------------------
  {
    for (s in 1:S) {
      int ni = n_i[s];
      real ll_s = 0;   // Accumulatore log-verosimiglianza per soggetto s

      for (j in 1:ni) {
        // Ricalcolo della media (identico al blocco model)
        real mu_sj = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));
        real sd_j  = sigma_t[j] + 1e-9;

        // Contributo scalare alla log-verosimiglianza
        ll_s += normal_lpdf(Y[start[s] + j - 1] | mu_sj, sd_j);

        // Campionamento dalla predittiva a posteriori per PPC
        Y_rep[start[s] + j - 1] = normal_rng(mu_sj, sd_j);
      }

      // Salva la log-verosimiglianza aggregata del soggetto
      // (usata da loo::loo() in R per il PSIS-LOO)
      log_lik[s] = ll_s;
    }
  }

  // ------------------------------------------------------------------
  // LOOP 2: Previsione su N_new nuove combinazioni di covariate
  //
  // Per ogni nuova unità:
  //   1. Calcola i predittori lineari di popolazione (effetti fissi).
  //   2. Campiona gli effetti random per propagare la variabilità
  //      inter-individuale: z_new ~ N(0,1).
  //   3. Costruisce i parametri individuali della nuova unità.
  //   4. Genera Y_new dal modello eteroschedastico.
  //
  // NOTA: b2 non ha effetto random, quindi NON si campiona z_b2_new.
  // ------------------------------------------------------------------
  matrix[N_new, 7] mu_new;    // Medie predette per ogni nuova unità x 7 time-point
  matrix[N_new, 7] Y_new;     // Campioni predittivi per ogni nuova unità x 7 time-point

  {
    int n_t = 7;   // Numero di time-point previsti (fisso per le previsioni)

    for (n in 1:N_new) {

      // --- Predittori lineari di popolazione per la nuova unità ---
      // (termini lineari + interazione + quadratici)
      real eta_Asym_new = beta0_Asym
                        + beta1_Asym * I_new[n]  + beta2_Asym * D_new[n]
                        + beta3_Asym * P_new[n]  + beta4_Asym * ID_new[n]
                        + beta5_Asym * I2_new[n]
                        + beta6_Asym * D2_new[n]
                        + beta7_Asym * P2_new[n];

      real eta_b2_new   = beta0_b2
                        + beta1_b2 * I_new[n]  + beta2_b2 * D_new[n]
                        + beta3_b2 * P_new[n]  + beta4_b2 * ID_new[n]
                        + beta5_b2 * I2_new[n]
                        + beta6_b2 * D2_new[n]
                        + beta7_b2 * P2_new[n];

      real eta_b3_new   = beta0_b3
                        + beta1_b3 * I_new[n]  + beta2_b3 * D_new[n]
                        + beta3_b3 * P_new[n]  + beta4_b3 * ID_new[n]
                        + beta5_b3 * I2_new[n]
                        + beta6_b3 * D2_new[n]
                        + beta7_b3 * P2_new[n];

      // --- Campionamento effetti random per la nuova unità ---
      // Propaga la variabilità inter-individuale nella predizione.
      // b2 NON ha effetto random: nessun campionamento di z_b2_new.
      real z_Asym_new = normal_rng(0, 1);
      real z_b3_new   = normal_rng(0, 1);

      // --- Parametri individuali della nuova unità ---
      real Asym_new_val = exp(eta_Asym_new + sigma_u_Asym * z_Asym_new);
      real b2_new_val   = exp(eta_b2_new);                                   // solo fissi
      real b3_new_val   = inv_logit(eta_b3_new + sigma_u_b3 * z_b3_new);

      // --- Generazione medie e campioni predittivi per i 7 time-point ---
      for (t in 1:n_t) {
        // Media non lineare
        real mu_t = Asym_new_val * exp(-b2_new_val * pow(b3_new_val, t));

        // SD eteroschedastica al time-point t (usa sigma_t già calcolato)
        real sd_t = sigma_t[t] + 1e-9;

        // Campionamento dalla predittiva (propagazione incertezza)
        mu_new[n, t] = mu_t;
        Y_new[n, t]  = normal_rng(mu_t, sd_t);
      }
    }
  }
}
