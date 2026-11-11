// =============================================================
//  Modello nonlineare misto con covarianza AR(1)
//  Y_i ~ N_{n_i}(mu_i(theta), V_i)
//  mu_ij = Asym_i * exp(-b2_i * b3_i^j)
//  Asym_i = eta_i + u_i,   u_i ~ N(0, sigma2_u)
//  b2_i   = eta_i
//  b3_i   = eta_i
//  eta_i  = beta0 + beta1*I_i + beta2*D_i + beta3*P_i + beta4*I_i*D_i
//  (V_i)_{jk} = sigma2 * rho^|j-k|
// =============================================================

// Definisce i dati di input che Stan riceve dall'esterno
data {
  int<lower=1> N;                     // numero totale osservazioni (= 53*7 - 2 = 369)
  int<lower=1> S;                     // numero soggetti (= 53)
  array[S] int<lower=1> n_i;          // n_i: 6 per i=12,15; 7 per gli altri
  array[N] int<lower=1> id_biomassa;  // indice soggetto per ogni osservazione
  array[N] int<lower=1> t_idx;        // indice temporale j (1..n_i) dentro ogni soggetto
  vector[N] Y;                        // variabile risposta (in ordine soggetto, poi tempo)

  // Covariate (una per soggetto)
  vector[S] I_cov;                     // covariata I
  vector[S] D_cov;                     // covariata D
  vector[S] P_cov;                     // covariata P
  
  // --- NUOVI DATI PER LA PREVISIONE ---
  int<lower=0> N_new;       // Sarà 27
  vector[N_new] I_new;
  vector[N_new] D_new;
  vector[N_new] P_new;
}

// viene eseguito una sola volta prima dell'inizio del campionamento MCMC
transformed data {
  // Vettori di lunghezza S che conterranno gli indici di inizio e
  // fine  di ogni soggetto all'interno del vettore globale Y.
  // Soggetto 1: da start[1] a end[1], Soggetto 2: da start[2] a
  // end[2], e così via.
  array[S] int start;
  array[S] int end;
  
  // Blocco locale: serve a usare la variabile temporanea 'pos' 
  // senza salvarla nell'output finale del modello
  {
    int pos = 1; // Punto di partenza nel vettore Y
    for (s in 1:S){
      start[s] = pos;
      end[s]   = pos + n_i[s] - 1;
      pos      = pos + n_i[s];
    }
  }
}

// dichiara i parametri latenti che l'algoritmo MCMC deve stimare.
// Ogni parametro della curva non lineare (Asym, b2, b3) ha i suoi
// 5 $\beta$ dedicati, quindi ci sono 15 $\beta$ in tutto.
// vettore u contiene le deviazioni individuali
// sigma, sigma_u e rho sono le varianze e le covarianza. Hanno
// supporto limitato.
parameters {
  // --- EFFETTI FISSI ----
  // Coefficienti per il predittore Asym (scala log)
  real beta0_Asym; // Intercetta
  real beta1_Asym; // Effetto covariata I
  real beta2_Asym; // Effetto covariata D
  real beta3_Asym; // Effetto covariata P
  real beta4_Asym; // Effetto interazione I * D

  // Coefficienti per il predittore del parametro di scala b2 (scala log)
  real beta0_b2;   // Intercetta
  real beta1_b2;   // Effetto covariata I
  real beta2_b2;   // Effetto covariata D
  real beta3_b2;   // Effetto covariata P
  real beta4_b2;   // Effetto interazione I * D

  // Coefficienti per il predittore del tasso di crescita b3 (scala logit, quindi b3 ∈ (0,1))
  real beta0_b3;   // Intercetta
  real beta1_b3;   // Effetto covariata I
  real beta2_b3;   // Effetto covariata D
  real beta3_b3;   // Effetto covariata P
  real beta4_b3;   // Effetto interazione I * D

  // --- EFFETTI CASUALI (Variazione individuale) ---
  // z[s] ~ N(0,1) sono i "raw effects"; u[s] = sigma_u * z[s]
  // è calcolato in transformed parameters.
  vector[S] z;                         // effetti random standardizzati

  // --- PARAMETRI DI VARIANZA E CORRELAZIONE ---
  real<lower=0> sigma;                 // deviazione standard errore (sigma)
  real<lower=0> sigma_u;              // deviazione standard effetto random
  // rho ∈ (-1,1): supporto già vincolato
  real<lower=-1, upper=1> rho;        // correlazione AR(1)
}

transformed parameters {
  
  // Predittori lineari della popolazione
  // (valori medi attesi in base alle covariate)
  vector[S] eta_Asym;
  vector[S] eta_b2;
  vector[S] eta_b3;

  // Parametri effettivi finali di ogni singolo soggetto
  // !! forzati positivi !!
  // (usati nella funzione non lineare)
  // possono includere una deviazione personale (l'effetto casuale)
  vector<lower=0>[S] Asym_i;
  vector<lower=0>[S] b2_i;
  vector<lower=0, upper=1>[S] b3_i;

  // Effetto random su scala originale (per uso in generated quantities
  // e per interpretazione diretta)
  vector[S] u; // u[s] = sigma_u * z[s]

  for (s in 1:S) { // ciclo sulle S biomasse
    real id = I_cov[s] * D_cov[s];   // interazione I*D

    // 1. Calcolo dei predittori lineari (effetti fissi della popolazione)
    eta_Asym[s] = beta0_Asym + beta1_Asym * I_cov[s] + beta2_Asym * D_cov[s]
                + beta3_Asym * P_cov[s] + beta4_Asym * id;

    eta_b2[s]   = beta0_b2 + beta1_b2 * I_cov[s] + beta2_b2 * D_cov[s]
                + beta3_b2 * P_cov[s] + beta4_b2 * id;

    eta_b3[s]   = beta0_b3 + beta1_b3 * I_cov[s] + beta2_b3 * D_cov[s]
                + beta3_b3 * P_cov[s] + beta4_b3 * id;

    // 2. Ricostruzione effetto random (da non centrato a centrato)
    u[s] = sigma_u * z[s];

    // 3. Parametri individuali
    Asym_i[s]   = exp(eta_Asym[s] + u[s]);   // scala log e effetto random
    b2_i[s]     = exp(eta_b2[s]); // scala log
    b3_i[s]     = inv_logit(eta_b3[s]); // scala logit
  }
}

model {
  real sigma2 = square(sigma); // sigma^2, usato nella costruzione V_i

  // -------------------------------------------------------
  // Prior
  // -------------------------------------------------------
  // Effetti fissi: prior debolmente informativi
  // I coefficienti di pendenza hanno prior N(0, 0.5): effetti moderati
  // sulla scala log/logit, che corrispondono a raddoppi/dimezzamenti
  // del parametro per una variazione di 1 unità della covariata.
  // -------------------------------------------------------
  // normal_lpdf: The log of the normal density of y given
    // location mu and scale sigma

  // su scala esponenziale
  target += normal_lpdf(beta0_Asym | 2, 1);
  target += normal_lpdf(beta1_Asym | 0.5, 0.5);
  target += normal_lpdf(beta2_Asym | 0, 0.5);
  target += normal_lpdf(beta3_Asym | 0, 0.5);
  target += normal_lpdf(beta4_Asym | 0, 0.5);

  // su scala esponenziale
  target += normal_lpdf(beta0_b2 | 1.61, 1.0);
  target += normal_lpdf(beta1_b2 | 0.5, 0.5);
  target += normal_lpdf(beta2_b2 | 0, 0.5);
  target += normal_lpdf(beta3_b2 | 0, 0.5);
  target += normal_lpdf(beta4_b2 | 0, 0.5);

  // su scala logit
  target += normal_lpdf(beta0_b3 | 2, 0.5);
  target += normal_lpdf(beta1_b3 | 0.5, 0.25);
  target += normal_lpdf(beta2_b3 | 0, 0.25);
  target += normal_lpdf(beta3_b3 | 0, 0.25);
  target += normal_lpdf(beta4_b3 | 0, 0.25);

    // -------------------------------------------------------
  // Prior — parametri di varianza
  // Sostituiamo exponential(1) con half-Normal(0, 0.3):
  // mette massa principalmente in [0, 0.6], riducendo la
  // sovra-stima della varianza pur restando debolmente informativa.
  // -------------------------------------------------------
  target += normal_lpdf(sigma       | 0, 0.3) - normal_lccdf(0 | 0, 0.3);
  target += normal_lpdf(sigma_u| 0, 0.3) - normal_lccdf(0 | 0, 0.3);
  
  // -------------------------------------------------------
  // PRIOR — correlazione rho [C3]
  //
  // Problema originale: uniform_lpdf(rho|-1,1) assegna massa uniforme
  // anche vicino a rho=±1, dove la matrice V_i diventa quasi-singolare
  // (autovalore minimo → 0), causando divergenze HMC.
  //
  // Soluzione: Beta(2,2) ricentrata su (-1,1).
  //   Se z = (rho+1)/2 ~ Beta(2,2), allora rho = 2z-1 ~ distribuzione
  //   simmetrica unimodale con moda in 0, che penalizza |rho| vicino a 1.
  //
  // Lo jacobiano della trasformazione z = (rho+1)/2 vale dz/drho = 1/2,
  // quindi: log p(rho) = log p_Beta(z) + log|dz/drho|
  //                     = beta_lpdf(z|2,2) + log(0.5)
  //
  // Stan gestisce automaticamente il jacobiano per il vincolo
  // <lower=-1, upper=1> (mappa interna a reali), quindi occorre
  // aggiungere solo il contributo della Beta via target +=.
  // La costante log(0.5) è ininfluente per il campionamento.
  // -------------------------------------------------------
  target += beta_lpdf((rho + 1) / 2 | 2, 2);

  // PRIOR — effetti random: parametrizzazione non centrata [C2]
  // z[s] ~ N(0,1) indipendenti: prior standard per il termine "raw".
  // L'effetto random effettivo è u[s] = sigma_u * z[s], quindi
  // implicitamente u[s] ~ N(0, sigma_u^2) come richiesto dal modello.
  // -------------------------------------------------------
  target += std_normal_lpdf(z);
  
  // -------------------------------------------------------
  // Verosimiglianza: una MVN per soggetto con covarianza AR(1)
  //
  // Per efficienza e stabilità numerica si usa multi_normal_cholesky_lpdf
  // che accetta direttamente L = chol(Sigma) evitando la decomposizione
  // interna di multi_normal_lpdf (che chiamerebbe cholesky_decompose
  // internamente ad ogni iterazione senza possibilità di controllo).
  //
  // La matrice di Cholesky di una AR(1) sigma^2 * R(rho) può essere
  // costruita direttamente sfruttando la struttura tridiagonale della
  // matrice di precisione, ma per semplicità e leggibilità si costruisce
  // Sigma_s esplicitamente e si chiama cholesky_decompose.
  // Se n_i è piccolo (6 o 7), il costo computazionale è trascurabile.
  //
  // NOTA: multi_normal_cholesky_lpdf vuole L t.c. Sigma = L * L^T.
  // -------------------------------------------------------
  for (s in 1:S) {
    int ni = n_i[s];

    // Costruisci il vettore mu_i e la matrice Sigma_i
    vector[ni] mu_s;
    matrix[ni, ni] Sigma_s;

    // Calcolo del vettore di media non lineare
    for (j in 1:ni) {
      mu_s[j] = Asym_i[s] * exp(-b2_i[s] * pow(b3_i[s], j));
    }

    // Costruzione della matrice di covarianza AR(1)
    // (V_i)_{jk} = sigma^2 * rho^|j-k|
    for (j in 1:ni) {
      // Aggiunta di 1e-9 (nugget) per stabilità numerica
      Sigma_s[j, j] = sigma2 + 1e-9; 
      for (k in (j+1):ni) {
        Sigma_s[j, k] = sigma2 * pow(rho, k - j);
        Sigma_s[k, j] = Sigma_s[j, k]; // Forza la simmetria perfetta
      }
    }

    // Decomposizione di Cholesky: L t.c. Sigma_s = L * L^T [C4]
    // cholesky_decompose fallisce se Sigma_s non è definita positiva;
    // la prior Beta(2,2) su rho riduce questo rischio.
    matrix[ni, ni] L_s = cholesky_decompose(Sigma_s);

    // Log-verosimiglianza MVN tramite fattore di Cholesky
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
  
  // ---- PREVISIONE SULLE NUOVE 27 COMBINAZIONI (TRAIETTORIE TEMPORALI) ----
  // OUTPUT SALVATI IN R: Matrici tridimensionali implicite [Iterazioni x 27 Combinazioni x 7 Settimane]
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

      // Parametri individuali della nuova unità
      real Asym_new_val = exp(eta_Asym_new + sigma_u * z_Asym_new);
      real b2_new_val   = exp(eta_b2_new);
      real b3_new_val   = inv_logit(eta_b3_new);

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
