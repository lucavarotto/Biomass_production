// =============================================================================
// gompertz_re_b3_het_pars
// Modello B — RE su b3, errore eteroschedastico i.i.d., prior parsimoniosa
//
// Modello:
//   OD_ij ~ Normal( mu_ij, sigma^2 * mu_ij^delta )
//   mu_ij = Asym_i * exp( -b2_i * b3_i^t_ij )
//
//   log(Asym_i)  = x_i' * beta_A
//   log(b2_i)    = x_i' * beta_B
//   logit(b3_i)  = x_i' * beta_C + v_i,   v_i ~ Normal(0, sigma_C^2)
//
// Prior semplificata:
//   Intercette: Normal(media_biologica, 0.5)
//   Tutti gli altri beta: Normal(0, 0.25)
//   delta ~ Normal(1, 0.5), delta > 0
// =============================================================================

functions {
  real pred_log_asym(row_vector x, vector beta) { return dot_product(x, beta); }
  real pred_log_b2(row_vector x, vector beta)   { return dot_product(x, beta); }
  real pred_logit_b3(row_vector x, vector beta) { return dot_product(x, beta); }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;

  vector[N] OD;
  vector[N] tempo;

  matrix[K, P] X;

  array[K] int<lower=1, upper=N> start_idx;
  array[K] int<lower=1, upper=N> end_idx;
}

parameters {
  vector[P] beta_A;
  vector[P] beta_B;
  vector[P] beta_C;

  vector[K] v_raw;
  real<lower=0> sigma_C;

  real<lower=0> sigma;
  real<lower=0> delta;
}

transformed parameters {
  vector[K] v = sigma_C * v_raw;
}

model {

  // --- Prior intercette ---
  beta_A[1] ~ normal(2.1, 0.2);   // prior più stretta sull'asintoto
  beta_B[1] ~ normal(0.7, 0.5);
  beta_C[1] ~ normal(0.5, 0.5);

  // --- Prior covariate: unica prior fissa per tutti ---
  beta_A[2:P] ~ normal(0, 0.25);
  beta_B[2:P] ~ normal(0, 0.25);
  beta_C[2:P] ~ normal(0, 0.25);

  // --- Prior RE ---
  v_raw   ~ normal(0, 1);
  sigma_C ~ exponential(1);

  // --- Prior strutturali ---
  sigma ~ exponential(1);
  delta ~ normal(1, 0.5);

  // --- Verosimiglianza eteroschedastica i.i.d. ---
  for (k in 1:K) {
    int n_k         = end_idx[k] - start_idx[k] + 1;
    vector[n_k] t_k = tempo[start_idx[k]:end_idx[k]];
    vector[n_k] y_k = OD[start_idx[k]:end_idx[k]];

    real Asym_k = exp(pred_log_asym(X[k], beta_A));
    real b2_k   = exp(pred_log_b2(X[k],   beta_B));
    real b3_k   = inv_logit(pred_logit_b3(X[k], beta_C) + v[k]);

    vector[n_k] mu_k;
    for (j in 1:n_k)
      mu_k[j] = Asym_k * exp(-b2_k * pow(b3_k, t_k[j]));

    vector[n_k] sd_k;
    for (j in 1:n_k)
      sd_k[j] = sigma * pow(mu_k[j], delta / 2.0);

    y_k ~ normal(mu_k, sd_k);
  }
}

generated quantities {
  vector[N] OD_rep;
  vector[K] log_lik;

  for (k in 1:K) {
    int n_k         = end_idx[k] - start_idx[k] + 1;
    vector[n_k] t_k = tempo[start_idx[k]:end_idx[k]];
    vector[n_k] y_k = OD[start_idx[k]:end_idx[k]];

    real Asym_k = exp(pred_log_asym(X[k], beta_A));
    real b2_k   = exp(pred_log_b2(X[k],   beta_B));
    real b3_k   = inv_logit(pred_logit_b3(X[k], beta_C) + v[k]);

    vector[n_k] mu_k;
    for (j in 1:n_k)
      mu_k[j] = Asym_k * exp(-b2_k * pow(b3_k, t_k[j]));

    vector[n_k] sd_k;
    for (j in 1:n_k)
      sd_k[j] = sigma * pow(mu_k[j], delta / 2.0);

    for (j in 1:n_k)
      OD_rep[start_idx[k] + j - 1] = normal_rng(mu_k[j], sd_k[j]);
    log_lik[k] = normal_lpdf(y_k | mu_k, sd_k);
  }
}
