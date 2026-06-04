// =============================================================================
// Nonlinear Ornstein-Uhlenbeck model with stochastic volatility and Student-t
// innovations, hierarchical (partial-pooling) priors, and a TMG interaction.
//
// SINGLE SOURCE OF TRUTH: this file is the canonical model definition. The R
// function ou_nonlinear_tmg_stan_code() reads it; fit_ou_nonlinear_tmg()
// compiles it directly. Do NOT keep a second copy of the model as an R string.
//
// Discretization: Euler-Maruyama with dt = 1. The discrete AR coefficient on the
// deviation z = Y - theta is (1 - kappa) for the linear part; interpret the
// half-life accordingly (see R documentation), not as log(2)/kappa.
//
// Likelihood window: the likelihood and log_lik are summed over t in 2..T_lik.
//   - T_lik = T_train  -> honest forecasting design (test block is held out).
//   - T_lik = T        -> full-information design (LOO valid over all t).
// This removes the train/test leakage of fitting on the full sample while
// labelling part of it "out-of-sample".
// =============================================================================

functions {
  // Partial likelihood for reduce_sum: accumulates the one-step-ahead Student-t
  // log density over the time indices in t_idx_slice that satisfy 1 < t <= T_lik.
  real ou_nl_partial_sum(array[] int t_idx_slice,
                         int start, int end,
                         int    T_lik,
                         matrix Yz,
                         matrix Xz,
                         matrix COM_ts,
                         vector zTMG_byK,
                         int    soft_wedge,
                         vector delta_z,
                         vector com_wmean_train,
                         vector com_wsd_train,
                         vector mu_xz,
                         vector theta_s,
                         vector kappa_s,
                         vector a3_s,
                         vector beta0_s,
                         real   beta1,
                         matrix h,
                         real   nu,
                         int    com_in_mean,
                         real   gamma) {
    real lp = 0;
    int S = cols(Yz);
    for (t_idx in t_idx_slice) {
      int t = t_idx;
      if (t <= 1) continue;
      if (t > T_lik) continue;
      real ztmg_eff = zTMG_byK[t];
      if (soft_wedge == 1) ztmg_eff += delta_z[t];
      ztmg_eff = fmin(fmax(ztmg_eff, -1e6), 1e6);
      for (s in 1:S) {
        real zlag   = Yz[t-1,s] - theta_s[s];
        real drift  = kappa_s[s] * (theta_s[s] - Yz[t-1,s] + a3_s[s] * zlag^3);
        real betaT  = beta0_s[s] + beta1 * ztmg_eff;

        real denom_sd = com_wsd_train[s];
        denom_sd = (denom_sd > 1e-12) ? denom_sd : 1.0;
        real com_std  = (COM_ts[t-1,s] - com_wmean_train[s]) / denom_sd;
        com_std = fmin(fmax(com_std, -1e6), 1e6);
        real com_term = (com_in_mean == 1) ? gamma * com_std : 0;

        real sd_safe = fmin(fmax(exp(0.5 * h[t,s]), 1e-8), 1e8);
        real mean_   = drift + betaT * (Xz[t-1,s] - mu_xz[s]) + com_term;
        real y_      = Yz[t,s] - Yz[t-1,s] - mean_;

        lp += student_t_lpdf(y_ | nu, 0, sd_safe);
      }
    }
    return lp;
  }
}

data {
  int<lower=2> T;
  int<lower=1> S;
  int<lower=2> T_train;
  int<lower=2, upper=T> T_lik;          // last time index included in likelihood
  matrix[T,S] Yz;
  matrix[T,S] Xz;
  vector[T] zTMG_byK;
  vector[T] zTMG_exo;
  int<lower=0,upper=1> soft_wedge;
  real<lower=0> sigma_delta_z;
  matrix[T,S] COM_ts;
  matrix[T,S] K_ts;
  int<lower=0,upper=1> com_in_mean;
  vector[S] mu_xz;

  // Configurable priors (defaults set on the R side; see fit_ou_nonlinear_tmg).
  real beta1_prior_mean;                // neutral default 0 (no sign baked in)
  real<lower=0> beta1_prior_sd;
  real<lower=0> nu_prior_shape;         // gamma() on nu_tilde; default (2, 0.1)
  real<lower=0> nu_prior_rate;          //   -> weakly-informative, mean nu ~ 22
  real rho_prior_mean;                  // SV persistence prior; default 0.7
  real<lower=0> rho_prior_sd;           //   sd default 0.2 (not the rigid 0.05)
}

transformed data {
  vector[S] com_wmean_train;
  vector[S] com_wsd_train;
  vector[S] COM_s;

  // Capital-weighted COM mean and sd in TRAIN
  for (s in 1:S) {
    real denom = 0;
    for (t in 1:T_train) denom += K_ts[t, s];
    if (denom <= 0) denom = 1;

    {
      real num = 0;
      for (t in 1:T_train) num += COM_ts[t, s] * (K_ts[t, s] / denom);
      com_wmean_train[s] = num;
    }
    {
      real v = 0;
      for (t in 1:T_train) {
        real wt = K_ts[t, s] / denom;
        v += wt * square(COM_ts[t, s] - com_wmean_train[s]);
      }
      com_wsd_train[s] = sqrt(fmax(v, 1e-16));
    }
  }

  // Cross-sectional standardization of the per-sector COM mean
  {
    real muS = mean(com_wmean_train);
    real sdS = sd(com_wmean_train);
    if (sdS <= 1e-8) sdS = 1.0;
    for (s in 1:S) COM_s[s] = (com_wmean_train[s] - muS) / sdS;
  }

  // reduce_sum index set: only the time points that enter the likelihood
  array[T_lik] int t_idx;
  for (t in 1:T_lik) t_idx[t] = t;
  int grainsize = 1;
}

parameters {
  // Cubic OU structure: NON-CENTERED hierarchical parameterization.
  // Standardized sector effects:
  vector[S] theta_z;
  vector[S] kappa_z;
  vector[S] a3_z;
  vector[S] beta0_z;

  // Hyperparameters (now actually used to build the sector parameters):
  real theta0;     real theta_COM;    real<lower=1e-6> sigma_theta;
  real kappa0;     real kappa_COM;    real<lower=1e-6> sigma_kappa;
  real a3_0;                          real<lower=1e-6> sigma_a3;
  real beta00;     real beta0_COM;    real<lower=1e-6> sigma_beta0;
  real beta1;

  // SV (stochastic volatility) and fat tails
  vector[S] alpha_s;
  vector<lower=-0.995, upper=0.995>[S] rho_s;
  vector<lower=1e-6>[S] sigma_eta_s;
  matrix[T,S] h_raw;
  real<lower=1e-6> nu_tilde;

  // TMG wedge (hard vs soft). Length 0 when hard (soft_wedge == 0): the wedge
  // is not used, so sampling 0 parameters avoids the degenerate near-zero-SD
  // geometry that throttled the NUTS step size.
  vector[soft_wedge == 1 ? T : 0] delta_z;

  // COM in mean
  real gamma;
}

transformed parameters {
  // Build sector parameters from the non-centered hierarchy.
  vector[S] theta_s;
  vector[S] kappa_tilde;
  vector[S] a3_tilde;
  vector[S] beta0_s;

  for (s in 1:S) {
    theta_s[s]     = theta0 + theta_COM * COM_s[s] + sigma_theta * theta_z[s];
    kappa_tilde[s] = kappa0 + kappa_COM * COM_s[s] + sigma_kappa * kappa_z[s];
    a3_tilde[s]    = a3_0                          + sigma_a3   * a3_z[s];
    beta0_s[s]     = beta00 + beta0_COM * COM_s[s] + sigma_beta0 * beta0_z[s];
  }

  vector<lower=0>[S] kappa_s = exp(kappa_tilde);
  vector<upper=0>[S] a3_s    = -exp(a3_tilde);   // restoring-force assumption a3<0
  real<lower=2> nu = 2 + nu_tilde;

  // Stationary non-centered AR(1) log-variance
  matrix[T,S] h;
  matrix[T,S] h_std;
  for (s in 1:S) {
    h_std[1,s] = h_raw[1,s] / sqrt(1 - square(rho_s[s]) + 1e-8);
    for (t in 2:T) {
      h_std[t,s] = rho_s[s] * h_std[t-1,s] + h_raw[t,s];
    }
  }
  for (t in 1:T) {
    for (s in 1:S) {
      h[t,s] = alpha_s[s] + sigma_eta_s[s] * h_std[t,s];
    }
  }
}

model {
  // ---- Hyperpriors (informative-but-weak; centers preserve the previous
  //      implied locations so results stay comparable to v0.1.3) ----
  theta0    ~ normal(0, 1);
  theta_COM ~ normal(0, 0.5);
  sigma_theta ~ normal(0, 1);

  kappa0    ~ normal(-1, 0.5);     // global log-kappa center (kappa ~ 0.37)
  kappa_COM ~ normal(0, 0.5);
  sigma_kappa ~ normal(0, 0.5);

  a3_0      ~ normal(log(0.05), 0.4);
  sigma_a3  ~ normal(0, 0.3);

  beta00    ~ normal(0, 0.5);
  beta0_COM ~ normal(0, 0.5);
  sigma_beta0 ~ normal(0, 0.5);

  // ---- Non-centered standardized sector effects ----
  theta_z ~ normal(0, 1);
  kappa_z ~ normal(0, 1);
  a3_z    ~ normal(0, 1);
  beta0_z ~ normal(0, 1);

  // ---- Key configurable priors ----
  beta1 ~ normal(beta1_prior_mean, beta1_prior_sd);   // neutral by default

  // ---- SV and tails ----
  alpha_s      ~ normal(0, 1);
  rho_s        ~ normal(rho_prior_mean, rho_prior_sd);
  sigma_eta_s  ~ normal(0, 0.5);
  to_vector(h_raw) ~ normal(0, 1);
  nu_tilde     ~ gamma(nu_prior_shape, nu_prior_rate);

  // ---- TMG wedge ----
  // When hard (soft_wedge == 0) delta_z has length 0 and needs no prior.
  if (soft_wedge == 1) delta_z ~ normal(zTMG_exo - zTMG_byK, sigma_delta_z);

  // ---- COM in mean ----
  gamma ~ normal(0, 0.5);

  // ---- Parallelized likelihood (reduce_sum), restricted to t <= T_lik ----
  target += reduce_sum(ou_nl_partial_sum, t_idx, grainsize,
                       T_lik, Yz, Xz, COM_ts, zTMG_byK, soft_wedge, delta_z,
                       com_wmean_train, com_wsd_train, mu_xz,
                       theta_s, kappa_s, a3_s, beta0_s, beta1,
                       h, nu, com_in_mean, gamma);
}

generated quantities {
  // Pointwise one-step-ahead log density, computed on EXACTLY the same window
  // used for fitting (2..T_lik). Zero elsewhere so the R side can subset
  // consistently for PSIS-LOO.
  matrix[T,S] log_lik;
  for (t in 1:T) for (s in 1:S) log_lik[t,s] = 0;

  for (t in 2:T_lik) {
    real ztmg_eff = zTMG_byK[t];
    if (soft_wedge == 1) ztmg_eff += delta_z[t];
    ztmg_eff = fmin(fmax(ztmg_eff, -1e6), 1e6);
    for (s in 1:S) {
      real zlag  = Yz[t-1,s] - theta_s[s];
      real drift = kappa_s[s] * (theta_s[s] - Yz[t-1,s] + a3_s[s] * zlag^3);
      real betaT = beta0_s[s] + beta1 * ztmg_eff;

      real denom_sd = com_wsd_train[s];
      denom_sd = (denom_sd > 1e-12) ? denom_sd : 1.0;
      real com_std = (COM_ts[t-1,s] - com_wmean_train[s]) / denom_sd;
      com_std = fmin(fmax(com_std, -1e6), 1e6);

      real mean_ = drift + betaT * (Xz[t-1,s] - mu_xz[s])
                   + (com_in_mean == 1 ? gamma * com_std : 0);
      real sd_   = fmin(fmax(exp(0.5 * h[t,s]), 1e-8), 1e8);

      log_lik[t,s] = student_t_lpdf( Yz[t,s] - Yz[t-1,s] - mean_ | nu, 0, sd_ );
    }
  }
}
