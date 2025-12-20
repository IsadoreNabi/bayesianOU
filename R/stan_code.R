#' Stan code for nonlinear OU model with SV and Student-t
#'
#' Returns the complete Stan code for the nonlinear Ornstein-Uhlenbeck model
#' with cubic drift, stochastic volatility, and Student-t innovations.
#' Includes numerical guardrails and parallel computation support.
#'
#' @return Character string containing Stan model code
#'
#' @details
#' The model implements:
#' \itemize{
#'   \item Cubic drift: \eqn{\kappa(\theta - Y + a_3(Y-\theta)^3)}
#'   \item Stochastic volatility with AR(1) log-variance
#'   \item Student-t innovations with estimated degrees of freedom
#'   \item Hierarchical priors for sector-specific parameters
#'   \item Optional soft constraint on TMG variable
#'   \item Parallel likelihood computation via reduce_sum
#' }
#'
#' @examples
#' code <- ou_nonlinear_tmg_stan_code()
#' cat(substr(code, 1, 500))
#'
#' @export
ou_nonlinear_tmg_stan_code <- function() {
  paste0(
'functions {
  // Partial likelihood for reduce_sum: sums over t in [start:end]
  real ou_nl_partial_sum(array[] int t_idx_slice,
                         int start, int end,
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
      real ztmg_eff = zTMG_byK[t] + (soft_wedge == 1 ? delta_z[t] : 0);
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
}

transformed data {
  vector[S] com_wmean_train;
  vector[S] com_wsd_train;
  vector[S] COM_s;

  // Weighted mean by CAPITAL_TOTAL in TRAIN
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

  // Cross-sectional standardization of COM_s
  {
    real muS = mean(com_wmean_train);
    real sdS = sd(com_wmean_train);
    if (sdS <= 1e-8) sdS = 1.0;
    for (s in 1:S) COM_s[s] = (com_wmean_train[s] - muS) / sdS;
  }

  // Indices for reduce_sum and grainsize
  array[T] int t_idx;
  for (t in 1:T) t_idx[t] = t;
  int grainsize = 4;
}

parameters {
  // Cubic OU structure
  vector[S] theta_s;
  vector[S] kappa_tilde;
  vector[S] a3_tilde;
  vector[S] beta0_s;
  real kappa0;     real kappa_COM;    real<lower=1e-6> sigma_kappa;
  real theta0;     real theta_COM;    real<lower=1e-6> sigma_theta;
  real a3_0;       real<lower=1e-6> sigma_a3;
  real beta00;     real beta0_COM;    real<lower=1e-6> sigma_beta0;
  real beta1;

  // SV (stochastic) and fat tails
  vector[S] alpha_s;
  vector<lower=-0.995, upper=0.995>[S] rho_s;
  vector<lower=1e-6>[S] sigma_eta_s;
  matrix[T,S] h_raw;
  real<lower=0> nu_tilde;

  // TMG wedge (hard vs soft)
  vector[T] delta_z;

  // COM in mean
  real gamma;
}

transformed parameters {
  vector<lower=0>[S] kappa_s = exp(kappa_tilde);
  vector[S] a3_s    = -exp(a3_tilde);
  real<lower=2> nu = 2 + nu_tilde;

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
  // Hyper and hierarchical priors
  sigma_kappa ~ normal(0,1);
  sigma_theta ~ normal(0,1);
  sigma_a3    ~ normal(0,1);
  sigma_beta0 ~ normal(0,1);

  kappa0   ~ normal(0,1);
  kappa_COM~ normal(0,1);
  theta0   ~ normal(0,1);
  theta_COM~ normal(0,1);
  a3_0     ~ normal(0,1);
  beta00   ~ normal(0,1);
  beta0_COM~ normal(0,1);

  theta_s      ~ normal(theta0 + theta_COM * COM_s, sigma_theta);
  kappa_tilde  ~ normal(-1, 0.5);
  a3_tilde     ~ normal(log(0.05), 0.4);
  beta0_s      ~ normal(0, 0.5);
  beta1        ~ normal(0.5, 0.25);

  // SV and tails
  alpha_s      ~ normal(0, 1);
  rho_s        ~ normal(0.90, 0.05);
  sigma_eta_s  ~ normal(0, 0.5);
  to_vector(h_raw) ~ normal(0,1);
  nu_tilde ~ exponential(3);

  // TMG wedge
  if (soft_wedge == 1) delta_z ~ normal(zTMG_exo - zTMG_byK, sigma_delta_z);
  else                 delta_z ~ normal(0, 1e-6);

  // COM in mean: soft prior
  gamma ~ normal(0, 0.5);

  // Parallelized likelihood (reduce_sum) with guardrails
  target += reduce_sum(ou_nl_partial_sum, t_idx, grainsize,
                       Yz, Xz, COM_ts, zTMG_byK, soft_wedge, delta_z,
                       com_wmean_train, com_wsd_train, mu_xz,
                       theta_s, kappa_s, a3_s, beta0_s, beta1,
                       h, nu, com_in_mean, gamma);
}

generated quantities {
  // log_lik only in TRAIN for PSIS-LOO
  matrix[T,S] log_lik;
  for (t in 1:T) for (s in 1:S) log_lik[t,s] = 0;

  {
    real nu_ = nu;
    for (t in 2:T_train) {
      real ztmg_eff = zTMG_byK[t] + (soft_wedge == 1 ? delta_z[t] : 0);
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
        real sd_   = exp(0.5 * h[t,s]);
        sd_ = fmin(fmax(sd_, 1e-8), 1e8);

        log_lik[t,s] = student_t_lpdf( Yz[t,s] - Yz[t-1,s] - mean_ | nu_, 0, sd_ );
      }
    }
  }
}
'
  )
}
