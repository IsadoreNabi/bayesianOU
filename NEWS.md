# bayesianOU 0.1.4

Robustness overhaul. Several changes alter results and are intentional; review
the methodology section of the README and re-run with your data.

## Correctness fixes

* **PSIS-LOO**: the log-likelihood was passed to `loo::loo()` as a 3-D array
  `[draws, time, sector]`, which `loo` interprets as `[iterations, chains,
  observations]` — it silently treated time as chains and only the sectors as
  observations. It is now reshaped to a proper `[draws x observations]` matrix
  over the fitted window, with a `chain_id` so `relative_eff()` is correct.
* **Train/test leakage**: the likelihood was summed over the full sample even
  when a train/test split was used, so the "out-of-sample" evaluation was
  contaminated. New `fit_window` argument (`"train"` default / `"full"`) keeps
  the likelihood, `log_lik`, and OOS evaluation coherent.
* **Real hierarchy**: `kappa0/kappa_COM/sigma_kappa`, `a3_0/sigma_a3`,
  `beta00/beta0_COM/sigma_beta0` were declared but unused — only `theta_s` was
  hierarchical. All four sector blocks now use a non-centered hierarchical
  (partial-pooling) parameterization that actually uses the hyperparameters.
* **`evaluate_oos` index bug**: for horizons longer than the test window the
  `(T_train+1):(Tn-hh+1)` sequence ran backwards and indexed out of range; it
  now guards the window and returns `NA`/`n_obs = 0` instead.
* **Factor leakage**: `use_train_loadings` now defaults to `TRUE` (loadings from
  training only), avoiding look-ahead in the orthogonalized TMG regressor.
* **`compare_models_loo`**: validates that both `loo` objects exist, are of
  class `loo`, and have matching observation counts; `deltaELPD` is now computed
  unambiguously as `elpd_new - elpd_base`.

## Priors (defaults changed — override via `priors`)

* `beta1 ~ Normal(0, 0.5)` (was `Normal(0.5, 0.25)`): neutral on the TMG-effect
  hypothesis (no sign baked into the prior).
* `nu_tilde ~ Gamma(2, 0.1)` (was `Exponential(3)`): weakly informative; the old
  prior forced `nu` into `(2,3)` (extreme heavy tails).
* `rho_s ~ Normal(0.7, 0.2)` (was `Normal(0.90, 0.05)`): less rigid SV
  persistence.
* New tunable prior entries: `beta1_mean`, `beta1_sd`, `nu_shape`, `nu_rate`,
  `rho_mean`, `rho_sd`.

## Diagnostics and robustness

* `validate_ou_fit()` is now a real validator: structural checks plus a Pareto-k
  summary with a warning when PSIS-LOO is unreliable; separates MCMC convergence
  from dynamic mean reversion.
* `extract_convergence_evidence()` no longer prints "CONVERGENCE GUARANTEED"; it
  reports dynamic mean-reversion evidence and is aliased by the clearer
  `kappa_stability_evidence()`.
* `count_divergences()` / LOO are surfaced with `loo_pareto_k` in `diagnostics`.

## Engineering

* The TMG wedge `delta_z` now has length 0 in the default hard case
  (`hard_sum_zero = TRUE`). Previously it was sampled as `T` parameters pinned by
  a `normal(0, 1e-6)` prior, whose near-zero scale throttled the NUTS step size
  and made sampling crawl. Removing them is both faster and better-mixing.
* The Stan model is now compiled into a writable user cache
  (`tools::R_user_dir`), not the (possibly read-only) package directory, and the
  compiled binary is no longer left in the source tree.
* Single source of truth for the Stan model: `inst/stan/ou_nonlinear_tmg.stan`
  is canonical and is read/compiled directly; the R string copy is gone.
* `fit_ou_nonlinear_tmg()` no longer calls `setwd()`; validates inputs
  (finiteness, `iter > warmup`, `train_frac`); caps threads to avoid
  oversubscription; exposes `train_frac`; `thin` default is now `1`; no longer
  duplicates the full draws array in the result (access via `stan_fit`).
* Declared `parallel` in Imports; removed stray `LazyData: true` (no `data/`);
  fixed repository URL capitalization.
* New tests: `evaluate_oos` index guard, PSIS-LOO reshape, and a Stan-based
  parameter-recovery test (skipped when no backend is available).
