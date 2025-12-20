# bayesianOU

Bayesian Nonlinear Ornstein-Uhlenbeck Models with Stochastic Volatility

## Overview

The `bayesianOU` package fits Bayesian nonlinear Ornstein-Uhlenbeck models 
with cubic drift, stochastic volatility (SV), and Student-t innovations. 
It implements hierarchical priors for sector-specific parameters and supports 
parallel MCMC sampling via Stan.

## Installation

```r
# Install from GitHub (development version)
# install.packages("remotes")
remotes::install_github("author/bayesianOU")

# For Stan backend, you need either cmdstanr or rstan
# cmdstanr (recommended):
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

# Or rstan:
install.packages("rstan")
```

## Quick Start

```r
library(bayesianOU)

# Prepare data
Y <- as.matrix(your_prices_data)
X <- as.matrix(your_production_prices_data)
TMG <- your_tmg_series
COM <- as.matrix(your_com_data)
K <- as.matrix(your_capital_data)

# Fit model
results <- fit_ou_nonlinear_tmg(
  results_robust = list(),
  Y = Y, X = X, TMG = TMG, COM = COM, CAPITAL_TOTAL = K,
  chains = 4, iter = 8000, warmup = 4000,
  verbose = TRUE
)

# Validate fit
validate_ou_fit(results)

# Extract convergence evidence
conv <- extract_convergence_evidence(results)

# Plot results
plot_beta_tmg(results)
plot_drift_curves(results)
```

## Model Specification

The model implements a nonlinear OU process with cubic drift:

$$dY_t = \kappa(\theta - Y_t + a_3 (Y_t - \theta)^3) dt + \sigma_t dW_t$$

where:
- $\kappa_s$ is the sector-specific mean reversion speed
- $\theta_s$ is the sector-specific equilibrium level  
- $a_3$ is the cubic nonlinearity coefficient
- $\sigma_t$ follows an AR(1) stochastic volatility process
- Innovations are Student-t distributed with estimated degrees of freedom

## Features

- Hierarchical priors for sector-specific parameters
- Stochastic volatility with AR(1) log-variance
- Student-t innovations for fat tails
- Parallel likelihood computation via Stan's reduce_sum
- PSIS-LOO cross-validation for model comparison
- Out-of-sample forecast evaluation

## Citation

If you use this package, please cite:

```
@software{bayesianOU,
  author = {Author Name},
  title = {bayesianOU: Bayesian Nonlinear Ornstein-Uhlenbeck Models},
  year = {2024},
  url = {https://github.com/author/bayesianOU}
}
```

## References

- Stan User's Guide (SV, HMM, parallelization)
- Vehtari, Gelman, Gabry (2017). Practical Bayesian model evaluation 
  using leave-one-out cross-validation. Statistics and Computing.
- Gelman (2006). Prior distributions for variance parameters in 
  hierarchical models. Bayesian Analysis.

## License

MIT License
