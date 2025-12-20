#' bayesianOU: Bayesian Nonlinear Ornstein-Uhlenbeck Models
#'
#' Fits Bayesian nonlinear Ornstein-Uhlenbeck models with cubic drift,
#' stochastic volatility (SV), and Student-t innovations. The package
#' implements hierarchical priors for sector-specific parameters and
#' supports parallel MCMC sampling via 'Stan'.
#'
#' @section Main Functions:
#' \itemize{
#'    \item \code{\link{fit_ou_nonlinear_tmg}}: Fit the main OU model
#'    \item \code{\link{extract_posterior_summary}}: Extract posterior summaries
#'    \item \code{\link{validate_ou_fit}}: Validate model fit
#'    \item \code{\link{compare_models_loo}}: Compare models via PSIS-LOO
#' }
#'
#' @section Model Specification:
#' The model implements a nonlinear OU process with cubic drift:
#' \deqn{dY_t = \kappa(\theta - Y_t + a_3 (Y_t - \theta)^3) dt + \sigma_t dW_t}
#'
#' @docType package
#' @name bayesianOU-package
#' @aliases bayesianOU
#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "bayesianOU: Bayesian Nonlinear OU Models\n",
    "For 'Stan' backend, ensure 'cmdstanr' or 'rstan' is installed."
  )
}