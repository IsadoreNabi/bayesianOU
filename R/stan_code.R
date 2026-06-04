#' Locate the canonical Stan model file
#'
#' Resolves the path to \code{inst/stan/ou_nonlinear_tmg.stan}, the single
#' source of truth for the model. Works both for the installed package
#' (via \code{system.file}) and during development (source-tree fallbacks).
#'
#' @return Character path to the \code{.stan} file.
#'
#' @keywords internal
#' @noRd
.stan_file_path <- function() {
  p <- system.file("stan", "ou_nonlinear_tmg.stan", package = "bayesianOU")
  if (nzchar(p) && file.exists(p)) {
    return(p)
  }

  # Development fallbacks (devtools::load_all, running tests from source, etc.)
  candidates <- c(
    file.path("inst", "stan", "ou_nonlinear_tmg.stan"),
    file.path("..", "inst", "stan", "ou_nonlinear_tmg.stan"),
    file.path("..", "..", "inst", "stan", "ou_nonlinear_tmg.stan")
  )
  for (cand in candidates) {
    if (file.exists(cand)) {
      return(normalizePath(cand))
    }
  }

  stop(
    "Could not locate 'ou_nonlinear_tmg.stan'. Reinstall the package or run ",
    "from the package root.",
    call. = FALSE
  )
}


#' Stan code for the nonlinear OU model with SV and Student-t
#'
#' Returns the complete Stan code for the nonlinear Ornstein-Uhlenbeck model
#' with cubic drift, stochastic volatility, Student-t innovations, and a
#' non-centered hierarchical (partial-pooling) parameterization for the
#' sector-specific parameters.
#'
#' The code is read from the canonical file \code{inst/stan/ou_nonlinear_tmg.stan}
#' (single source of truth); this function does not embed a duplicate copy.
#'
#' @return Character string containing the Stan model code.
#'
#' @details
#' The model implements:
#' \itemize{
#'   \item Cubic drift: \eqn{\kappa(\theta - Y + a_3(Y-\theta)^3)} with the
#'     restoring-force assumption \eqn{a_3 < 0} (enforced via \code{-exp()}).
#'   \item Stochastic volatility with a stationary AR(1) log-variance.
#'   \item Student-t innovations with estimated degrees of freedom (\eqn{\nu>2}).
#'   \item Real non-centered hierarchical priors for \code{theta_s},
#'     \code{kappa_s}, \code{a3_s} and \code{beta0_s} (the hyperparameters are
#'     actually used to build the sector parameters).
#'   \item Configurable priors for \code{beta1}, \code{nu} and \code{rho_s}.
#'   \item A likelihood/log-lik window controlled by \code{T_lik} to keep the
#'     training/test split coherent (no leakage).
#'   \item Parallel likelihood computation via \code{reduce_sum}.
#' }
#'
#' @examples
#' code <- ou_nonlinear_tmg_stan_code()
#' cat(substr(code, 1, 300))
#'
#' @export
ou_nonlinear_tmg_stan_code <- function() {
  paste(readLines(.stan_file_path(), warn = FALSE), collapse = "\n")
}
