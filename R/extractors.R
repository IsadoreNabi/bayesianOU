#' Extract posterior summary from fitted model
#'
#' Extracts median point estimates and credible intervals for all
#' model parameters from a fitted Stan model.
#'
#' @param fit Fitted Stan model object (CmdStanMCMC or stanfit)
#'
#' @return List with components:
#'   \describe{
#'     \item{beta1}{Median of global TMG effect}
#'     \item{beta0_s}{Vector of sector-specific intercepts}
#'     \item{kappa_s}{Vector of mean reversion speeds}
#'     \item{a3_s}{Vector of cubic drift coefficients}
#'     \item{theta_s}{Vector of equilibrium levels}
#'     \item{rho_s}{Vector of SV persistence parameters}
#'     \item{alpha_s}{Vector of SV level parameters}
#'     \item{sigma_eta_s}{Vector of SV volatility parameters}
#'     \item{nu}{Degrees of freedom for Student-t}
#'     \item{gamma}{COM effect in mean}
#'     \item{rhat}{R-hat convergence diagnostics}
#'     \item{ess}{Effective sample sizes}
#'   }
#'
#' @examples
#' \donttest{
#' # 1. Create a mock CmdStanMCMC object
#' # We simulate a posterior distribution for 2 sectors
#' S <- 2
#' n_draws <- 100
#' 
#' # Helper to generate random draws
#' mock_draws <- function(name, n_cols=1) {
#'   m <- matrix(rnorm(n_draws * n_cols), nrow = n_draws, ncol = n_cols)
#'   if (n_cols > 1) {
#'     colnames(m) <- paste0(name, "[", 1:n_cols, "]")
#'   } else {
#'     colnames(m) <- name
#'   }
#'   as.data.frame(m)
#' }
#' 
#' # Combine draws into one data frame
#' df_draws <- cbind(
#'   mock_draws("beta1", 1),
#'   mock_draws("beta0_s", S),
#'   mock_draws("kappa_tilde", S), # Note: function expects log-scale kappa
#'   mock_draws("a3_tilde", S),    # Note: function expects log-scale a3
#'   mock_draws("theta_s", S),
#'   mock_draws("rho_s", S),
#'   mock_draws("alpha_s", S),
#'   mock_draws("sigma_eta_s", S),
#'   mock_draws("nu_tilde", 1),
#'   mock_draws("gamma", 1)
#' )
#' 
#' # Mock fit object
#' mock_fit <- structure(list(
#'   draws = function(vars, format="df") {
#'      # Simple regex matching for the mock
#'      if (length(vars) == 1) {
#'        # Check if it's a scalar or vector parameter request
#'        if (vars %in% names(df_draws)) return(df_draws[vars])
#'        # Pattern match for vectors like "beta0_s" -> "beta0_s[1]", "beta0_s[2]"
#'        cols <- grep(paste0("^", vars, "\\["), names(df_draws), value = TRUE)
#'        if (length(cols) > 0) return(df_draws[cols])
#'      }
#'      return(df_draws) 
#'   },
#'   summary = function() {
#'     data.frame(
#'       variable = names(df_draws),
#'       rhat = rep(1.0, ncol(df_draws)),
#'       ess_bulk = rep(400, ncol(df_draws))
#'     )
#'   }
#' ), class = "CmdStanMCMC")
#' 
#' # 2. Run extraction
#' summ <- extract_posterior_summary(mock_fit)
#' print(summ$kappa_s)
#' }
#'
#' @export
extract_posterior_summary <- function(fit) {
  
  draw_mat <- function(p) {
    if (inherits(fit, "CmdStanMCMC")) {
      df <- fit$draws(p, format = "df")
      keep <- grep(sprintf("^%s\\[\\d+\\]$", p), names(df), value = TRUE)
      if (length(keep) > 0) {
        return(as.matrix(df[, keep, drop = FALSE]))
      }
      if (p %in% names(df)) {
        return(as.matrix(df[, p, drop = FALSE]))
      }
      warning(sprintf("Parameter '%s' not found in draws", p))
      return(matrix(NA_real_, nrow = nrow(df), ncol = 1))
    } else {
      rstan::extract(fit, pars = p)[[1]]
    }
  }
  
  med <- function(M) {
    M <- if (is.null(dim(M))) matrix(M, ncol = 1) else M
    apply(M, 2, stats::median)
  }
  
  beta1 <- stats::median(draw_mat("beta1"))
  beta0_s <- med(draw_mat("beta0_s"))
  kappa_s <- med(exp(draw_mat("kappa_tilde")))
  a3_s <- -med(exp(draw_mat("a3_tilde")))
  theta_s <- med(draw_mat("theta_s"))
  rho_s <- med(draw_mat("rho_s"))
  alpha_s <- med(draw_mat("alpha_s"))
  sigma_eta_s <- med(draw_mat("sigma_eta_s"))
  nu <- stats::median(2 + draw_mat("nu_tilde"))
  gamma <- stats::median(draw_mat("gamma"))
  
  if (inherits(fit, "CmdStanMCMC")) {
    sf <- fit$summary()
    rhat <- sf$rhat
    ess <- sf$ess_bulk
  } else {
    sm <- rstan::summary(fit)$summary
    rhat <- sm[, "Rhat"]
    ess <- sm[, "n_eff"]
  }
  
  list(
    beta1 = beta1,
    beta0_s = beta0_s,
    kappa_s = kappa_s,
    a3_s = a3_s,
    theta_s = theta_s,
    rho_s = rho_s,
    alpha_s = alpha_s,
    sigma_eta_s = sigma_eta_s,
    nu = nu,
    gamma = gamma,
    rhat = rhat,
    ess = ess
  )
}


#' Extract convergence evidence for kappa parameters
#'
#' Computes 95 percent credible intervals for each kappa_s (mean reversion speed)
#' and verifies formal convergence conditions.
#'
#' @param fit_res List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param verbose Logical. Print summary to console. Default TRUE.
#'
#' @return List with components:
#'   \describe{
#'     \item{kappa_ic95}{Matrix (S x 3) with columns q2.5, median, q97.5}
#'     \item{convergence}{Logical indicating if all kappa in (0,1)}
#'     \item{prob_convergence}{Posterior probability of joint convergence}
#'   }
#'
#' @examples
#' \donttest{
#' # 1. Create a mock fit object with kappa draws
#' # kappa_tilde is log(kappa), so we use log(0.5) roughly -0.69
#' n_draws <- 100
#' S <- 2
#' kappa_tilde_draws <- matrix(rnorm(n_draws * S, mean = -0.7, sd = 0.1), 
#'                             nrow = n_draws, ncol = S)
#' colnames(kappa_tilde_draws) <- c("kappa_tilde[1]", "kappa_tilde[2]")
#' 
#' mock_fit <- structure(list(
#'   draws = function(vars, format="matrix") {
#'     if (vars == "kappa_tilde") return(kappa_tilde_draws)
#'     return(NULL)
#'   }
#' ), class = "CmdStanMCMC")
#' 
#' # 2. Wrap in the results list structure
#' results_mock <- list(
#'   factor_ou = list(
#'     stan_fit = mock_fit
#'   )
#' )
#' 
#' # 3. Extract evidence
#' conv <- extract_convergence_evidence(results_mock)
#' print(conv$prob_convergence)
#' }
#'
#' @export
extract_convergence_evidence <- function(fit_res, verbose = TRUE) {
  fit <- fit_res$factor_ou$stan_fit
  
  if (inherits(fit, "CmdStanMCMC")) {
    kappa_draws <- exp(fit$draws("kappa_tilde", format = "matrix"))
  } else {
    kappa_draws <- exp(rstan::extract(fit, "kappa_tilde")[[1]])
  }
  
  kappa_ic95 <- t(apply(
    kappa_draws, 2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  ))
  colnames(kappa_ic95) <- c("q2.5", "median", "q97.5")
  
  all_positive <- all(kappa_ic95[, "q2.5"] > 0)
  all_less_one <- all(kappa_ic95[, "q97.5"] < 1)
  convergence <- all_positive && all_less_one
  
  prob_conv <- mean(apply(kappa_draws > 0 & kappa_draws < 1, 1, all))
  
  if (verbose) {
    message("\n=== FORMAL CONVERGENCE EVIDENCE (NONLINEAR) ===")
    message(paste(rep("=", 50), collapse = ""))
    message("95% CI for kappa_s (mean reversion speed) by sector:\n")
    print(utils::head(kappa_ic95, 10))
    message("\nConvergence verification:")
    message(sprintf("- All kappa > 0: %s", all_positive))
    message(sprintf("- All kappa < 1: %s", all_less_one))
    message(sprintf("- CONVERGENCE GUARANTEED: %s", convergence))
    message(sprintf("\nP(convergence | data) = %.3f", prob_conv))
  }
  
  list(
    kappa_ic95 = kappa_ic95,
    convergence = convergence,
    prob_convergence = prob_conv
  )
}


#' Build beta(TMG_t) table by sector and time
#'
#' Constructs the time-varying beta matrix using posterior medians.
#'
#' @param fit Fitted Stan model object
#' @param zTMG_use Numeric vector. Standardized TMG series used in fitting.
#'
#' @return List with components:
#'   \describe{
#'     \item{beta_point}{Matrix (T x S) of beta values}
#'     \item{meta}{List with description metadata}
#'   }
#'
#' @export
build_beta_tmg_table <- function(fit, zTMG_use) {
  summ <- extract_posterior_summary(fit)
  
  beta_ts <- outer(zTMG_use, rep(1, length(summ$beta0_s))) * summ$beta1 +
    matrix(
      rep(summ$beta0_s, each = length(zTMG_use)),
      nrow = length(zTMG_use)
    )
  
  list(
    beta_point = beta_ts,
    meta = list(description = "beta_s(TMG_t) with median beta1 and beta0_s")
  )
}


#' Summarize stochastic volatility sigmas
#'
#' Extracts median volatility paths from the SV component.
#'
#' @param fit Fitted Stan model object
#'
#' @return List with component:
#'   \describe{
#'     \item{sigma_t}{Matrix (T x S) of median volatilities}
#'   }
#'
#' @export
summarize_sv_sigmas <- function(fit) {
  if (inherits(fit, "CmdStanMCMC")) {
    H <- fit$draws("h", format = "matrix")
    medH <- apply(H, 2, stats::median)
    dfh <- fit$draws("h", format = "df")
    ts_labels <- colnames(dfh)
    ts_labels <- ts_labels[grepl("^h\\[", ts_labels)]
    Tguess <- max(as.integer(gsub("^h\\[(\\d+),.*$", "\\1", ts_labels)))
    Sguess <- max(as.integer(gsub("^h\\[\\d+,(\\d+)\\]$", "\\1", ts_labels)))
    matrix_sigma <- matrix(exp(0.5 * medH), nrow = Tguess, ncol = Sguess)
    list(sigma_t = matrix_sigma)
  } else {
    H <- rstan::extract(fit, "h")[[1]]
    Hmed <- apply(H, c(2, 3), stats::median)
    list(sigma_t = exp(0.5 * Hmed))
  }
}


#' Drift decomposition over grid
#'
#' Computes the cubic OU drift function over a grid of z values.
#'
#' @param fit Fitted Stan model object (reserved for future use)
#' @param summ List. Output from \code{\link{extract_posterior_summary}}
#' @param z_grid Numeric vector. Grid of z values to evaluate drift.
#'
#' @return List with components:
#'   \describe{
#'     \item{z}{The input z grid}
#'     \item{drift}{Matrix (length(z) x S) of drift values by sector}
#'   }
#'
#' @export
drift_decomposition_grid <- function(fit, summ,
                                     z_grid = seq(-2.5, 2.5, length.out = 101)) {
  m <- sapply(seq_along(summ$kappa_s), function(s) {
    summ$kappa_s[s] * (-z_grid + summ$a3_s[s] * z_grid^3)
  })
  
  list(z = z_grid, drift = m)
}


#' Build accounting block for TMG
#'
#' Creates accounting information for TMG transformations.
#'
#' @param TMG_raw Numeric vector. Original TMG series.
#' @param zTMG_exo Numeric vector. Exogenous z-scored TMG.
#' @param zTMG_use Numeric vector. TMG used in model (possibly orthogonalized).
#' @param mu_tmg Numeric. Training mean of TMG.
#' @param sd_tmg Numeric. Training SD of TMG.
#' @param hard Logical. Whether hard sum-to-zero constraint was used.
#' @param sigma_delta Numeric. Prior SD for wedge in original units.
#'
#' @return List with components:
#'   \describe{
#'     \item{tmg_byK}{Back-transformed TMG used in model}
#'     \item{tmg_exo}{Back-transformed exogenous TMG}
#'     \item{wedge_delta}{Difference (zero if hard=TRUE)}
#'     \item{sigma_delta_prior}{Prior SD for wedge}
#'     \item{note}{Description of constraint type}
#'   }
#'
#' @export
build_accounting_block <- function(TMG_raw, zTMG_exo, zTMG_use,
                                   mu_tmg, sd_tmg, hard, sigma_delta) {
  tmg_byK <- zTMG_use * sd_tmg + mu_tmg
  tmg_exo <- zTMG_exo * sd_tmg + mu_tmg
  delta <- tmg_byK - tmg_exo
  
  list(
    tmg_byK = tmg_byK,
    tmg_exo = tmg_exo,
    wedge_delta = if (hard) rep(0, length(delta)) else delta,
    sigma_delta_prior = sigma_delta,
    note = if (hard) {
      "hard=TRUE: TMG identical to byK"
    } else {
      "hard=FALSE: delta_t reported in original scale"
    }
  )
}
