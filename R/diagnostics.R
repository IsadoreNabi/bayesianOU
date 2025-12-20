#' Count HMC divergences
#'
#' Extracts the number of divergent transitions from a fitted Stan model.
#'
#' @param fit Fitted Stan model object (CmdStanMCMC or stanfit)
#'
#' @return Integer. Number of divergent transitions (post-warmup).
#'
#' @examples
#' \donttest{
#' # Create a "mock" CmdStanMCMC object for demonstration
#' # (This simulates a model with 0 divergences)
#' mock_fit <- structure(list(
#'   sampler_diagnostics = function() {
#'     # Return a 3D array: [iterations, chains, variables]
#'     # Variable 1 is usually accept_stat__, let's say var 2 is divergent__
#'     ar <- array(0, dim = c(100, 4, 6)) 
#'     dimnames(ar)[[3]] <- c("accept_stat__", "divergent__", "energy__", 
#'                            "n_leapfrog__", "stepsize__", "treedepth__")
#'     return(ar)
#'   }
#' ), class = "CmdStanMCMC")
#' 
#' # Now the example can run without errors:
#' n_div <- count_divergences(mock_fit)
#' print(n_div)
#' }
#'
#' @export
count_divergences <- function(fit) {
  tryCatch({
    if (inherits(fit, "CmdStanMCMC")) {
      sd <- fit$sampler_diagnostics()
      vnames <- dimnames(sd)[[3]]
      if (!("divergent__" %in% vnames)) {
        return(NA_integer_)
      }
      i <- match("divergent__", vnames)
      return(sum(sd[, , i], na.rm = TRUE))
    } else {
      sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      return(sum(vapply(sp, function(m) sum(m[, "divergent__"]), numeric(1))))
    }
  }, error = function(e) {
    if (inherits(fit, "CmdStanMCMC")) {
      csvs <- tryCatch(fit$output_files(), error = function(e) NULL)
      if (is.null(csvs) || length(csvs) == 0) {
        return(NA_integer_)
      }
      sd2 <- cmdstanr::read_cmdstan_csv(csvs)$sampler_diagnostics
      if (length(dim(sd2)) == 3 && "divergent__" %in% dimnames(sd2)[[3]]) {
        return(sum(sd2[, , "divergent__"], na.rm = TRUE))
      }
    }
    NA_integer_
  })
}


#' Evaluate out-of-sample forecast metrics
#'
#' Computes RMSE and MAE for multiple forecast horizons.
#'
#' @param summ List. Posterior summary from \code{\link{extract_posterior_summary}}
#' @param Yz Numeric matrix. Standardized Y values (T x S)
#' @param Xz Numeric matrix. Standardized X values (T x S)
#' @param zTMG Numeric vector. Standardized TMG series
#' @param T_train Integer. End of training period
#' @param COM_ts Numeric matrix. COM values by time and sector (T x S)
#' @param K_ts Numeric matrix. Capital values by time and sector (T x S)
#' @param com_in_mean Logical. Whether COM is included in mean equation
#' @param horizons Integer vector. Forecast horizons to evaluate
#'
#' @return Named list with one element per horizon...
#'
#' @examples
#' # 1. Generate dummy data for testing
#' T_obs <- 20
#' S <- 2
#' Yz <- matrix(rnorm(T_obs * S), nrow = T_obs, ncol = S)
#' Xz <- matrix(rnorm(T_obs * S), nrow = T_obs, ncol = S)
#' COM_ts <- matrix(abs(rnorm(T_obs * S)), nrow = T_obs, ncol = S)
#' K_ts <- matrix(abs(rnorm(T_obs * S)) + 1, nrow = T_obs, ncol = S)
#' zTMG <- rnorm(T_obs)
#' 
#' # 2. Create a dummy summary list (mimicking extract_posterior_summary)
#' summ <- list(
#'   theta_s = runif(S),
#'   kappa_s = runif(S),
#'   a3_s = runif(S),
#'   beta0_s = runif(S),
#'   beta1 = 0.5,
#'   gamma = 0.1
#' )
#' 
#' # 3. Run the function
#' metrics <- evaluate_oos(summ, Yz, Xz, zTMG, T_train = 15, 
#'                         COM_ts, K_ts, horizons = c(1, 2))
#' print(metrics)
#'
#' @export
evaluate_oos <- function(summ, Yz, Xz, zTMG, T_train,
                         COM_ts, K_ts, com_in_mean = FALSE,
                         horizons = c(1, 4, 8)) {
  
  Tn <- nrow(Yz)
  S <- ncol(Yz)
  
  com_wmean <- numeric(S)
  com_wsd <- numeric(S)
  
  for (s in seq_len(S)) {
    denom <- sum(K_ts[seq_len(T_train), s], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    w <- K_ts[seq_len(T_train), s] / denom
    com_wmean[s] <- sum(COM_ts[seq_len(T_train), s] * w, na.rm = TRUE)
    v <- sum(w * (COM_ts[seq_len(T_train), s] - com_wmean[s])^2, na.rm = TRUE)
    com_wsd[s] <- sqrt(max(v, 1e-16))
  }
  
  res <- lapply(horizons, function(hh) {
    errs <- c()
    
    for (t in (T_train + 1):(Tn - hh + 1)) {
      if (t - 1 < 1) next
      
      y_pred <- Yz[t - 1, ]
      
      for (h in seq_len(hh)) {
        t_pred <- t - 1 + h
        if (t_pred > Tn) break
        
        ztmg <- if (t_pred <= Tn) zTMG[t_pred] else zTMG[Tn]
        
        for (s in seq_len(S)) {
          zlag <- y_pred[s] - summ$theta_s[s]
          drift <- summ$kappa_s[s] * (
            summ$theta_s[s] - y_pred[s] + summ$a3_s[s] * zlag^3
          )
          xterm <- (summ$beta0_s[s] + summ$beta1 * ztmg) *
            (if (t_pred <= Tn) Xz[t_pred - 1, s] else Xz[Tn, s])
          
          com_std_denom <- if (com_wsd[s] > 0) com_wsd[s] else 1
          com_std <- (COM_ts[min(t_pred - 1, Tn), s] - com_wmean[s]) / com_std_denom
          
          mu <- drift + xterm + (if (com_in_mean) summ$gamma * com_std else 0)
          y_pred[s] <- y_pred[s] + mu
        }
      }
      
      errs <- c(errs, Yz[t + hh - 1, ] - y_pred)
    }
    
    list(
      h = hh,
      RMSE = sqrt(mean(errs^2, na.rm = TRUE)),
      MAE = mean(abs(errs), na.rm = TRUE)
    )
  })
  
  names(res) <- paste0("h", horizons)
  res
}


#' Validate OU model fit
#'
#' Prints diagnostic summaries and hypothesis tests for a fitted model.
#'
#' @param fit_res List returned by \code{\link{fit_ou_nonlinear_tmg}}
#' @param verbose Logical. Print detailed output. Default TRUE.
#'
#' @return Invisibly returns the diagnostics list
#'
#' @examples
#' \donttest{
#' # Create a dummy results list that mimics the output of fit_ou_nonlinear_tmg
#' dummy_results <- list(
#'   diagnostics = list(
#'     rhat = c(alpha = 1.01, beta = 1.00),
#'     ess = c(alpha = 400, beta = 350),
#'     loo = list(estimates = matrix(c(1, 0.1), ncol=2, 
#'                dimnames=list("elpd_loo", c("Estimate", "SE")))),
#'     oos = list(h1 = list(RMSE = 0.5))
#'   ),
#'   factor_ou = list(beta1 = 0.3),
#'   nonlinear = list(a3 = -0.5),
#'   sv = list(rho_s = 0.2)
#' )
#' 
#' # Run validation on the dummy object
#' validate_ou_fit(dummy_results)
#' }
#'
#' @export
validate_ou_fit <- function(fit_res, verbose = TRUE) {
  dg <- fit_res$diagnostics
  
  if (verbose) {
    message("\n==== NONLINEAR OU MODEL VALIDATION ====")
    message("R-hat (summary):")
    print(utils::head(dg$rhat))
    message("\nESS (summary):")
    print(utils::head(dg$ess))
    
    if (!is.null(dg$loo)) {
      message("\nPSIS-LOO:")
      print(dg$loo)
    }
    
    message("\nOOS metrics:")
    print(dg$oos)
    
    beta1 <- fit_res$factor_ou$beta1
    message(sprintf(
      "\nH1: beta1 sign - median point estimate: %s",
      if (beta1 > 0) ">0" else "<=0"
    ))
    
    a3_med <- fit_res$nonlinear$a3
    message(sprintf(
      "H4: a3<0 (increasing restoring force): %.2f proportion of sectors",
      mean(a3_med < 0, na.rm = TRUE)
    ))
    
    rho <- fit_res$sv$rho_s
    message(sprintf("H6: median rho_s: %.3f", stats::median(rho, na.rm = TRUE)))
  }
  
  invisible(dg)
}


#' Compare models using PSIS-LOO
#'
#' Compares two fitted models using PSIS-LOO cross-validation.
#'
#' @param results_new List. Results from new model.
#' @param results_base List. Results from base model.
#'
#' @return List with components:
#'   \describe{
#'     \item{loo_table}{Comparison table from loo::loo_compare}
#'     \item{deltaELPD}{Numeric difference in ELPD}
#'   }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   # 1. Create mock 'loo' objects manually to avoid computation errors
#'   # Structure required by loo_compare: list with 'estimates' and class 'psis_loo'
#'   
#'   # Mock Model 1: ELPD = -100
#'   est1 <- matrix(c(-100, 5), nrow = 1, dimnames = list("elpd_loo", c("Estimate", "SE")))
#'   # Pointwise data (required for diff calculation). 10 observations.
#'   pw1 <- matrix(c(rep(-10, 10)), ncol = 1, dimnames = list(NULL, "elpd_loo"))
#'   
#'   loo_obj1 <- list(estimates = est1, pointwise = pw1)
#'   class(loo_obj1) <- c("psis_loo", "loo")
#'   
#'   # Mock Model 2: ELPD = -102 (worse)
#'   est2 <- matrix(c(-102, 5), nrow = 1, dimnames = list("elpd_loo", c("Estimate", "SE")))
#'   pw2 <- matrix(c(rep(-10.2, 10)), ncol = 1, dimnames = list(NULL, "elpd_loo"))
#'   
#'   loo_obj2 <- list(estimates = est2, pointwise = pw2)
#'   class(loo_obj2) <- c("psis_loo", "loo")
#'   
#'   # 2. Wrap in the structure expected by your package
#'   res_new <- list(diagnostics = list(loo = loo_obj1))
#'   res_base <- list(diagnostics = list(loo = loo_obj2))
#'   
#'   # 3. Compare (This will now run cleanly without warnings)
#'   cmp <- compare_models_loo(res_new, res_base)
#'   print(cmp)
#' }
#' }
#'
#' @export
compare_models_loo <- function(results_new, results_base) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required for model comparison.")
  }
  
  loo_new <- results_new$diagnostics$loo
  loo_base <- results_base$diagnostics$loo %||% results_base$loo
  
  cmp <- loo::loo_compare(loo_new, loo_base)
  
  list(
    loo_table = cmp,
    deltaELPD = as.numeric(loo::loo_compare(loo_new, loo_base)[1, "elpd_diff"])
  )
}
