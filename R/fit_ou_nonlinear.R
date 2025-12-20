#' Fit Bayesian nonlinear OU model with TMG effect and SV
#'
#' Fits a Bayesian nonlinear Ornstein-Uhlenbeck model with cubic drift,
#' stochastic volatility, and Student-t innovations using Stan.
#'
#' @param results_robust List. Previous results object to extend (can be empty list).
#' @param Y Numeric matrix (T x S). Dependent variable (prices/values by sector).
#' @param X Numeric matrix (T x S). Independent variable (production prices).
#' @param TMG Numeric vector (length T). Aggregate TMG series.
#' @param COM Numeric matrix (T x S). Composition of capital by sector.
#' @param CAPITAL_TOTAL Numeric matrix (T x S). Total capital by sector.
#' @param model Character. Model type. Currently only "base" supported.
#' @param priors List. Prior specifications. Currently supports sigma_delta.
#' @param com_in_mean Logical. Include COM effect in mean equation. Default TRUE.
#' @param chains Integer. Number of MCMC chains. Default 6.
#' @param iter Integer. Total iterations per chain. Default 12000.
#' @param warmup Integer. Warmup iterations. Default 6000.
#' @param thin Integer. Thinning interval. Default 2.
#' @param cores Integer. Number of cores for parallel chains.
#' @param threads_per_chain Integer. Threads per chain for within-chain parallelism.
#' @param hard_sum_zero Logical. If TRUE, TMG wedge is fixed at zero. Default TRUE.
#' @param orthogonalize_tmg Logical. Orthogonalize TMG w.r.t. common factor. Default TRUE.
#' @param factor_from Character. Source for common factor: "X" or "Y". Default "X".
#' @param use_train_loadings Logical. Compute factor loadings from training only. Default FALSE.
#' @param adapt_delta Numeric. Target acceptance rate (0-1). Default 0.97.
#' @param max_treedepth Integer. Maximum tree depth for NUTS. Default 12.
#' @param seed Integer. Random seed for reproducibility.
#' @param init Numeric or function. Initial values for parameters.
#' @param moment_match Logical. Use moment matching for LOO. Default NULL.
#' @param verbose Logical. Print progress messages. Default FALSE.
#'
#' @return List containing:
#'   \describe{
#'     \item{factor_ou}{Model results including draws and parameter estimates}
#'     \item{beta_tmg}{Time-varying beta estimates}
#'     \item{sv}{Stochastic volatility summaries}
#'     \item{nonlinear}{Nonlinearity diagnostics}
#'     \item{accounting}{TMG accounting block}
#'     \item{diagnostics}{MCMC diagnostics, LOO, and OOS metrics}
#'   }
#'
#' @details
#' The model uses hierarchical priors for sector-specific parameters.
#' Training period is set to 70 percent of observations by default.
#' All data standardization uses training period statistics only.
#'
#' @examples
#' \donttest{
#' # 1. Prepare dummy data
#' T_obs <- 20
#' S_sectors <- 2
#' Y <- matrix(rnorm(T_obs * S_sectors), nrow = T_obs, ncol = S_sectors)
#' X <- matrix(rnorm(T_obs * S_sectors), nrow = T_obs, ncol = S_sectors)
#' TMG <- rnorm(T_obs)
#' COM <- matrix(runif(T_obs * S_sectors), nrow = T_obs, ncol = S_sectors)
#' K <- matrix(runif(T_obs * S_sectors, 100, 1000), nrow = T_obs, ncol = S_sectors)
#'
#' # 2. Run model (conditional on Stan backend availability)
#' # We use very short chains just to demonstrate execution
#' if (requireNamespace("cmdstanr", quietly = TRUE) || 
#'     requireNamespace("rstan", quietly = TRUE)) {
#'   
#'   # Wrap in try to avoid failure if Stan is not configured locally
#'   try({
#'     results <- fit_ou_nonlinear_tmg(
#'       results_robust = list(),
#'       Y = Y, X = X, TMG = TMG, COM = COM, CAPITAL_TOTAL = K,
#'       chains = 1, iter = 100, warmup = 50, # Short run for example
#'       verbose = FALSE
#'     )
#'   }, silent = TRUE)
#' }
#' }
#'
#' @export
fit_ou_nonlinear_tmg <- function(
    results_robust,
    Y, X,
    TMG,
    COM,
    CAPITAL_TOTAL,
    model = c("base"),
    priors = list(sigma_delta = 0.002),
    com_in_mean = TRUE,
    chains = 6,
    iter = 12000,
    warmup = 6000,
    thin = 2,
    cores = max(1, parallel::detectCores() - 1),
    threads_per_chain = 2,
    hard_sum_zero = TRUE,
    orthogonalize_tmg = TRUE,
    factor_from = c("X", "Y"),
    use_train_loadings = FALSE,
    adapt_delta = 0.97,
    max_treedepth = 12,
    seed = 1234,
    init = NULL,
    moment_match = NULL,
    verbose = FALSE
) {
  
  old_wd <- getwd()
  temp_dir <- tempdir()
  setwd(temp_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  factor_from <- match.arg(factor_from)
  model <- match.arg(model)
  
  stopifnot(is.matrix(Y) || is.data.frame(Y))
  stopifnot(is.matrix(X) || is.data.frame(X))
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  stopifnot(nrow(Y) == nrow(X), ncol(Y) == ncol(X))
  stopifnot(length(TMG) == nrow(Y))
  
  stopifnot(is.matrix(COM) || is.data.frame(COM))
  stopifnot(is.matrix(CAPITAL_TOTAL) || is.data.frame(CAPITAL_TOTAL))
  COM_ts <- as.matrix(COM)
  K_ts <- as.matrix(CAPITAL_TOTAL)
  
  Tn <- nrow(Y)
  S <- ncol(Y)
  T_train <- max(2L, floor(Tn * 0.70))
  
  vmsg(sprintf("Data dimensions: T=%d, S=%d, T_train=%d", Tn, S, T_train), verbose)
  
  if (!is.null(colnames(Y)) && !is.null(colnames(COM_ts))) {
    common <- intersect(colnames(Y), colnames(COM_ts))
    if (length(common) != S) {
      missing <- setdiff(colnames(Y), colnames(COM_ts))
      extra <- setdiff(colnames(COM_ts), colnames(Y))
      stop(sprintf(
        "Column mismatch COM vs Y. Missing in COM: %s. Extra in COM: %s",
        if (length(missing) == 0) "(none)" else paste(missing, collapse = ", "),
        if (length(extra) == 0) "(none)" else paste(extra, collapse = ", ")
      ))
    }
    COM_ts <- COM_ts[, colnames(Y), drop = FALSE]
  } else if (ncol(COM_ts) != S) {
    stop("Dimension mismatch COM vs Y.")
  }
  
  if (!is.null(colnames(Y)) && !is.null(colnames(K_ts))) {
    commonK <- intersect(colnames(Y), colnames(K_ts))
    if (length(commonK) != S) {
      missing <- setdiff(colnames(Y), colnames(K_ts))
      extra <- setdiff(colnames(K_ts), colnames(Y))
      stop(sprintf(
        "Column mismatch CAPITAL_TOTAL vs Y. Missing in K: %s. Extra in K: %s",
        if (length(missing) == 0) "(none)" else paste(missing, collapse = ", "),
        if (length(extra) == 0) "(none)" else paste(extra, collapse = ", ")
      ))
    }
    K_ts <- K_ts[, colnames(Y), drop = FALSE]
  } else if (ncol(K_ts) != S) {
    stop("Dimension mismatch CAPITAL_TOTAL vs Y.")
  }
  
  vmsg("Standardizing data using training period statistics", verbose)
  zY <- zscore_train(Y, T_train)
  zX <- zscore_train(X, T_train)
  
  mu_tmg <- mean(TMG[seq_len(T_train)], na.rm = TRUE)
  sd_tmg <- stats::sd(TMG[seq_len(T_train)], na.rm = TRUE)
  if (!is.finite(sd_tmg) || sd_tmg < 1e-8) sd_tmg <- 1
  zTMG <- (TMG - mu_tmg) / sd_tmg
  
  vmsg(sprintf("Computing common factor from %s", factor_from), verbose)
  Mz_factor <- if (factor_from == "X") zX$Mz else zY$Mz
  Ft <- compute_common_factor(Mz_factor, T_train, use_train_loadings, verbose)
  
  if (orthogonalize_tmg) {
    vmsg("Orthogonalizing TMG with respect to common factor", verbose)
    fit_t <- stats::lm(zTMG[seq_len(T_train)] ~ Ft[seq_len(T_train)])
    zTMG_ortho <- zTMG - cbind(1, Ft) %*% stats::coef(fit_t)
    zTMG_use <- as.numeric(zTMG_ortho)
  } else {
    zTMG_use <- zTMG
  }
  
  sigma_delta_z <- priors$sigma_delta / sd_tmg
  soft_wedge <- as.integer(!hard_sum_zero)
  
  stan_dat <- list(
    T = Tn,
    S = S,
    T_train = T_train,
    Yz = zY$Mz,
    Xz = zX$Mz,
    zTMG_byK = as.vector(zTMG_use),
    zTMG_exo = as.vector(zTMG),
    soft_wedge = soft_wedge,
    sigma_delta_z = sigma_delta_z,
    COM_ts = COM_ts,
    K_ts = K_ts,
    com_in_mean = as.integer(isTRUE(com_in_mean)),
    mu_xz = rep(0.0, S)
  )
  
  stan_src <- ou_nonlinear_tmg_stan_code()
  fit <- NULL
  backend <- check_stan_backend(verbose)
  
  if (backend == "none") {
    stop("Stan backend required. Please install cmdstanr or rstan.")
  }
  
  if (backend == "cmdstanr") {
    vmsg("Compiling Stan model with cmdstanr", verbose)
    tf <- cmdstanr::write_stan_file(stan_src)
    mod <- cmdstanr::cmdstan_model(
      tf,
      pedantic = FALSE,
      cpp_options = list(stan_threads = TRUE)
    )
    
    vmsg("Running MCMC sampling", verbose)
    fit <- mod$sample(
      data = stan_dat,
      chains = chains,
      parallel_chains = min(chains, cores),
      iter_warmup = warmup,
      iter_sampling = iter - warmup,
      thin = thin,
      seed = seed,
      refresh = if (verbose) 200 else 0,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      threads_per_chain = threads_per_chain,
      init = init
    )
  } else {
    vmsg("Compiling Stan model with rstan", verbose)
    sm <- rstan::stan_model(model_code = stan_src)
    
    vmsg("Running MCMC sampling", verbose)
    fit <- rstan::sampling(
      sm,
      data = stan_dat,
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      seed = seed,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      refresh = if (verbose) 200 else 0,
      init = init
    )
  }
  
  vmsg("Extracting posterior summaries", verbose)
  posterior <- if (inherits(fit, "CmdStanMCMC")) {
    fit$draws()
  } else {
    rstan::As.mcmc.list(fit)
  }
  
  summ <- extract_posterior_summary(fit)
  rhat_vec <- as.numeric(summ$rhat)
  rhat_max <- max(rhat_vec, na.rm = TRUE)
  rhat_share <- mean(rhat_vec > 1.01, na.rm = TRUE)
  
  vmsg("Computing PSIS-LOO", verbose)
  if (inherits(fit, "CmdStanMCMC")) {
    log_lik <- fit$draws("log_lik", format = "matrix")
    loglik_arr <- array(log_lik, dim = c(nrow(log_lik), Tn, S))
  } else {
    loglik_arr <- rstan::extract(fit, pars = "log_lik")[[1]]
  }
  loglik_train <- loglik_arr[, 2:T_train, , drop = FALSE]
  
  loo_res <- NULL
  if (requireNamespace("loo", quietly = TRUE)) {
    args <- list(loglik_train)
    if (!is.null(moment_match)) args$moment_match <- moment_match
    loo_res <- do.call(loo::loo, args)
  }
  
  vmsg("Computing out-of-sample metrics", verbose)
  oos <- evaluate_oos(
    summ, zY$Mz, zX$Mz, zTMG_use, T_train,
    COM_ts = COM_ts,
    K_ts = K_ts,
    com_in_mean = isTRUE(com_in_mean),
    horizons = c(1, 4, 8)
  )
  
  vmsg("Computing divergence count", verbose)
  dv <- count_divergences(fit)
  
  out <- results_robust
  out$factor_ou <- c(
    out$factor_ou %||% list(),
    list(
      model = "ou_nonlinear_tmg",
      draws = posterior,
      stan_fit = fit,
      beta1 = summ$beta1,
      beta0_s = summ$beta0_s,
      kappa_s = summ$kappa_s,
      a3_s = summ$a3_s,
      theta_s = summ$theta_s,
      sv = list(
        alpha = summ$alpha_s,
        rho = summ$rho_s,
        sigma_eta = summ$sigma_eta_s
      ),
      nu = summ$nu,
      gamma = summ$gamma,
      factor_ou_info = list(
        T_train = T_train,
        com_in_mean = isTRUE(com_in_mean),
        factor_from = factor_from,
        use_train_loadings = isTRUE(use_train_loadings)
      )
    )
  )
  
  out$beta_tmg <- build_beta_tmg_table(fit, zTMG_use)
  out$sv <- list(h_summary = summarize_sv_sigmas(fit), rho_s = summ$rho_s)
  out$nonlinear <- list(
    a3 = summ$a3_s,
    drift_decomp = drift_decomposition_grid(fit, summ)
  )
  out$accounting <- build_accounting_block(
    TMG, zTMG, zTMG_use, mu_tmg, sd_tmg,
    hard_sum_zero, priors$sigma_delta
  )
  out$diagnostics <- list(
    rhat = summ$rhat,
    ess = summ$ess,
    rhat_max = rhat_max,
    rhat_share = rhat_share,
    divergences = dv,
    loo = loo_res,
    oos = oos
  )
  
  class(out$factor_ou) <- c("ou_nonlinear_tmg", "list")
  
  vmsg("Model fitting complete", verbose)
  out
}
