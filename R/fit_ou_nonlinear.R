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
#' @param priors List. Prior specifications (partial override allowed; missing
#'   entries fall back to robust defaults). Supported names: \code{sigma_delta}
#'   (wedge SD, original units, default 0.002); \code{beta1_mean}/\code{beta1_sd}
#'   (prior for the global TMG effect, default 0 / 0.5 -- neutral, no sign baked
#'   in); \code{nu_shape}/\code{nu_rate} (gamma prior for the Student-t degrees of
#'   freedom shift, default 2 / 0.1 -> weakly informative, prior mean nu ~ 22);
#'   \code{rho_mean}/\code{rho_sd} (SV persistence prior, default 0.7 / 0.2).
#' @param com_in_mean Logical. Include COM effect in mean equation. Default TRUE.
#' @param train_frac Numeric in (0,1). Fraction of observations used for the
#'   training window. Default 0.70.
#' @param fit_window Character. \code{"train"} (default) fits the likelihood and
#'   computes \code{log_lik} ONLY on the training window (2:T_train), so the
#'   test window is genuinely held out and \code{evaluate_oos} is a real
#'   out-of-sample evaluation. \code{"full"} fits on all observations
#'   (full-information / rstanarm-style); then PSIS-LOO is valid over all points
#'   but \code{evaluate_oos} becomes in-sample. See the README methodology note.
#' @param chains Integer. Number of MCMC chains. Default 6.
#' @param iter Integer. Total iterations per chain (must exceed warmup). Default 12000.
#' @param warmup Integer. Warmup iterations. Default 6000.
#' @param thin Integer. Thinning interval. Default 1 (thinning is statistically
#'   wasteful in HMC; increase only to save memory).
#' @param cores Integer. Number of cores for parallel chains.
#' @param threads_per_chain Integer. Threads per chain for within-chain
#'   parallelism. Capped internally so parallel_chains * threads <= cores.
#' @param hard_sum_zero Logical. If TRUE, TMG wedge is fixed at zero. Default TRUE.
#' @param orthogonalize_tmg Logical. Orthogonalize TMG w.r.t. common factor. Default TRUE.
#' @param factor_from Character. Source for common factor: "X" or "Y". Default "X".
#' @param use_train_loadings Logical. Compute factor loadings from training only
#'   (avoids look-ahead leakage). Default TRUE.
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
#' The model uses a non-centered hierarchical (partial-pooling) parameterization
#' for the sector-specific parameters \code{theta_s}, \code{kappa_s},
#' \code{a3_s} and \code{beta0_s}: each is built as
#' \code{hyper_intercept + hyper_slope * COM_s + hyper_sd * z}, with
#' \code{z ~ N(0,1)}. The training window is \code{floor(T * train_frac)}.
#' All data standardization, the COM weighting, and (by default) the common
#' factor loadings are computed from the training window only. The likelihood
#' window is controlled by \code{fit_window} to keep train/test coherent.
#'
#' The cubic drift enforces \code{a3 < 0} (a restoring force that strengthens
#' with the deviation); this is a stability \emph{assumption}, not an estimated
#' result, and it precludes detecting locally expansive (self-amplifying)
#' regimes. The Euler discretization uses dt = 1, so the discrete persistence of
#' the linear part is \code{1 - kappa}; interpret half-lives accordingly.
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
    priors = list(),
    com_in_mean = TRUE,
    train_frac = 0.70,
    fit_window = c("train", "full"),
    chains = 6,
    iter = 12000,
    warmup = 6000,
    thin = 1,
    cores = max(1, parallel::detectCores() - 1),
    threads_per_chain = 2,
    hard_sum_zero = TRUE,
    orthogonalize_tmg = TRUE,
    factor_from = c("X", "Y"),
    use_train_loadings = TRUE,
    adapt_delta = 0.97,
    max_treedepth = 12,
    seed = 1234,
    init = NULL,
    moment_match = NULL,
    verbose = FALSE
) {

  factor_from <- match.arg(factor_from)
  model <- match.arg(model)
  fit_window <- match.arg(fit_window)

  # Merge user priors over robust defaults (partial override allowed).
  prior_defaults <- list(
    sigma_delta = 0.002,
    beta1_mean  = 0,      beta1_sd = 0.5,   # neutral on the TMG effect H1
    nu_shape    = 2,      nu_rate  = 0.1,   # weakly-informative df (mean ~22)
    rho_mean    = 0.7,    rho_sd   = 0.2    # not the rigid N(0.90, 0.05)
  )
  priors <- utils::modifyList(prior_defaults, priors %||% list())

  stopifnot(is.matrix(Y) || is.data.frame(Y))
  stopifnot(is.matrix(X) || is.data.frame(X))
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  stopifnot(nrow(Y) == nrow(X), ncol(Y) == ncol(X))
  stopifnot(length(TMG) == nrow(Y))

  # ---- Hard input validation (fail early with a clear message) ----
  if (!all(is.finite(Y))) stop("`Y` contains non-finite values (NA/NaN/Inf).", call. = FALSE)
  if (!all(is.finite(X))) stop("`X` contains non-finite values (NA/NaN/Inf).", call. = FALSE)
  if (!all(is.finite(TMG))) stop("`TMG` contains non-finite values.", call. = FALSE)
  if (!is.numeric(train_frac) || train_frac <= 0 || train_frac >= 1) {
    stop("`train_frac` must be in (0, 1).", call. = FALSE)
  }
  if (warmup < 1L || iter <= warmup) {
    stop(sprintf("`iter` (%d) must be strictly greater than `warmup` (%d).",
                 as.integer(iter), as.integer(warmup)), call. = FALSE)
  }
  if (chains < 1L) stop("`chains` must be >= 1.", call. = FALSE)
  if (thin < 1L) stop("`thin` must be >= 1.", call. = FALSE)
  
  stopifnot(is.matrix(COM) || is.data.frame(COM))
  stopifnot(is.matrix(CAPITAL_TOTAL) || is.data.frame(CAPITAL_TOTAL))
  COM_ts <- as.matrix(COM)
  K_ts <- as.matrix(CAPITAL_TOTAL)
  
  Tn <- nrow(Y)
  S <- ncol(Y)
  T_train <- max(2L, floor(Tn * train_frac))
  if (T_train >= Tn) {
    stop("`train_frac` leaves no observations for the test window; lower it.",
         call. = FALSE)
  }
  # Likelihood/log-lik window: train-only (honest forecasting) vs full sample.
  T_lik <- if (fit_window == "train") T_train else Tn

  vmsg(sprintf("Data dimensions: T=%d, S=%d, T_train=%d, fit_window=%s (T_lik=%d)",
               Tn, S, T_train, fit_window, T_lik), verbose)
  
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
    mu_xz = rep(0.0, S),
    T_lik = as.integer(T_lik),
    beta1_prior_mean = as.numeric(priors$beta1_mean),
    beta1_prior_sd = as.numeric(priors$beta1_sd),
    nu_prior_shape = as.numeric(priors$nu_shape),
    nu_prior_rate = as.numeric(priors$nu_rate),
    rho_prior_mean = as.numeric(priors$rho_mean),
    rho_prior_sd = as.numeric(priors$rho_sd)
  )

  fit <- NULL
  backend <- check_stan_backend(verbose)

  if (backend == "none") {
    stop("Stan backend required. Please install cmdstanr or rstan.")
  }

  # Avoid thread oversubscription: parallel_chains * threads_per_chain <= cores.
  par_chains <- max(1L, min(as.integer(chains), as.integer(cores)))
  thr_per_chain <- max(1L, min(as.integer(threads_per_chain),
                               as.integer(floor(cores / par_chains))))

  if (backend == "cmdstanr") {
    vmsg("Compiling Stan model with cmdstanr (from canonical .stan file)", verbose)
    # Compile into a writable user cache (NOT the package/library directory,
    # which may be read-only, and NOT the source tree). cmdstanr reuses the
    # cached binary across sessions when the model is unchanged.
    cache_dir <- tryCatch({
      d <- tools::R_user_dir("bayesianOU", "cache")
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(d)) d else tempdir()
    }, error = function(e) tempdir())
    mod <- cmdstanr::cmdstan_model(
      .stan_file_path(),
      dir = cache_dir,
      pedantic = FALSE,
      cpp_options = list(stan_threads = TRUE)
    )

    vmsg("Running MCMC sampling", verbose)
    fit <- mod$sample(
      data = stan_dat,
      chains = chains,
      parallel_chains = par_chains,
      iter_warmup = warmup,
      iter_sampling = iter - warmup,
      thin = thin,
      seed = seed,
      refresh = if (verbose) 200 else 0,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      threads_per_chain = thr_per_chain,
      init = init
    )
  } else {
    vmsg("Compiling Stan model with rstan", verbose)
    sm <- rstan::stan_model(model_code = ou_nonlinear_tmg_stan_code())
    
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
  summ <- extract_posterior_summary(fit)
  rhat_vec <- as.numeric(summ$rhat)
  rhat_max <- max(rhat_vec, na.rm = TRUE)
  rhat_share <- mean(rhat_vec > 1.01, na.rm = TRUE)

  if (isTRUE(moment_match)) {
    vmsg(paste("Note: moment_match is not applied to the array-based PSIS-LOO.",
               "Use loo::loo_moment_match() on the returned stan_fit if needed."),
         verbose)
  }

  vmsg("Computing PSIS-LOO over the fitted window (2:T_lik)", verbose)
  loo_res <- NULL
  loo_pareto_k_summary <- NULL
  if (requireNamespace("loo", quietly = TRUE)) {
    loo_res <- tryCatch(
      .compute_loo(fit, T_lik, S),
      error = function(e) {
        warning("PSIS-LOO could not be computed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(loo_res)) loo_pareto_k_summary <- .summarize_pareto_k(loo_res)
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
  
  out <- results_robust %||% list()
  # Direct assignment (do NOT c()-append: that produced duplicated names on
  # repeated calls). The full draws array is intentionally NOT duplicated here;
  # access draws via the returned stan_fit. To persist across sessions with
  # cmdstanr, call results$factor_ou$stan_fit$save_object(file).
  out$factor_ou <- list(
    model = "ou_nonlinear_tmg",
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
      T_lik = T_lik,
      fit_window = fit_window,
      com_in_mean = isTRUE(com_in_mean),
      factor_from = factor_from,
      use_train_loadings = isTRUE(use_train_loadings)
    )
  )

  out$beta_tmg <- build_beta_tmg_table(fit, zTMG_use, summ = summ)
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
    loo_pareto_k = loo_pareto_k_summary,
    oos = oos
  )
  
  class(out$factor_ou) <- c("ou_nonlinear_tmg", "list")
  
  vmsg("Model fitting complete", verbose)
  out
}
