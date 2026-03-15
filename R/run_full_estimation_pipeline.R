# =============================================================================
# run_full_estimation_pipeline.R
# -----------------------------------------------------------------------------
# Full estimation pipeline for Lin & Wang (2022):
#
#   1. bandwidth selection for mu(t)        [optional if h_mu provided]
#   2. mean estimation                      -> compute_mu_hat()
#   3. mean on observed times               -> compute_mu_on_observed_times()
#   4. squared residuals                    -> compute_residuals()
#   5. bandwidth selection for varsigma^2   [optional if h_sigma provided]
#   6. varsigma^2 estimation                -> compute_varsigma2_hat()
#   7. h0 rule                              -> compute_h0()
#   8. sigma0^2 estimation                  -> compute_sigma0()
#   9. sigmaX^2 estimation                  -> compute_sigmaX2_hat()
#  10. covariance estimation                -> estimate_covariance()
#
# This wrapper is designed to maximize user control while keeping the pipeline
# reproducible and transparent.
#
# Dependencies:
#   - bandwidth_cv_mu()
#   - compute_mu_hat()
#   - compute_mu_on_observed_times()
#   - compute_residuals()
#   - bandwidth_cv_sigma()
#   - compute_varsigma2_hat()
#   - compute_h0()
#   - compute_sigma0()
#   - compute_sigmaX2_hat()
#   - estimate_covariance()
# =============================================================================


#' Run the full estimation pipeline
#'
#' @param T_list         List of observation times per subject.
#' @param Y_list         List of observed values per subject.
#' @param kernel         Character. Kernel type. Default "epanechnikov".
#' @param scheme         Character. Weighting scheme. Default "OBS".
#' @param kappa          Integer >= 2. Number of CV folds. Default 5.
#' @param n_h            Integer >= 2. Number of bandwidth candidates. Default 20.
#' @param n_grid         Integer >= 2. Number of grid points for function estimates.
#'                       Default 100.
#' @param model          Character. Correlation model for covariance estimation.
#'                       One of "power_exponential", "rational_quadratic", "matern".
#' @param h_mu           Optional numeric scalar. If provided, skip CV for mu.
#' @param h_sigma        Optional numeric scalar. If provided, skip CV for varsigma^2.
#' @param theta_init     Optional numeric vector of length 2. Initial value for theta.
#' @param s_grid         Optional numeric vector. Grid for covariance surface.
#' @param t_grid         Optional numeric vector. Grid for covariance surface.
#' @param optim_method   Character. Optimization method for estimate_theta().
#'                       Default "BFGS".
#' @param optim_control  Optional list of control parameters for optim().
#' @param seed_global    Optional integer. Global seed set at the start.
#' @param seed_cv_mu     Optional integer. Seed for CV of mu. If NULL, uses
#'                       seed_global when available.
#' @param seed_cv_sigma  Optional integer. Seed for CV of varsigma^2. If NULL,
#'                       uses seed_global + 1 when available.
#' @param return_all     Logical. If TRUE, return the full pipeline output.
#'                       If FALSE, return a compact output.
#' @param run_checks     Logical. If TRUE, run data checks and reorder each
#'                       subject by increasing time if needed.
#' @param verbose        Logical. If TRUE, print progress information.
#'
#' @return
#' If return_all = TRUE, a list with:
#'   - inputs
#'   - seeds
#'   - cv_mu
#'   - h_mu
#'   - res_mu
#'   - mu_obs_list
#'   - Z_list
#'   - cv_sigma
#'   - h_sigma
#'   - res_varsigma
#'   - h0_res
#'   - sigma0_res
#'   - res_sigmaX2
#'   - covariance_fit
#'   - diagnostics
#'
#' If return_all = FALSE, a compact list with:
#'   - h_mu
#'   - h_sigma
#'   - sigma0_hat
#'   - theta_hat
#'   - covariance_hat
#'   - diagnostics
#'
#' @examples
#' # fit <- run_full_estimation_pipeline(
#' #   T_list = T_list,
#' #   Y_list = Y_list,
#' #   model  = "power_exponential",
#' #   seed_global = 123,
#' #   verbose = TRUE
#' # )

run_full_estimation_pipeline <- function(
    T_list,
    Y_list,
    kernel = "epanechnikov",
    scheme = "OBS",
    kappa = 5,
    n_h = 20,
    n_grid = 100,
    model = "power_exponential",
    h_mu = NULL,
    h_sigma = NULL,
    theta_init = NULL,
    s_grid = NULL,
    t_grid = NULL,
    optim_method = "BFGS",
    optim_control = NULL,
    seed_global = NULL,
    seed_cv_mu = NULL,
    seed_cv_sigma = NULL,
    return_all = TRUE,
    run_checks = TRUE,
    verbose = FALSE
) {
  
  # ---------------------------------------------------------------------------
  # helper
  # ---------------------------------------------------------------------------
  say <- function(...) {
    if (isTRUE(verbose)) cat(...)
  }
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list)) {
    stop("'T_list' and 'Y_list' must both be lists.")
  }
  
  if (length(T_list) != length(Y_list)) {
    stop("'T_list' and 'Y_list' must have the same length.")
  }
  
  if (length(T_list) == 0) {
    stop("'T_list' and 'Y_list' must not be empty.")
  }
  
  if (!is.character(kernel) || length(kernel) != 1) {
    stop("'kernel' must be a single character string.")
  }
  
  if (!is.character(scheme) || length(scheme) != 1) {
    stop("'scheme' must be a single character string.")
  }
  
  if (!is.numeric(kappa) || length(kappa) != 1 || kappa < 2) {
    stop("'kappa' must be a single number >= 2.")
  }
  
  if (!is.numeric(n_h) || length(n_h) != 1 || n_h < 2) {
    stop("'n_h' must be a single number >= 2.")
  }
  
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid < 2) {
    stop("'n_grid' must be a single number >= 2.")
  }
  
  valid_models <- c("power_exponential", "rational_quadratic", "matern")
  if (!is.character(model) || length(model) != 1 || !model %in% valid_models) {
    stop(sprintf(
      "'model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }
  
  if (!is.null(h_mu) && (!is.numeric(h_mu) || length(h_mu) != 1 || h_mu <= 0)) {
    stop("'h_mu' must be NULL or a strictly positive scalar.")
  }
  
  if (!is.null(h_sigma) && (!is.numeric(h_sigma) || length(h_sigma) != 1 || h_sigma <= 0)) {
    stop("'h_sigma' must be NULL or a strictly positive scalar.")
  }
  
  if (!is.null(theta_init)) {
    if (!is.numeric(theta_init) || length(theta_init) != 2 || any(is.na(theta_init))) {
      stop("'theta_init' must be NULL or a numeric vector of length 2.")
    }
  }
  
  if (!is.null(seed_global) && (!is.numeric(seed_global) || length(seed_global) != 1)) {
    stop("'seed_global' must be NULL or a single numeric value.")
  }
  if (!is.null(seed_cv_mu) && (!is.numeric(seed_cv_mu) || length(seed_cv_mu) != 1)) {
    stop("'seed_cv_mu' must be NULL or a single numeric value.")
  }
  if (!is.null(seed_cv_sigma) && (!is.numeric(seed_cv_sigma) || length(seed_cv_sigma) != 1)) {
    stop("'seed_cv_sigma' must be NULL or a single numeric value.")
  }
  
  if (!is.logical(return_all) || length(return_all) != 1) {
    stop("'return_all' must be TRUE or FALSE.")
  }
  
  if (!is.logical(run_checks) || length(run_checks) != 1) {
    stop("'run_checks' must be TRUE or FALSE.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. seed logic
  # ---------------------------------------------------------------------------
  if (!is.null(seed_global)) {
    set.seed(seed_global)
  }
  
  if (is.null(seed_cv_mu) && !is.null(seed_global)) {
    seed_cv_mu <- as.integer(seed_global)
  }
  
  if (is.null(seed_cv_sigma) && !is.null(seed_global)) {
    seed_cv_sigma <- as.integer(seed_global) + 1L
  }
  
  # ---------------------------------------------------------------------------
  # 3. data checks / ordering
  # ---------------------------------------------------------------------------
  T_in <- T_list
  Y_in <- Y_list
  
  if (run_checks) {
    for (i in seq_along(T_in)) {
      Ti <- T_in[[i]]
      Yi <- Y_in[[i]]
      
      if (!is.numeric(Ti) || !is.numeric(Yi)) {
        stop(sprintf("T_list[[%d]] and Y_list[[%d]] must be numeric vectors.", i, i))
      }
      if (length(Ti) != length(Yi)) {
        stop(sprintf("T_list[[%d]] and Y_list[[%d]] must have the same length.", i, i))
      }
      if (length(Ti) == 0) {
        stop(sprintf("Subject %d has zero observations.", i))
      }
      if (any(is.na(Ti)) || any(is.na(Yi))) {
        stop(sprintf("T_list[[%d]] and Y_list[[%d]] must not contain NA values.", i, i))
      }
      
      # reorder by increasing time if needed
      if (length(Ti) > 1 && any(diff(Ti) < 0)) {
        ord <- order(Ti)
        T_in[[i]] <- Ti[ord]
        Y_in[[i]] <- Yi[ord]
        say(sprintf("Subject %d reordered by increasing time.\n", i))
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # 4. bandwidth for mu
  # ---------------------------------------------------------------------------
  cv_mu <- NULL
  
  if (is.null(h_mu)) {
    say("--------------------------------------------------\n")
    say("Step 1/10 - Selecting bandwidth for mu(t)\n")
    say("--------------------------------------------------\n")
    
    if (is.null(seed_cv_mu)) {
      stop("No seed available for bandwidth_cv_mu(). Provide 'seed_cv_mu' or 'seed_global'.")
    }
    
    cv_mu <- bandwidth_cv_mu(
      T_list = T_in,
      Y_list = Y_in,
      kappa  = kappa,
      n_h    = n_h,
      kernel = kernel,
      scheme = scheme,
      seed   = seed_cv_mu
    )
    
    h_mu <- cv_mu$h_opt
  } else {
    say("Skipping CV for mu(t): using provided h_mu.\n")
  }
  
  # ---------------------------------------------------------------------------
  # 5. mean estimation
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 2/10 - Estimating mu(t)\n")
  say("--------------------------------------------------\n")
  
  res_mu <- compute_mu_hat(
    T_list = T_in,
    Y_list = Y_in,
    h_mu   = h_mu,
    kernel = kernel,
    scheme = scheme,
    n_t    = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # 6. mean on observed times
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 3/10 - Evaluating mu_hat on observed times\n")
  say("--------------------------------------------------\n")
  
  mu_obs_list <- compute_mu_on_observed_times(
    T_list = T_in,
    Y_list = Y_in,
    h_mu   = h_mu,
    kernel = kernel,
    scheme = scheme
  )
  
  # ---------------------------------------------------------------------------
  # 7. squared residuals
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 4/10 - Computing squared residuals\n")
  say("--------------------------------------------------\n")
  
  Z_list <- compute_residuals(
    T_list = T_in,
    Y_list = Y_in,
    h_mu   = h_mu,
    kernel = kernel,
    scheme = scheme
  )
  
  # ---------------------------------------------------------------------------
  # 8. bandwidth for varsigma^2
  # ---------------------------------------------------------------------------
  cv_sigma <- NULL
  
  if (is.null(h_sigma)) {
    say("--------------------------------------------------\n")
    say("Step 5/10 - Selecting bandwidth for varsigma^2(t)\n")
    say("--------------------------------------------------\n")
    
    if (is.null(seed_cv_sigma)) {
      stop("No seed available for bandwidth_cv_sigma(). Provide 'seed_cv_sigma' or 'seed_global'.")
    }
    
    cv_sigma <- bandwidth_cv_sigma(
      T_list = T_in,
      Z_list = Z_list,
      kappa  = kappa,
      n_h    = n_h,
      kernel = kernel,
      scheme = scheme,
      seed   = seed_cv_sigma
    )
    
    h_sigma <- cv_sigma$h_opt
  } else {
    say("Skipping CV for varsigma^2(t): using provided h_sigma.\n")
  }
  
  # ---------------------------------------------------------------------------
  # 9. varsigma^2 estimation
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 6/10 - Estimating varsigma^2(t)\n")
  say("--------------------------------------------------\n")
  
  res_varsigma <- compute_varsigma2_hat(
    T_list  = T_in,
    Z_list  = Z_list,
    h_sigma = h_sigma,
    kernel  = kernel,
    scheme  = scheme,
    n_t     = n_grid
  )
  
  # ---------------------------------------------------------------------------
  # 10. h0 rule
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 7/10 - Computing h0\n")
  say("--------------------------------------------------\n")
  
  h0_res <- compute_h0(
    T_list       = T_in,
    res_varsigma = res_varsigma
  )
  
  # ---------------------------------------------------------------------------
  # 11. sigma0^2 estimation
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 8/10 - Estimating sigma0^2\n")
  say("--------------------------------------------------\n")
  
  sigma0_res <- compute_sigma0(
    T_list = T_in,
    Y_list = Y_in,
    h0     = h0_res$h0
  )
  
  # ---------------------------------------------------------------------------
  # 12. sigmaX^2 estimation
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 9/10 - Estimating sigmaX^2(t)\n")
  say("--------------------------------------------------\n")
  
  res_sigmaX2 <- compute_sigmaX2_hat(
    res_varsigma = res_varsigma,
    sigma0_hat   = sigma0_res$sigma0_hat
  )
  
  # ---------------------------------------------------------------------------
  # 13. covariance estimation
  # ---------------------------------------------------------------------------
  say("--------------------------------------------------\n")
  say("Step 10/10 - Estimating covariance surface\n")
  say("--------------------------------------------------\n")
  
  covariance_fit <- estimate_covariance(
    T_list       = T_in,
    Y_list       = Y_in,
    mu_obs_list  = mu_obs_list,
    res_sigmaX2  = res_sigmaX2,
    model        = model,
    theta_init   = theta_init,
    s_grid       = s_grid,
    t_grid       = t_grid,
    rule         = 2,
    optim_method = optim_method,
    optim_control = optim_control,
    return_optim = TRUE,
    verbose      = verbose
  )
  
  # ---------------------------------------------------------------------------
  # 14. diagnostics
  # ---------------------------------------------------------------------------
  diag_matches_sigmaX2 <- NULL
  is_cov_symmetric <- NULL
  
  cov_hat <- covariance_fit$covariance_hat
  
  if (is.matrix(cov_hat) &&
      length(covariance_fit$s_grid) == length(covariance_fit$t_grid) &&
      isTRUE(all.equal(covariance_fit$s_grid, covariance_fit$t_grid, tolerance = 1e-12))) {
    
    is_cov_symmetric <- isTRUE(all.equal(cov_hat, t(cov_hat), tolerance = 1e-8))
    diag_matches_sigmaX2 <- isTRUE(all.equal(
      diag(cov_hat),
      approx(
        x    = res_sigmaX2$t_grid,
        y    = res_sigmaX2$sigmaX2_hat,
        xout = covariance_fit$s_grid,
        rule = 2
      )$y,
      tolerance = 1e-8
    ))
  }
  
  diagnostics <- list(
    sigma0_finite         = is.finite(sigma0_res$sigma0_hat),
    theta_finite          = all(is.finite(covariance_fit$theta_hat)),
    covariance_finite     = all(is.finite(cov_hat)),
    n_truncated_sigmaX2   = res_sigmaX2$n_truncated,
    is_cov_symmetric      = is_cov_symmetric,
    diag_matches_sigmaX2  = diag_matches_sigmaX2
  )
  
  # ---------------------------------------------------------------------------
  # 15. output
  # ---------------------------------------------------------------------------
  if (!return_all) {
    return(list(
      h_mu           = h_mu,
      h_sigma        = h_sigma,
      sigma0_hat     = sigma0_res$sigma0_hat,
      theta_hat      = covariance_fit$theta_hat,
      covariance_hat = covariance_fit$covariance_hat,
      diagnostics    = diagnostics
    ))
  }
  
  list(
    inputs = list(
      kernel        = kernel,
      scheme        = scheme,
      kappa         = kappa,
      n_h           = n_h,
      n_grid        = n_grid,
      model         = model,
      h_mu          = h_mu,
      h_sigma       = h_sigma,
      theta_init    = theta_init,
      s_grid        = s_grid,
      t_grid        = t_grid,
      optim_method  = optim_method,
      optim_control = optim_control
    ),
    seeds = list(
      seed_global   = seed_global,
      seed_cv_mu    = seed_cv_mu,
      seed_cv_sigma = seed_cv_sigma
    ),
    cv_mu          = cv_mu,
    h_mu           = h_mu,
    res_mu         = res_mu,
    mu_obs_list    = mu_obs_list,
    Z_list         = Z_list,
    cv_sigma       = cv_sigma,
    h_sigma        = h_sigma,
    res_varsigma   = res_varsigma,
    h0_res         = h0_res,
    sigma0_res     = sigma0_res,
    res_sigmaX2    = res_sigmaX2,
    covariance_fit = covariance_fit,
    diagnostics    = diagnostics
  )
}
