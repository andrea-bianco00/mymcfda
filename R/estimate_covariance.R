# =============================================================================
# estimate_covariance.R
# -----------------------------------------------------------------------------
# Full wrapper for Block 5: covariance estimation.
#
# Pipeline:
#   1. evaluate_sigmaX_on_observed_times()
#   2. compute_raw_covariances()
#   3. estimate_theta()
#   4. compute_covariance_surface()
#
# Final output:
#   C_hat(s,t) = sigmaX_hat(s) * rho_{theta_hat}(s,t) * sigmaX_hat(t)
#
# Required inputs:
#   - T_list
#   - Y_list
#   - mu_obs_list
#   - res_sigmaX2
#   - model
#
# Dependencies:
#   - evaluate_sigmaX_on_observed_times()
#   - compute_raw_covariances()
#   - estimate_theta()
#   - compute_covariance_surface()
# =============================================================================


#' Estimate the covariance surface (full Block 5 wrapper)
#'
#' @param T_list        List of length n. Observation times per subject.
#' @param Y_list        List of length n. Observed values per subject.
#' @param mu_obs_list   List of length n. mu_hat(T_ij) for each subject.
#' @param res_sigmaX2   List returned by compute_sigmaX2_hat().
#' @param model         Character string. One of:
#'                      - "power_exponential"
#'                      - "rational_quadratic"
#'                      - "matern"
#' @param theta_init    Optional numeric vector of length 2. Initial value for theta.
#' @param s_grid        Optional numeric vector of s-values for covariance surface.
#'                      If NULL, use res_sigmaX2$t_grid.
#' @param t_grid        Optional numeric vector of t-values for covariance surface.
#'                      If NULL, use res_sigmaX2$t_grid.
#' @param rule          Integer passed to approx() in interpolation.
#'                      Default is 2.
#' @param optim_method  Optimization method passed to estimate_theta().
#'                      Default is "BFGS".
#' @param optim_control Optional list of control parameters for optim().
#' @param return_optim  Logical. If TRUE, include full optim() output.
#' @param verbose       Logical. If TRUE, print progress information.
#'
#' @return A list with:
#'   - model
#'   - theta_hat
#'   - Qn_hat
#'   - sigmaX_obs_list
#'   - raw_cov_df
#'   - fit_theta
#'   - covariance_surface
#'   - s_grid
#'   - t_grid
#'   - covariance_hat
#'
#' @examples
#' # cov_fit <- estimate_covariance(
#' #   T_list      = T_list,
#' #   Y_list      = Y_list,
#' #   mu_obs_list = mu_obs_list,
#' #   res_sigmaX2 = res_sigmaX2,
#' #   model       = "power_exponential"
#' # )
#' # cov_fit$theta_hat
#' # cov_fit$covariance_hat

#' @export
estimate_covariance <- function(T_list,
                                Y_list,
                                mu_obs_list,
                                res_sigmaX2,
                                model,
                                theta_init   = NULL,
                                s_grid       = NULL,
                                t_grid       = NULL,
                                rule         = 2,
                                optim_method = "BFGS",
                                optim_control = NULL,
                                return_optim = TRUE,
                                verbose      = FALSE) {
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list) || !is.list(mu_obs_list)) {
    stop("'T_list', 'Y_list', and 'mu_obs_list' must all be lists.")
  }
  
  if (length(T_list) != length(Y_list) || length(T_list) != length(mu_obs_list)) {
    stop("'T_list', 'Y_list', and 'mu_obs_list' must all have the same length.")
  }
  
  n <- length(T_list)
  
  for (i in seq_len(n)) {
    Ti  <- T_list[[i]]
    Yi  <- Y_list[[i]]
    mui <- mu_obs_list[[i]]
    
    if (!is.numeric(Ti) || !is.numeric(Yi) || !is.numeric(mui)) {
      stop(sprintf(
        "T_list[[%d]], Y_list[[%d]], and mu_obs_list[[%d]] must be numeric vectors.",
        i, i, i
      ))
    }
    
    if (!(length(Ti) == length(Yi) && length(Ti) == length(mui))) {
      stop(sprintf(
        "T_list[[%d]], Y_list[[%d]], and mu_obs_list[[%d]] must have the same length.",
        i, i, i
      ))
    }
    
    if (any(is.na(Ti)) || any(is.na(Yi)) || any(is.na(mui))) {
      stop(sprintf(
        "T_list[[%d]], Y_list[[%d]], and mu_obs_list[[%d]] must not contain NA values.",
        i, i, i
      ))
    }
  }
  
  if (!is.list(res_sigmaX2)) {
    stop("'res_sigmaX2' must be a list.")
  }
  
  required_res_names <- c("t_grid", "sigmaX2_hat")
  if (!all(required_res_names %in% names(res_sigmaX2))) {
    stop("'res_sigmaX2' must contain at least 't_grid' and 'sigmaX2_hat'.")
  }
  
  if (!is.character(model) || length(model) != 1) {
    stop("'model' must be a single character string.")
  }
  
  valid_models <- c("power_exponential", "rational_quadratic", "matern")
  if (!model %in% valid_models) {
    stop(sprintf(
      "'model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }
  
  # ---------------------------------------------------------------------------
  # 2. sigmaX_hat evaluated at observed times
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("--------------------------------------------------\n")
    cat("Step 1/4 - Evaluating sigmaX_hat on observed times\n")
    cat("--------------------------------------------------\n")
  }
  
  sigmaX_obs_list <- evaluate_sigmaX_on_observed_times(
    T_list      = T_list,
    res_sigmaX2 = res_sigmaX2,
    rule        = rule
  )
  
  # ---------------------------------------------------------------------------
  # 3. raw off-diagonal covariance pseudo-observations
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("--------------------------------------------------\n")
    cat("Step 2/4 - Computing raw covariance pseudo-observations\n")
    cat("--------------------------------------------------\n")
  }
  
  raw_cov_df <- compute_raw_covariances(
    T_list      = T_list,
    Y_list      = Y_list,
    mu_obs_list = mu_obs_list
  )
  
  # ---------------------------------------------------------------------------
  # 4. estimate theta
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("--------------------------------------------------\n")
    cat("Step 3/4 - Estimating theta\n")
    cat("--------------------------------------------------\n")
  }
  
  fit_theta <- estimate_theta(
    raw_cov_df      = raw_cov_df,
    sigmaX_obs_list = sigmaX_obs_list,
    model           = model,
    theta_init      = theta_init,
    method          = optim_method,
    control         = optim_control,
    return_optim    = return_optim,
    verbose         = verbose
  )
  
  # ---------------------------------------------------------------------------
  # 5. covariance surface
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("--------------------------------------------------\n")
    cat("Step 4/4 - Computing covariance surface\n")
    cat("--------------------------------------------------\n")
  }
  
  cov_surface <- compute_covariance_surface(
    res_sigmaX2 = res_sigmaX2,
    theta       = fit_theta$theta_hat,
    model       = model,
    s_grid      = s_grid,
    t_grid      = t_grid,
    rule        = rule
  )
  
  # ---------------------------------------------------------------------------
  # 6. output
  # ---------------------------------------------------------------------------
  out <- list(
    model              = model,
    theta_hat          = fit_theta$theta_hat,
    Qn_hat             = fit_theta$Qn_hat,
    sigmaX_obs_list    = sigmaX_obs_list,
    raw_cov_df         = raw_cov_df,
    fit_theta          = fit_theta,
    covariance_surface = cov_surface,
    s_grid             = cov_surface$s_grid,
    t_grid             = cov_surface$t_grid,
    covariance_hat     = cov_surface$covariance_hat
  )
  
  return(out)
}
