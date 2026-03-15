# =============================================================================
# estimate_theta.R
# -----------------------------------------------------------------------------
# Estimate the parameter vector theta by minimizing Q_n(theta) using optim().
#
# Supported models:
#   - "power_exponential"
#   - "rational_quadratic"
#   - "matern"
#
# This function uses an unconstrained reparameterization eta -> theta so that
# optim() can run safely in R^2 while respecting the model constraints:
#
#   power_exponential:
#     theta1 in (0, 2], theta2 > 0
#
#   rational_quadratic:
#     theta1 > 0, theta2 > 0
#
#   matern:
#     theta1 > 0, theta2 > 0
#
# Dependencies:
#   - compute_Qn_theta()  defined in compute_Qn_theta.R
# =============================================================================


#' Estimate theta by minimizing Q_n(theta)
#'
#' @param raw_cov_df       Data frame returned by compute_raw_covariances().
#' @param sigmaX_obs_list  List of sigmaX_hat values at observed times.
#' @param model            Character string. One of:
#'                         - "power_exponential"
#'                         - "rational_quadratic"
#'                         - "matern"
#' @param theta_init       Optional numeric vector of length 2 giving the
#'                         initial value for theta. If NULL, a default is used.
#' @param method           Optimization method passed to optim().
#'                         Default: "BFGS".
#' @param control          Optional list of control parameters for optim().
#' @param return_optim     Logical. If TRUE, return the full optim() object.
#' @param verbose          Logical. If TRUE, print progress information.
#'
#' @return A list with:
#'   - theta_hat        Estimated parameter vector c(theta1, theta2)
#'   - Qn_hat           Minimum value of Q_n(theta)
#'   - model            Model name
#'   - theta_init       Initial theta used
#'   - eta_init         Initial unconstrained parameter used by optim()
#'   - eta_hat          Final unconstrained parameter
#'   - convergence      optim() convergence code
#'   - message          optim() message, if available
#'   - counts           Function/gradient evaluation counts, if available
#'   - optim_result     Full optim() output if return_optim = TRUE, else NULL

estimate_theta <- function(raw_cov_df,
                           sigmaX_obs_list,
                           model,
                           theta_init   = NULL,
                           method       = "BFGS",
                           control      = NULL,
                           return_optim = TRUE,
                           verbose      = FALSE) {
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.data.frame(raw_cov_df)) {
    stop("'raw_cov_df' must be a data.frame.")
  }
  
  if (!is.list(sigmaX_obs_list)) {
    stop("'sigmaX_obs_list' must be a list.")
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
  
  if (!is.null(theta_init)) {
    if (!is.numeric(theta_init) || length(theta_init) != 2 || any(is.na(theta_init))) {
      stop("'theta_init' must be NULL or a numeric vector of length 2.")
    }
  }
  
  if (!is.character(method) || length(method) != 1) {
    stop("'method' must be a single character string.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. default initialization
  # ---------------------------------------------------------------------------
  if (is.null(theta_init)) {
    theta_init <- switch(
      model,
      "power_exponential"  = c(1.0, 0.2),
      "rational_quadratic" = c(1.0, 0.2),
      "matern"             = c(1.0, 0.2)
    )
  }
  
  # ---------------------------------------------------------------------------
  # 3. validate theta_init according to model
  # ---------------------------------------------------------------------------
  if (model == "power_exponential") {
    if (theta_init[1] <= 0 || theta_init[1] > 2) {
      stop("For 'power_exponential', theta_init[1] must be in (0, 2].")
    }
    if (theta_init[2] <= 0) {
      stop("For 'power_exponential', theta_init[2] must be > 0.")
    }
  } else {
    if (any(theta_init <= 0)) {
      stop(sprintf("For '%s', both entries of theta_init must be > 0.", model))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 4. transformation helpers: eta <-> theta
  # ---------------------------------------------------------------------------
  logistic <- function(x) 1 / (1 + exp(-x))
  
  theta_from_eta <- function(eta, model) {
    if (!is.numeric(eta) || length(eta) != 2 || any(!is.finite(eta))) {
      stop("'eta' must be a finite numeric vector of length 2.")
    }
    
    if (model == "power_exponential") {
      theta1 <- 2 * logistic(eta[1])   # in (0, 2)
      theta2 <- exp(eta[2])            # > 0
    } else {
      theta1 <- exp(eta[1])            # > 0
      theta2 <- exp(eta[2])            # > 0
    }
    
    c(theta1, theta2)
  }
  
  eta_from_theta <- function(theta, model) {
    if (!is.numeric(theta) || length(theta) != 2 || any(!is.finite(theta))) {
      stop("'theta' must be a finite numeric vector of length 2.")
    }
    
    if (model == "power_exponential") {
      # theta1 = 2 * logistic(eta1)
      # => p = theta1 / 2 in (0,1), eta1 = log(p / (1-p))
      p <- theta[1] / 2
      
      # numerical safeguard
      eps <- 1e-8
      p <- min(max(p, eps), 1 - eps)
      
      eta1 <- log(p / (1 - p))
      eta2 <- log(theta[2])
    } else {
      eta1 <- log(theta[1])
      eta2 <- log(theta[2])
    }
    
    c(eta1, eta2)
  }
  
  eta_init <- eta_from_theta(theta_init, model = model)
  
  # ---------------------------------------------------------------------------
  # 5. objective in eta-space
  # ---------------------------------------------------------------------------
  objective_eta <- function(eta) {
    theta <- theta_from_eta(eta, model = model)
    
    qn_val <- compute_Qn_theta(
      raw_cov_df      = raw_cov_df,
      sigmaX_obs_list = sigmaX_obs_list,
      theta           = theta,
      model           = model,
      return_details  = FALSE
    )
    
    if (!is.finite(qn_val)) {
      return(.Machine$double.xmax / 1000)
    }
    
    qn_val
  }
  
  if (verbose) {
    cat("--------------------------------------------------\n")
    cat("Estimating theta\n")
    cat("Model       :", model, "\n")
    cat("theta_init  :", paste(signif(theta_init, 6), collapse = ", "), "\n")
    cat("eta_init    :", paste(signif(eta_init,   6), collapse = ", "), "\n")
    cat("Qn(theta_init) =", objective_eta(eta_init), "\n")
    cat("--------------------------------------------------\n")
  }
  
  # ---------------------------------------------------------------------------
  # 6. optimization
  # ---------------------------------------------------------------------------
  optim_res <- optim(
    par     = eta_init,
    fn      = objective_eta,
    method  = method,
    control = control
  )
  
  eta_hat   <- optim_res$par
  theta_hat <- theta_from_eta(eta_hat, model = model)
  Qn_hat    <- optim_res$value
  
  if (verbose) {
    cat("Optimization completed.\n")
    cat("theta_hat   :", paste(signif(theta_hat, 6), collapse = ", "), "\n")
    cat("Qn_hat      :", signif(Qn_hat, 6), "\n")
    cat("convergence :", optim_res$convergence, "\n")
    if (!is.null(optim_res$message)) {
      cat("message     :", optim_res$message, "\n")
    }
  }
  
  # ---------------------------------------------------------------------------
  # 7. output
  # ---------------------------------------------------------------------------
  out <- list(
    theta_hat    = theta_hat,
    Qn_hat       = Qn_hat,
    model        = model,
    theta_init   = theta_init,
    eta_init     = eta_init,
    eta_hat      = eta_hat,
    convergence  = optim_res$convergence,
    message      = if (!is.null(optim_res$message)) optim_res$message else NULL,
    counts       = if (!is.null(optim_res$counts)) optim_res$counts else NULL,
    optim_result = if (isTRUE(return_optim)) optim_res else NULL
  )
  
  return(out)
}
