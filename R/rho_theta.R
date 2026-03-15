# =============================================================================
# rho_theta.R
# -----------------------------------------------------------------------------
# Unified wrapper for the parametric correlation family rho_theta(s, t)
# used in the covariance estimation block of Lin & Wang (2022).
#
# This function does NOT re-implement the correlation formulas. Instead, it
# dispatches to the correlation families already defined in:
#
#   R/functions/data_generation/
#     - matern_correlation.R
#     - power_exponential_correlation.R
#     - rational_quadratic_correlation.R
#
# Supported models:
#   - "matern"
#   - "power_exponential"
#   - "rational_quadratic"
#
# Convention:
#   theta = c(theta1, theta2)
#
# Output:
#   A matrix of dimension length(s) x length(t), whose (a,b)-entry is
#   rho_theta(s[a], t[b]).
# =============================================================================


#' Unified parametric correlation function rho_theta(s, t)
#'
#' @param s      Numeric vector. First set of time points.
#' @param t      Numeric vector. Second set of time points.
#' @param theta  Numeric vector of length 2: c(theta1, theta2).
#' @param model  Character string. One of:
#'               - "matern"
#'               - "power_exponential"
#'               - "rational_quadratic"
#'
#' @return Numeric matrix of dimension length(s) x length(t).
#'
#' @details
#' This is a wrapper/dispatcher used in the covariance estimation block.
#' The actual formulas are delegated to the corresponding functions already
#' implemented in the data_generation folder.
#'
#' @examples
#' s <- c(0.1, 0.3, 0.5)
#' t <- c(0.2, 0.4)
#'
#' rho_theta(s, t, theta = c(1, 0.5), model = "power_exponential")
#' rho_theta(s, t, theta = c(1, 0.5), model = "rational_quadratic")
#' rho_theta(s, t, theta = c(0.5, 1), model = "matern")

#' @export
rho_theta <- function(s, t, theta, model) {
  
  # ---------------------------------------------------------------------------
  # input checks
  # ---------------------------------------------------------------------------
  if (!is.numeric(s) || length(s) < 1 || any(is.na(s))) {
    stop("'s' must be a non-empty numeric vector without NA values.")
  }
  
  if (!is.numeric(t) || length(t) < 1 || any(is.na(t))) {
    stop("'t' must be a non-empty numeric vector without NA values.")
  }
  
  if (!is.numeric(theta) || length(theta) != 2 || any(is.na(theta))) {
    stop("'theta' must be a numeric vector of length 2: c(theta1, theta2).")
  }
  
  if (!is.character(model) || length(model) != 1) {
    stop("'model' must be a single character string.")
  }
  
  valid_models <- c("matern", "power_exponential", "rational_quadratic")
  if (!model %in% valid_models) {
    stop(sprintf(
      "'model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }
  
  theta1 <- theta[1]
  theta2 <- theta[2]
  
  if (!is.finite(theta1) || !is.finite(theta2)) {
    stop("'theta1' and 'theta2' must be finite.")
  }
  
  # ---------------------------------------------------------------------------
  # check required underlying function exists
  # ---------------------------------------------------------------------------
  required_fun <- switch(
    model,
    "matern"              = "matern_correlation",
    "power_exponential"   = "power_exponential_correlation",
    "rational_quadratic"  = "rational_quadratic_correlation"
  )
  
  if (!exists(required_fun, mode = "function")) {
    stop(sprintf(
      paste0(
        "Required function '%s' was not found in the current environment.\n",
        "Please source the corresponding file from:\n",
        "  R/functions/data_generation/\n",
        "before calling rho_theta()."
      ),
      required_fun
    ))
  }
  
  # ---------------------------------------------------------------------------
  # dispatch to the chosen family
  # ---------------------------------------------------------------------------
  rho_mat <- switch(
    model,
    
    "matern" = matern_correlation(
      s = s, t = t,
      theta1 = theta1, theta2 = theta2
    ),
    
    "power_exponential" = power_exponential_correlation(
      s = s, t = t,
      theta1 = theta1, theta2 = theta2
    ),
    
    "rational_quadratic" = rational_quadratic_correlation(
      s = s, t = t,
      theta1 = theta1, theta2 = theta2
    )
  )
  
  # ---------------------------------------------------------------------------
  # final sanity checks
  # ---------------------------------------------------------------------------
  if (!is.matrix(rho_mat)) {
    stop("Internal error: rho_theta() expected a matrix output.")
  }
  
  if (!all(dim(rho_mat) == c(length(s), length(t)))) {
    stop("Internal error: returned matrix has incorrect dimensions.")
  }
  
  if (any(!is.finite(rho_mat))) {
    stop("Internal error: rho_theta() produced non-finite values.")
  }
  
  return(rho_mat)
}
