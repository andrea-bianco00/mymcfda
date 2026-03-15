# =============================================================================
# power_exponential_correlation.R
# -----------------------------------------------------------------------------
# Implements the Power Exponential correlation function as defined in
# Lin & Wang (2022), Section 2:
#
#   rho_theta(s, t) = exp{ -|s-t|^theta1 / theta2^theta1 }
#
# with 0 < theta1 <= 2, theta2 > 0.
#
# Special cases:
#   theta1 = 1  ->  exponential:        rho(s,t) = exp(-|s-t| / theta2)
#   theta1 = 2  ->  squared exponential (Gaussian): rho(s,t) = exp(-|s-t|^2 / theta2^2)
#
# Note: theta1 controls smoothness of the process.
#   theta1 < 2  -> process is NOT mean-square differentiable (rough)
#   theta1 = 2  -> process is infinitely mean-square differentiable (smooth)
#
# Dependencies: base R only
# =============================================================================


#' Power Exponential correlation function (Lin & Wang 2022)
#'
#' @param s       Numeric vector. First set of time points.
#' @param t       Numeric vector. Second set of time points.
#'                The function computes the full outer matrix rho(s_i, t_j).
#' @param theta1  Numeric in (0, 2). Shape parameter.
#'                theta1 = 1  -> exponential
#'                theta1 = 2  -> squared exponential (Gaussian)
#' @param theta2  Numeric > 0. Scale (range) parameter.
#'
#' @return Numeric matrix of dimension length(s) x length(t).
#'         Entry (i,j) is rho_theta(s_i, t_j).
#'
#' @examples
#' s <- seq(0, 1, length.out = 5)
#' t <- seq(0, 1, length.out = 5)
#'
#' # Exponential (theta1=1)
#' rho <- power_exponential_correlation(s, t, theta1 = 1, theta2 = 0.5)
#'
#' # Squared exponential (theta1=2)
#' rho <- power_exponential_correlation(s, t, theta1 = 2, theta2 = 0.5)

#' @export
power_exponential_correlation <- function(s, t, theta1, theta2) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(s) || !is.numeric(t))
    stop("s and t must be numeric vectors.")
  if (!is.numeric(theta1) || length(theta1) != 1 || theta1 <= 0 || theta1 > 2)
    stop("theta1 must be a single number in (0, 2].")
  if (!is.numeric(theta2) || length(theta2) != 1 || theta2 <= 0)
    stop("theta2 must be a single positive number.")
  
  # --- compute distance matrix ------------------------------------------------
  d <- abs(outer(s, t, "-"))   # |s_i - t_j|
  
  # --- Power Exponential formula ----------------------------------------------
  rho <- exp(-(d^theta1) / (theta2^theta1))
  
  return(rho)
}
