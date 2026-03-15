# =============================================================================
# matern_correlation.R
# -----------------------------------------------------------------------------
# Implements the Matérn correlation function as defined in equation (4) of
# Lin & Wang (2022):
#
#   rho_theta(s, t) = 1 / (Gamma(theta1) * 2^{theta1-1}) *
#                     (sqrt(2*theta1) * |s-t| / theta2)^theta1 *
#                     B_{theta1}(sqrt(2*theta1) * |s-t| / theta2)
#
# where B_{theta1}(.) is the modified Bessel function of the second kind
# of order theta1 (implemented in R as besselK(..., nu = theta1)).
#
# Special case: rho_theta(t, t) = 1  (by continuity, as |s-t| -> 0)
#
# For Cov I in the simulation study, the paper uses theta = (0.5, 1).
# With theta1 = 0.5 this reduces to the exponential correlation:
#   rho(s,t) = exp(-|s-t| / theta2)
#
# Dependencies: base R only (besselK is in base R)
# =============================================================================


#' Matérn correlation function (Lin & Wang 2022, equation 4)
#'
#' @param s       Numeric vector. First set of time points.
#' @param t       Numeric vector. Second set of time points.
#'                The function computes the full outer matrix rho(s_i, t_j).
#' @param theta1  Numeric > 0. Smoothness parameter.
#'                theta1 = 0.5 -> exponential (rough process)
#'                theta1 = 1.5 -> once differentiable
#'                theta1 = 2.5 -> twice differentiable
#'                theta1 -> Inf -> squared exponential (infinitely smooth)
#' @param theta2  Numeric > 0. Range (scale) parameter. Larger theta2 means
#'                stronger correlation at the same distance.
#'
#' @return Numeric matrix of dimension length(s) x length(t).
#'         Entry (i,j) is rho_theta(s_i, t_j).
#'
#' @examples
#' # Cov I from Lin & Wang (2022): theta = (0.5, 1)
#' s <- seq(0, 1, length.out = 5)
#' t <- seq(0, 1, length.out = 5)
#' rho <- matern_correlation(s, t, theta1 = 0.5, theta2 = 1)
#'
#' # Check diagonal = 1
#' all(diag(rho) == 1)  # TRUE
#'
#' # With theta1=0.5 reduces to exp(-|s-t|/theta2)
#' d <- abs(outer(s, t, "-"))
#' max(abs(rho - exp(-d)))  # ~0 (up to numerical precision)

matern_correlation <- function(s, t, theta1, theta2) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(s) || !is.numeric(t))
    stop("s and t must be numeric vectors.")
  if (!is.numeric(theta1) || length(theta1) != 1 || theta1 <= 0)
    stop("theta1 must be a single positive number.")
  if (!is.numeric(theta2) || length(theta2) != 1 || theta2 <= 0)
    stop("theta2 must be a single positive number.")
  
  # --- compute distance matrix ------------------------------------------------
  d <- abs(outer(s, t, "-"))   # |s_i - t_j|
  
  # --- scaled distance --------------------------------------------------------
  # u = sqrt(2*theta1) * |s-t| / theta2
  u <- sqrt(2 * theta1) * d / theta2
  
  # --- Matérn formula (eq. 4) -------------------------------------------------
  # At u=0 (i.e. s=t): rho = 1 by continuity
  # At u>0: rho = 1/(Gamma(theta1)*2^{theta1-1}) * u^theta1 * B_{theta1}(u)
  rho <- ifelse(
    u == 0,
    1,
    (1 / (gamma(theta1) * 2^(theta1 - 1))) *
      u^theta1 *
      besselK(u, nu = theta1)
  )
  
  return(rho)
}
