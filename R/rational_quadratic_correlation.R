# =============================================================================
# rational_quadratic_correlation.R
# -----------------------------------------------------------------------------
# Implements the Rational Quadratic (Cauchy) correlation function as defined
# in Lin & Wang (2022), Section 2:
#
#   rho_theta(s, t) = { 1 + |s-t|^2 / theta2^2 }^{-theta1}
#
# with theta1, theta2 > 0.
#
# Interpretation:
#   - theta2  controls the range (scale): larger theta2 -> slower decay
#   - theta1  controls the tail heaviness: larger theta1 -> faster decay
#              (as theta1 -> Inf, approaches the squared exponential)
#
# Note: unlike the Power Exponential and Matérn families, the Rational
# Quadratic has polynomial (heavy) tails in the distance, meaning
# correlations decay more slowly for large |s-t|.
#
# Dependencies: base R only
# =============================================================================


#' Rational Quadratic (Cauchy) correlation function (Lin & Wang 2022)
#'
#' @param s       Numeric vector. First set of time points.
#' @param t       Numeric vector. Second set of time points.
#'                The function computes the full outer matrix rho(s_i, t_j).
#' @param theta1  Numeric > 0. Tail decay parameter. Larger values give
#'                faster decay and lighter tails.
#' @param theta2  Numeric > 0. Scale (range) parameter.
#'
#' @return Numeric matrix of dimension length(s) x length(t).
#'         Entry (i,j) is rho_theta(s_i, t_j).
#'
#' @examples
#' s <- seq(0, 1, length.out = 5)
#' t <- seq(0, 1, length.out = 5)
#' rho <- rational_quadratic_correlation(s, t, theta1 = 1, theta2 = 0.5)

#' @export
rational_quadratic_correlation <- function(s, t, theta1, theta2) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(s) || !is.numeric(t))
    stop("s and t must be numeric vectors.")
  if (!is.numeric(theta1) || length(theta1) != 1 || theta1 <= 0)
    stop("theta1 must be a single positive number.")
  if (!is.numeric(theta2) || length(theta2) != 1 || theta2 <= 0)
    stop("theta2 must be a single positive number.")
  
  # --- compute distance matrix ------------------------------------------------
  d <- abs(outer(s, t, "-"))   # |s_i - t_j|
  
  # --- Rational Quadratic formula ---------------------------------------------
  rho <- (1 + d^2 / theta2^2)^(-theta1)
  
  return(rho)
}
