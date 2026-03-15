# =============================================================================
# density_fzero_truncexp.R
# -----------------------------------------------------------------------------
# Defines the Truncated Exponential density for the measurement times f_0,
# as described in Section 1.6 of Lin & Wang (2022), Assumption 3:
#
#   f_0(s) = lambda * exp(-lambda * s) / (1 - exp(-lambda * delta)),
#                                               s in [0, delta]
#   f_0(s) = 0,                                 otherwise
#
# ROLE OF f_0 IN THE MODEL (Section 1.6, Lin & Wang 2022):
#   f_0 is defined on the reference interval [0, delta]. The conditional
#   density of T_ij given O_i = u is obtained by shifting f_0 to the
#   subject's window [u - delta/2, u + delta/2]:
#
#     f_{T|O}(t | u) = f_0(t - u + delta/2)
#
#   So f_0 describes HOW measurements are spread within a subject's window.
#   With the Truncated Exponential, measurements tend to cluster towards
#   the beginning of each subject's window, with probability decaying
#   exponentially towards the end.
#
# STRICT POSITIVITY (Assumption 3, Lin & Wang 2022):
#   exp(-lambda * s) > 0 for all s, so f_0(s) > 0 on all of [0, delta].
#   Satisfied for any lambda > 0.
#
# LIPSCHITZ CONTINUITY OF f_0' (Assumption 3, Lin & Wang 2022):
#   f_0'(s) = -lambda^2 * exp(-lambda*s) / (1 - exp(-lambda*delta)),
#   which is bounded and continuous on [0, delta], hence Lipschitz.
#
# DEFAULT PARAMETER:
#   lambda = 1 / delta
#   Moderate decay proportional to the support length, consistent across
#   all values of delta. As lambda -> 0, f_0 converges to the Uniform.
#
# RELATIONSHIP WITH density_fO_truncexp.R:
#   Same mathematical form, different support:
#     density_fO_truncexp.R  : support [delta/2, 1-delta/2], shift lo=delta/2
#     density_fzero_truncexp.R: support [0, delta],          shift lo=0
#
# This function is called by:
#   eval_fzero.R   — to evaluate f_0(s) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Truncated Exponential density for measurement times f_0
#'
#' @param s      Numeric vector. Points at which to evaluate f_0.
#'               Values outside [0, delta] return 0.
#' @param delta  Numeric in (0, 1). Snippet length. Determines the support
#'               [0, delta].
#' @param lambda Numeric > 0. Decay rate parameter.
#'                 small lambda  -> nearly uniform (slow decay)
#'                 large lambda  -> concentrated near s = 0 (fast decay)
#'               Default: 1 / delta, i.e. the reciprocal of the support
#'               length.
#'
#' @return Numeric vector of the same length as s with values of f_0(s).
#'
#' @examples
#' s <- seq(-0.1, 0.6, length.out = 500)
#'
#' # Default lambda = 1/delta = 2.5
#' f <- fzero_truncexp(s, delta = 0.4)
#' plot(s, f, type = "l", main = "f_0 truncexp, delta = 0.4")
#'
#' # Fast decay: concentrated near s = 0
#' f2 <- fzero_truncexp(s, delta = 0.4, lambda = 10)
#' lines(s, f2, col = "red")
#'
#' # Slow decay: nearly uniform
#' f3 <- fzero_truncexp(s, delta = 0.4, lambda = 0.1)
#' lines(s, f3, col = "blue")
#'
#' # Verify it integrates to 1
#' s_fine <- seq(0, 0.4, length.out = 10000)
#' sum(fzero_truncexp(s_fine, delta = 0.4)) * (0.4 / 10000)  # approximately 1

fzero_truncexp <- function(s, delta, lambda = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s))
    stop("s must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  
  # ---------------------------------------------------------------------------
  # DEFAULT PARAMETER
  # lambda = 1 / delta — moderate decay proportional to support length
  # ---------------------------------------------------------------------------
  if (is.null(lambda)) lambda <- 1 / delta
  
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("lambda must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # NORMALISATION CONSTANT
  # Z = integral_0^delta exp(-lambda*s) ds = (1 - exp(-lambda*delta)) / lambda
  # ---------------------------------------------------------------------------
  Z <- (1 - exp(-lambda * delta)) / lambda
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # f_0(s) = exp(-lambda * s) / Z   inside [0, delta]
  # f_0(s) = 0                      outside
  # ---------------------------------------------------------------------------
  inside <- (s >= 0 & s <= delta)
  
  f <- numeric(length(s))
  f[inside] <- exp(-lambda * s[inside]) / Z
  
  return(f)
}
