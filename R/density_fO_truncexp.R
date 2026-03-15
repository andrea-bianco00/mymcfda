# =============================================================================
# density_fO_truncexp.R
# -----------------------------------------------------------------------------
# Defines the Truncated Exponential density for the reference times f_O,
# as described in Section 1.5 of Lin & Wang (2022):
#
#   "Truncated Exponential. f_O(u) proportional to exp(-lambda*(u - delta/2))
#    for u in [delta/2, 1 - delta/2], with lambda > 0. This models a higher
#    recruitment rate early in the study window, decaying over time."
#
# Explicitly:
#
#   f_O(u) = lambda * exp(-lambda*(u - lo)) / (1 - exp(-lambda*(hi - lo))),
#                                               u in [lo, hi]
#   f_O(u) = 0,                                 otherwise
#
# where lo = delta/2, hi = 1 - delta/2.
#
# INTERPRETATION:
#   lambda > 0 controls the decay rate of recruitment over time:
#     small lambda  -> nearly uniform (slow decay, spread recruitment)
#     large lambda  -> concentrated near lo (fast decay, early recruitment)
#   As lambda -> 0, f_O converges to the Uniform density (fO_uniform).
#
# STRICT POSITIVITY (Assumption 2, Lin & Wang 2022):
#   exp(-lambda*(u - lo)) > 0 for all u, so f_O(u) > 0 on [lo, hi].
#   Satisfied for any lambda > 0.
#
# DEFAULT PARAMETER:
#   lambda = 1 / (hi - lo) = 1 / (1 - delta)
#   This gives a moderate decay proportional to the support length,
#   shifting the mean slightly towards lo relative to the Uniform.
#   Consistent across all values of delta in (0, 1).
#
# This function is called by:
#   eval_fO.R    — to evaluate f_O(u) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Truncated Exponential density for reference times f_O
#'
#' @param u      Numeric vector. Points at which to evaluate f_O.
#'               Values outside (delta/2, 1 - delta/2) return 0.
#' @param delta  Numeric in (0, 1). Snippet length. Determines the support
#'               (delta/2, 1 - delta/2).
#' @param lambda Numeric > 0. Decay rate parameter. Controls how quickly
#'               the recruitment probability decreases over time:
#'                 small lambda  -> nearly uniform
#'                 large lambda  -> concentrated near delta/2
#'               Default: 1 / (1 - delta), i.e. the reciprocal of the
#'               support length.
#'
#' @return Numeric vector of the same length as u with values of f_O(u).
#'
#' @examples
#' u <- seq(0, 1, length.out = 500)
#'
#' # Default lambda
#' f <- fO_truncexp(u, delta = 0.4)
#' plot(u, f, type = "l", main = "f_O truncexp, delta = 0.4")
#'
#' # Fast decay: heavy concentration near lo = 0.2
#' f2 <- fO_truncexp(u, delta = 0.4, lambda = 10)
#' lines(u, f2, col = "red")
#'
#' # Slow decay: nearly uniform
#' f3 <- fO_truncexp(u, delta = 0.4, lambda = 0.1)
#' lines(u, f3, col = "blue")
#'
#' # Verify it integrates to 1
#' sum(f) * (1 / 500)   # approximately 1

#' @export
fO_truncexp <- function(u, delta, lambda = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(u))
    stop("u must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  
  # ---------------------------------------------------------------------------
  # SUPPORT
  # ---------------------------------------------------------------------------
  lo <- delta / 2
  hi <- 1 - delta / 2
  
  # ---------------------------------------------------------------------------
  # DEFAULT PARAMETER
  # lambda = 1 / (hi - lo) = 1 / (1 - delta)
  # Moderate decay, consistent across all values of delta.
  # ---------------------------------------------------------------------------
  if (is.null(lambda)) lambda <- 1 / (hi - lo)
  
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("lambda must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # NORMALISATION CONSTANT
  # Z = integral_{lo}^{hi} exp(-lambda*(u-lo)) du = (1 - exp(-lambda*(hi-lo))) / lambda
  # ---------------------------------------------------------------------------
  Z <- (1 - exp(-lambda * (hi - lo))) / lambda
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # f_O(u) = exp(-lambda*(u - lo)) / Z   inside support
  # f_O(u) = 0                           outside support
  # ---------------------------------------------------------------------------
  inside <- (u >= lo & u <= hi)
  
  f <- numeric(length(u))
  f[inside] <- exp(-lambda * (u[inside] - lo)) / Z
  
  return(f)
}
