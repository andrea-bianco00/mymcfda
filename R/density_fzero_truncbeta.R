# =============================================================================
# density_fzero_truncbeta.R
# -----------------------------------------------------------------------------
# Defines the Rescaled Beta density for the measurement times f_0, as
# described in Section 1.6 of Lin & Wang (2022), Assumption 3:
#
#   If X ~ Beta(a_0, b_0) on [0, 1], then S = delta * X has density:
#
#   f_0(s) = dbeta(s / delta, a_0, b_0) / delta,   s in [0, delta]
#   f_0(s) = 0,                                     otherwise
#
# ROLE OF f_0 IN THE MODEL (Section 1.6, Lin & Wang 2022):
#   f_0 is defined on the reference interval [0, delta]. The conditional
#   density of T_ij given O_i = u is obtained by shifting f_0 to the
#   subject's window [u - delta/2, u + delta/2]:
#
#     f_{T|O}(t | u) = f_0(t - u + delta/2)
#
#   So f_0 describes HOW measurements are spread within a subject's window.
#   With the Rescaled Beta, measurements can be concentrated towards the
#   beginning, end, or centre of each subject's window.
#
# SHAPE INTERPRETATION:
#   a_0 = b_0        -> symmetric around delta/2
#   a_0 > b_0        -> concentrated towards the start of the window
#   a_0 < b_0        -> concentrated towards the end of the window
#   a_0 = b_0 = 1    -> reduces to the Uniform density (fzero_uniform)
#
# STRICT POSITIVITY (Assumption 3, Lin & Wang 2022):
#   Beta(a_0, b_0) is strictly positive on (0,1) only when a_0 > 1 AND
#   b_0 > 1. A warning is raised if this condition is violated.
#
# LIPSCHITZ CONTINUITY OF f_0' (Assumption 3, Lin & Wang 2022):
#   The derivative of the rescaled Beta density is Lipschitz continuous
#   on [0, delta] when a_0 > 1 and b_0 > 1, since the Beta density is
#   then continuously differentiable on the open interval (0, 1).
#
# DEFAULT PARAMETERS:
#   a_0 = b_0 = 2  -> symmetric, bell-shaped, strictly positive at
#                     boundaries. Simplest choice satisfying a_0 > 1, b_0 > 1.
#
# RELATIONSHIP WITH density_fO_truncbeta.R:
#   Same mathematical form, different support:
#     density_fO_truncbeta.R  : support [delta/2, 1-delta/2]
#     density_fzero_truncbeta.R: support [0, delta]
#
# This function is called by:
#   eval_fzero.R   — to evaluate f_0(s) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Rescaled Beta density for measurement times f_0
#'
#' @param s     Numeric vector. Points at which to evaluate f_0.
#'              Values outside [0, delta] return 0.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [0, delta].
#' @param a_0   Numeric > 0. First shape parameter of Beta(a_0, b_0).
#'              Must be > 1 for strict positivity at s = 0.
#'              Default: 2.
#' @param b_0   Numeric > 0. Second shape parameter of Beta(a_0, b_0).
#'              Must be > 1 for strict positivity at s = delta.
#'              Default: 2.
#'
#' @return Numeric vector of the same length as s with values of f_0(s).
#'
#' @examples
#' s <- seq(-0.1, 0.6, length.out = 500)
#'
#' # Default: symmetric Beta(2,2)
#' f <- fzero_truncbeta(s, delta = 0.4)
#' plot(s, f, type = "l", main = "f_0 truncbeta, delta = 0.4")
#'
#' # Concentrated towards the start of the window
#' f2 <- fzero_truncbeta(s, delta = 0.4, a_0 = 3, b_0 = 2)
#' lines(s, f2, col = "red")
#'
#' # Concentrated towards the end of the window
#' f3 <- fzero_truncbeta(s, delta = 0.4, a_0 = 2, b_0 = 3)
#' lines(s, f3, col = "blue")
#'
#' # Verify it integrates to 1
#' s_fine <- seq(0, 0.4, length.out = 10000)
#' sum(fzero_truncbeta(s_fine, delta = 0.4)) * (0.4 / 10000)  # approximately 1

fzero_truncbeta <- function(s, delta, a_0 = 2, b_0 = 2) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s))
    stop("s must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  if (!is.numeric(a_0) || length(a_0) != 1 || a_0 <= 0)
    stop("a_0 must be a single positive number.")
  if (!is.numeric(b_0) || length(b_0) != 1 || b_0 <= 0)
    stop("b_0 must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # STRICT POSITIVITY WARNING
  # Required by Assumption 3 of Lin & Wang (2022).
  # ---------------------------------------------------------------------------
  if (a_0 <= 1 || b_0 <= 1)
    warning("a_0 <= 1 or b_0 <= 1: f_0 may be zero or infinite at the ",
            "boundaries of [0, delta], violating the strict positivity ",
            "condition of Assumption 3 in Lin & Wang (2022). ",
            "Use a_0 > 1 and b_0 > 1.")
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # Rescaling: if X ~ Beta(a_0, b_0) on [0,1] and S = delta * X,
  # then f_S(s) = f_X(s/delta) / delta = dbeta(s/delta, a_0, b_0) / delta
  # ---------------------------------------------------------------------------
  inside <- (s >= 0 & s <= delta)
  
  f <- numeric(length(s))
  s_scaled  <- s[inside] / delta       # rescale to [0, 1]
  f[inside] <- dbeta(s_scaled, shape1 = a_0, shape2 = b_0) / delta
  
  return(f)
}
