# =============================================================================
# density_fO_truncbeta.R
# -----------------------------------------------------------------------------
# Defines the Rescaled Beta density for the reference times f_O, as
# described in Section 1.5 of Lin & Wang (2022):
#
#   If X ~ Beta(a_O, b_O) on [0, 1], then U = lo + (hi - lo) * X
#   has density on [lo, hi] = [delta/2, 1 - delta/2]:
#
#   f_O(u) = dbeta((u - lo) / (hi - lo), a_O, b_O) / (hi - lo),
#                                               u in [delta/2, 1 - delta/2]
#   f_O(u) = 0,                                 otherwise
#
# INTERPRETATION (Section 1.5, Lin & Wang 2022):
#   "Allows flexible asymmetric recruitment patterns."
#   The shape of f_O is controlled by a_O and b_O:
#     a_O = b_O        -> symmetric around the midpoint of the support
#     a_O > b_O        -> left-skewed  (more subjects enter early)
#     a_O < b_O        -> right-skewed (more subjects enter late)
#     a_O = b_O = 1    -> reduces to the Uniform density (fO_uniform)
#     a_O, b_O -> Inf  -> concentrates around the midpoint (like truncnorm)
#
# STRICT POSITIVITY (Assumption 2, Lin & Wang 2022):
#   The Beta(a_O, b_O) density is strictly positive on (0,1) if and only if
#   a_O > 1 AND b_O > 1. At the boundaries:
#     - if a_O <= 1: dbeta(0, a_O, b_O) = 0 or Inf
#     - if b_O <= 1: dbeta(1, a_O, b_O) = 0 or Inf
#   Therefore a_O > 1 and b_O > 1 is required to satisfy Assumption 2.
#   A warning is raised if this condition is violated.
#
# DEFAULT PARAMETERS:
#   a_O = b_O = 2  -> symmetric, bell-shaped, strictly positive at boundaries.
#   This is the simplest choice satisfying Assumption 2 with a_O > 1, b_O > 1.
#
# This function is called by:
#   eval_fO.R    — to evaluate f_O(u) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Rescaled Beta density for reference times f_O
#'
#' @param u     Numeric vector. Points at which to evaluate f_O.
#'              Values outside [delta/2, 1 - delta/2] return 0.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [delta/2, 1 - delta/2].
#' @param a_O   Numeric > 0. First shape parameter of Beta(a_O, b_O).
#'              Must be > 1 for strict positivity at the left boundary.
#'              Default: 2.
#' @param b_O   Numeric > 0. Second shape parameter of Beta(a_O, b_O).
#'              Must be > 1 for strict positivity at the right boundary.
#'              Default: 2.
#'
#' @return Numeric vector of the same length as u with values of f_O(u).
#'
#' @examples
#' u <- seq(0, 1, length.out = 500)
#'
#' # Symmetric Beta(2,2) — default
#' f <- fO_truncbeta(u, delta = 0.4)
#' plot(u, f, type = "l", main = "f_O truncbeta, delta = 0.4")
#'
#' # Left-skewed Beta(3,2): more subjects enter early
#' f2 <- fO_truncbeta(u, delta = 0.4, a_O = 3, b_O = 2)
#' lines(u, f2, col = "red")
#'
#' # Right-skewed Beta(2,3): more subjects enter late
#' f3 <- fO_truncbeta(u, delta = 0.4, a_O = 2, b_O = 3)
#' lines(u, f3, col = "blue")
#'
#' # Verify it integrates to 1
#' sum(f) * (1 / 500)   # approximately 1

fO_truncbeta <- function(u, delta, a_O = 2, b_O = 2) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(u))
    stop("u must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  if (!is.numeric(a_O) || length(a_O) != 1 || a_O <= 0)
    stop("a_O must be a single positive number.")
  if (!is.numeric(b_O) || length(b_O) != 1 || b_O <= 0)
    stop("b_O must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # STRICT POSITIVITY WARNING
  # Required by Assumption 2 of Lin & Wang (2022): f_O(u) > 0 on full support.
  # Beta(a_O, b_O) is strictly positive on (0,1) only when a_O > 1 AND b_O > 1.
  # ---------------------------------------------------------------------------
  if (a_O <= 1 || b_O <= 1)
    warning("a_O <= 1 or b_O <= 1: f_O may be zero or infinite at the ",
            "boundaries of [delta/2, 1-delta/2], violating the strict ",
            "positivity condition of Assumption 2 in Lin & Wang (2022). ",
            "Use a_O > 1 and b_O > 1 to ensure identifiability of mu(t) ",
            "at all t in [0, 1].")
  
  # ---------------------------------------------------------------------------
  # SUPPORT
  # ---------------------------------------------------------------------------
  lo <- delta / 2
  hi <- 1 - delta / 2
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # Rescaling: if X ~ Beta(a_O, b_O) on [0,1] and U = lo + (hi-lo)*X,
  # then f_U(u) = f_X((u-lo)/(hi-lo)) / (hi-lo) = dbeta(...) / (hi-lo)
  # ---------------------------------------------------------------------------
  inside <- (u >= lo & u <= hi)
  
  f <- numeric(length(u))
  
  u_scaled  <- (u[inside] - lo) / (hi - lo)   # rescale to [0, 1]
  f[inside] <- dbeta(u_scaled, shape1 = a_O, shape2 = b_O) / (hi - lo)
  
  return(f)
}
