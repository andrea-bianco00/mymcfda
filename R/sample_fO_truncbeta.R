# =============================================================================
# sample_fO_truncbeta.R
# -----------------------------------------------------------------------------
# Draws n i.i.d. samples from the Rescaled Beta density of the reference
# times f_O, as described in Section 1.5 of Lin & Wang (2022):
#
#   If X ~ Beta(a_O, b_O) on [0, 1], then
#   O_i = lo + (hi - lo) * X ~ Rescaled Beta on [lo, hi]
#
# where lo = delta/2, hi = 1 - delta/2.
#
# SAMPLING METHOD — direct rescaling:
#   1. Draw X ~ Beta(a_O, b_O) via rbeta()
#   2. Rescale: O_i = lo + (hi - lo) * X
#   No approximation involved.
#
# RELATIONSHIP WITH density_fO_truncbeta.R:
#   This file is the sampling companion of density_fO_truncbeta.R:
#     density_fO_truncbeta.R  -> fO_truncbeta()        evaluates f_O(u)
#     sample_fO_truncbeta.R   -> sample_fO_truncbeta()  draws O_i ~ f_O
#
# This function is called by:
#   sample_fO.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw n i.i.d. samples from the Rescaled Beta density of f_O
#'
#' @param n     Integer >= 1. Number of samples to draw.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [delta/2, 1 - delta/2].
#' @param a_O   Numeric > 0. First shape parameter of Beta(a_O, b_O).
#'              Must be > 1 for strict positivity at the left boundary.
#'              Default: 2.
#' @param b_O   Numeric > 0. Second shape parameter of Beta(a_O, b_O).
#'              Must be > 1 for strict positivity at the right boundary.
#'              Default: 2.
#'
#' @return Numeric vector of length n with samples O_1, ..., O_n,
#'         all in [delta/2, 1 - delta/2].
#'
#' @examples
#' # Default: symmetric Beta(2,2)
#' O <- sample_fO_truncbeta(n = 200, delta = 0.4)
#' hist(O, breaks = 20, main = "f_O truncbeta Beta(2,2), delta = 0.4")
#'
#' # Left-skewed Beta(3,2): more subjects enter early
#' O2 <- sample_fO_truncbeta(n = 500, delta = 0.4, a_O = 3, b_O = 2)
#'
#' # Right-skewed Beta(2,3): more subjects enter late
#' O3 <- sample_fO_truncbeta(n = 500, delta = 0.4, a_O = 2, b_O = 3)

sample_fO_truncbeta <- function(n, delta, a_O = 2, b_O = 2) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n))
    stop("n must be a positive integer.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  if (!is.numeric(a_O) || length(a_O) != 1 || a_O <= 0)
    stop("a_O must be a single positive number.")
  if (!is.numeric(b_O) || length(b_O) != 1 || b_O <= 0)
    stop("b_O must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # STRICT POSITIVITY WARNING
  # Required by Assumption 2 of Lin & Wang (2022).
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
  # SAMPLING — direct rescaling
  # X ~ Beta(a_O, b_O) on [0,1]
  # O_i = lo + (hi - lo) * X  lies in [lo, hi]
  # ---------------------------------------------------------------------------
  X       <- rbeta(n, shape1 = a_O, shape2 = b_O)
  samples <- lo + (hi - lo) * X
  
  return(samples)
}
