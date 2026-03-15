# =============================================================================
# sample_fzero_truncbeta.R
# -----------------------------------------------------------------------------
# Draws m i.i.d. samples from the Rescaled Beta density of the measurement
# times f_0, as described in Section 1.6 of Lin & Wang (2022):
#
#   If X ~ Beta(a_0, b_0) on [0, 1], then S = delta * X has density:
#
#   f_0(s) = dbeta(s / delta, a_0, b_0) / delta,   s in [0, delta]
#
# The actual measurement times are then obtained by shifting to the
# subject's window [O_i - delta/2, O_i + delta/2]:
#
#   T_ij = S_j + O_i - delta/2
#
# This shift is applied in draw_measurement_times.R, NOT here.
# This function only samples S_j from f_0 on [0, delta].
#
# SAMPLING METHOD — direct rescaling:
#   1. Draw X ~ Beta(a_0, b_0) via rbeta()
#   2. Rescale: S_j = delta * X
#   No approximation involved.
#
# DEFAULT PARAMETERS:
#   a_0 = b_0 = 2  -> symmetric, bell-shaped, strictly positive at
#                     boundaries. Simplest choice satisfying a_0 > 1, b_0 > 1,
#                     as required by Assumption 3 of Lin & Wang (2022).
#
# RELATIONSHIP WITH density_fzero_truncbeta.R:
#   density_fzero_truncbeta.R  -> fzero_truncbeta()        evaluates f_0(s)
#   sample_fzero_truncbeta.R   -> sample_fzero_truncbeta()  draws S_j ~ f_0
#
# This function is called by:
#   sample_fzero.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw m i.i.d. samples from the Rescaled Beta density of f_0 on (0, delta)
#'
#' @param m     Integer >= 1. Number of samples to draw (measurements per
#'              subject).
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              (0, delta).
#' @param a_0   Numeric > 0. First shape parameter of Beta(a_0, b_0).
#'              Must be > 1 for strict positivity at s = 0.
#'              Default: 2.
#' @param b_0   Numeric > 0. Second shape parameter of Beta(a_0, b_0).
#'              Must be > 1 for strict positivity at s = delta.
#'              Default: 2.
#'
#' @return Numeric vector of length m with samples S_1, ..., S_m,
#'         all in (0, delta). These are RELATIVE positions within the
#'         window — add O_i - delta/2 to obtain the actual T_ij.
#'
#' @examples
#' # Default: symmetric Beta(2,2)
#' S <- sample_fzero_truncbeta(m = 5, delta = 0.4)
#' range(S)   # inside (0, 0.4)
#'
#' # Actual measurement times for subject i with O_i = 0.6
#' O_i <- 0.6; delta <- 0.4
#' T_ij <- S + O_i - delta / 2   # shift to (O_i - 0.2, O_i + 0.2) = (0.4, 0.8)
#'
#' # Concentrated towards the start of the window: Beta(3,2)
#' S2 <- sample_fzero_truncbeta(m = 1000, delta = 0.4, a_0 = 3, b_0 = 2)
#' hist(S2, breaks = 30, main = "f_0 truncbeta Beta(3,2)")
#'
#' # Concentrated towards the end of the window: Beta(2,3)
#' S3 <- sample_fzero_truncbeta(m = 1000, delta = 0.4, a_0 = 2, b_0 = 3)

#' @export
sample_fzero_truncbeta <- function(m, delta, a_0 = 2, b_0 = 2) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(m) || length(m) != 1 || m < 1 || m != round(m))
    stop("m must be a positive integer.")
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
  # SAMPLING — direct rescaling
  # X ~ Beta(a_0, b_0) on [0, 1]
  # S_j = delta * X  lies in [0, delta]
  # ---------------------------------------------------------------------------
  X       <- rbeta(m, shape1 = a_0, shape2 = b_0)
  samples <- delta * X
  
  return(samples)
}
