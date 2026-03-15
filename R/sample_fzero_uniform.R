# =============================================================================
# sample_fzero_uniform.R
# -----------------------------------------------------------------------------
# Draws m_i i.i.d. samples from the Uniform density of the measurement
# times f_0, as described in Section 1.6 of Lin & Wang (2022):
#
#   S_j ~ Uniform(0, delta),   j = 1, ..., m_i
#
# The actual measurement times are then obtained by shifting to the
# subject's window [O_i - delta/2, O_i + delta/2]:
#
#   T_ij = S_j + O_i - delta/2
#
# This shift is applied in draw_measurement_times.R, NOT here.
# This function only samples S_j from f_0 on [0, delta].
#
# SAMPLING METHOD:
#   Direct simulation via runif(). No approximation involved.
#
# RELATIONSHIP WITH density_fzero_uniform.R:
#   density_fzero_uniform.R  -> fzero_uniform()        evaluates f_0(s)
#   sample_fzero_uniform.R   -> sample_fzero_uniform()  draws S_j ~ f_0
#
# This function is called by:
#   sample_fzero.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw m i.i.d. samples from the Uniform density of f_0 on [0, delta]
#'
#' @param m     Integer >= 1. Number of samples to draw (measurements per
#'              subject).
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [0, delta].
#'
#' @return Numeric vector of length m with samples S_1, ..., S_m,
#'         all in [0, delta]. These are RELATIVE positions within the
#'         window — add O_i - delta/2 to obtain the actual T_ij.
#'
#' @examples
#' # Sample 5 relative positions within a window of length 0.4
#' S <- sample_fzero_uniform(m = 5, delta = 0.4)
#' range(S)   # inside [0, 0.4]
#'
#' # Actual measurement times for subject i with O_i = 0.6
#' O_i <- 0.6; delta <- 0.4
#' T_ij <- S + O_i - delta / 2   # shift to [O_i - 0.2, O_i + 0.2] = [0.4, 0.8]

sample_fzero_uniform <- function(m, delta) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(m) || length(m) != 1 || m < 1 || m != round(m))
    stop("m must be a positive integer.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  
  # ---------------------------------------------------------------------------
  # SAMPLING
  # S_j ~ Uniform(0, delta) via direct simulation
  # ---------------------------------------------------------------------------
  samples <- runif(m, min = 0, max = delta)
  
  return(samples)
}
