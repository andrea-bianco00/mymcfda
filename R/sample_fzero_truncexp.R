# =============================================================================
# sample_fzero_truncexp.R
# -----------------------------------------------------------------------------
# Draws m i.i.d. samples from the Truncated Exponential density of the
# measurement times f_0, as described in Section 1.6 of Lin & Wang (2022):
#
#   S_j ~ TruncatedExponential(lambda) on [0, delta],   j = 1, ..., m_i
#
# The actual measurement times are then obtained by shifting to the
# subject's window [O_i - delta/2, O_i + delta/2]:
#
#   T_ij = S_j + O_i - delta/2
#
# This shift is applied in draw_measurement_times.R, NOT here.
# This function only samples S_j from f_0 on [0, delta].
#
# SAMPLING METHOD — inverse-CDF:
#   The CDF of the Truncated Exponential on [0, delta] is:
#
#     F(s) = (1 - exp(-lambda * s)) / (1 - exp(-lambda * delta))
#
#   Inverting: if U ~ Uniform(0, 1), then
#
#     S_j = -(1/lambda) * log(1 - U * (1 - exp(-lambda * delta)))
#
#   No extra packages needed — only base R runif, log, exp.
#
# DEFAULT PARAMETER:
#   lambda = 1 / delta — moderate decay proportional to the support length,
#   consistent across all values of delta. As lambda -> 0, f_0 converges
#   to the Uniform density.
#
# RELATIONSHIP WITH density_fzero_truncexp.R:
#   density_fzero_truncexp.R  -> fzero_truncexp()        evaluates f_0(s)
#   sample_fzero_truncexp.R   -> sample_fzero_truncexp()  draws S_j ~ f_0
#
# This function is called by:
#   sample_fzero.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw m i.i.d. samples from the Truncated Exponential density of f_0 on (0, delta)
#'
#' @param m      Integer >= 1. Number of samples to draw (measurements per
#'               subject).
#' @param delta  Numeric in (0, 1). Snippet length. Determines the support
#'               (0, delta).
#' @param lambda Numeric > 0. Decay rate parameter.
#'                 small lambda  -> nearly uniform (slow decay)
#'                 large lambda  -> concentrated near s = 0 (fast decay)
#'               Default: 1 / delta, i.e. the reciprocal of the support
#'               length.
#'
#' @return Numeric vector of length m with samples S_1, ..., S_m,
#'         all in (0, delta). These are RELATIVE positions within the
#'         window — add O_i - delta/2 to obtain the actual T_ij.
#'
#' @examples
#' # Default lambda = 1/delta = 2.5
#' S <- sample_fzero_truncexp(m = 5, delta = 0.4)
#' range(S)   # inside (0, 0.4)
#'
#' # Actual measurement times for subject i with O_i = 0.6
#' O_i <- 0.6; delta <- 0.4
#' T_ij <- S + O_i - delta / 2   # shift to (O_i - 0.2, O_i + 0.2) = (0.4, 0.8)
#'
#' # Fast decay: concentrated near s = 0
#' S2 <- sample_fzero_truncexp(m = 1000, delta = 0.4, lambda = 10)
#' hist(S2, breaks = 30, main = "f_0 truncexp fast decay")
#'
#' # Slow decay: nearly uniform
#' S3 <- sample_fzero_truncexp(m = 1000, delta = 0.4, lambda = 0.01)

#' @export
sample_fzero_truncexp <- function(m, delta, lambda = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(m) || length(m) != 1 || m < 1 || m != round(m))
    stop("m must be a positive integer.")
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
  # SAMPLING — inverse-CDF method
  #
  # CDF: F(s) = (1 - exp(-lambda * s)) / (1 - exp(-lambda * delta))
  #
  # Inverse: F^{-1}(p) = -(1/lambda) * log(1 - p * (1 - exp(-lambda * delta)))
  #
  # So: draw U ~ Uniform(0, 1), then set
  #   S_j = -(1/lambda) * log(1 - U * (1 - exp(-lambda * delta)))
  # ---------------------------------------------------------------------------
  U       <- runif(m)
  C       <- 1 - exp(-lambda * delta)    # normalisation constant
  samples <- -(1 / lambda) * log(1 - U * C)
  
  # clip any tiny numerical overflows to the boundary
  samples <- pmax(pmin(samples, delta), 0)
  
  return(samples)
}
