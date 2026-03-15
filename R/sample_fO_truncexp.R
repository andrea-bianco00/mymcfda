# =============================================================================
# sample_fO_truncexp.R
# -----------------------------------------------------------------------------
# Draws n i.i.d. samples from the Truncated Exponential density of the
# reference times f_O, as described in Section 1.5 of Lin & Wang (2022):
#
#   O_i ~ TruncatedExponential(lambda) on [delta/2, 1 - delta/2]
#
# SAMPLING METHOD — inverse-CDF:
#   The CDF of the Truncated Exponential on [lo, hi] is:
#
#     F(u) = (1 - exp(-lambda*(u - lo))) / (1 - exp(-lambda*(hi - lo)))
#
#   Inverting: if U ~ Uniform(0, 1), then
#
#     O_i = lo - (1/lambda) * log(1 - U * (1 - exp(-lambda*(hi - lo))))
#
#   No extra packages needed — only base R runif, log, exp.
#
# RELATIONSHIP WITH density_fO_truncexp.R:
#   This file is the sampling companion of density_fO_truncexp.R:
#     density_fO_truncexp.R  -> fO_truncexp()        evaluates f_O(u)
#     sample_fO_truncexp.R   -> sample_fO_truncexp()  draws O_i ~ f_O
#
# This function is called by:
#   sample_fO.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw n i.i.d. samples from the Truncated Exponential density of f_O
#'
#' @param n      Integer >= 1. Number of samples to draw.
#' @param delta  Numeric in (0, 1). Snippet length. Determines the support
#'               (delta/2, 1 - delta/2).
#' @param lambda Numeric > 0. Decay rate parameter.
#'                 small lambda  -> nearly uniform (slow decay)
#'                 large lambda  -> concentrated near delta/2 (fast decay)
#'               Default: 1 / (1 - delta), i.e. the reciprocal of the
#'               support length.
#'
#' @return Numeric vector of length n with samples O_1, ..., O_n,
#'         all in (delta/2, 1 - delta/2).
#'
#' @examples
#' # Default lambda
#' O <- sample_fO_truncexp(n = 200, delta = 0.4)
#' hist(O, breaks = 20, main = "f_O truncexp, delta = 0.4")
#'
#' # Fast decay: heavy concentration near lo = 0.2
#' O2 <- sample_fO_truncexp(n = 500, delta = 0.4, lambda = 10)
#'
#' # Slow decay: nearly uniform
#' O3 <- sample_fO_truncexp(n = 500, delta = 0.4, lambda = 0.01)

#' @export
sample_fO_truncexp <- function(n, delta, lambda = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n))
    stop("n must be a positive integer.")
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
  # ---------------------------------------------------------------------------
  if (is.null(lambda)) lambda <- 1 / (hi - lo)
  
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("lambda must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # SAMPLING — inverse-CDF method
  #
  # CDF: F(u) = (1 - exp(-lambda*(u-lo))) / (1 - exp(-lambda*(hi-lo)))
  #
  # Inverse: F^{-1}(p) = lo - (1/lambda) * log(1 - p*(1-exp(-lambda*(hi-lo))))
  #
  # So: draw U ~ Uniform(0,1), then set
  #   O_i = lo - (1/lambda) * log(1 - U * (1 - exp(-lambda*(hi-lo))))
  # ---------------------------------------------------------------------------
  U       <- runif(n)
  C       <- 1 - exp(-lambda * (hi - lo))   # normalisation constant
  samples <- lo - (1 / lambda) * log(1 - U * C)
  
  # clip any tiny numerical overflows to the boundary
  samples <- pmax(pmin(samples, hi), lo)
  
  return(samples)
}
