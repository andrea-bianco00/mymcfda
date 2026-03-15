# =============================================================================
# sample_fO_uniform.R
# -----------------------------------------------------------------------------
# Draws n i.i.d. samples from the Uniform density of the reference times f_O,
# as described in Section 1.5 of Lin & Wang (2022):
#
#   O_i ~ Uniform(delta/2, 1 - delta/2)
#
# This is the DEFAULT choice in Lin & Wang (2022):
#   "The simplest choice: subjects enter the study completely at random,
#    with no preferred entry time. It is the default assumption in many
#    simulation studies, including Lin & Wang (2022)."
#
# SAMPLING METHOD:
#   Direct simulation via runif(). No approximation involved.
#
# RELATIONSHIP WITH density_fO_uniform.R:
#   This file is the sampling companion of density_fO_uniform.R:
#     density_fO_uniform.R  -> fO_uniform()        evaluates f_O(u)
#     sample_fO_uniform.R   -> sample_fO_uniform()  draws O_i ~ f_O
#
# This function is called by:
#   sample_fO.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw n i.i.d. samples from the Uniform density of f_O
#'
#' @param n     Integer >= 1. Number of samples to draw.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              (delta/2, 1 - delta/2).
#'
#' @return Numeric vector of length n with samples O_1, ..., O_n,
#'         all in (delta/2, 1 - delta/2).
#'
#' @examples
#' O <- sample_fO_uniform(n = 200, delta = 0.4)
#' hist(O, breaks = 20, main = "f_O uniform, delta = 0.4")
#' range(O)   # should be inside (0.2, 0.8)

#' @export
sample_fO_uniform <- function(n, delta) {
  
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
  # SAMPLING
  # O_i ~ Uniform(lo, hi) via direct simulation
  # ---------------------------------------------------------------------------
  samples <- runif(n, min = lo, max = hi)
  
  return(samples)
}
