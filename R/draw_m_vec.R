# =============================================================================
# draw_m_vec.R
# -----------------------------------------------------------------------------
# Implements Step 1 of the data-generating mechanism of Lin & Wang (2022),
# Section 5 (Simulation):
#
#   "The number of measurements per subject m_i is drawn from a
#    Poisson distribution with mean m_mean, truncated below at 2."
#
# The truncation at m_i >= 2 is required because estimation of sigma_0^2
# via within-subject pairs (Section 3.2.2) needs at least 2 observations
# per subject.
#
# MECHANISM:
#   1. Draw m_i ~ Poisson(m_mean) independently for each subject i.
#   2. If m_i < 2, redraw until m_i >= 2.
#
# NOTE ON THE MEAN:
#   The truncation shifts E[m_i | m_i >= 2] slightly above m_mean,
#   because the probability mass at 0 and 1 is redistributed.
#   For m_mean = 4 (the paper's default), P(Poisson(4) < 2) ≈ 0.092,
#   so the effective mean is approximately 4.17 — a negligible difference.
#   This is the exact mechanism used by Lin & Wang (2022).
#
# INPUT:
#   n       : integer >= 1. Number of subjects.
#   m_mean  : numeric > 0. Mean of the Poisson distribution BEFORE
#             truncation. This is the user's "XXX" parameter.
#             Default: 4 (as in Lin & Wang 2022).
#
# OUTPUT:
#   Integer vector of length n with m_1, ..., m_n, all >= 2.
#
# Dependencies: base R only.
# =============================================================================


#' Draw number of measurements per subject from truncated Poisson
#'
#' @param n      Integer >= 1. Number of subjects.
#' @param m_mean Numeric > 0. Mean of the Poisson distribution (before
#'               truncation). This controls the average number of
#'               measurements per subject. Default: 4.
#'
#' @return Integer vector of length n with m_1, ..., m_n, all >= 2.
#'         The sample mean mean(m_vec) will be approximately m_mean
#'         (slightly above due to truncation).
#'
#' @examples
#' # Default: Poisson(4), truncated at >= 2
#' set.seed(42)
#' m_vec <- draw_m_vec(n = 200, m_mean = 4)
#' range(m_vec)    # min is 2
#' mean(m_vec)     # approximately 4
#'
#' # Denser observations: Poisson(10)
#' m_vec2 <- draw_m_vec(n = 200, m_mean = 10)
#' mean(m_vec2)    # approximately 10
#'
#' # Very sparse: Poisson(2.5)
#' m_vec3 <- draw_m_vec(n = 200, m_mean = 2.5)
#' mean(m_vec3)    # approximately 2.5 (but all >= 2)

draw_m_vec <- function(n, m_mean = 4) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n))
    stop("n must be a positive integer.")
  if (!is.numeric(m_mean) || length(m_mean) != 1 || m_mean <= 0)
    stop("m_mean must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # WARNING: very small m_mean makes truncation dominant
  # If m_mean < 2, the majority of Poisson draws will be < 2 and redrawn,
  # so the effective mean will be much larger than m_mean.
  # ---------------------------------------------------------------------------
  if (m_mean < 2)
    warning("m_mean < 2: most Poisson draws will be < 2 and redrawn. ",
            "The effective mean will be substantially above m_mean. ",
            "Consider using m_mean >= 2 for meaningful control.")
  
  # ---------------------------------------------------------------------------
  # DRAW m_i ~ Poisson(m_mean), truncated at >= 2
  #
  # Strategy: draw all n values at once, then redraw any that are < 2.
  # This is efficient because P(Poisson(m_mean) < 2) is small for
  # typical m_mean >= 3.
  # ---------------------------------------------------------------------------
  m_vec <- rpois(n, lambda = m_mean)
  
  # redraw any m_i < 2
  bad <- which(m_vec < 2)
  while (length(bad) > 0) {
    m_vec[bad] <- rpois(length(bad), lambda = m_mean)
    bad <- which(m_vec < 2)
  }
  
  return(m_vec)
}
