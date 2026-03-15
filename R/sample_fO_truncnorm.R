# =============================================================================
# sample_fO_truncnorm.R
# -----------------------------------------------------------------------------
# Draws n i.i.d. samples from the Truncated Normal density of the reference
# times f_O, as described in Section 1.5 of Lin & Wang (2022):
#
#   O_i ~ TruncatedNormal(mu_O, sd_O) on [delta/2, 1 - delta/2]
#
# SAMPLING METHOD — inverse-CDF:
#   If U ~ Uniform(Phi(lo; mu_O, sd_O), Phi(hi; mu_O, sd_O))
#   then Phi^{-1}(U; mu_O, sd_O) ~ TruncNorm(mu_O, sd_O) on [lo, hi].
#
#   Explicitly:
#     p_lo <- pnorm(lo, mean = mu_O, sd = sd_O)
#     p_hi <- pnorm(hi, mean = mu_O, sd = sd_O)
#     U    <- runif(n, min = p_lo, max = p_hi)
#     O_i  <- qnorm(U, mean = mu_O, sd = sd_O)
#
#   No extra packages needed — only base R pnorm, qnorm, runif.
#
# RELATIONSHIP WITH density_fO_truncnorm.R:
#   This file is the sampling companion of density_fO_truncnorm.R:
#     density_fO_truncnorm.R  -> fO_truncnorm()        evaluates f_O(u)
#     sample_fO_truncnorm.R   -> sample_fO_truncnorm()  draws O_i ~ f_O
#
# This function is called by:
#   sample_fO.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw n i.i.d. samples from the Truncated Normal density of f_O
#'
#' @param n     Integer >= 1. Number of samples to draw.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              (delta/2, 1 - delta/2).
#' @param mu_O  Numeric. Centre of the normal distribution. Default is the
#'              midpoint of the support: 0.5.
#' @param sd_O  Numeric > 0. Scale of the normal distribution. Default is
#'              (1 - delta) / 4, placing ~95% of the normal mass inside
#'              the support.
#'
#' @return Numeric vector of length n with samples O_1, ..., O_n,
#'         all in (delta/2, 1 - delta/2).
#'
#' @examples
#' # Default: centred at 0.5
#' O <- sample_fO_truncnorm(n = 200, delta = 0.4)
#' hist(O, breaks = 20, main = "f_O truncnorm, delta = 0.4")
#'
#' # Shifted towards 0.7: more subjects enter late
#' O2 <- sample_fO_truncnorm(n = 500, delta = 0.4, mu_O = 0.7, sd_O = 0.1)

#' @export
sample_fO_truncnorm <- function(n, delta, mu_O = NULL, sd_O = NULL) {
  
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
  # DEFAULT PARAMETERS
  # mu_O = midpoint of support = 0.5
  # sd_O = (hi - lo) / 4  so that mu_O +/- 2*sd_O coincides with [lo, hi]
  #        -> ~95.45% of the normal mass falls inside the support
  # ---------------------------------------------------------------------------
  if (is.null(mu_O)) mu_O <- (lo + hi) / 2
  if (is.null(sd_O)) sd_O <- (hi - lo) / 4
  
  if (!is.numeric(mu_O) || length(mu_O) != 1)
    stop("mu_O must be a single number.")
  if (!is.numeric(sd_O) || length(sd_O) != 1 || sd_O <= 0)
    stop("sd_O must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # CHECK: sufficient mass inside the support
  # ---------------------------------------------------------------------------
  p_lo <- pnorm(lo, mean = mu_O, sd = sd_O)
  p_hi <- pnorm(hi, mean = mu_O, sd = sd_O)
  
  if (p_hi - p_lo < 1e-10)
    stop("Truncated normal has negligible mass on [delta/2, 1-delta/2]. ",
         "Increase sd_O or move mu_O closer to the centre of the support.")
  
  # ---------------------------------------------------------------------------
  # SAMPLING — inverse-CDF method
  # U ~ Uniform(Phi(lo), Phi(hi))  =>  Phi^{-1}(U) ~ TruncNorm on [lo, hi]
  # ---------------------------------------------------------------------------
  U       <- runif(n, min = p_lo, max = p_hi)
  samples <- qnorm(U, mean = mu_O, sd = sd_O)
  
  # clip any tiny numerical overflows to the boundary
  samples <- pmax(pmin(samples, hi), lo)
  
  return(samples)
}
