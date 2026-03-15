# =============================================================================
# sample_fzero_truncnorm.R
# -----------------------------------------------------------------------------
# Draws m i.i.d. samples from the Truncated Normal density of the measurement
# times f_0, as described in Section 1.6 of Lin & Wang (2022):
#
#   S_j ~ TruncatedNormal(mu_0, sd_0) on [0, delta],   j = 1, ..., m_i
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
#   If U ~ Uniform(Phi(lo; mu_0, sd_0), Phi(hi; mu_0, sd_0))
#   then Phi^{-1}(U; mu_0, sd_0) ~ TruncNorm(mu_0, sd_0) on [lo, hi].
#
#   Explicitly:
#     p_lo <- pnorm(lo, mean = mu_0, sd = sd_0)   with lo = 0
#     p_hi <- pnorm(hi, mean = mu_0, sd = sd_0)   with hi = delta
#     U    <- runif(m, min = p_lo, max = p_hi)
#     S_j  <- qnorm(U, mean = mu_0, sd = sd_0)
#
#   No extra packages needed — only base R pnorm, qnorm, runif.
#
# DEFAULT PARAMETERS:
#   mu_0 = delta / 2   midpoint of [0, delta] -> measurements cluster
#                      around the centre of each subject's window
#   sd_0 = delta / 4   so that mu_0 +/- 2*sd_0 = [0, delta] ->
#                      ~95.45% of the normal mass falls inside [0, delta]
#
# RELATIONSHIP WITH density_fzero_truncnorm.R:
#   density_fzero_truncnorm.R  -> fzero_truncnorm()        evaluates f_0(s)
#   sample_fzero_truncnorm.R   -> sample_fzero_truncnorm()  draws S_j ~ f_0
#
# This function is called by:
#   sample_fzero.R   — dispatcher that selects the sampling function by type
#
# Dependencies: base R only.
# =============================================================================


#' Draw m i.i.d. samples from the Truncated Normal density of f_0 on [0, delta]
#'
#' @param m     Integer >= 1. Number of samples to draw (measurements per
#'              subject).
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [0, delta].
#' @param mu_0  Numeric. Centre of the normal distribution on [0, delta].
#'              Default: delta/2 (midpoint of the window).
#' @param sd_0  Numeric > 0. Scale of the normal distribution.
#'              Default: delta/4, placing ~95% of the mass inside [0, delta].
#'
#' @return Numeric vector of length m with samples S_1, ..., S_m,
#'         all in [0, delta]. These are RELATIVE positions within the
#'         window — add O_i - delta/2 to obtain the actual T_ij.
#'
#' @examples
#' # Default: centred at delta/2 = 0.2
#' S <- sample_fzero_truncnorm(m = 5, delta = 0.4)
#' range(S)   # inside [0, 0.4]
#'
#' # Actual measurement times for subject i with O_i = 0.6
#' O_i <- 0.6; delta <- 0.4
#' T_ij <- S + O_i - delta / 2   # shift to [O_i - 0.2, O_i + 0.2] = [0.4, 0.8]
#'
#' # Shifted towards the end of the window
#' S2 <- sample_fzero_truncnorm(m = 1000, delta = 0.4, mu_0 = 0.3, sd_0 = 0.08)
#' hist(S2, breaks = 30, main = "f_0 truncnorm shifted")

sample_fzero_truncnorm <- function(m, delta, mu_0 = NULL, sd_0 = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(m) || length(m) != 1 || m < 1 || m != round(m))
    stop("m must be a positive integer.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  
  # ---------------------------------------------------------------------------
  # SUPPORT: [0, delta]
  # ---------------------------------------------------------------------------
  lo <- 0
  hi <- delta
  
  # ---------------------------------------------------------------------------
  # DEFAULT PARAMETERS
  # mu_0 = delta/2  — midpoint of [0, delta]
  # sd_0 = delta/4  — so that mu_0 +/- 2*sd_0 coincides with [0, delta]
  #                   -> ~95.45% of the normal mass inside the support
  # ---------------------------------------------------------------------------
  if (is.null(mu_0)) mu_0 <- delta / 2
  if (is.null(sd_0)) sd_0 <- delta / 4
  
  if (!is.numeric(mu_0) || length(mu_0) != 1)
    stop("mu_0 must be a single number.")
  if (!is.numeric(sd_0) || length(sd_0) != 1 || sd_0 <= 0)
    stop("sd_0 must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # CHECK: sufficient mass inside the support
  # ---------------------------------------------------------------------------
  p_lo <- pnorm(lo, mean = mu_0, sd = sd_0)
  p_hi <- pnorm(hi, mean = mu_0, sd = sd_0)
  
  if (p_hi - p_lo < 1e-10)
    stop("Truncated normal has negligible mass on [0, delta]. ",
         "Increase sd_0 or move mu_0 closer to the centre delta/2.")
  
  # ---------------------------------------------------------------------------
  # SAMPLING — inverse-CDF method
  # U ~ Uniform(Phi(0), Phi(delta))  =>  Phi^{-1}(U) ~ TruncNorm on [0, delta]
  # ---------------------------------------------------------------------------
  U       <- runif(m, min = p_lo, max = p_hi)
  samples <- qnorm(U, mean = mu_0, sd = sd_0)
  
  # clip any tiny numerical overflows to the boundary
  samples <- pmax(pmin(samples, hi), lo)
  
  return(samples)
}
