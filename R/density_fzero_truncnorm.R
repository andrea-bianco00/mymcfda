# =============================================================================
# density_fzero_truncnorm.R
# -----------------------------------------------------------------------------
# Defines the Truncated Normal density for the measurement times f_0, as
# described in Section 1.6 of Lin & Wang (2022), Assumption 3:
#
#   f_0(s) = phi((s - mu_0) / sd_0) / (sd_0 * Z),   s in [0, delta]
#   f_0(s) = 0,                                       otherwise
#
# where:
#   phi(.) = standard normal pdf
#   Z      = Phi((delta - mu_0)/sd_0) - Phi((0 - mu_0)/sd_0)  normalisation
#
# ROLE OF f_0 IN THE MODEL (Section 1.6, Lin & Wang 2022):
#   f_0 is defined on the reference interval [0, delta]. The conditional
#   density of T_ij given O_i = u is obtained by shifting f_0 to the
#   subject's window [u - delta/2, u + delta/2]:
#
#     f_{T|O}(t | u) = f_0(t - u + delta/2)
#
#   So f_0 describes HOW measurements are spread within a subject's window.
#   With the Truncated Normal, measurements tend to cluster around the
#   relative position mu_0 within the window.
#
# STRICT POSITIVITY (Assumption 3, Lin & Wang 2022):
#   phi(.) > 0 everywhere, so f_0(s) > 0 for all s in [0, delta],
#   provided Z > 0.
#
# LIPSCHITZ CONTINUITY OF f_0' (Assumption 3, Lin & Wang 2022):
#   The derivative of the truncated normal density is Lipschitz continuous
#   on [0, delta] since phi' is bounded and continuous.
#
# DEFAULT PARAMETERS:
#   mu_0 = delta / 2   midpoint of [0, delta] -> measurements cluster
#                      around the centre of each subject's window
#   sd_0 = delta / 4   so that mu_0 +/- 2*sd_0 = [0, delta] ->
#                      ~95.45% of the normal mass falls inside [0, delta]
#
# RELATIONSHIP WITH density_fO_truncnorm.R:
#   Same mathematical form, different support:
#     density_fO_truncnorm.R  : support [delta/2, 1-delta/2]
#     density_fzero_truncnorm.R: support [0, delta]
#
# This function is called by:
#   eval_fzero.R   — to evaluate f_0(s) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Truncated Normal density for measurement times f_0
#'
#' @param s     Numeric vector. Points at which to evaluate f_0.
#'              Values outside [0, delta] return 0.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [0, delta].
#' @param mu_0  Numeric. Centre of the normal distribution on [0, delta].
#'              Default: delta/2 (midpoint of the window).
#' @param sd_0  Numeric > 0. Scale of the normal distribution.
#'              Default: delta/4, placing ~95% of the mass inside [0, delta].
#'
#' @return Numeric vector of the same length as s with values of f_0(s).
#'
#' @examples
#' s <- seq(-0.1, 0.6, length.out = 500)
#'
#' # Default: centred at delta/2 = 0.2
#' f <- fzero_truncnorm(s, delta = 0.4)
#' plot(s, f, type = "l", main = "f_0 truncnorm, delta = 0.4")
#'
#' # Shifted towards the end of the window
#' f2 <- fzero_truncnorm(s, delta = 0.4, mu_0 = 0.3, sd_0 = 0.08)
#' lines(s, f2, col = "red")
#'
#' # Verify it integrates to 1
#' s_fine <- seq(0, 0.4, length.out = 10000)
#' sum(fzero_truncnorm(s_fine, delta = 0.4)) * (0.4 / 10000)  # approximately 1

fzero_truncnorm <- function(s, delta, mu_0 = NULL, sd_0 = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s))
    stop("s must be a numeric vector.")
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
  # NORMALISATION CONSTANT
  # Z = P(0 <= X <= delta) where X ~ Normal(mu_0, sd_0)
  # ---------------------------------------------------------------------------
  Z <- pnorm(hi, mean = mu_0, sd = sd_0) - pnorm(lo, mean = mu_0, sd = sd_0)
  
  if (Z < 1e-10)
    stop("Truncated normal has negligible mass on [0, delta]. ",
         "Increase sd_0 or move mu_0 closer to the centre delta/2.")
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # f_0(s) = phi((s - mu_0) / sd_0) / (sd_0 * Z)   inside [0, delta]
  # f_0(s) = 0                                       outside
  # ---------------------------------------------------------------------------
  inside <- (s >= lo & s <= hi)
  
  f <- numeric(length(s))
  f[inside] <- dnorm(s[inside], mean = mu_0, sd = sd_0) / Z
  
  return(f)
}
