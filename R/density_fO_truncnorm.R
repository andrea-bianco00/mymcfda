# =============================================================================
# density_fO_truncnorm.R
# -----------------------------------------------------------------------------
# Defines the Truncated Normal density for the reference times f_O, as
# described in Section 1.5 of Lin & Wang (2022):
#
#   f_O(u) = phi((u - mu_O) / sd_O) / (sd_O * Z),   u in [delta/2, 1-delta/2]
#   f_O(u) = 0,                                       otherwise
#
# where:
#   phi(.)  = standard normal pdf
#   Z       = Phi((hi - mu_O)/sd_O) - Phi((lo - mu_O)/sd_O)  normalisation
#   lo      = delta/2
#   hi      = 1 - delta/2
#
# INTERPRETATION (Section 1.5, Lin & Wang 2022):
#   Models a tendency for subjects to be recruited around a preferred time
#   mu_O, while still covering the full domain. sd_O controls how spread
#   out the recruitment is around mu_O.
#
# STRICT POSITIVITY (Assumption 2, Lin & Wang 2022):
#   phi(.) > 0 everywhere, so f_O(u) > 0 for all u in [lo, hi], provided
#   Z > 0. Z is negligible only if mu_O is extremely far from [lo, hi]
#   relative to sd_O — this is caught by a runtime check.
#
# This function is called by:
#   eval_fO.R    — to evaluate f_O(u) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Truncated Normal density for reference times f_O
#'
#' @param u     Numeric vector. Points at which to evaluate f_O.
#'              Values outside (delta/2, 1 - delta/2) return 0.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              (delta/2, 1 - delta/2).
#' @param mu_O  Numeric. Centre of the normal distribution. Default is the
#'              midpoint of the support: (delta/2 + 1 - delta/2) / 2 = 0.5.
#' @param sd_O  Numeric > 0. Scale of the normal distribution. Default is
#'              (1 - delta) / 4, i.e. one quarter of the support length,
#'              which places approximately 95% of the mass inside the support.
#'
#' @return Numeric vector of the same length as u with values of f_O(u).
#'
#' @examples
#' u <- seq(0, 1, length.out = 500)
#'
#' # Centred at 0.5
#' f <- fO_truncnorm(u, delta = 0.4, mu_O = 0.5, sd_O = 0.15)
#' plot(u, f, type = "l", main = "f_O truncnorm, delta = 0.4")
#'
#' # Shifted towards 0.7 (more subjects enter late)
#' f2 <- fO_truncnorm(u, delta = 0.4, mu_O = 0.7, sd_O = 0.1)
#' lines(u, f2, col = "red")
#'
#' # Verify it integrates to 1
#' sum(f) * (1 / 500)   # approximately 1

#' @export
fO_truncnorm <- function(u, delta, mu_O = NULL, sd_O = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(u))
    stop("u must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  
  # ---------------------------------------------------------------------------
  # SUPPORT
  # ---------------------------------------------------------------------------
  lo <- delta / 2
  hi <- 1 - delta / 2
  
  # ---------------------------------------------------------------------------
  # DEFAULT PARAMETERS
  # ---------------------------------------------------------------------------
  if (is.null(mu_O)) mu_O <- (lo + hi) / 2   # midpoint of the support
  if (is.null(sd_O)) sd_O <- (hi - lo) / 4   # quarter of the support length
  
  if (!is.numeric(mu_O) || length(mu_O) != 1)
    stop("mu_O must be a single number.")
  if (!is.numeric(sd_O) || length(sd_O) != 1 || sd_O <= 0)
    stop("sd_O must be a single positive number.")
  
  # ---------------------------------------------------------------------------
  # NORMALISATION CONSTANT
  # Z = P(lo <= X <= hi) where X ~ Normal(mu_O, sd_O)
  # ---------------------------------------------------------------------------
  Z <- pnorm(hi, mean = mu_O, sd = sd_O) - pnorm(lo, mean = mu_O, sd = sd_O)
  
  if (Z < 1e-10)
    stop("Truncated normal has negligible mass on [delta/2, 1-delta/2]. ",
         "Increase sd_O or move mu_O closer to the centre of the support.")
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # f_O(u) = phi((u - mu_O) / sd_O) / (sd_O * Z)   inside support
  # f_O(u) = 0                                       outside support
  # ---------------------------------------------------------------------------
  inside <- (u >= lo & u <= hi)
  
  f <- numeric(length(u))
  f[inside] <- dnorm(u[inside], mean = mu_O, sd = sd_O) / Z
  
  return(f)
}
