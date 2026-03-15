# =============================================================================
# density_fzero_uniform.R
# -----------------------------------------------------------------------------
# Defines the Uniform density for the measurement times f_0, as described
# in Section 1.6 of Lin & Wang (2022), Assumption 3:
#
#   f_0(s) = 1 / delta,   s in [0, delta]
#   f_0(s) = 0,           otherwise
#
# ROLE OF f_0 IN THE MODEL (Section 1.6, Lin & Wang 2022):
#   f_0 is defined on the reference interval [0, delta]. The conditional
#   density of T_ij given O_i = u is obtained by shifting f_0 to the
#   subject's window [u - delta/2, u + delta/2]:
#
#     f_{T|O}(t | u) = f_0(t - u + delta/2)
#
#   So f_0 describes HOW measurements are spread within a subject's window,
#   while O_i determines WHERE the window is located on [0, 1].
#
# THIS CHOICE (Uniform):
#   f_0(s) = 1/delta on [0, delta] means that measurement times are spread
#   uniformly within each subject's window — no preferred measurement time
#   within the window. This is the default choice in Lin & Wang (2022).
#
# STRICT POSITIVITY (Assumption 3, Lin & Wang 2022):
#   f_0(s) = 1/delta > 0 for all s in [0, delta]. Satisfied trivially
#   since delta in (0, 1) implies 1/delta > 0.
#
# LIPSCHITZ CONTINUITY OF f_0' (Assumption 3, Lin & Wang 2022):
#   f_0 is constant on (0, delta), so f_0'(s) = 0 everywhere, which is
#   trivially Lipschitz continuous.
#
# RELATIONSHIP WITH density_fO_uniform.R:
#   Same mathematical form (uniform density), different support:
#     density_fO_uniform.R  : support [delta/2, 1-delta/2], value 1/(1-delta)
#     density_fzero_uniform.R: support [0, delta],          value 1/delta
#
# This function is called by:
#   eval_fzero.R   — to evaluate f_0(s) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Uniform density for measurement times f_0
#'
#' @param s     Numeric vector. Points at which to evaluate f_0.
#'              Values outside [0, delta] return 0.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              [0, delta] and the constant density value 1/delta.
#'
#' @return Numeric vector of the same length as s with values of f_0(s).
#'
#' @examples
#' s <- seq(-0.1, 0.6, length.out = 500)
#' f <- fzero_uniform(s, delta = 0.4)
#' plot(s, f, type = "l", main = "f_0 uniform, delta = 0.4")
#' # f_0 = 1/0.4 = 2.5 on [0, 0.4], zero elsewhere
#'
#' # Verify it integrates to 1
#' s_fine <- seq(0, 0.4, length.out = 10000)
#' sum(fzero_uniform(s_fine, delta = 0.4)) * (0.4 / 10000)  # approximately 1

fzero_uniform <- function(s, delta) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s))
    stop("s must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  
  # ---------------------------------------------------------------------------
  # DENSITY
  # f_0(s) = 1 / delta   for s in [0, delta]
  # f_0(s) = 0           otherwise
  # ---------------------------------------------------------------------------
  f <- ifelse(s >= 0 & s <= delta, 1 / delta, 0)
  
  return(f)
}
