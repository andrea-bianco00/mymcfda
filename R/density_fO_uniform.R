# =============================================================================
# density_fO_uniform.R
# -----------------------------------------------------------------------------
# Defines the Uniform density for the reference times f_O, as described in
# Section 1.5 of Lin & Wang (2022):
#
#   f_O(u) = 1 / (1 - delta),   u in [delta/2, 1 - delta/2]
#   f_O(u) = 0,                  otherwise
#
# This is the DEFAULT choice in Lin & Wang (2022):
#   "The simplest choice: subjects enter the study completely at random,
#    with no preferred entry time. It is the default assumption in many
#    simulation studies, including Lin & Wang (2022)."
#
# STRICT POSITIVITY (Assumption 2, Lin & Wang 2022):
#   f_O(u) = 1/(1-delta) > 0 for all u in [delta/2, 1-delta/2]. Satisfied
#   trivially since delta in (0,1) implies 1-delta > 0.
#
# This function is called by:
#   eval_fO.R    — to evaluate f_O(u) at a vector of points
#
# Dependencies: base R only.
# =============================================================================


#' Uniform density for reference times f_O
#'
#' @param u     Numeric vector. Points at which to evaluate f_O.
#'              Values outside (delta/2, 1 - delta/2) return 0.
#' @param delta Numeric in (0, 1). Snippet length. Determines the support
#'              (delta/2, 1 - delta/2) and the constant density value
#'              1 / (1 - delta).
#'
#' @return Numeric vector of the same length as u with values of f_O(u).
#'
#' @examples
#' u <- seq(0, 1, length.out = 500)
#' f <- fO_uniform(u, delta = 0.4)
#' plot(u, f, type = "l", main = "f_O uniform, delta = 0.4")
#' # f_O = 1/(1-0.4) = 1.667 on (0.2, 0.8), zero elsewhere
#'
#' # Verify it integrates to 1
#' sum(f) * (1 / 500)   # approximately 1

#' @export
fO_uniform <- function(u, delta) {
  
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
  # DENSITY
  # f_O(u) = 1 / (hi - lo) = 1 / (1 - delta)   inside support
  # f_O(u) = 0                                   outside support
  # ---------------------------------------------------------------------------
  f <- ifelse(u >= lo & u <= hi, 1 / (hi - lo), 0)
  
  return(f)
}
