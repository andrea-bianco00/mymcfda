# =============================================================================
# eval_fO.R
# -----------------------------------------------------------------------------
# Dispatcher: evaluates the density f_O of the reference times O_i at one
# or more points u, by calling the appropriate density file based on `type`.
#
# This function is the single entry point for evaluating f_O. It does NOT
# contain any mathematical logic itself — it delegates entirely to:
#
#   type = "uniform"   -> density_fO_uniform.R   -> fO_uniform(u, delta)
#   type = "truncnorm" -> density_fO_truncnorm.R -> fO_truncnorm(u, delta, mu_O, sd_O)
#   type = "truncbeta" -> density_fO_truncbeta.R -> fO_truncbeta(u, delta, a_O, b_O)
#   type = "truncexp"  -> density_fO_truncexp.R  -> fO_truncexp(u, delta, lambda)
#
# REQUIRED SOURCES:
#   source("density_fO_uniform.R")
#   source("density_fO_truncnorm.R")
#   source("density_fO_truncbeta.R")
#   source("density_fO_truncexp.R")
#
# For the theoretical background on f_O and its role in the data-generating
# mechanism, see Section 1.5 of Lin & Wang (2022) and the individual
# density files above.
#
# Dependencies: base R only (via the sourced density files).
# =============================================================================



#' Evaluate the density f_O of reference times at u
#'
#' @param u      Numeric vector. Points at which to evaluate f_O.
#'               Values outside (delta/2, 1 - delta/2) return 0.
#' @param delta  Numeric in (0, 1). Snippet length.
#' @param type   Character. Family of f_O. One of:
#'                 "uniform"   (default, Lin & Wang 2022)
#'                 "truncnorm"
#'                 "truncbeta"
#'                 "truncexp"
#' @param params Named list of parameters for the chosen family:
#'                 "uniform"  : no parameters needed
#'                 "truncnorm": list(mu_O, sd_O)
#'                 "truncbeta": list(a_O, b_O)
#'                 "truncexp" : list(lambda)
#'               If params is empty or a parameter is missing, the default
#'               value defined in each density file is used.
#'
#' @return Numeric vector of the same length as u with values of f_O(u).
#'
#' @examples
#' u <- seq(0, 1, length.out = 500)
#'
#' # Uniform (default)
#' f1 <- eval_fO(u, delta = 0.4)
#'
#' # Truncated Normal
#' f2 <- eval_fO(u, delta = 0.4, type = "truncnorm",
#'               params = list(mu_O = 0.5, sd_O = 0.15))
#'
#' # Rescaled Beta
#' f3 <- eval_fO(u, delta = 0.4, type = "truncbeta",
#'               params = list(a_O = 2, b_O = 3))
#'
#' # Truncated Exponential
#' f4 <- eval_fO(u, delta = 0.4, type = "truncexp",
#'               params = list(lambda = 2))
#'
#' # Plot all four
#' plot(u, f1, type = "l", ylim = c(0, 4), main = "f_O families")
#' lines(u, f2, col = "red")
#' lines(u, f3, col = "blue")
#' lines(u, f4, col = "green")
#' legend("topright", legend = c("uniform","truncnorm","truncbeta","truncexp"),
#'        col = c("black","red","blue","green"), lty = 1)

#' @export
eval_fO <- function(u, delta, type = "uniform", params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(u))
    stop("u must be a numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  if (!type %in% c("uniform", "truncnorm", "truncbeta", "truncexp"))
    stop("type must be one of: 'uniform', 'truncnorm', 'truncbeta', 'truncexp'.")
  if (!is.list(params))
    stop("params must be a named list.")
  
  # ---------------------------------------------------------------------------
  # DISPATCH to the appropriate density function
  # ---------------------------------------------------------------------------
  f <- switch(type,
              
              "uniform" = {
                fO_uniform(u, delta)
              },
              
              "truncnorm" = {
                fO_truncnorm(u, delta,
                             mu_O = params$mu_O,
                             sd_O = params$sd_O)
              },
              
              "truncbeta" = {
                fO_truncbeta(u, delta,
                             a_O = if (!is.null(params$a_O)) params$a_O else 2,
                             b_O = if (!is.null(params$b_O)) params$b_O else 2)
              },
              
              "truncexp" = {
                fO_truncexp(u, delta,
                            lambda = params$lambda)
              }
  )
  
  return(f)
}
