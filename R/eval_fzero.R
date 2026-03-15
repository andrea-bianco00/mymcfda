# =============================================================================
# eval_fzero.R
# -----------------------------------------------------------------------------
# Dispatcher: evaluates the density f_0 of the measurement times T_ij at one
# or more points s, by calling the appropriate density file based on `type`.
#
# This function is the single entry point for evaluating f_0. It does NOT
# contain any mathematical logic itself — it delegates entirely to:
#
#   type = "uniform"   -> density_fzero_uniform.R   -> fzero_uniform(s, delta)
#   type = "truncnorm" -> density_fzero_truncnorm.R -> fzero_truncnorm(s, delta, mu_0, sd_0)
#   type = "truncbeta" -> density_fzero_truncbeta.R -> fzero_truncbeta(s, delta, a_0, b_0)
#   type = "truncexp"  -> density_fzero_truncexp.R  -> fzero_truncexp(s, delta, lambda)
#
# REQUIRED SOURCES:
#   source("density_fzero_uniform.R")
#   source("density_fzero_truncnorm.R")
#   source("density_fzero_truncbeta.R")
#   source("density_fzero_truncexp.R")
#
# RELATIONSHIP WITH eval_fO.R:
#   eval_fO.R    — evaluates f_O(u) on [delta/2, 1-delta/2]  (reference times)
#   eval_fzero.R — evaluates f_0(s) on [0, delta]            (measurement times)
#   Identical interface, different support.
#
# Dependencies: base R only (via the sourced density files).
# =============================================================================



#' Evaluate the density f_0 of measurement times at s
#'
#' @param s      Numeric vector. Points at which to evaluate f_0.
#'               Values outside [0, delta] return 0.
#' @param delta  Numeric in (0, 1). Snippet length.
#' @param type   Character. Family of f_0. One of:
#'                 "uniform"   (default, Lin & Wang 2022)
#'                 "truncnorm"
#'                 "truncbeta"
#'                 "truncexp"
#' @param params Named list of parameters for the chosen family:
#'                 "uniform"  : no parameters needed
#'                 "truncnorm": list(mu_0, sd_0)
#'                 "truncbeta": list(a_0, b_0)
#'                 "truncexp" : list(lambda)
#'               If params is empty or a parameter is missing, the default
#'               value defined in each density file is used.
#'
#' @return Numeric vector of the same length as s with values of f_0(s).
#'
#' @examples
#' s <- seq(-0.1, 0.6, length.out = 500)
#'
#' # Uniform (default)
#' f1 <- eval_fzero(s, delta = 0.4)
#'
#' # Truncated Normal
#' f2 <- eval_fzero(s, delta = 0.4, type = "truncnorm",
#'                  params = list(mu_0 = 0.2, sd_0 = 0.08))
#'
#' # Rescaled Beta
#' f3 <- eval_fzero(s, delta = 0.4, type = "truncbeta",
#'                  params = list(a_0 = 2, b_0 = 3))
#'
#' # Truncated Exponential
#' f4 <- eval_fzero(s, delta = 0.4, type = "truncexp",
#'                  params = list(lambda = 3))
#'
#' # Plot all four
#' plot(s, f1, type = "l", ylim = c(0, 5), main = "f_0 families")
#' lines(s, f2, col = "red")
#' lines(s, f3, col = "blue")
#' lines(s, f4, col = "green")
#' legend("topright", legend = c("uniform","truncnorm","truncbeta","truncexp"),
#'        col = c("black","red","blue","green"), lty = 1)

eval_fzero <- function(s, delta, type = "uniform", params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s))
    stop("s must be a numeric vector.")
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
                fzero_uniform(s, delta)
              },
              
              "truncnorm" = {
                fzero_truncnorm(s, delta,
                                mu_0 = params$mu_0,
                                sd_0 = params$sd_0)
              },
              
              "truncbeta" = {
                fzero_truncbeta(s, delta,
                                a_0 = if (!is.null(params$a_0)) params$a_0 else 2,
                                b_0 = if (!is.null(params$b_0)) params$b_0 else 2)
              },
              
              "truncexp" = {
                fzero_truncexp(s, delta,
                               lambda = params$lambda)
              }
  )
  
  return(f)
}
