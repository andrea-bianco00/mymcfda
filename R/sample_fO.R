# =============================================================================
# sample_fO.R
# -----------------------------------------------------------------------------
# Dispatcher: draws n i.i.d. samples from the density f_O of the reference
# times O_i, by calling the appropriate sampling file based on `type`.
#
# This function is the single entry point for sampling from f_O. It does NOT
# contain any mathematical logic itself — it delegates entirely to:
#
#   type = "uniform"   -> sample_fO_uniform.R   -> sample_fO_uniform(n, delta)
#   type = "truncnorm" -> sample_fO_truncnorm.R -> sample_fO_truncnorm(n, delta, mu_O, sd_O)
#   type = "truncbeta" -> sample_fO_truncbeta.R -> sample_fO_truncbeta(n, delta, a_O, b_O)
#   type = "truncexp"  -> sample_fO_truncexp.R  -> sample_fO_truncexp(n, delta, lambda)
#
# REQUIRED SOURCES:
#   source("sample_fO_uniform.R")
#   source("sample_fO_truncnorm.R")
#   source("sample_fO_truncbeta.R")
#   source("sample_fO_truncexp.R")
#
# RELATIONSHIP WITH eval_fO.R:
#   eval_fO.R    — evaluates f_O(u) at a vector of points  (density values)
#   sample_fO.R  — draws O_1, ..., O_n from f_O            (random samples)
#   The two dispatchers are perfectly symmetric and share the same interface.
#
# For the theoretical background on f_O and its role in the data-generating
# mechanism, see Section 1.5 of Lin & Wang (2022) and the individual
# sampling files above.
#
# Dependencies: base R only (via the sourced sampling files).
# =============================================================================



#' Draw n i.i.d. reference times from f_O on (delta/2, 1 - delta/2)
#'
#' @param n      Integer >= 1. Number of samples to draw.
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
#'               value defined in each sampling file is used.
#'
#' @return Numeric vector of length n with samples O_1, ..., O_n,
#'         all in (delta/2, 1 - delta/2).
#'
#' @examples
#' # Uniform (default in Lin & Wang 2022)
#' O <- sample_fO(n = 200, delta = 0.4)
#'
#' # Truncated Normal
#' O <- sample_fO(n = 200, delta = 0.4, type = "truncnorm",
#'                params = list(mu_O = 0.5, sd_O = 0.15))
#'
#' # Rescaled Beta
#' O <- sample_fO(n = 200, delta = 0.4, type = "truncbeta",
#'                params = list(a_O = 2, b_O = 3))
#'
#' # Truncated Exponential
#' O <- sample_fO(n = 200, delta = 0.4, type = "truncexp",
#'                params = list(lambda = 2))

#' @export
sample_fO <- function(n, delta, type = "uniform", params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n))
    stop("n must be a positive integer.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  if (!type %in% c("uniform", "truncnorm", "truncbeta", "truncexp"))
    stop("type must be one of: 'uniform', 'truncnorm', 'truncbeta', 'truncexp'.")
  if (!is.list(params))
    stop("params must be a named list.")
  
  # ---------------------------------------------------------------------------
  # DISPATCH to the appropriate sampling function
  # ---------------------------------------------------------------------------
  samples <- switch(type,
                    
                    "uniform" = {
                      sample_fO_uniform(n, delta)
                    },
                    
                    "truncnorm" = {
                      sample_fO_truncnorm(n, delta,
                                          mu_O = params$mu_O,
                                          sd_O = params$sd_O)
                    },
                    
                    "truncbeta" = {
                      sample_fO_truncbeta(n, delta,
                                          a_O = if (!is.null(params$a_O)) params$a_O else 2,
                                          b_O = if (!is.null(params$b_O)) params$b_O else 2)
                    },
                    
                    "truncexp" = {
                      sample_fO_truncexp(n, delta,
                                         lambda = params$lambda)
                    }
  )
  
  return(samples)
}
