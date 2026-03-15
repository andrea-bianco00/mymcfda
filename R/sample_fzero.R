# =============================================================================
# sample_fzero.R
# -----------------------------------------------------------------------------
# Dispatcher: draws m i.i.d. samples from the density f_0 of the measurement
# times S_j, by calling the appropriate sampling file based on `type`.
#
# This function is the single entry point for sampling from f_0. It does NOT
# contain any mathematical logic itself — it delegates entirely to:
#
#   type = "uniform"   -> sample_fzero_uniform.R   -> sample_fzero_uniform(m, delta)
#   type = "truncnorm" -> sample_fzero_truncnorm.R -> sample_fzero_truncnorm(m, delta, mu_0, sd_0)
#   type = "truncbeta" -> sample_fzero_truncbeta.R -> sample_fzero_truncbeta(m, delta, a_0, b_0)
#   type = "truncexp"  -> sample_fzero_truncexp.R  -> sample_fzero_truncexp(m, delta, lambda)
#
# REQUIRED SOURCES:
#   source("sample_fzero_uniform.R")
#   source("sample_fzero_truncnorm.R")
#   source("sample_fzero_truncbeta.R")
#   source("sample_fzero_truncexp.R")
#
# RELATIONSHIP WITH eval_fzero.R:
#   eval_fzero.R    — evaluates f_0(s) at a vector of points  (density values)
#   sample_fzero.R  — draws S_1, ..., S_m from f_0             (random samples)
#   The two dispatchers are perfectly symmetric and share the same interface.
#
# RELATIONSHIP WITH sample_fO.R:
#   sample_fO.R    — draws O_1, ..., O_n from f_O on [delta/2, 1 - delta/2]
#   sample_fzero.R — draws S_1, ..., S_m from f_0 on [0, delta]
#   Both dispatchers share the same structure and type/params interface.
#
# For the theoretical background on f_0 and its role in the data-generating
# mechanism, see Section 1.6 of Lin & Wang (2022) and the individual
# sampling files above.
#
# Dependencies: base R only (via the sourced sampling files).
# =============================================================================



#' Draw m i.i.d. measurement times from f_0 on (0, delta)
#'
#' @param m      Integer >= 1. Number of samples to draw (measurements per
#'               subject).
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
#'               value defined in each sampling file is used.
#'
#' @return Numeric vector of length m with samples S_1, ..., S_m,
#'         all in (0, delta). These are RELATIVE positions within the
#'         window — add O_i - delta/2 to obtain the actual T_ij.
#'
#' @examples
#' # Uniform (default in Lin & Wang 2022)
#' S <- sample_fzero(m = 5, delta = 0.4)
#'
#' # Truncated Normal
#' S <- sample_fzero(m = 5, delta = 0.4, type = "truncnorm",
#'                   params = list(mu_0 = 0.2, sd_0 = 0.08))
#'
#' # Rescaled Beta
#' S <- sample_fzero(m = 5, delta = 0.4, type = "truncbeta",
#'                   params = list(a_0 = 2, b_0 = 3))
#'
#' # Truncated Exponential
#' S <- sample_fzero(m = 5, delta = 0.4, type = "truncexp",
#'                   params = list(lambda = 2))

#' @export
sample_fzero <- function(m, delta, type = "uniform", params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(m) || length(m) != 1 || m < 1 || m != round(m))
    stop("m must be a positive integer.")
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
                      sample_fzero_uniform(m, delta)
                    },
                    
                    "truncnorm" = {
                      sample_fzero_truncnorm(m, delta,
                                             mu_0 = params$mu_0,
                                             sd_0 = params$sd_0)
                    },
                    
                    "truncbeta" = {
                      sample_fzero_truncbeta(m, delta,
                                             a_0 = if (!is.null(params$a_0)) params$a_0 else 2,
                                             b_0 = if (!is.null(params$b_0)) params$b_0 else 2)
                    },
                    
                    "truncexp" = {
                      sample_fzero_truncexp(m, delta,
                                            lambda = params$lambda)
                    }
  )
  
  return(samples)
}
