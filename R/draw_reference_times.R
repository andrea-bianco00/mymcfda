# =============================================================================
# draw_reference_times.R
# -----------------------------------------------------------------------------
# Implements Step 2 of the data-generating mechanism of Lin & Wang (2022),
# Section 1.7:
#
#   "Draw O_1, ..., O_n i.i.d. from a density f_O satisfying
#    f_O(u) > 0 for all u in [delta/2, 1 - delta/2].
#    The observation window of subject i is then [O_i - delta/2, O_i + delta/2].
#    The reference times are latent: they are never observed by the analyst."
#
# This function does NOT contain any mathematical or sampling logic.
# It delegates sampling entirely to sample_fO.R, which in turn calls
# the appropriate density-specific file.
#
# OUTPUT:
#   A list with two components:
#     - O_vec    : numeric vector of length n with the reference times
#                  O_1, ..., O_n, all in [delta/2, 1 - delta/2]
#     - windows  : data frame with columns A_i = O_i - delta/2 and
#                  B_i = O_i + delta/2, the observation window boundaries
#                  for each subject. These are LATENT in real data —
#                  returned here only for simulation diagnostics.
#
# REQUIRED SOURCES:
#   source("sample_fO.R")   which sources all sample_fO_*.R files
#
# Dependencies: base R only (via the sourced files).
# =============================================================================



#' Draw reference times O_1, ..., O_n (Step 2 of Lin & Wang 2022)
#'
#' @param n        Integer >= 1. Number of subjects.
#' @param delta    Numeric in (0, 1). Snippet length. Each subject's
#'                 observation window (O_i - delta/2, O_i + delta/2)
#'                 has length exactly delta.
#' @param type     Character. Distribution of reference times. One of:
#'                   "uniform"   (default, as in Lin & Wang 2022)
#'                   "truncnorm"
#'                   "truncbeta"
#'                   "truncexp"
#' @param params   Named list of extra parameters for the chosen density.
#'                 Passed directly to sample_fO(). See sample_fO.R for
#'                 the full parameter documentation of each family.
#'
#' @return A list with components:
#'   - O_vec   : numeric vector of length n with O_1, ..., O_n in
#'               (delta/2, 1 - delta/2). LATENT in real data.
#'   - windows : data frame with n rows and columns:
#'       - A : left endpoint  A_i = O_i - delta/2
#'       - B : right endpoint B_i = O_i + delta/2
#'     All windows satisfy (A_i, B_i) subset (0, 1). LATENT in real data.
#'
#' @examples
#' # Default: uniform (as in Lin & Wang 2022)
#' result <- draw_reference_times(n = 100, delta = 0.4)
#' range(result$O_vec)          # inside (0.2, 0.8)
#' range(result$windows$A)      # inside (0.0, 0.6)
#' range(result$windows$B)      # inside (0.4, 1.0)
#'
#' # Truncated Normal
#' result <- draw_reference_times(n = 100, delta = 0.4,
#'                                type   = "truncnorm",
#'                                params = list(mu_O = 0.5, sd_O = 0.1))
#'
#' # Truncated Exponential
#' result <- draw_reference_times(n = 100, delta = 0.4,
#'                                type   = "truncexp",
#'                                params = list(lambda = 2))

#' @export
draw_reference_times <- function(n,
                                 delta,
                                 type   = "uniform",
                                 params = list()) {
  
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
  # STEP 2 (Lin & Wang 2022, Section 1.7):
  # Draw O_1, ..., O_n i.i.d. from f_O on [delta/2, 1 - delta/2]
  # ---------------------------------------------------------------------------
  O_vec <- sample_fO(n      = n,
                     delta  = delta,
                     type   = type,
                     params = params)
  
  # ---------------------------------------------------------------------------
  # COMPUTE OBSERVATION WINDOWS [A_i, B_i] = [O_i - delta/2, O_i + delta/2]
  # These are latent in real data — returned for simulation diagnostics only.
  # ---------------------------------------------------------------------------
  windows <- data.frame(
    A = O_vec - delta / 2,
    B = O_vec + delta / 2
  )
  
  # ---------------------------------------------------------------------------
  # SANITY CHECK: all windows must lie inside [0, 1]
  # ---------------------------------------------------------------------------
  if (any(windows$A < -1e-10) || any(windows$B > 1 + 1e-10))
    stop("Internal error: some observation windows fall outside [0, 1].")
  
  return(list(
    O_vec   = O_vec,
    windows = windows
  ))
}
