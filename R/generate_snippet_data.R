# =============================================================================
# generate_snippet_data.R
# -----------------------------------------------------------------------------
# Master function for the data-generating mechanism of Lin & Wang (2022).
# Orchestrates the three steps in a single call:
#
#   Step 1: Draw m_1, ..., m_n ~ Poisson(m_mean) truncated at >= 2
#           (via draw_m_vec.R)
#
#   Step 2: Draw O_1, ..., O_n i.i.d. from f_O on [delta/2, 1 - delta/2]
#           (via draw_reference_times.R)
#
#   Step 3: For each subject i, draw T_i1, ..., T_im_i from f_{T|O}
#           by sampling from f_0 on [0, delta] and shifting by O_i - delta/2
#           (via draw_measurement_times.R)
#
# The user controls:
#   - n         : number of subjects
#   - m_mean    : average number of measurements per snippet (XXX)
#   - delta     : snippet length
#   - seed      : reproducibility
#   - fO_type   : density family for reference times O_i
#   - fO_params : parameters for the chosen f_O family
#   - f0_type   : density family for measurement times T_ij
#   - f0_params : parameters for the chosen f_0 family
#
# OUTPUT:
#   A list with all generated data AND all input parameters stored,
#   so that the simulation is fully reproducible and self-documenting.
#
# REQUIRED SOURCES:
#   source("draw_m_vec.R")
#   source("draw_reference_times.R")    (which sources sample_fO.R -> ...)
#   source("draw_measurement_times.R")  (which sources sample_fzero.R -> ...)
#
# Dependencies: base R only (via the sourced files).
# =============================================================================



#' Generate complete snippet observation structure
#'
#' @param n         Integer >= 1. Number of subjects.
#' @param m_mean    Numeric > 0. Average number of measurements per snippet.
#'                  Passed as the Poisson mean (before truncation at >= 2).
#'                  Default: 4 (as in Lin & Wang 2022).
#' @param delta     Numeric in (0, 1). Snippet length.
#'                  Default: 0.25 (as in Lin & Wang 2022).
#' @param seed      Integer or NULL. Random seed for reproducibility.
#'                  If NULL, no seed is set. Default: NULL.
#' @param fO_type   Character. Density family for reference times O_i.
#'                  One of: "uniform" (default), "truncnorm", "truncbeta",
#'                  "truncexp".
#' @param fO_params Named list. Parameters for the chosen f_O family.
#'                  See sample_fO.R for documentation of each family.
#'                  Default: list() (uses each family's defaults).
#' @param f0_type   Character. Density family for measurement times T_ij.
#'                  One of: "uniform" (default), "truncnorm", "truncbeta",
#'                  "truncexp".
#' @param f0_params Named list. Parameters for the chosen f_0 family.
#'                  See sample_fzero.R for documentation of each family.
#'                  Default: list() (uses each family's defaults).
#'
#' @return A list with components:
#'   - T_list    : list of n numeric vectors. T_list[[i]] contains the m_i
#'                 measurement times for subject i, sorted in increasing
#'                 order, all in [0, 1].
#'   - m_vec     : integer vector of length n. Number of measurements per
#'                 subject (m_1, ..., m_n), all >= 2.
#'   - O_vec     : numeric vector of length n. Reference times O_1, ..., O_n
#'                 in [delta/2, 1 - delta/2]. LATENT in real data.
#'   - windows   : data frame with columns A (left endpoint) and B (right
#'                 endpoint) of each subject's observation window. LATENT.
#'   - config    : list containing all input parameters (n, m_mean, delta,
#'                 seed, fO_type, fO_params, f0_type, f0_params) for full
#'                 reproducibility.
#'   - summary   : list with quick diagnostics:
#'       - n              : number of subjects
#'       - m_mean_actual  : actual sample mean of m_vec
#'       - m_range        : range of m_vec (min, max)
#'       - total_obs      : total number of observations sum(m_vec)
#'       - snippet_widths : summary stats of snippet widths (min, mean, max)
#'
#' @examples
#' # === Example 1: Default (as in Lin & Wang 2022) ===
#' dat <- generate_snippet_data(n = 200, m_mean = 4, delta = 0.25, seed = 42)
#' str(dat$T_list[1:3])        # first 3 subjects' measurement times
#' dat$summary                  # quick diagnostics
#'
#' # === Example 2: Truncated Normal for both O_i and T_ij ===
#' dat2 <- generate_snippet_data(
#'   n         = 100,
#'   m_mean    = 6,
#'   delta     = 0.30,
#'   seed      = 123,
#'   fO_type   = "truncnorm",
#'   fO_params = list(mu_O = 0.5, sd_O = 0.15),
#'   f0_type   = "truncnorm",
#'   f0_params = list(mu_0 = 0.15, sd_0 = 0.06)
#' )
#'
#' # === Example 3: Beta reference times, Exponential measurement times ===
#' dat3 <- generate_snippet_data(
#'   n         = 50,
#'   m_mean    = 3,
#'   delta     = 0.40,
#'   seed      = 7,
#'   fO_type   = "truncbeta",
#'   fO_params = list(a_O = 2, b_O = 5),
#'   f0_type   = "truncexp",
#'   f0_params = list(lambda = 3)
#' )

generate_snippet_data <- function(n,
                                  m_mean    = 4,
                                  delta     = 0.25,
                                  seed      = NULL,
                                  fO_type   = "uniform",
                                  fO_params = list(),
                                  f0_type   = "uniform",
                                  f0_params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n))
    stop("n must be a positive integer.")
  if (!is.numeric(m_mean) || length(m_mean) != 1 || m_mean <= 0)
    stop("m_mean must be a single positive number.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1))
    stop("seed must be a single integer or NULL.")
  
  valid_types <- c("uniform", "truncnorm", "truncbeta", "truncexp")
  if (!fO_type %in% valid_types)
    stop("fO_type must be one of: ", paste(valid_types, collapse = ", "))
  if (!f0_type %in% valid_types)
    stop("f0_type must be one of: ", paste(valid_types, collapse = ", "))
  if (!is.list(fO_params))
    stop("fO_params must be a named list.")
  if (!is.list(f0_params))
    stop("f0_params must be a named list.")
  
  # ---------------------------------------------------------------------------
  # SET SEED (if provided)
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  # ---------------------------------------------------------------------------
  # STEP 1: Draw m_1, ..., m_n ~ Poisson(m_mean) truncated at >= 2
  # ---------------------------------------------------------------------------
  m_vec <- draw_m_vec(n = n, m_mean = m_mean)
  
  # ---------------------------------------------------------------------------
  # STEP 2: Draw O_1, ..., O_n i.i.d. from f_O
  # ---------------------------------------------------------------------------
  ref <- draw_reference_times(n      = n,
                              delta  = delta,
                              type   = fO_type,
                              params = fO_params)
  
  # ---------------------------------------------------------------------------
  # STEP 3: Draw T_i1, ..., T_im_i for each subject from f_0
  # ---------------------------------------------------------------------------
  T_list <- draw_measurement_times(O_vec  = ref$O_vec,
                                   m_vec  = m_vec,
                                   delta  = delta,
                                   type   = f0_type,
                                   params = f0_params)
  
  # sort each subject's times in increasing order
  T_list <- lapply(T_list, sort)
  
  # ---------------------------------------------------------------------------
  # SUMMARY DIAGNOSTICS
  # ---------------------------------------------------------------------------
  snippet_widths <- sapply(T_list, function(t_i) max(t_i) - min(t_i))
  
  summary_info <- list(
    n              = n,
    m_mean_actual  = mean(m_vec),
    m_range        = range(m_vec),
    total_obs      = sum(m_vec),
    snippet_widths = c(min  = min(snippet_widths),
                       mean = mean(snippet_widths),
                       max  = max(snippet_widths))
  )
  
  # ---------------------------------------------------------------------------
  # STORE CONFIG for full reproducibility
  # ---------------------------------------------------------------------------
  config <- list(
    n         = n,
    m_mean    = m_mean,
    delta     = delta,
    seed      = seed,
    fO_type   = fO_type,
    fO_params = fO_params,
    f0_type   = f0_type,
    f0_params = f0_params
  )
  
  # ---------------------------------------------------------------------------
  # RETURN
  # ---------------------------------------------------------------------------
  return(list(
    T_list  = T_list,
    m_vec   = m_vec,
    O_vec   = ref$O_vec,
    windows = ref$windows,
    config  = config,
    summary = summary_info
  ))
}
