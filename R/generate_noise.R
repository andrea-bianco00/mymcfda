# =============================================================================
# generate_noise.R
# -----------------------------------------------------------------------------
# Generates measurement noise epsilon_{ij} for each subject i and each
# observation j, as described in Step 5 of Lin & Wang (2022):
#
#   E[epsilon_{ij}] = 0,  Var(epsilon_{ij}) = sigma_0^2,  E[epsilon_{ij}^4] < Inf
#
# This function supports two modes:
#
#   1. HOMOSCEDASTIC (default, as in the paper):
#      sigma0 is a single positive number. All observations have the same
#      noise variance sigma_0^2, regardless of time or subject.
#
#   2. HETEROSCEDASTIC (time-dependent, extension in Remark of Section 3):
#      sigma0 is a function sigma0(t) : [0,1] -> R+. The noise variance
#      at observation T_ij is sigma_0^2(T_ij), which varies with time.
#
# The distribution of the noise is controlled by noise_type and noise_params,
# which are passed to build_noise_fn() to construct a standardised noise
# function (mean 0, variance 1). The function then scales by sigma_0
# (or sigma_0(T_ij) in the heteroscedastic case).
#
# The parameter noise_mean is explicitly declared to make the zero-mean
# assumption visible and transparent. It must always be 0, as required
# by the observation model Y_ij = X_i(T_ij) + epsilon_ij: a nonzero
# mean would be confounded with the process mean mu(t) and could not
# be identified from the data.
#
# ROLE IN THE PIPELINE:
#   build_noise_fn         -->  noise_fn (standardised, called internally)
#   generate_snippet_data  -->  T_list
#   generate_noise         -->  noise_list  (THIS FILE)
#   assemble_observations  -->  Y_list = mu_values + Z_list + noise_list
#
# REQUIRED SOURCES:
#   source("build_noise_fn.R")
#
# INPUT:
#   T_list       : list of n numeric vectors. T_list[[i]] contains the m_i
#                  measurement times for subject i, all in [0, 1].
#   sigma0       : either a single non-negative number (homoscedastic case)
#                  or a vectorised function sigma0(t) : [0,1] -> R+ returning
#                  the noise standard deviation at time t (heteroscedastic).
#   noise_mean   : numeric. Mean of the noise distribution. Must be 0.
#                  This parameter is declared explicitly for transparency:
#                  the zero-mean assumption is a fundamental requirement
#                  of the model, not a hidden default. Default: 0.
#   noise_type   : character. Distribution family for the noise.
#                  One of: "gaussian" (default), "t", "uniform", "laplace",
#                  "custom". Passed to build_noise_fn().
#   noise_params : named list of parameters for the chosen distribution.
#                  Passed to build_noise_fn(). See build_noise_fn.R for
#                  documentation of each family's parameters.
#                  Default: list() (uses family defaults).
#
# OUTPUT:
#   A list of n numeric vectors. noise_list[[i]] is a vector of length m_i
#   containing epsilon_{i1}, ..., epsilon_{im_i}.
#
# EXAMPLES:
#   T_list <- list(c(0.1, 0.2, 0.3), c(0.5, 0.6))
#
#   # Homoscedastic Gaussian noise (default, as in the paper)
#   noise <- generate_noise(T_list, sigma0 = 0.5, noise_mean = 0)
#
#   # Homoscedastic t-distributed noise with df = 5
#   noise_t <- generate_noise(T_list, sigma0 = 0.5, noise_mean = 0,
#                             noise_type = "t", noise_params = list(df = 5))
#
#   # Heteroscedastic Gaussian noise (variance increases with time)
#   sigma0_fn <- function(t) 0.3 + 0.4 * t
#   noise_het <- generate_noise(T_list, sigma0 = sigma0_fn, noise_mean = 0)
#
# Dependencies: base R only (via build_noise_fn.R).
# =============================================================================



#' Generate measurement noise for each subject
#'
#' @param T_list       List of n numeric vectors of measurement times.
#' @param sigma0       Either a single non-negative number (homoscedastic)
#'                     or a vectorised function sigma0(t) -> R+
#'                     (heteroscedastic).
#' @param noise_mean   Numeric. Mean of the noise. Must be 0. Declared
#'                     explicitly for transparency. Default: 0.
#' @param noise_type   Character. Distribution family. Default: "gaussian".
#' @param noise_params Named list of parameters for the distribution.
#'                     Default: list().
#'
#' @return A list of n numeric vectors. noise_list((i)) contains the m_i
#'         noise values for subject i.

#' @export
generate_noise <- function(T_list,
                           sigma0,
                           noise_mean   = 0,
                           noise_type   = "gaussian",
                           noise_params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || length(T_list) < 1)
    stop("T_list must be a non-empty list of numeric vectors.")
  
  # --- noise_mean must be 0 ---
  if (!is.numeric(noise_mean) || length(noise_mean) != 1)
    stop("noise_mean must be a single number.")
  
  if (noise_mean != 0)
    stop("noise_mean must be 0. A nonzero noise mean would be confounded ",
         "with the process mean mu(t) and cannot be identified from the ",
         "data. See Step 5 of Lin & Wang (2022).")
  
  # ---------------------------------------------------------------------------
  # BUILD THE STANDARDISED NOISE FUNCTION via build_noise_fn
  # ---------------------------------------------------------------------------
  noise_fn <- build_noise_fn(type = noise_type, params = noise_params)
  
  # ---------------------------------------------------------------------------
  # DETERMINE MODE: homoscedastic (number) or heteroscedastic (function)
  # ---------------------------------------------------------------------------
  if (is.numeric(sigma0) && length(sigma0) == 1) {
    
    # --- HOMOSCEDASTIC MODE ---
    if (sigma0 < 0)
      stop("sigma0 must be non-negative (sigma0 = 0 means no noise).")
    
    mode <- "homoscedastic"
    
  } else if (is.function(sigma0)) {
    
    # --- HETEROSCEDASTIC MODE ---
    test_input  <- c(0.0, 0.5, 1.0)
    test_output <- tryCatch(sigma0(test_input),
                            error = function(e) {
                              stop("sigma0 function failed on test input. ",
                                   "Error: ", e$message)
                            })
    
    if (!is.numeric(test_output) || length(test_output) != length(test_input))
      stop("sigma0 function must be vectorised: sigma0(t) must return a ",
           "vector of the same length as t.")
    
    if (any(test_output < 0))
      stop("sigma0 function must return non-negative values for all t.")
    
    mode <- "heteroscedastic"
    
  } else {
    stop("sigma0 must be either a single non-negative number ",
         "(homoscedastic) or a vectorised function of t (heteroscedastic).")
  }
  
  # ---------------------------------------------------------------------------
  # GENERATE NOISE
  # ---------------------------------------------------------------------------
  n <- length(T_list)
  noise_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    t_i <- T_list[[i]]
    m_i <- length(t_i)
    
    if (!is.numeric(t_i) || m_i < 1)
      stop("T_list[[", i, "]] must be a non-empty numeric vector.")
    
    # generate m_i standardised noise values (mean 0, variance 1)
    eta_i <- noise_fn(m_i)
    
    # scale by sigma0
    if (mode == "homoscedastic") {
      noise_list[[i]] <- sigma0 * eta_i
    } else {
      sigma0_i <- sigma0(t_i)
      noise_list[[i]] <- sigma0_i * eta_i
    }
  }
  
  return(noise_list)
}
