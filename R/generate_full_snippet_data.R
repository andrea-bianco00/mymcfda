# =============================================================================
# generate_full_snippet_data.R
# -----------------------------------------------------------------------------
# Full data-generating pipeline for functional snippets, with optional
# block-specific seeds for advanced reproducibility control.
#
# FINAL MODEL:
#
#   Y_ij = mu(T_ij) + Z_i(T_ij) + epsilon_ij
#
# BLOCKS:
#   1. design   : snippet structure (m_i, O_i, T_ij)
#   2. process  : centred process Z_i(T_ij)
#   3. noise    : measurement noise epsilon_ij
#
# REPRODUCIBILITY MODES:
#   - seed      : simple global seed
#   - seed_list : advanced block-specific seeds
#
# If seed_list is provided, it has priority over seed.
#
# NOISE MODES:
#   - noise_structure = "homoscedastic"
#       -> sigma0 must be a single non-negative numeric value
#   - noise_structure = "heteroscedastic"
#       -> sigma0 must be a vectorised function sigma0(t)
#
# REQUIRED SOURCES:
#   source("generate_snippet_data.R")
#   source("eval_mean_at_times.R")
#   source("generate_centered_process.R")
#   source("generate_noise.R")
#   source("assemble_observations.R")
# =============================================================================



#' Generate full functional snippet data
#'
#' @param n               Integer >= 1. Number of subjects.
#' @param m_mean          Numeric > 0. Mean of truncated Poisson for m_i.
#'                        Default: 4.
#' @param delta           Numeric in (0,1). Snippet length. Default: 0.25.
#' @param seed            Integer or NULL. Simple global seed. Default: NULL.
#' @param seed_list       Named list of optional block-specific seeds. Possible
#'                        names:
#'                          - design
#'                          - process
#'                          - noise
#'                        If provided, it has priority over seed.
#'                        Default: NULL.
#'
#' @param fO_type         Character. Density family for O_i. Default: "uniform".
#' @param fO_params       Named list. Parameters for f_O. Default: list().
#' @param f0_type         Character. Density family for within-window times.
#'                        Default: "uniform".
#' @param f0_params       Named list. Parameters for f_0. Default: list().
#'
#' @param mu_fn           Vectorised mean function mu(t).
#'                        Default: 2*t^2*cos(2*pi*t).
#'
#' @param road            Character. Must be "kl" or "decomposition".
#'
#' @param eigenfn_list    List of eigenfunctions, used only if road = "kl".
#' @param eigenval_vec    Numeric vector of eigenvalues, used only if road = "kl".
#' @param score_type      Character. Score distribution for KL road.
#'                        Default: "gaussian".
#' @param score_params    Named list. Score parameters. Default: list().
#'
#' @param sigmaX_fn       Function sigma_X(t), used only if
#'                        road = "decomposition".
#' @param corr_type       Character. Correlation family, used only if
#'                        road = "decomposition".
#' @param corr_params     Named list. Correlation parameters, used only if
#'                        road = "decomposition". Default: list().
#' @param jitter          Small non-negative diagonal inflation for decomposition
#'                        road. Default: 1e-10.
#'
#' @param noise_structure Character. Must be either "homoscedastic" or
#'                        "heteroscedastic".
#'                        - "homoscedastic"  -> sigma0 must be numeric
#'                        - "heteroscedastic"-> sigma0 must be a function
#'                        Default: "homoscedastic".
#' @param sigma0          If noise_structure = "homoscedastic", sigma0 must be
#'                        a single non-negative numeric value.
#'                        If noise_structure = "heteroscedastic", sigma0 must be
#'                        a vectorised function sigma0(t).
#'                        Default: 0.
#' @param noise_mean      Numeric. Must be 0. Default: 0.
#' @param noise_type      Character. Noise distribution. Default: "gaussian".
#' @param noise_params    Named list. Noise parameters. Default: list().
#'
#' @return A list with:
#'   - Y_list
#'   - T_list
#'   - mu_values
#'   - Z_list
#'   - noise_list
#'   - m_vec
#'   - O_vec
#'   - windows
#'   - process_info
#'   - design_info
#'   - config
#'   - summary

generate_full_snippet_data <- function(
    n,
    m_mean          = 4,
    delta           = 0.25,
    seed            = NULL,
    seed_list       = NULL,
    fO_type         = "uniform",
    fO_params       = list(),
    f0_type         = "uniform",
    f0_params       = list(),
    mu_fn           = function(t) 2 * t^2 * cos(2 * pi * t),
    road,
    eigenfn_list    = NULL,
    eigenval_vec    = NULL,
    score_type      = "gaussian",
    score_params    = list(),
    sigmaX_fn       = NULL,
    corr_type       = NULL,
    corr_params     = list(),
    jitter          = 1e-10,
    noise_structure = c("homoscedastic", "heteroscedastic"),
    sigma0          = 0,
    noise_mean      = 0,
    noise_type      = "gaussian",
    noise_params    = list()
) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (missing(road)) {
    stop("You must specify road = 'kl' or road = 'decomposition'.")
  }
  
  if (!is.character(road) || length(road) != 1 ||
      !road %in% c("kl", "decomposition")) {
    stop("road must be exactly one of: 'kl', 'decomposition'.")
  }
  
  if (!is.function(mu_fn)) {
    stop("mu_fn must be a function.")
  }
  
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1)) {
    stop("seed must be a single number or NULL.")
  }
  
  if (!is.null(seed_list)) {
    if (!is.list(seed_list)) {
      stop("seed_list must be a named list or NULL.")
    }
    
    allowed_seed_names <- c("design", "process", "noise")
    bad_names <- setdiff(names(seed_list), allowed_seed_names)
    
    if (length(bad_names) > 0) {
      stop("seed_list contains invalid names: ",
           paste(bad_names, collapse = ", "),
           ". Allowed names are: design, process, noise.")
    }
    
    for (nm in names(seed_list)) {
      if (!is.null(seed_list[[nm]]) &&
          (!is.numeric(seed_list[[nm]]) || length(seed_list[[nm]]) != 1)) {
        stop("seed_list$", nm, " must be a single number or NULL.")
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # NOISE STRUCTURE VALIDATION
  # ---------------------------------------------------------------------------
  noise_structure <- match.arg(noise_structure)
  
  if (noise_structure == "homoscedastic") {
    if (!is.numeric(sigma0) || length(sigma0) != 1 || sigma0 < 0) {
      stop("If noise_structure = 'homoscedastic', sigma0 must be a single ",
           "non-negative numeric value.")
    }
  }
  
  if (noise_structure == "heteroscedastic") {
    if (!is.function(sigma0)) {
      stop("If noise_structure = 'heteroscedastic', sigma0 must be a ",
           "vectorised function sigma0(t).")
    }
    
    test_input  <- c(0, 0.5, 1)
    test_output <- tryCatch(
      sigma0(test_input),
      error = function(e) {
        stop("The sigma0 function failed on test input c(0, 0.5, 1). ",
             "Error: ", e$message)
      }
    )
    
    if (!is.numeric(test_output) || length(test_output) != length(test_input)) {
      stop("If noise_structure = 'heteroscedastic', sigma0 must be vectorised: ",
           "sigma0(t) must return a numeric vector of the same length as t.")
    }
    
    if (any(!is.finite(test_output)) || any(test_output < 0)) {
      stop("The sigma0 function must return finite non-negative values.")
    }
  }
  
  if (!is.numeric(noise_mean) || length(noise_mean) != 1 || noise_mean != 0) {
    stop("noise_mean must be exactly 0.")
  }
  
  # ---------------------------------------------------------------------------
  # SEED RESOLUTION
  # ---------------------------------------------------------------------------
  # Priority:
  #   1. seed_list (advanced control)
  #   2. seed      (simple global seed)
  #   3. NULL      (no seed)
  # ---------------------------------------------------------------------------
  resolved_seeds <- list(
    design  = NULL,
    process = NULL,
    noise   = NULL
  )
  
  if (!is.null(seed_list)) {
    for (nm in names(seed_list)) {
      resolved_seeds[[nm]] <- seed_list[[nm]]
    }
  } else if (!is.null(seed)) {
    resolved_seeds$design  <- seed
    resolved_seeds$process <- seed + 1
    resolved_seeds$noise   <- seed + 2
  }
  
  # ---------------------------------------------------------------------------
  # STEP 1: GENERATE SNIPPET DESIGN
  # ---------------------------------------------------------------------------
  if (!is.null(resolved_seeds$design)) {
    set.seed(resolved_seeds$design)
  }
  
  design_info <- generate_snippet_data(
    n         = n,
    m_mean    = m_mean,
    delta     = delta,
    seed      = NULL,
    fO_type   = fO_type,
    fO_params = fO_params,
    f0_type   = f0_type,
    f0_params = f0_params
  )
  
  T_list  <- design_info$T_list
  m_vec   <- design_info$m_vec
  O_vec   <- design_info$O_vec
  windows <- design_info$windows
  
  # ---------------------------------------------------------------------------
  # STEP 2: EVALUATE MEAN
  # ---------------------------------------------------------------------------
  mu_values <- eval_mean_at_times(
    T_list = T_list,
    mu_fn  = mu_fn
  )
  
  # ---------------------------------------------------------------------------
  # STEP 3: GENERATE CENTRED PROCESS
  # ---------------------------------------------------------------------------
  if (!is.null(resolved_seeds$process)) {
    set.seed(resolved_seeds$process)
  }
  
  process_info <- generate_centered_process(
    T_list        = T_list,
    road          = road,
    eigenfn_list  = eigenfn_list,
    eigenval_vec  = eigenval_vec,
    score_type    = score_type,
    score_params  = score_params,
    sigmaX_fn     = sigmaX_fn,
    corr_type     = corr_type,
    corr_params   = corr_params,
    jitter        = jitter,
    seed          = NULL
  )
  
  Z_list <- process_info$Z_list
  
  # ---------------------------------------------------------------------------
  # STEP 4: GENERATE NOISE
  # ---------------------------------------------------------------------------
  if (!is.null(resolved_seeds$noise)) {
    set.seed(resolved_seeds$noise)
  }
  
  noise_list <- generate_noise(
    T_list       = T_list,
    sigma0       = sigma0,
    noise_mean   = noise_mean,
    noise_type   = noise_type,
    noise_params = noise_params
  )
  
  # ---------------------------------------------------------------------------
  # STEP 5: ASSEMBLE FINAL OBSERVATIONS
  # ---------------------------------------------------------------------------
  obs_info <- assemble_observations(
    mu_values  = mu_values,
    Z_list     = Z_list,
    noise_list = noise_list
  )
  
  Y_list <- obs_info$Y_list
  
  # ---------------------------------------------------------------------------
  # SUMMARY
  # ---------------------------------------------------------------------------
  summary_info <- list(
    n_subjects           = length(T_list),
    total_observations   = sum(lengths(T_list)),
    mean_obs_per_subject = mean(lengths(T_list)),
    road                 = road,
    noise_structure      = noise_structure,
    noise_type           = noise_type,
    has_noise            = !(noise_structure == "homoscedastic" &&
                               is.numeric(sigma0) &&
                               length(sigma0) == 1 &&
                               sigma0 == 0),
    seeds_used           = resolved_seeds
  )
  
  # ---------------------------------------------------------------------------
  # CONFIG
  # ---------------------------------------------------------------------------
  config <- list(
    n              = n,
    m_mean         = m_mean,
    delta          = delta,
    seed           = seed,
    seed_list      = seed_list,
    resolved_seeds = resolved_seeds,
    fO_type        = fO_type,
    fO_params      = fO_params,
    f0_type        = f0_type,
    f0_params      = f0_params,
    mu_fn          = mu_fn,
    road           = road,
    eigenfn_list   = eigenfn_list,
    eigenval_vec   = eigenval_vec,
    score_type     = score_type,
    score_params   = score_params,
    sigmaX_fn      = sigmaX_fn,
    corr_type      = corr_type,
    corr_params    = corr_params,
    jitter         = jitter,
    noise_structure = noise_structure,
    sigma0         = sigma0,
    noise_mean     = noise_mean,
    noise_type     = noise_type,
    noise_params   = noise_params
  )
  
  # ---------------------------------------------------------------------------
  # RETURN
  # ---------------------------------------------------------------------------
  return(list(
    Y_list       = Y_list,
    T_list       = T_list,
    mu_values    = mu_values,
    Z_list       = Z_list,
    noise_list   = noise_list,
    m_vec        = m_vec,
    O_vec        = O_vec,
    windows      = windows,
    process_info = process_info,
    design_info  = design_info,
    config       = config,
    summary      = summary_info
  ))
}
