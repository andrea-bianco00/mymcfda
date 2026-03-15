# =============================================================================
# generate_decomposition_process.R
# -----------------------------------------------------------------------------
# Generates the centred process values Z_i(T_i) using the DECOMPOSITION road:
#
#   C(s,t) = sigma_X(s) * rho_theta(s,t) * sigma_X(t)
#
# For each subject i with measurement times T_i = (T_i1, ..., T_im_i),
# this function:
#
#   1. builds the subject-specific covariance matrix
#        Sigma_i = [ C(T_ij, T_il) ]_{j,l=1}^{m_i}
#
#   2. generates
#        Z_i(T_i) ~ N(0, Sigma_i)
#
# This function is the decomposition counterpart of generate_KL_process(),
# which instead generates Z_i(T_i) directly from the KL representation.
#
# REQUIRED SOURCES:
#   source("compute_cov_from_dec.R")
#
# REQUIRED PACKAGE:
#   MASS::mvrnorm
#
# Dependencies: base R + MASS
# =============================================================================



#' Generate centred process via variance-correlation decomposition
#'
#' @param T_list       List of n numeric vectors. T_list[[i]] contains the m_i
#'                     measurement times for subject i.
#' @param sigmaX_fn    Function. Standard deviation function sigma_X(t).
#'                     Must be vectorised and non-negative.
#' @param corr_type    Character. Correlation family. One of:
#'                       "matern"
#'                       "power_exponential"
#'                       "rational_quadratic"
#' @param corr_params  Named list of parameters for the chosen correlation.
#' @param seed         Integer or NULL. Optional seed for reproducibility.
#' @param return_Sigma Logical. If TRUE, also returns Sigma_list.
#'                     Default: TRUE.
#' @param jitter       Small non-negative number added to the diagonal of
#'                     Sigma_i for numerical stability. Default: 1e-10.
#'
#' @return A list with:
#'   - Z_list     : list of n numeric vectors. Z_list[[i]] contains
#'                  Z_i(T_i1), ..., Z_i(T_im_i).
#'   - Sigma_list : optional list of subject-specific covariance matrices.
#'   - config     : list of arguments used for generation.
#'
#' @examples
#' T_list <- list(c(0.1, 0.2, 0.3), c(0.55, 0.60, 0.72, 0.80))
#'
#' sigmaX_fn <- function(t) {
#'   sqrt( sqrt(t) * sqrt(1 - t) / 10 + 1 )
#' }
#'
#' out <- generate_decomposition_process(
#'   T_list = T_list,
#'   sigmaX_fn = sigmaX_fn,
#'   corr_type = "matern",
#'   corr_params = list(theta1 = 0.5, theta2 = 1),
#'   seed = 123
#' )
#'
#' str(out$Z_list)
#' str(out$Sigma_list)

generate_decomposition_process <- function(T_list,
                                           sigmaX_fn,
                                           corr_type,
                                           corr_params = list(),
                                           seed = NULL,
                                           return_Sigma = TRUE,
                                           jitter = 1e-10) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || length(T_list) < 1)
    stop("T_list must be a non-empty list of numeric vectors.")
  
  if (!is.function(sigmaX_fn))
    stop("sigmaX_fn must be a function.")
  
  valid_types <- c("matern", "power_exponential", "rational_quadratic")
  if (!is.character(corr_type) || length(corr_type) != 1 || !corr_type %in% valid_types)
    stop("corr_type must be one of: ", paste(valid_types, collapse = ", "), ".")
  
  if (!is.list(corr_params))
    stop("corr_params must be a named list.")
  
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1))
    stop("seed must be a single number or NULL.")
  
  if (!is.logical(return_Sigma) || length(return_Sigma) != 1)
    stop("return_Sigma must be TRUE or FALSE.")
  
  if (!is.numeric(jitter) || length(jitter) != 1 || jitter < 0)
    stop("jitter must be a single non-negative number.")
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install it with install.packages('MASS').")
  }
  
  # Validate T_list
  for (i in seq_along(T_list)) {
    if (!is.numeric(T_list[[i]]) || length(T_list[[i]]) < 1)
      stop("T_list[[", i, "]] must be a non-empty numeric vector.")
    
    if (any(!is.finite(T_list[[i]])))
      stop("T_list[[", i, "]] contains non-finite values.")
    
    if (any(T_list[[i]] < 0 | T_list[[i]] > 1))
      stop("T_list[[", i, "]] contains values outside [0,1].")
  }
  
  # ---------------------------------------------------------------------------
  # SET SEED
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  # ---------------------------------------------------------------------------
  # GENERATE SUBJECT-SPECIFIC PROCESS VALUES
  # ---------------------------------------------------------------------------
  n <- length(T_list)
  Z_list <- vector("list", n)
  
  if (return_Sigma) {
    Sigma_list <- vector("list", n)
  }
  
  for (i in seq_len(n)) {
    
    T_i <- T_list[[i]]
    m_i <- length(T_i)
    
    # -------------------------------------------------------------------------
    # Build Sigma_i = C(T_i, T_i) using decomposition only
    # -------------------------------------------------------------------------
    Sigma_i <- compute_cov_from_dec(
      s = T_i,
      t = T_i,
      sigmaX_fn = sigmaX_fn,
      corr_type = corr_type,
      corr_params = corr_params
    )
    
    if (!is.matrix(Sigma_i) || nrow(Sigma_i) != m_i || ncol(Sigma_i) != m_i) {
      stop("Internal error: Sigma_i for subject ", i,
           " is not a square matrix of size m_i x m_i.")
    }
    
    # Numerical symmetrisation
    Sigma_i <- (Sigma_i + t(Sigma_i)) / 2
    
    # Small diagonal inflation for stability
    if (jitter > 0) {
      diag(Sigma_i) <- diag(Sigma_i) + jitter
    }
    
    if (return_Sigma) {
      Sigma_list[[i]] <- Sigma_i
    }
    
    # -------------------------------------------------------------------------
    # Generate Z_i(T_i) ~ N(0, Sigma_i)
    # -------------------------------------------------------------------------
    Z_i <- tryCatch(
      MASS::mvrnorm(
        n = 1,
        mu = rep(0, m_i),
        Sigma = Sigma_i
      ),
      error = function(e) {
        stop("Gaussian simulation failed for subject ", i, ". ",
             "Possible reason: Sigma_i is not numerically positive semidefinite. ",
             "Original error: ", e$message)
      }
    )
    
    Z_list[[i]] <- as.numeric(Z_i)
  }
  
  # ---------------------------------------------------------------------------
  # RETURN
  # ---------------------------------------------------------------------------
  config <- list(
    seed         = seed,
    return_Sigma = return_Sigma,
    jitter       = jitter,
    corr_type    = corr_type,
    corr_params  = corr_params,
    sigmaX_fn    = sigmaX_fn
  )
  
  if (return_Sigma) {
    return(list(
      Z_list     = Z_list,
      Sigma_list = Sigma_list,
      config     = config
    ))
  } else {
    return(list(
      Z_list = Z_list,
      config = config
    ))
  }
}
