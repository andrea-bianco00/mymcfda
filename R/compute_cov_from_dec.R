# =============================================================================
# compute_cov_from_dec.R
# -----------------------------------------------------------------------------
# Computes the covariance function C(s,t) through the decomposition used in
# Lin & Wang (2022):
#
#   C(s,t) = sigma_X(s) * rho_theta(s,t) * sigma_X(t)
#
# where:
#   - sigma_X(.) is a user-supplied standard deviation function
#   - rho_theta(s,t) is a user-chosen correlation family
#
# This function is the decomposition analogue of compute_cov_from_KL.R:
#
#   - compute_cov_from_KL.R   : builds C(s,t) from eigenfunctions/eigenvalues
#   - compute_cov_from_dec.R  : builds C(s,t) from sigma_X and rho_theta
#
# IMPORTANT:
#   sigmaX_fn must be the STANDARD DEVIATION function sigma_X(t),
#   NOT the variance function sigma_X^2(t).
#
#   If the user starts from a variance function sigma_X^2(t), they must pass:
#
#     sigmaX_fn <- function(t) sqrt(sigmaX2_fn(t))
#
# CORRELATION FAMILIES SUPPORTED:
#   - "matern"
#   - "power_exponential"
#   - "rational_quadratic"
#
# REQUIRED SOURCES:
#   source("matern_correlation.R")
#   source("power_exponential_correlation.R")
#   source("rational_quadratic_correlation.R")
#
# OUTPUT:
#   Numeric matrix of dimension length(s) x length(t), exactly like
#   compute_cov_from_KL.R
#
# Dependencies: base R only (through the sourced correlation files).
# =============================================================================



#' Compute covariance function from variance-correlation decomposition
#'
#' @param s            Numeric vector. First set of time points.
#' @param t            Numeric vector. Second set of time points.
#' @param sigmaX_fn    Function. User-supplied standard deviation function
#'                     sigma_X(t). Must be vectorised and return a numeric
#'                     vector of the same length as input, with non-negative
#'                     values.
#' @param corr_type    Character. Correlation family to use. One of:
#'                       "matern"
#'                       "power_exponential"
#'                       "rational_quadratic"
#' @param corr_params  Named list of parameters for the chosen correlation:
#'                     - matern:             list(theta1, theta2)
#'                     - power_exponential:  list(theta1, theta2)
#'                     - rational_quadratic: list(theta1, theta2)
#'
#' @return Numeric matrix of dimension length(s) x length(t).
#'         Entry (i,j) is:
#'           C(s_i, t_j) = sigma_X(s_i) * rho(s_i, t_j) * sigma_X(t_j)
#'
#' @examples
#' # ---------------------------------------------------------------------------
#' # Example 1: Cov I from Lin & Wang (2022)
#' # sigma_X^2(t) = sqrt(t)*(1-t)^(1/2)/10 + 1
#' # rho = Matern(theta1 = 0.5, theta2 = 1)
#' # ---------------------------------------------------------------------------
#' sigmaX_fn <- function(t) {
#'   sqrt( sqrt(t) * sqrt(1 - t) / 10 + 1 )
#' }
#'
#' s <- seq(0, 1, length.out = 100)
#' t <- seq(0, 1, length.out = 100)
#'
#' Cmat <- compute_cov_from_dec(
#'   s = s,
#'   t = t,
#'   sigmaX_fn = sigmaX_fn,
#'   corr_type = "matern",
#'   corr_params = list(theta1 = 0.5, theta2 = 1)
#' )
#'
#' dim(Cmat)   # 100 x 100
#'
#' # ---------------------------------------------------------------------------
#' # Example 2: Power Exponential correlation
#' # ---------------------------------------------------------------------------
#' sigmaX_fn2 <- function(t) 1 + 0.5 * sin(2 * pi * t)
#'
#' Cmat2 <- compute_cov_from_dec(
#'   s = s,
#'   t = t,
#'   sigmaX_fn = sigmaX_fn2,
#'   corr_type = "power_exponential",
#'   corr_params = list(theta1 = 1.5, theta2 = 0.3)
#' )
#'
#' # ---------------------------------------------------------------------------
#' # Example 3: Rational Quadratic correlation
#' # ---------------------------------------------------------------------------
#' Cmat3 <- compute_cov_from_dec(
#'   s = s,
#'   t = t,
#'   sigmaX_fn = sigmaX_fn2,
#'   corr_type = "rational_quadratic",
#'   corr_params = list(theta1 = 1, theta2 = 0.25)
#' )

compute_cov_from_dec <- function(s,
                                 t,
                                 sigmaX_fn,
                                 corr_type,
                                 corr_params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s) || length(s) < 1)
    stop("s must be a non-empty numeric vector.")
  
  if (!is.numeric(t) || length(t) < 1)
    stop("t must be a non-empty numeric vector.")
  
  if (!is.function(sigmaX_fn))
    stop("sigmaX_fn must be a function.")
  
  valid_types <- c("matern", "power_exponential", "rational_quadratic")
  if (!is.character(corr_type) || length(corr_type) != 1 || !corr_type %in% valid_types)
    stop("corr_type must be one of: ",
         paste(valid_types, collapse = ", "), ".")
  
  if (!is.list(corr_params))
    stop("corr_params must be a named list.")
  
  # ---------------------------------------------------------------------------
  # CHECK sigmaX_fn VECTORISATION
  # ---------------------------------------------------------------------------
  test_input <- c(0, 0.5, 1)
  
  test_output <- tryCatch(
    sigmaX_fn(test_input),
    error = function(e) {
      stop("sigmaX_fn failed when applied to a test vector c(0, 0.5, 1). ",
           "Error: ", e$message)
    }
  )
  
  if (!is.numeric(test_output))
    stop("sigmaX_fn must return a numeric vector.")
  
  if (length(test_output) != length(test_input))
    stop("sigmaX_fn must be vectorised: sigmaX_fn(x) must return a vector ",
         "of the same length as x.")
  
  if (any(!is.finite(test_output)))
    stop("sigmaX_fn returned non-finite values on the test input.")
  
  if (any(test_output < 0))
    stop("sigmaX_fn must return non-negative values.")
  
  # ---------------------------------------------------------------------------
  # EVALUATE sigma_X(s) AND sigma_X(t)
  # ---------------------------------------------------------------------------
  sigma_s <- sigmaX_fn(s)
  sigma_t <- sigmaX_fn(t)
  
  if (!is.numeric(sigma_s) || length(sigma_s) != length(s))
    stop("sigmaX_fn(s) must return a numeric vector of length length(s).")
  
  if (!is.numeric(sigma_t) || length(sigma_t) != length(t))
    stop("sigmaX_fn(t) must return a numeric vector of length length(t).")
  
  if (any(!is.finite(sigma_s)) || any(!is.finite(sigma_t)))
    stop("sigmaX_fn returned non-finite values on s or t.")
  
  if (any(sigma_s < 0) || any(sigma_t < 0))
    stop("sigmaX_fn returned negative values on s or t.")
  
  # ---------------------------------------------------------------------------
  # BUILD CORRELATION MATRIX rho(s,t)
  # ---------------------------------------------------------------------------
  rho_mat <- switch(
    corr_type,
    
    "matern" = {
      if (is.null(corr_params$theta1) || is.null(corr_params$theta2))
        stop("For corr_type = 'matern', corr_params must contain ",
             "'theta1' and 'theta2'.")
      
      matern_correlation(
        s = s,
        t = t,
        theta1 = corr_params$theta1,
        theta2 = corr_params$theta2
      )
    },
    
    "power_exponential" = {
      if (is.null(corr_params$theta1) || is.null(corr_params$theta2))
        stop("For corr_type = 'power_exponential', corr_params must contain ",
             "'theta1' and 'theta2'.")
      
      power_exponential_correlation(
        s = s,
        t = t,
        theta1 = corr_params$theta1,
        theta2 = corr_params$theta2
      )
    },
    
    "rational_quadratic" = {
      if (is.null(corr_params$theta1) || is.null(corr_params$theta2))
        stop("For corr_type = 'rational_quadratic', corr_params must contain ",
             "'theta1' and 'theta2'.")
      
      rational_quadratic_correlation(
        s = s,
        t = t,
        theta1 = corr_params$theta1,
        theta2 = corr_params$theta2
      )
    }
  )
  
  # ---------------------------------------------------------------------------
  # CHECK CORRELATION MATRIX
  # ---------------------------------------------------------------------------
  if (!is.matrix(rho_mat))
    stop("The selected correlation function did not return a matrix.")
  
  if (nrow(rho_mat) != length(s) || ncol(rho_mat) != length(t))
    stop("The selected correlation function must return a matrix of ",
         "dimension length(s) x length(t).")
  
  if (any(!is.finite(rho_mat)))
    stop("The correlation matrix contains non-finite values.")
  
  # ---------------------------------------------------------------------------
  # COMPUTE C(s,t) = sigma_X(s) * rho(s,t) * sigma_X(t)
  # ---------------------------------------------------------------------------
  C_mat <- outer(sigma_s, sigma_t) * rho_mat
  
  return(C_mat)
}
