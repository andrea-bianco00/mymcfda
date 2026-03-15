# =============================================================================
# evaluate_estimation_performance.R
# -----------------------------------------------------------------------------
# Wrapper to evaluate Monte Carlo performance of the estimation procedure.
#
# It can compute:
#   - RMSE  for sigma0^2
#   - RMISE for sigmaX^2(t)
#   - RMISE for C(s,t)
#
# depending on which inputs are provided.
#
# Dependencies:
#   - compute_rmse()
#   - compute_rmise_variance()
#   - compute_rmise_covariance()
# =============================================================================


#' Evaluate estimation performance over Monte Carlo replicates
#'
#' @param sigma0_estimates           Optional numeric vector. Estimated sigma0^2
#'                                   values over Monte Carlo replicates.
#' @param sigma2_0_true             Optional numeric scalar. True sigma0^2.
#'
#' @param variance_estimates_list   Optional list. Estimated variance functions
#'                                   over Monte Carlo replicates.
#' @param t_grid                    Optional numeric vector. Common grid for the
#'                                   variance function.
#' @param sigmaX2_true             Optional numeric vector. True sigmaX^2(t) on t_grid.
#' @param sigmaX2_true_fun         Optional function. True sigmaX^2(t).
#'
#' @param covariance_estimates_list Optional list. Estimated covariance surfaces
#'                                   over Monte Carlo replicates.
#' @param s_grid                    Optional numeric vector. Grid for first argument
#'                                   of covariance.
#' @param t_grid_cov                Optional numeric vector. Grid for second argument
#'                                   of covariance.
#' @param C_true_mat               Optional matrix. True covariance on s_grid x t_grid_cov.
#' @param C_true_fun               Optional function. True covariance function.
#'
#' @param na_rm_sigma0             Logical. Passed to compute_rmse().
#'
#' @return A list with components:
#'   - rmse_sigma0
#'   - rmise_variance
#'   - rmise_covariance
#'
#' Components are NULL when the required inputs were not provided.
#'
#' @examples
#' # perf <- evaluate_estimation_performance(
#' #   sigma0_estimates         = c(0.11, 0.09, 0.12),
#' #   sigma2_0_true           = 0.10,
#' #   variance_estimates_list = list(v1, v2, v3),
#' #   t_grid                  = t_grid,
#' #   sigmaX2_true_fun        = sigmaX2_true_fun,
#' #   covariance_estimates_list = list(C1, C2, C3),
#' #   s_grid                  = s_grid,
#' #   t_grid_cov              = t_grid_cov,
#' #   C_true_fun              = C_true_fun
#' # )

evaluate_estimation_performance <- function(
    sigma0_estimates = NULL,
    sigma2_0_true = NULL,
    variance_estimates_list = NULL,
    t_grid = NULL,
    sigmaX2_true = NULL,
    sigmaX2_true_fun = NULL,
    covariance_estimates_list = NULL,
    s_grid = NULL,
    t_grid_cov = NULL,
    C_true_mat = NULL,
    C_true_fun = NULL,
    na_rm_sigma0 = FALSE
) {
  
  # ---------------------------------------------------------------------------
  # initialize outputs
  # ---------------------------------------------------------------------------
  rmse_sigma0      <- NULL
  rmise_variance   <- NULL
  rmise_covariance <- NULL
  
  # ---------------------------------------------------------------------------
  # 1. RMSE for sigma0^2
  # ---------------------------------------------------------------------------
  if (!is.null(sigma0_estimates) || !is.null(sigma2_0_true)) {
    
    if (is.null(sigma0_estimates) || is.null(sigma2_0_true)) {
      stop("To compute RMSE for sigma0^2, provide both 'sigma0_estimates' and 'sigma2_0_true'.")
    }
    
    rmse_sigma0 <- compute_rmse(
      estimates      = sigma0_estimates,
      sigma2_0_true  = sigma2_0_true,
      na_rm          = na_rm_sigma0
    )
  }
  
  # ---------------------------------------------------------------------------
  # 2. RMISE for sigmaX^2(t)
  # ---------------------------------------------------------------------------
  if (!is.null(variance_estimates_list) ||
      !is.null(t_grid) ||
      !is.null(sigmaX2_true) ||
      !is.null(sigmaX2_true_fun)) {
    
    if (is.null(variance_estimates_list) || is.null(t_grid)) {
      stop("To compute RMISE for variance, provide at least 'variance_estimates_list' and 't_grid'.")
    }
    
    rmise_variance <- compute_rmise_variance(
      estimates_list   = variance_estimates_list,
      t_grid           = t_grid,
      sigmaX2_true     = sigmaX2_true,
      sigmaX2_true_fun = sigmaX2_true_fun
    )
  }
  
  # ---------------------------------------------------------------------------
  # 3. RMISE for covariance C(s,t)
  # ---------------------------------------------------------------------------
  if (!is.null(covariance_estimates_list) ||
      !is.null(s_grid) ||
      !is.null(t_grid_cov) ||
      !is.null(C_true_mat) ||
      !is.null(C_true_fun)) {
    
    if (is.null(covariance_estimates_list) || is.null(s_grid) || is.null(t_grid_cov)) {
      stop(
        paste(
          "To compute RMISE for covariance, provide at least",
          "'covariance_estimates_list', 's_grid', and 't_grid_cov'."
        )
      )
    }
    
    rmise_covariance <- compute_rmise_covariance(
      estimates_list = covariance_estimates_list,
      s_grid         = s_grid,
      t_grid         = t_grid_cov,
      C_true_mat     = C_true_mat,
      C_true_fun     = C_true_fun
    )
  }
  
  # ---------------------------------------------------------------------------
  # output
  # ---------------------------------------------------------------------------
  list(
    rmse_sigma0      = rmse_sigma0,
    rmise_variance   = rmise_variance,
    rmise_covariance = rmise_covariance
  )
}
