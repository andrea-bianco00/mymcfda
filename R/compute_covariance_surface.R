# =============================================================================
# compute_covariance_surface.R
# -----------------------------------------------------------------------------
# Compute the final covariance surface:
#
#   C_hat(s, t) = sigmaX_hat(s) * rho_theta(s, t) * sigmaX_hat(t)
#
# using:
#   - res_sigmaX2 from compute_sigmaX2_hat()
#   - theta from estimate_theta()
#   - rho_theta() for the chosen model
#
# sigmaX_hat is obtained from sigmaX2_hat via square root, after interpolation
# on the requested grids.
#
# Dependencies:
#   - rho_theta()  defined in rho_theta.R
# =============================================================================


#' Compute the estimated covariance surface C_hat(s, t)
#'
#' @param res_sigmaX2  List returned by compute_sigmaX2_hat(), containing:
#'                     - t_grid
#'                     - sigmaX2_hat
#' @param theta        Numeric vector of length 2. Estimated parameter vector.
#' @param model        Character string. Correlation family passed to rho_theta().
#' @param s_grid       Optional numeric vector of s-values. If NULL, use
#'                     res_sigmaX2$t_grid.
#' @param t_grid       Optional numeric vector of t-values. If NULL, use
#'                     res_sigmaX2$t_grid.
#' @param rule         Integer passed to approx() for extrapolation behavior.
#'                     Default is 2 (constant extrapolation at boundaries).
#'
#' @return A list with:
#'   - s_grid
#'   - t_grid
#'   - sigma_s
#'   - sigma_t
#'   - rho_mat
#'   - covariance_hat
#'
#' @examples
#' # cov_res <- compute_covariance_surface(
#' #   res_sigmaX2 = res_sigmaX2,
#' #   theta       = fit_pe$theta_hat,
#' #   model       = "power_exponential"
#' # )

compute_covariance_surface <- function(res_sigmaX2,
                                       theta,
                                       model,
                                       s_grid = NULL,
                                       t_grid = NULL,
                                       rule = 2) {
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.list(res_sigmaX2)) {
    stop("'res_sigmaX2' must be a list.")
  }
  
  required_names <- c("t_grid", "sigmaX2_hat")
  if (!all(required_names %in% names(res_sigmaX2))) {
    stop("'res_sigmaX2' must contain 't_grid' and 'sigmaX2_hat'.")
  }
  
  base_grid <- res_sigmaX2$t_grid
  sigmaX2_hat <- res_sigmaX2$sigmaX2_hat
  
  if (!is.numeric(base_grid) || !is.numeric(sigmaX2_hat)) {
    stop("'res_sigmaX2$t_grid' and 'res_sigmaX2$sigmaX2_hat' must be numeric.")
  }
  
  if (length(base_grid) != length(sigmaX2_hat)) {
    stop("'t_grid' and 'sigmaX2_hat' must have the same length.")
  }
  
  if (length(base_grid) < 2) {
    stop("'t_grid' must contain at least two points.")
  }
  
  if (any(is.na(base_grid)) || any(is.na(sigmaX2_hat))) {
    stop("'t_grid' and 'sigmaX2_hat' must not contain NA values.")
  }
  
  if (any(diff(base_grid) <= 0)) {
    stop("'t_grid' must be strictly increasing.")
  }
  
  if (any(sigmaX2_hat < 0)) {
    stop("'sigmaX2_hat' must be nonnegative.")
  }
  
  if (!is.numeric(theta) || length(theta) != 2 || any(is.na(theta))) {
    stop("'theta' must be a numeric vector of length 2.")
  }
  
  if (!is.character(model) || length(model) != 1) {
    stop("'model' must be a single character string.")
  }
  
  if (!is.numeric(rule) || length(rule) != 1 || !(rule %in% c(1, 2))) {
    stop("'rule' must be either 1 or 2, as in approx().")
  }
  
  # ---------------------------------------------------------------------------
  # 2. default grids
  # ---------------------------------------------------------------------------
  if (is.null(s_grid)) s_grid <- base_grid
  if (is.null(t_grid)) t_grid <- base_grid
  
  if (!is.numeric(s_grid) || length(s_grid) < 1 || any(is.na(s_grid))) {
    stop("'s_grid' must be a non-empty numeric vector without NA.")
  }
  
  if (!is.numeric(t_grid) || length(t_grid) < 1 || any(is.na(t_grid))) {
    stop("'t_grid' must be a non-empty numeric vector without NA.")
  }
  
  # ---------------------------------------------------------------------------
  # 3. interpolate sigmaX2_hat on s_grid and t_grid
  # ---------------------------------------------------------------------------
  sigmaX2_s <- approx(
    x    = base_grid,
    y    = sigmaX2_hat,
    xout = s_grid,
    rule = rule
  )$y
  
  sigmaX2_t <- approx(
    x    = base_grid,
    y    = sigmaX2_hat,
    xout = t_grid,
    rule = rule
  )$y
  
  # numerical safeguard
  sigmaX2_s <- pmax(sigmaX2_s, 0)
  sigmaX2_t <- pmax(sigmaX2_t, 0)
  
  sigma_s <- sqrt(sigmaX2_s)
  sigma_t <- sqrt(sigmaX2_t)
  
  # ---------------------------------------------------------------------------
  # 4. correlation matrix rho_theta(s, t)
  # ---------------------------------------------------------------------------
  rho_mat <- rho_theta(
    s     = s_grid,
    t     = t_grid,
    theta = theta,
    model = model
  )
  
  if (!is.matrix(rho_mat)) {
    stop("rho_theta() must return a matrix.")
  }
  
  if (!all(dim(rho_mat) == c(length(s_grid), length(t_grid)))) {
    stop("rho_theta() returned a matrix with incorrect dimensions.")
  }
  
  # ---------------------------------------------------------------------------
  # 5. covariance surface
  # ---------------------------------------------------------------------------
  covariance_hat <- outer(sigma_s, sigma_t) * rho_mat
  
  if (any(!is.finite(covariance_hat))) {
    stop("Non-finite values produced in covariance_hat.")
  }
  
  # ---------------------------------------------------------------------------
  # 6. output
  # ---------------------------------------------------------------------------
  return(list(
    s_grid         = s_grid,
    t_grid         = t_grid,
    sigma_s        = sigma_s,
    sigma_t        = sigma_t,
    rho_mat        = rho_mat,
    covariance_hat = covariance_hat
  ))
}
