# =============================================================================
# evaluate_sigmaX_on_observed_times.R
# -----------------------------------------------------------------------------
# Evaluates the estimated process standard deviation sigma_X(t) at the observed
# times T_ij for each subject.
#
# Input:
#   - T_list: list of observed times for each subject
#   - res_sigmaX2: output of compute_sigmaX2_hat(), containing:
#         * t_grid
#         * sigmaX2_hat
#
# The function linearly interpolates sigmaX2_hat(t) from t_grid to the observed
# times T_ij, then returns:
#
#   sigmaX_obs_list[[i]] = ( sigmaX_hat(T_i1), ..., sigmaX_hat(T_i m_i) )
#
# where:
#
#   sigmaX_hat(t) = sqrt( sigmaX2_hat(t) )
#
# Dependencies:
#   - none beyond base R
# =============================================================================


#' Evaluate sigma_X(t) at observed times
#'
#' @param T_list       List of length n. T_list((i)) contains the observation
#'                     times for subject i.
#' @param res_sigmaX2  List. Output of compute_sigmaX2_hat(), containing:
#'                     - t_grid
#'                     - sigmaX2_hat
#' @param rule         Integer. Passed to approx() for extrapolation behavior.
#'                     Default is 2, meaning constant extrapolation at boundaries.
#'
#' @return A list of length n. For each subject i:
#'         - sigmaX_obs_list((i)) is a numeric vector of length m_i
#'           containing sigmaX_hat(T_ij) = sqrt(sigmaX2_hat(T_ij))
#'
#' @examples
#' # res_sigmaX2 <- compute_sigmaX2_hat(res_varsigma, sigma0_hat)
#' # sigmaX_obs_list <- evaluate_sigmaX_on_observed_times(T_list, res_sigmaX2)

#' @export
evaluate_sigmaX_on_observed_times <- function(T_list, res_sigmaX2, rule = 2) {
  
  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!is.list(T_list)) {
    stop("T_list must be a list.")
  }
  
  if (!is.list(res_sigmaX2)) {
    stop("res_sigmaX2 must be a list.")
  }
  
  required_names <- c("t_grid", "sigmaX2_hat")
  if (!all(required_names %in% names(res_sigmaX2))) {
    stop("res_sigmaX2 must contain 't_grid' and 'sigmaX2_hat'.")
  }
  
  t_grid <- res_sigmaX2$t_grid
  sigmaX2_hat <- res_sigmaX2$sigmaX2_hat
  
  if (!is.numeric(t_grid) || !is.numeric(sigmaX2_hat)) {
    stop("'t_grid' and 'sigmaX2_hat' must be numeric vectors.")
  }
  
  if (length(t_grid) != length(sigmaX2_hat)) {
    stop("'t_grid' and 'sigmaX2_hat' must have the same length.")
  }
  
  if (length(t_grid) < 2) {
    stop("'t_grid' must contain at least two points for interpolation.")
  }
  
  if (any(is.na(t_grid)) || any(is.na(sigmaX2_hat))) {
    stop("'t_grid' and 'sigmaX2_hat' must not contain NA values.")
  }
  
  if (any(diff(t_grid) <= 0)) {
    stop("'t_grid' must be strictly increasing.")
  }
  
  if (any(sigmaX2_hat < 0)) {
    stop("'sigmaX2_hat' must be nonnegative.")
  }
  
  if (!is.numeric(rule) || length(rule) != 1 || !(rule %in% c(1, 2))) {
    stop("'rule' must be either 1 or 2, as in approx().")
  }
  
  n <- length(T_list)
  
  for (i in seq_len(n)) {
    if (!is.numeric(T_list[[i]])) {
      stop(sprintf("T_list[[%d]] must be a numeric vector.", i))
    }
    if (any(is.na(T_list[[i]]))) {
      stop(sprintf("T_list[[%d]] contains NA values.", i))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Interpolate sigmaX2_hat on observed times and take square root
  # ---------------------------------------------------------------------------
  sigmaX_obs_list <- vector("list", length = n)
  
  for (i in seq_len(n)) {
    Ti <- T_list[[i]]
    
    if (length(Ti) == 0) {
      sigmaX_obs_list[[i]] <- numeric(0)
      next
    }
    
    sigmaX2_interp <- approx(
      x    = t_grid,
      y    = sigmaX2_hat,
      xout = Ti,
      rule = rule
    )$y
    
    # numerical safety
    sigmaX2_interp <- pmax(sigmaX2_interp, 0)
    
    sigmaX_obs_list[[i]] <- sqrt(sigmaX2_interp)
  }
  
  return(sigmaX_obs_list)
}
