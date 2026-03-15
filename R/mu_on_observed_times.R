# =============================================================================
# mu_on_observed_times.R
# -----------------------------------------------------------------------------
# Evaluates the estimated mean function mu_hat(t) at the observed times T_ij
# for each subject.
#
# For each subject i, this function computes:
#
#   mu_obs_list[[i]] = ( mu_hat(T_i1), ..., mu_hat(T_i m_i) )
#
# using the already implemented function compute_mu_hat_t().
#
# This is a utility / bridge function that avoids recomputing mu_hat(T_ij)
# repeatedly in later blocks (for example, raw covariance construction).
#
# Dependencies:
#   - compute_mu_hat_t() defined in mu_hat_t.R
# =============================================================================


#' Evaluate mu_hat at the observed times of each subject
#'
#' @param T_list  List of length n. T_list((i)) contains the observation times
#'                for subject i.
#' @param Y_list  List of length n. Y_list((i)) contains the observed values
#'                for subject i.
#' @param h_mu    Numeric scalar. Bandwidth used for mean estimation.
#' @param kernel  Character string. Kernel type passed to compute_mu_hat_t().
#' @param scheme  Character string. Weighting scheme passed to compute_mu_hat_t().
#'
#' @return A list of length n. For each subject i:
#'   - mu_obs_list((i)) is a numeric vector of length m_i
#'     containing mu_hat(T_ij) for j = 1, ..., m_i.
#'
#' @examples
#' # mu_obs_list <- compute_mu_on_observed_times(
#' #   T_list = T_list,
#' #   Y_list = Y_list,
#' #   h_mu   = h_mu_opt,
#' #   kernel = "epanechnikov",
#' #   scheme = "OBS"
#' # )

#' @export
compute_mu_on_observed_times <- function(T_list, Y_list, h_mu, kernel, scheme) {
  
  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list)) {
    stop("T_list and Y_list must be lists.")
  }
  
  if (length(T_list) != length(Y_list)) {
    stop("T_list and Y_list must have the same length.")
  }
  
  if (!is.numeric(h_mu) || length(h_mu) != 1 || is.na(h_mu) || h_mu <= 0) {
    stop("h_mu must be a single strictly positive numeric value.")
  }
  
  n <- length(T_list)
  
  for (i in seq_len(n)) {
    Ti <- T_list[[i]]
    Yi <- Y_list[[i]]
    
    if (!is.numeric(Ti) || !is.numeric(Yi)) {
      stop(sprintf("T_list[[%d]] and Y_list[[%d]] must both be numeric vectors.", i, i))
    }
    
    if (length(Ti) != length(Yi)) {
      stop(sprintf("T_list[[%d]] and Y_list[[%d]] must have the same length.", i, i))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Evaluate mu_hat(T_ij) subject by subject
  # ---------------------------------------------------------------------------
  mu_obs_list <- vector("list", length = n)
  
  for (i in seq_len(n)) {
    Ti <- T_list[[i]]
    
    mu_obs_list[[i]] <- compute_mu_hat_t(
      t_vec  = Ti,
      T_list = T_list,
      Y_list = Y_list,
      h_mu   = h_mu,
      kernel = kernel,
      scheme = scheme
    )
  }
  
  return(mu_obs_list)
}
