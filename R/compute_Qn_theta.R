# =============================================================================
# compute_Qn_theta.R
# -----------------------------------------------------------------------------
# Optimized computation of the empirical criterion Q_n(theta) used to estimate
# the parametric correlation rho_theta(s, t) in the covariance block.
#
# Paper criterion:
#
#   Q_n(theta) =
#     sum_{i=1}^n [ 1 / {m_i (m_i - 1)} ] *
#       sum_{j != l}
#         [ sigmaX_hat(T_ij) * rho_theta(T_ij, T_il) * sigmaX_hat(T_il)
#           - C_ijl ]^2
#
# where:
#   C_ijl = {Y_ij - mu_hat(T_ij)} {Y_il - mu_hat(T_il)}
#
# IMPORTANT:
#   This implementation is optimized:
#   - no rowwise mapply over rho_theta()
#   - no N x N matrix construction
#   - rho_theta is evaluated directly on the vector of pairwise distances
#
# Supported models:
#   - "power_exponential"
#   - "rational_quadratic"
#   - "matern"
#
# Dependencies:
#   - raw_cov_df from compute_raw_covariances()
#   - sigmaX_obs_list from evaluate_sigmaX_on_observed_times()
# =============================================================================


#' Compute the empirical criterion Q_n(theta) (optimized version)
#'
#' @param raw_cov_df       Data frame returned by compute_raw_covariances().
#'                         Must contain columns:
#'                         - id
#'                         - j
#'                         - l
#'                         - t1
#'                         - t2
#'                         - raw_cov
#' @param sigmaX_obs_list  List of length n. sigmaX_obs_list[[i]][j] is
#'                         sigmaX_hat(T_ij).
#' @param theta            Numeric vector of length 2: c(theta1, theta2).
#' @param model            Character string. One of:
#'                         - "power_exponential"
#'                         - "rational_quadratic"
#'                         - "matern"
#' @param return_details   Logical. If TRUE, return diagnostics in addition
#'                         to Q_n(theta).
#'
#' @return
#' If return_details = FALSE:
#'   Numeric scalar Q_n(theta).
#'
#' If return_details = TRUE:
#'   A list with components:
#'   - Qn
#'   - n_pairs
#'   - rho
#'   - sigma1
#'   - sigma2
#'   - fitted_cov
#'   - residuals
#'   - weights
#'
#' @examples
#' # Qn_val <- compute_Qn_theta(
#' #   raw_cov_df      = raw_cov_df,
#' #   sigmaX_obs_list = sigmaX_obs_list,
#' #   theta           = c(1, 0.5),
#' #   model           = "power_exponential"
#' # )

compute_Qn_theta <- function(raw_cov_df,
                             sigmaX_obs_list,
                             theta,
                             model,
                             return_details = FALSE) {
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.data.frame(raw_cov_df)) {
    stop("'raw_cov_df' must be a data.frame.")
  }
  
  required_cols <- c("id", "j", "l", "t1", "t2", "raw_cov")
  missing_cols <- setdiff(required_cols, names(raw_cov_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "'raw_cov_df' is missing required column(s): %s.",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  if (!is.list(sigmaX_obs_list)) {
    stop("'sigmaX_obs_list' must be a list.")
  }
  
  if (!is.numeric(theta) || length(theta) != 2 || any(is.na(theta))) {
    stop("'theta' must be a numeric vector of length 2: c(theta1, theta2).")
  }
  
  if (!is.character(model) || length(model) != 1) {
    stop("'model' must be a single character string.")
  }
  
  valid_models <- c("power_exponential", "rational_quadratic", "matern")
  if (!model %in% valid_models) {
    stop(sprintf(
      "'model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }
  
  theta1 <- theta[1]
  theta2 <- theta[2]
  
  if (!is.finite(theta1) || !is.finite(theta2)) {
    stop("'theta1' and 'theta2' must be finite.")
  }
  
  if (model == "power_exponential") {
    if (theta1 <= 0 || theta1 > 2) {
      stop("For 'power_exponential', theta1 must be in (0, 2].")
    }
    if (theta2 <= 0) {
      stop("For 'power_exponential', theta2 must be > 0.")
    }
  } else {
    if (theta1 <= 0 || theta2 <= 0) {
      stop(sprintf("For '%s', theta1 and theta2 must both be > 0.", model))
    }
  }
  
  # Empty case
  if (nrow(raw_cov_df) == 0) {
    warning("raw_cov_df has 0 rows. Returning NA.")
    if (!return_details) return(NA_real_)
    return(list(
      Qn         = NA_real_,
      n_pairs    = 0L,
      rho        = numeric(0),
      sigma1     = numeric(0),
      sigma2     = numeric(0),
      fitted_cov = numeric(0),
      residuals  = numeric(0),
      weights    = numeric(0)
    ))
  }
  
  # ---------------------------------------------------------------------------
  # 2. validate indices and extract m_i
  # ---------------------------------------------------------------------------
  n <- length(sigmaX_obs_list)
  m_vec <- sapply(sigmaX_obs_list, length)
  
  if (any(raw_cov_df$id < 1) || any(raw_cov_df$id > n)) {
    stop("'raw_cov_df$id' contains invalid subject indices.")
  }
  
  for (r in seq_len(nrow(raw_cov_df))) {
    i <- raw_cov_df$id[r]
    j <- raw_cov_df$j[r]
    l <- raw_cov_df$l[r]
    
    if (!is.numeric(sigmaX_obs_list[[i]]) || any(is.na(sigmaX_obs_list[[i]]))) {
      stop(sprintf("sigmaX_obs_list[[%d]] must be a numeric vector without NA.", i))
    }
    
    if (j < 1 || j > m_vec[i]) {
      stop(sprintf("Row %d: j=%d is out of bounds for subject %d.", r, j, i))
    }
    if (l < 1 || l > m_vec[i]) {
      stop(sprintf("Row %d: l=%d is out of bounds for subject %d.", r, l, i))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 3. extract sigmaX_hat(T_ij) and sigmaX_hat(T_il)
  # ---------------------------------------------------------------------------
  sigma1 <- mapply(
    function(i, j) sigmaX_obs_list[[i]][j],
    i = raw_cov_df$id,
    j = raw_cov_df$j
  )
  
  sigma2 <- mapply(
    function(i, l) sigmaX_obs_list[[i]][l],
    i = raw_cov_df$id,
    l = raw_cov_df$l
  )
  
  if (any(!is.finite(sigma1)) || any(!is.finite(sigma2))) {
    stop("Non-finite values found in sigma1 or sigma2.")
  }
  
  if (any(sigma1 < 0) || any(sigma2 < 0)) {
    stop("sigmaX values must be nonnegative.")
  }
  
  # ---------------------------------------------------------------------------
  # 4. optimized computation of rho for pairwise rows
  # ---------------------------------------------------------------------------
  d <- abs(raw_cov_df$t1 - raw_cov_df$t2)
  
  rho_vec <- switch(
    model,
    
    "power_exponential" = {
      exp(-(d^theta1) / (theta2^theta1))
    },
    
    "rational_quadratic" = {
      (1 + d^2 / theta2^2)^(-theta1)
    },
    
    "matern" = {
      u <- sqrt(2 * theta1) * d / theta2
      
      rho <- numeric(length(u))
      zero_idx <- (u == 0)
      rho[zero_idx] <- 1
      
      if (any(!zero_idx)) {
        u_nz <- u[!zero_idx]
        rho[!zero_idx] <- (1 / (gamma(theta1) * 2^(theta1 - 1))) *
          (u_nz^theta1) *
          besselK(u_nz, nu = theta1)
      }
      
      rho
    }
  )
  
  if (any(!is.finite(rho_vec))) {
    stop("rho_theta produced non-finite values.")
  }
  
  # numerical safeguard
  rho_vec <- pmin(pmax(rho_vec, 0), 1)
  
  # ---------------------------------------------------------------------------
  # 5. modeled covariance and residuals
  # ---------------------------------------------------------------------------
  fitted_cov <- sigma1 * rho_vec * sigma2
  residuals  <- fitted_cov - raw_cov_df$raw_cov
  
  # ---------------------------------------------------------------------------
  # 6. subject weights from the paper:
  #    weight(row from subject i) = 1 / {m_i (m_i - 1)}
  # ---------------------------------------------------------------------------
  weights <- 1 / (m_vec[raw_cov_df$id] * (m_vec[raw_cov_df$id] - 1))
  
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Invalid subject weights encountered in Q_n(theta).")
  }
  
  Qn <- sum(weights * residuals^2)
  
  if (!is.finite(Qn)) {
    stop("Q_n(theta) is not finite.")
  }
  
  # ---------------------------------------------------------------------------
  # 7. output
  # ---------------------------------------------------------------------------
  if (!return_details) {
    return(Qn)
  }
  
  return(list(
    Qn         = Qn,
    n_pairs    = nrow(raw_cov_df),
    rho        = rho_vec,
    sigma1     = sigma1,
    sigma2     = sigma2,
    fitted_cov = fitted_cov,
    residuals  = residuals,
    weights    = weights
  ))
}
