# =============================================================================
# raw_covariances.R
# -----------------------------------------------------------------------------
# Computes raw off-diagonal covariance pseudo-observations for the covariance
# estimation block, using precomputed values of mu_hat(T_ij).
#
# For each subject i and each ordered pair (j, l) with j != l, this function
# computes:
#
#   centered_ij = Y_ij - mu_hat(T_ij)
#   centered_il = Y_il - mu_hat(T_il)
#   raw_cov_ijl = centered_ij * centered_il
#
# IMPORTANT:
#   - Only within-subject pairs are used.
#   - Only off-diagonal pairs are used (j != l).
#   - Ordered pairs are kept, so each subject contributes m_i (m_i - 1) rows.
#
# Dependencies:
#   - none beyond the supplied inputs
# =============================================================================


#' Compute raw off-diagonal covariance pseudo-observations
#'
#' @param T_list       List of length n. T_list[[i]] contains the observation
#'                     times for subject i.
#' @param Y_list       List of length n. Y_list[[i]] contains the observed
#'                     values for subject i.
#' @param mu_obs_list  List of length n. mu_obs_list[[i]] contains the values
#'                     mu_hat(T_ij) for j = 1, ..., m_i.
#'
#' @return A data.frame with one row for each ordered off-diagonal pair (j,l)
#'         within the same subject. Columns:
#'         - id
#'         - j
#'         - l
#'         - t1
#'         - t2
#'         - y1
#'         - y2
#'         - mu1
#'         - mu2
#'         - centered1
#'         - centered2
#'         - raw_cov
#'
#' @examples
#' # mu_obs_list <- compute_mu_on_observed_times(T_list, Y_list, h_mu, kernel, scheme)
#' # raw_cov_df  <- compute_raw_covariances(T_list, Y_list, mu_obs_list)

compute_raw_covariances <- function(T_list, Y_list, mu_obs_list) {
  
  # ---------------------------------------------------------------------------
  # Input checks
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list) || !is.list(mu_obs_list)) {
    stop("T_list, Y_list, and mu_obs_list must all be lists.")
  }
  
  if (length(T_list) != length(Y_list) || length(T_list) != length(mu_obs_list)) {
    stop("T_list, Y_list, and mu_obs_list must all have the same length.")
  }
  
  n <- length(T_list)
  
  for (i in seq_len(n)) {
    Ti  <- T_list[[i]]
    Yi  <- Y_list[[i]]
    mui <- mu_obs_list[[i]]
    
    if (!is.numeric(Ti) || !is.numeric(Yi) || !is.numeric(mui)) {
      stop(sprintf(
        "T_list[[%d]], Y_list[[%d]], and mu_obs_list[[%d]] must all be numeric vectors.",
        i, i, i
      ))
    }
    
    if (!(length(Ti) == length(Yi) && length(Ti) == length(mui))) {
      stop(sprintf(
        "T_list[[%d]], Y_list[[%d]], and mu_obs_list[[%d]] must have the same length.",
        i, i, i
      ))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Build output subject by subject
  # ---------------------------------------------------------------------------
  out_list <- vector("list", length = n)
  
  for (i in seq_len(n)) {
    Ti  <- T_list[[i]]
    Yi  <- Y_list[[i]]
    mui <- mu_obs_list[[i]]
    
    mi <- length(Ti)
    
    # No off-diagonal pairs if fewer than 2 observations
    if (mi < 2) {
      out_list[[i]] <- NULL
      next
    }
    
    centered_i <- Yi - mui
    
    n_pairs <- mi * (mi - 1)
    
    # Preallocate
    id_vec        <- integer(n_pairs)
    j_vec         <- integer(n_pairs)
    l_vec         <- integer(n_pairs)
    t1_vec        <- numeric(n_pairs)
    t2_vec        <- numeric(n_pairs)
    y1_vec        <- numeric(n_pairs)
    y2_vec        <- numeric(n_pairs)
    mu1_vec       <- numeric(n_pairs)
    mu2_vec       <- numeric(n_pairs)
    cent1_vec     <- numeric(n_pairs)
    cent2_vec     <- numeric(n_pairs)
    raw_cov_vec   <- numeric(n_pairs)
    
    idx <- 1L
    
    for (j in seq_len(mi)) {
      for (l in seq_len(mi)) {
        
        if (j == l) next
        
        id_vec[idx]      <- i
        j_vec[idx]       <- j
        l_vec[idx]       <- l
        
        t1_vec[idx]      <- Ti[j]
        t2_vec[idx]      <- Ti[l]
        
        y1_vec[idx]      <- Yi[j]
        y2_vec[idx]      <- Yi[l]
        
        mu1_vec[idx]     <- mui[j]
        mu2_vec[idx]     <- mui[l]
        
        cent1_vec[idx]   <- centered_i[j]
        cent2_vec[idx]   <- centered_i[l]
        
        raw_cov_vec[idx] <- centered_i[j] * centered_i[l]
        
        idx <- idx + 1L
      }
    }
    
    out_list[[i]] <- data.frame(
      id        = id_vec,
      j         = j_vec,
      l         = l_vec,
      t1        = t1_vec,
      t2        = t2_vec,
      y1        = y1_vec,
      y2        = y2_vec,
      mu1       = mu1_vec,
      mu2       = mu2_vec,
      centered1 = cent1_vec,
      centered2 = cent2_vec,
      raw_cov   = raw_cov_vec
    )
  }
  
  # ---------------------------------------------------------------------------
  # Combine all subjects
  # ---------------------------------------------------------------------------
  out_df <- do.call(rbind, out_list)
  
  # If all subjects had mi < 2, return an empty data.frame with correct columns
  if (is.null(out_df)) {
    out_df <- data.frame(
      id        = integer(0),
      j         = integer(0),
      l         = integer(0),
      t1        = numeric(0),
      t2        = numeric(0),
      y1        = numeric(0),
      y2        = numeric(0),
      mu1       = numeric(0),
      mu2       = numeric(0),
      centered1 = numeric(0),
      centered2 = numeric(0),
      raw_cov   = numeric(0)
    )
  }
  
  rownames(out_df) <- NULL
  
  return(out_df)
}
