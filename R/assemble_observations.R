# =============================================================================
# assemble_observations.R
# -----------------------------------------------------------------------------
# Assembles the final observed data:
#
#   Y_ij = mu(T_ij) + Z_i(T_ij) + epsilon_ij
#
# from the three building blocks:
#   - mu_values  : deterministic mean evaluated at the observed times
#   - Z_list     : centred process values at the observed times
#   - noise_list : measurement noise values
#
# INPUT:
#   mu_values   : list of n numeric vectors. mu_values[[i]] contains
#                 mu(T_i1), ..., mu(T_im_i)
#   Z_list      : list of n numeric vectors. Z_list[[i]] contains
#                 Z_i(T_i1), ..., Z_i(T_im_i)
#   noise_list  : list of n numeric vectors. noise_list[[i]] contains
#                 epsilon_i1, ..., epsilon_im_i
#
# OUTPUT:
#   A list with:
#   - Y_list  : list of n numeric vectors. Y_list[[i]] contains the final
#               observed values Y_i1, ..., Y_im_i
#   - config  : basic metadata about the assembly
#
# Dependencies: base R only.
# =============================================================================


#' Assemble final observed data from mean, centred process, and noise
#'
#' @param mu_values   List of n numeric vectors. Mean values at observed times.
#' @param Z_list      List of n numeric vectors. Centred process values.
#' @param noise_list  List of n numeric vectors. Measurement noise values.
#'
#' @return A list with:
#'   - Y_list  : list of n numeric vectors with final observed data
#'   - config  : list with assembly metadata
#'
#' @examples
#' mu_values  <- list(c(1, 2, 3), c(0.5, 0.7))
#' Z_list     <- list(c(0.1, -0.2, 0.3), c(-0.1, 0.2))
#' noise_list <- list(c(0.01, -0.03, 0.02), c(0.05, -0.02))
#'
#' out <- assemble_observations(mu_values, Z_list, noise_list)
#' out$Y_list

assemble_observations <- function(mu_values, Z_list, noise_list) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(mu_values) || length(mu_values) < 1) {
    stop("mu_values must be a non-empty list of numeric vectors.")
  }
  
  if (!is.list(Z_list) || length(Z_list) < 1) {
    stop("Z_list must be a non-empty list of numeric vectors.")
  }
  
  if (!is.list(noise_list) || length(noise_list) < 1) {
    stop("noise_list must be a non-empty list of numeric vectors.")
  }
  
  n <- length(mu_values)
  
  if (length(Z_list) != n || length(noise_list) != n) {
    stop("mu_values, Z_list, and noise_list must have the same length.")
  }
  
  # ---------------------------------------------------------------------------
  # ASSEMBLE Y_list
  # ---------------------------------------------------------------------------
  Y_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    mu_i <- mu_values[[i]]
    Z_i  <- Z_list[[i]]
    e_i  <- noise_list[[i]]
    
    if (!is.numeric(mu_i) || length(mu_i) < 1) {
      stop("mu_values[[", i, "]] must be a non-empty numeric vector.")
    }
    
    if (!is.numeric(Z_i) || length(Z_i) < 1) {
      stop("Z_list[[", i, "]] must be a non-empty numeric vector.")
    }
    
    if (!is.numeric(e_i) || length(e_i) < 1) {
      stop("noise_list[[", i, "]] must be a non-empty numeric vector.")
    }
    
    if (length(mu_i) != length(Z_i) || length(mu_i) != length(e_i)) {
      stop("For subject ", i, ", mu_values[[i]], Z_list[[i]], and ",
           "noise_list[[i]] must have the same length.")
    }
    
    if (any(!is.finite(mu_i)) || any(!is.finite(Z_i)) || any(!is.finite(e_i))) {
      stop("Non-finite values detected for subject ", i, ".")
    }
    
    Y_list[[i]] <- mu_i + Z_i + e_i
  }
  
  # ---------------------------------------------------------------------------
  # RETURN
  # ---------------------------------------------------------------------------
  config <- list(
    n_subjects = n,
    assembled_from = c("mu_values", "Z_list", "noise_list")
  )
  
  return(list(
    Y_list = Y_list,
    config = config
  ))
}
