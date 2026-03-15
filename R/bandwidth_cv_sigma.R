# =============================================================================
# bandwidth_cv_sigma.R
# -----------------------------------------------------------------------------
# Selects the bandwidth h_sigma for varsigma^2(t) estimation via kappa-fold
# cross-validation, as described in Section 2.2 of Lin & Wang (2022).
#
# CV criterion (analogue of equation 3, applied to squared residuals Z_ij):
#
#   CV(h) = sum_{k=1}^{kappa} sum_{i in P_k} sum_{j=1}^{m_i}
#             {Z_ij - varsigma2_hat_{h,-k}(T_ij)}^2
#
# where varsigma2_hat_{h,-k} is estimated with bandwidth h excluding fold k.
#
# OPTIMIZED IMPLEMENTATION:
#   The innermost loop over observations j is eliminated. For each (h, k) pair,
#   all test observations are flattened into a single vector and
#   compute_varsigma2_hat_t() is called once per (h, k) instead of once per
#   observation. This reduces n_h * kappa * N calls to n_h * kappa calls.
#
# Dependencies:
#   - compute_varsigma2_hat_t()  defined in varsigma2_hat_t.R
#   - compute_Sr_sigma()         defined in S_r_sigma.R
#   - compute_Rr_sigma()         defined in R_r_sigma.R
#   - compute_weights()          defined in weights.R
#   - kernel_h_sigma()           defined in kernel_h_sigma.R
#   - kernel_fun()               defined in kernel.R
# =============================================================================


#' Select bandwidth h_sigma via kappa-fold cross-validation
#'
#' @param T_list  List of length n.
#' @param Z_list  List of length n. Z_list[[i]] = {Y_ij - mu_hat(T_ij)}^2.
#' @param kappa   Integer >= 2. Number of folds. Default 5.
#' @param n_h     Integer >= 2. Number of bandwidth candidates. Default 40.
#' @param h_min   Numeric > 0. Lower bound of grid. If NULL, set automatically.
#' @param h_max   Numeric > 0. Upper bound of grid. If NULL, set automatically.
#' @param kernel  Character. Kernel type.
#' @param scheme  Character. Weighting scheme: "OBS", "SUBJ", or "OPTIMAL".
#' @param seed    Integer. Random seed for fold assignment (required).
#'
#' @return A list: $h_opt, $h_grid, $cv_vals.

bandwidth_cv_sigma <- function(T_list, Z_list, kappa = 5, n_h = 40,
                               h_min = NULL, h_max = NULL,
                               kernel, scheme, seed) {
  
  # --- input checks -----------------------------------------------------------
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1)
    stop("'seed' must be specified as a single integer.\n  Example: seed = 123")
  if (!is.list(T_list) || !is.list(Z_list))
    stop("T_list and Z_list must be lists.")
  if (length(T_list) != length(Z_list))
    stop("T_list and Z_list must have the same length n.")
  if (!is.numeric(kappa) || kappa < 2)
    stop("kappa must be an integer >= 2.")
  if (!is.numeric(n_h) || n_h < 2)
    stop("n_h must be an integer >= 2.")
  if (!is.null(h_min) && (!is.numeric(h_min) || h_min <= 0))
    stop("h_min must be a strictly positive number.")
  if (!is.null(h_max) && (!is.numeric(h_max) || h_max <= 0))
    stop("h_max must be a strictly positive number.")
  if (!is.null(h_min) && !is.null(h_max) && h_min >= h_max)
    stop("h_min must be strictly less than h_max.")
  
  # --- basic quantities -------------------------------------------------------
  m_vec    <- sapply(T_list, length)
  n        <- length(m_vec)
  m_bar    <- mean(m_vec)
  h_center <- (n * m_bar)^(-1/5)
  
  if (is.null(h_min)) h_min <- h_center / 10
  if (is.null(h_max)) h_max <- 3 * h_center
  
  h_grid <- exp(seq(log(h_min), log(h_max), length.out = n_h))
  
  cat(sprintf("Bandwidth grid: %d values in [%.4f, %.4f] (h_center = %.4f)\n",
              n_h, h_min, h_max, h_center))
  
  # --- construct kappa folds --------------------------------------------------
  set.seed(seed)
  fold_ids <- sample(rep(seq_len(kappa), length.out = n))
  cat(sprintf("Fold sizes: %s\n", paste(tabulate(fold_ids), collapse = " ")))
  
  # --- cross-validation loop --------------------------------------------------
  cv_vals <- numeric(n_h)
  
  for (h_idx in seq_len(n_h)) {
    h    <- h_grid[h_idx]
    cv_h <- 0
    
    for (k in seq_len(kappa)) {
      
      # training set
      train_idx <- which(fold_ids != k)
      T_train   <- T_list[train_idx]
      Z_train   <- Z_list[train_idx]
      
      # test set: flatten all test observations into one vector
      test_idx <- which(fold_ids == k)
      T_test   <- unlist(T_list[test_idx])
      Z_test   <- unlist(Z_list[test_idx])
      
      # evaluate varsigma2_hat_{h,-k} at ALL test points in one call
      varsig_pred <- compute_varsigma2_hat_t(t_vec   = T_test,
                                             T_list  = T_train,
                                             Z_list  = Z_train,
                                             h_sigma = h,
                                             kernel  = kernel,
                                             scheme  = scheme)
      
      cv_h <- cv_h + sum((Z_test - varsig_pred)^2)
    }
    
    cv_vals[h_idx] <- cv_h
  }
  
  # --- select optimal bandwidth -----------------------------------------------
  h_opt   <- h_grid[which.min(cv_vals)]
  opt_idx <- which.min(cv_vals)
  
  if (opt_idx == 1)
    warning(sprintf(
      "h_opt = %.4f is at the LEFT boundary (h_min = %.4f). Consider smaller h_min.",
      h_opt, h_min))
  if (opt_idx == n_h)
    warning(sprintf(
      "h_opt = %.4f is at the RIGHT boundary (h_max = %.4f). Consider larger h_max.",
      h_opt, h_max))
  
  cat(sprintf("Optimal bandwidth: h_opt = %.4f  (CV = %.6f)\n",
              h_opt, min(cv_vals)))
  
  return(list(h_opt = h_opt, h_grid = h_grid, cv_vals = cv_vals))
}
