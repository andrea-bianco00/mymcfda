# =============================================================================
# bandwidth_cv_mu.R
# -----------------------------------------------------------------------------
# Selects the bandwidth h_mu for mean function estimation via kappa-fold
# cross-validation, as described in Section 2.1 of Lin & Wang (2022).
#
# Cross-validation criterion (equation 3, Lin & Wang 2022):
#
#   CV(h) = sum_{k=1}^{kappa} sum_{i in P_k} sum_{j=1}^{m_i}
#             {Y_ij - mu_hat_{h,-k}(T_ij)}^2
#
# where mu_hat_{h,-k} is the local linear estimator computed with bandwidth h
# and with all subjects in fold P_k excluded from estimation.
#
# The candidate bandwidth grid H is constructed automatically as a log-spaced
# grid of n_h values centered around the theoretically optimal rate
# h_center = (n * m_bar)^{-1/5} from condition (H4) of Lin & Wang (2022).
#
# OPTIMIZED IMPLEMENTATION:
#   The innermost loop over observations j is eliminated. For each (h, k) pair,
#   all test observations are flattened into a single vector and
#   compute_mu_hat_t() is called once per (h, k) instead of once per
#   observation. This reduces n_h * kappa * N calls to n_h * kappa calls.
#
# Dependencies:
#   - compute_mu_hat_t()  defined in mu_hat_t.R
#   - compute_Sr()        defined in S_r.R
#   - compute_Rr()        defined in R_r.R
#   - compute_weights()   defined in weights.R
#   - kernel_h_mu()       defined in kernel_h_mu.R
#   - kernel_fun()        defined in kernel.R
# =============================================================================


#' Select bandwidth h_mu via kappa-fold cross-validation
#'
#' @param T_list  List of length n.
#' @param Y_list  List of length n.
#' @param kappa   Integer >= 2. Number of folds. Default 5.
#' @param n_h     Integer >= 2. Number of bandwidth candidates. Default 40.
#' @param h_min   Numeric > 0. Lower bound of grid. If NULL, set automatically.
#' @param h_max   Numeric > 0. Upper bound of grid. If NULL, set automatically.
#' @param kernel  Character. Kernel type.
#' @param scheme  Character. Weighting scheme: "OBS", "SUBJ", or "OPTIMAL".
#' @param seed    Integer. Random seed for fold assignment (required).
#'
#' @return A list: $h_opt, $h_grid, $cv_vals.

bandwidth_cv_mu <- function(T_list, Y_list, kappa = 5, n_h = 40,
                            h_min = NULL, h_max = NULL,
                            kernel, scheme, seed) {
  
  # --- input checks -----------------------------------------------------------
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1)
    stop("'seed' must be specified as a single integer.\n  Example: seed = 123")
  if (!is.list(T_list) || !is.list(Y_list))
    stop("T_list and Y_list must be lists.")
  if (length(T_list) != length(Y_list))
    stop("T_list and Y_list must have the same length n.")
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
      Y_train   <- Y_list[train_idx]
      
      # test set: flatten all test observations into one vector
      test_idx  <- which(fold_ids == k)
      T_test    <- unlist(T_list[test_idx])
      Y_test    <- unlist(Y_list[test_idx])
      
      # evaluate mu_hat_{h,-k} at ALL test points in one call
      mu_pred <- compute_mu_hat_t(t_vec  = T_test,
                                  T_list = T_train,
                                  Y_list = Y_train,
                                  h_mu   = h,
                                  kernel = kernel,
                                  scheme = scheme)
      
      cv_h <- cv_h + sum((Y_test - mu_pred)^2)
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
