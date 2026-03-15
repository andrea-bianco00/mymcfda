# =============================================================================
# generate_KL_process.R
# -----------------------------------------------------------------------------
# Generates the centred (zero-mean) random process Z_i(t) = X_i(t) - mu(t)
# via truncated Karhunen-Loève expansion:
#
#   Z_i(t) = sum_{k=1}^{K} xi_{ik} * phi_k(t)
#
# where:
#   - phi_k(t) are the eigenfunctions (deterministic, user-supplied)
#   - lambda_k are the eigenvalues  (deterministic, user-supplied)
#   - xi_{ik} are the scores: random, with E[xi_{ik}] = 0, Var(xi_{ik}) = lambda_k
#
# The scores are generated as xi_{ik} = sqrt(lambda_k) * eta_{ik}, where
# eta_{ik} ~ score_fn with mean 0 and variance 1. The score_fn is built
# internally by calling build_score_fn() with the user's choice of
# score_type and score_params.
#
# By default, score_type = "gaussian", giving a Gaussian process as in
# Lin & Wang (2022), Section 5.
#
# ROLE IN THE PIPELINE:
#   generate_snippet_data  -->  T_list
#   eval_mean_at_times     -->  mu_values
#   build_score_fn         -->  score_fn (called internally)
#   generate_KL_process    -->  Z_list + scores  (THIS FILE)
#   generate_noise         -->  noise_list
#   assemble_observations  -->  Y_list = mu_values + Z_list + noise_list
#
# This function does NOT add the mean mu(t) — it returns only the centred
# part Z_i(T_ij). The mean is added later by assemble_observations.
#
# REQUIRED SOURCES:
#   source("build_score_fn.R")
#
# INPUT:
#   T_list        : list of n numeric vectors. T_list[[i]] contains the m_i
#                   measurement times for subject i, all in [0, 1].
#   eigenfn_list  : list of K functions. eigenfn_list[[k]] is a function
#                   phi_k : numeric vector -> numeric vector of same length.
#                   Each must be vectorised.
#   eigenval_vec  : numeric vector of length K with eigenvalues lambda_1, ...,
#                   lambda_K, all strictly positive.
#   score_type    : character. Distribution family for the scores.
#                   One of: "gaussian" (default), "t", "uniform", "laplace",
#                   "custom". Passed to build_score_fn().
#   score_params  : named list of parameters for the chosen distribution.
#                   Passed to build_score_fn(). See build_score_fn.R for
#                   documentation of each family's parameters.
#                   Default: list() (uses family defaults).
#
# OUTPUT:
#   A list with two components:
#   - Z_list : list of n numeric vectors. Z_list[[i]] is a vector of length
#              m_i containing Z_i(T_i1), ..., Z_i(T_im_i).
#   - scores : numeric matrix of dimension n x K. scores[i, k] = xi_{ik},
#              the score of subject i for component k. Stored for diagnostics
#              and debugging.
#
# EXAMPLES:
#   # --- Cov II from Lin & Wang (2022) ---
#   K <- 50
#   eigenfn_list <- lapply(1:K, function(k) {
#     function(t) sqrt(2) * sin(2 * k * pi * t)
#   })
#   eigenval_vec <- 2 * (1:K)^(-2)
#
#   T_list <- list(c(0.1, 0.2, 0.3), c(0.5, 0.6))
#   result <- generate_KL_process(T_list, eigenfn_list, eigenval_vec)
#   str(result$Z_list)   # list of 2 vectors
#   dim(result$scores)   # 2 x 50
#
#   # --- Non-Gaussian scores (t-distribution with df=5) ---
#   result2 <- generate_KL_process(T_list, eigenfn_list, eigenval_vec,
#                                  score_type = "t",
#                                  score_params = list(df = 5))
#
# Dependencies: base R only (via build_score_fn.R).
# =============================================================================



#' Generate centred process via Karhunen-Loève expansion
#'
#' @param T_list        List of n numeric vectors of measurement times.
#' @param eigenfn_list  List of K vectorised functions phi_k(t).
#' @param eigenval_vec  Numeric vector of K positive eigenvalues.
#' @param score_type    Character. Distribution family for scores.
#'                      Default: "gaussian".
#' @param score_params  Named list of parameters for the distribution.
#'                      Default: list().
#'
#' @return A list with components:
#'   - Z_list : list of n numeric vectors with Z_i(T_ij) values
#'   - scores : n x K matrix of scores xi_{ik}

generate_KL_process <- function(T_list,
                                eigenfn_list,
                                eigenval_vec,
                                score_type   = "gaussian",
                                score_params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || length(T_list) < 1)
    stop("T_list must be a non-empty list of numeric vectors.")
  
  if (!is.list(eigenfn_list) || length(eigenfn_list) < 1)
    stop("eigenfn_list must be a non-empty list of functions.")
  
  K <- length(eigenfn_list)
  
  if (!is.numeric(eigenval_vec) || length(eigenval_vec) != K)
    stop("eigenval_vec must be a numeric vector of the same length as ",
         "eigenfn_list (K = ", K, ").")
  
  if (any(eigenval_vec <= 0))
    stop("All eigenvalues must be strictly positive.")
  
  # ---------------------------------------------------------------------------
  # BUILD THE STANDARDISED SCORE FUNCTION via build_score_fn
  # ---------------------------------------------------------------------------
  score_fn <- build_score_fn(type = score_type, params = score_params)
  
  # ---------------------------------------------------------------------------
  # VALIDATE EIGENFUNCTIONS (vectorisation check)
  # ---------------------------------------------------------------------------
  test_input <- c(0.0, 0.5, 1.0)
  
  for (k in seq_len(K)) {
    if (!is.function(eigenfn_list[[k]]))
      stop("eigenfn_list[[", k, "]] is not a function.")
    
    test_out <- tryCatch(eigenfn_list[[k]](test_input),
                         error = function(e) {
                           stop("eigenfn_list[[", k, "]] failed on test input. ",
                                "Error: ", e$message)
                         })
    
    if (!is.numeric(test_out) || length(test_out) != length(test_input))
      stop("eigenfn_list[[", k, "]] must be vectorised: it must return a ",
           "numeric vector of the same length as its input. ",
           "Input length: ", length(test_input),
           ", output length: ", length(test_out), ".")
  }
  
  # ---------------------------------------------------------------------------
  # GENERATE SCORES
  # ---------------------------------------------------------------------------
  # For each subject i and each component k, generate:
  #   xi_{ik} = sqrt(lambda_k) * eta_{ik},   eta_{ik} ~ score_fn (mean 0, var 1)
  #
  # We generate all scores at once as an n x K matrix for efficiency.
  # ---------------------------------------------------------------------------
  n <- length(T_list)
  sqrt_eigenvals <- sqrt(eigenval_vec)
  
  scores <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    scores[, k] <- score_fn(n) * sqrt_eigenvals[k]
  }
  
  # ---------------------------------------------------------------------------
  # EVALUATE THE KL EXPANSION AT MEASUREMENT TIMES
  # ---------------------------------------------------------------------------
  # For each subject i:
  #   Z_i(T_ij) = sum_{k=1}^{K} xi_{ik} * phi_k(T_ij)
  #
  # Strategy: for each subject, build an m_i x K matrix of eigenfunction
  # values, then multiply by the score vector (length K) to get m_i values.
  # ---------------------------------------------------------------------------
  Z_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    t_i  <- T_list[[i]]
    m_i  <- length(t_i)
    
    if (!is.numeric(t_i) || m_i < 1)
      stop("T_list[[", i, "]] must be a non-empty numeric vector.")
    
    # build m_i x K matrix: phi_mat[j, k] = phi_k(T_ij)
    phi_mat <- matrix(NA_real_, nrow = m_i, ncol = K)
    for (k in seq_len(K)) {
      phi_mat[, k] <- eigenfn_list[[k]](t_i)
    }
    
    # Z_i = phi_mat %*% xi_i  (m_i x K) %*% (K x 1) = (m_i x 1)
    Z_list[[i]] <- as.numeric(phi_mat %*% scores[i, ])
  }
  
  # ---------------------------------------------------------------------------
  # RETURN
  # ---------------------------------------------------------------------------
  return(list(
    Z_list = Z_list,
    scores = scores
  ))
}
