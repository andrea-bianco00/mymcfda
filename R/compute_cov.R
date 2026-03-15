# =============================================================================
# compute_cov.R
# -----------------------------------------------------------------------------
# Unified dispatcher for the TRUE covariance function.
#
# The user MUST always specify which road to use:
#
#   road = "decomposition"   -> calls compute_cov_from_dec()
#   road = "kl"              -> calls compute_cov_from_KL()
#
# If route-specific arguments are omitted, the function uses sensible defaults:
#
#   DECOMPOSITION DEFAULTS  (Cov I, Lin & Wang 2022):
#     sigma_X^2(t) = sqrt(t) * sqrt(1 - t) / 10 + 1
#     rho(s,t)     = Matern(theta1 = 0.5, theta2 = 1)
#
#   KL DEFAULTS             (Cov II, Lin & Wang 2022):
#     C(s,t) = sum_{k=1}^{50} 2*k^{-2} * phi_k(s) * phi_k(t)
#     phi_k(t) = sqrt(2) * sin(2*k*pi*t)
#
# OUTPUT:
#   Numeric matrix of dimension length(s) x length(t), regardless of road.
#
# REQUIRED SOURCES:
#   source("compute_cov_from_dec.R")
#   source("compute_cov_from_KL.R")
#
# Dependencies: base R only (through the sourced files).
# =============================================================================



#' Compute true covariance using either decomposition or KL road
#'
#' @param s            Numeric vector. First set of time points.
#' @param t            Numeric vector. Second set of time points.
#' @param road         Character. Must be either "decomposition" or "kl".
#'
#' @param sigmaX_fn    Function. Used only if road = "decomposition".
#'                     Standard deviation function sigma_X(t), not variance.
#'                     If NULL, defaults to Cov I from Lin & Wang (2022).
#' @param corr_type    Character. Used only if road = "decomposition".
#'                     One of "matern", "power_exponential",
#'                     "rational_quadratic".
#'                     If NULL, defaults to "matern".
#' @param corr_params  Named list. Used only if road = "decomposition".
#'                     If NULL, defaults to list(theta1 = 0.5, theta2 = 1).
#'
#' @param eigenfn_list List of functions. Used only if road = "kl".
#'                     If NULL, defaults to Cov II from Lin & Wang (2022).
#' @param eigenval_vec Numeric vector. Used only if road = "kl".
#'                     If NULL, defaults to Cov II from Lin & Wang (2022).
#'
#' @return Numeric matrix of dimension length(s) x length(t).
#'
#' @examples
#' grid <- seq(0, 1, length.out = 100)
#'
#' # ---------------------------------------------------------------------------
#' # 1) Decomposition road with defaults (Cov I)
#' # ---------------------------------------------------------------------------
#' C1 <- compute_cov(s = grid, t = grid, road = "decomposition")
#'
#' # ---------------------------------------------------------------------------
#' # 2) Decomposition road with custom correlation
#' # ---------------------------------------------------------------------------
#' my_sigmaX <- function(t) 1 + 0.5 * sin(2 * pi * t)
#'
#' C2 <- compute_cov(
#'   s = grid,
#'   t = grid,
#'   road = "decomposition",
#'   sigmaX_fn = my_sigmaX,
#'   corr_type = "power_exponential",
#'   corr_params = list(theta1 = 1.5, theta2 = 0.3)
#' )
#'
#' # ---------------------------------------------------------------------------
#' # 3) KL road with defaults (Cov II)
#' # ---------------------------------------------------------------------------
#' C3 <- compute_cov(s = grid, t = grid, road = "kl")
#'
#' # ---------------------------------------------------------------------------
#' # 4) KL road with custom eigenstructure
#' # ---------------------------------------------------------------------------
#' K <- 10
#' my_eigenfn_list <- lapply(1:K, function(k) {
#'   function(t) sqrt(2) * sin(2 * pi * k * t)
#' })
#' my_eigenval_vec <- 1 / (1:K)^2
#'
#' C4 <- compute_cov(
#'   s = grid,
#'   t = grid,
#'   road = "kl",
#'   eigenfn_list = my_eigenfn_list,
#'   eigenval_vec = my_eigenval_vec
#' )

#' @export
compute_cov <- function(s,
                        t,
                        road,
                        sigmaX_fn    = NULL,
                        corr_type    = NULL,
                        corr_params  = NULL,
                        eigenfn_list = NULL,
                        eigenval_vec = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s) || length(s) < 1)
    stop("s must be a non-empty numeric vector.")
  
  if (!is.numeric(t) || length(t) < 1)
    stop("t must be a non-empty numeric vector.")
  
  if (missing(road))
    stop("You must specify road = 'decomposition' or road = 'kl'.")
  
  if (!is.character(road) || length(road) != 1 ||
      !road %in% c("decomposition", "kl")) {
    stop("road must be exactly one of: 'decomposition', 'kl'.")
  }
  
  # ---------------------------------------------------------------------------
  # ROAD 1: DECOMPOSITION
  # ---------------------------------------------------------------------------
  if (road == "decomposition") {
    
    # ----- default sigma_X(t): Cov I from Lin & Wang (2022) -------------------
    if (is.null(sigmaX_fn)) {
      sigmaX2_default <- function(t) {
        sqrt(t) * sqrt(1 - t) / 10 + 1
      }
      sigmaX_fn <- function(t) sqrt(sigmaX2_default(t))
    }
    
    # ----- default correlation family -----------------------------------------
    if (is.null(corr_type)) {
      corr_type <- "matern"
    }
    
    # ----- default correlation parameters: Cov I ------------------------------
    default_corr_params <- list(theta1 = 0.5, theta2 = 1)
    
    if (is.null(corr_params)) {
      corr_params <- default_corr_params
    } else {
      # merge user params over defaults
      corr_params <- modifyList(default_corr_params, corr_params)
    }
    
    # ----- compute covariance via decomposition -------------------------------
    C_mat <- compute_cov_from_dec(
      s = s,
      t = t,
      sigmaX_fn = sigmaX_fn,
      corr_type = corr_type,
      corr_params = corr_params
    )
    
    return(C_mat)
  }
  
  # ---------------------------------------------------------------------------
  # ROAD 2: KL
  # ---------------------------------------------------------------------------
  if (road == "kl") {
    
    # ----- defaults: Cov II from Lin & Wang (2022) ----------------------------
    if (is.null(eigenfn_list) && is.null(eigenval_vec)) {
      K_default <- 50
      
      eigenfn_list <- lapply(seq_len(K_default), function(k) {
        force(k)
        function(t) sqrt(2) * sin(2 * k * pi * t)
      })
      
      eigenval_vec <- 2 * (seq_len(K_default)^(-2))
      
    } else if (is.null(eigenfn_list) || is.null(eigenval_vec)) {
      stop("For road = 'kl', either provide BOTH eigenfn_list and ",
           "eigenval_vec, or provide neither and use the defaults.")
    }
    
    # ----- compute covariance via KL ------------------------------------------
    C_mat <- compute_cov_from_KL(
      s = s,
      t = t,
      eigenfn_list = eigenfn_list,
      eigenval_vec = eigenval_vec
    )
    
    return(C_mat)
  }
}
