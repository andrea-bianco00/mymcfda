#' Empirical criterion Q_n(theta) for covariance parameter estimation
#'
#' Computes the empirical objective function used to estimate the parameters
#' of a parametric correlation model \eqn{\rho_\theta(s,t)} in the covariance
#' structure of sparse functional data.
#'
#' The criterion corresponds to:
#'
#' \deqn{
#' Q_n(\theta) =
#' \sum_{i=1}^{n} \frac{1}{m_i(m_i-1)}
#' \sum_{j \neq l}
#' \left[
#' \sigma_X(T_{ij}) \rho_\theta(T_{ij},T_{il}) \sigma_X(T_{il})
#' - C_{ijl}
#' \right]^2
#' }
#'
#' where
#'
#' \deqn{
#' C_{ijl} = (Y_{ij} - \hat{\mu}(T_{ij})) (Y_{il} - \hat{\mu}(T_{il}))
#' }
#'
#' The function is designed to be used inside an optimizer (typically
#' \code{optim()}) to estimate the parameter vector \eqn{\theta}.
#'
#' Supported correlation models:
#'
#' * `"power_exponential"`
#' * `"rational_quadratic"`
#' * `"matern"`
#'
#' Numerical safeguards are implemented to ensure robustness during
#' optimization.
#'
#' @param raw_cov_df Data frame produced by \code{compute_raw_covariances()}.
#' Must contain columns:
#' \itemize{
#' \item id
#' \item j
#' \item l
#' \item t1
#' \item t2
#' \item raw_cov
#' }
#'
#' @param sigmaX_obs_list List where element \code{i} contains the estimated
#' values of \eqn{\sigma_X(T_{ij})} evaluated at the observed time points of
#' subject \code{i}.
#'
#' @param theta Numeric vector of length 2 containing the parameters of the
#' correlation model.
#'
#' @param model Character string specifying the correlation model.
#' Allowed values are:
#'
#' \itemize{
#' \item `"power_exponential"`
#' \item `"rational_quadratic"`
#' \item `"matern"`
#' }
#'
#' @param return_details Logical. If \code{FALSE} (default), the function
#' returns only the value of the criterion \eqn{Q_n(\theta)}.
#' If \code{TRUE}, diagnostic quantities are also returned.
#'
#' @return
#' If \code{return_details = FALSE}:
#'
#' A numeric scalar equal to \eqn{Q_n(\theta)}.
#'
#' If \code{return_details = TRUE}:
#'
#' A list containing
#'
#' \itemize{
#' \item Qn
#' \item n_pairs
#' \item rho
#' \item sigma1
#' \item sigma2
#' \item fitted_cov
#' \item residuals
#' \item weights
#' }
#'
#' @examples
#' # Example usage (assuming required objects exist)
#' # Qn_val <- compute_Qn_theta(
#' #   raw_cov_df      = raw_cov_df,
#' #   sigmaX_obs_list = sigmaX_obs_list,
#' #   theta           = c(1, 0.2),
#' #   model           = "matern"
#' # )
#'
#' @export
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
    stop("'theta' must be a numeric vector of length 2.")
  }

  if (!is.character(model) || length(model) != 1) {
    stop("'model' must be a character string.")
  }

  valid_models <- c("power_exponential", "rational_quadratic", "matern")

  if (!model %in% valid_models) {
    stop("Invalid model.")
  }

  theta1 <- theta[1]
  theta2 <- theta[2]

  if (!is.finite(theta1) || !is.finite(theta2)) {
    if (!return_details) return(1e12)
  }

  # ---------------------------------------------------------------------------
  # 2. subject sizes
  # ---------------------------------------------------------------------------
  n <- length(sigmaX_obs_list)
  m_vec <- sapply(sigmaX_obs_list, length)

  # ---------------------------------------------------------------------------
  # 3. sigmaX extraction
  # ---------------------------------------------------------------------------
  sigma1 <- mapply(
    function(i,j) sigmaX_obs_list[[i]][j],
    raw_cov_df$id,
    raw_cov_df$j
  )

  sigma2 <- mapply(
    function(i,l) sigmaX_obs_list[[i]][l],
    raw_cov_df$id,
    raw_cov_df$l
  )

  if (any(!is.finite(sigma1)) || any(!is.finite(sigma2))) {
    if (!return_details) return(1e12)
  }

  # ---------------------------------------------------------------------------
  # 4. compute pairwise distances
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

        # Numerical safeguard
        u_nz <- pmax(u[!zero_idx], 1e-12)

        const <- 1 / (gamma(theta1) * 2^(theta1 - 1))

        rho_nz <- const *
          (u_nz^theta1) *
          besselK(u_nz, nu = theta1)

        rho[!zero_idx] <- rho_nz
      }

      rho
    }
  )

  # ---------------------------------------------------------------------------
  # 5. robustness safeguard
  # ---------------------------------------------------------------------------

  if (any(!is.finite(rho_vec))) {
    if (!return_details) return(1e12)

    return(list(
      Qn         = 1e12,
      n_pairs    = nrow(raw_cov_df),
      rho        = rho_vec,
      sigma1     = sigma1,
      sigma2     = sigma2,
      fitted_cov = rep(NA_real_, length(rho_vec)),
      residuals  = rep(NA_real_, length(rho_vec)),
      weights    = rep(NA_real_, length(rho_vec))
    ))
  }

  # keep correlations inside [0,1]
  rho_vec <- pmin(pmax(rho_vec, 0), 1)

  # ---------------------------------------------------------------------------
  # 6. covariance model
  # ---------------------------------------------------------------------------

  fitted_cov <- sigma1 * rho_vec * sigma2
  residuals  <- fitted_cov - raw_cov_df$raw_cov

  weights <- 1 / (m_vec[raw_cov_df$id] * (m_vec[raw_cov_df$id] - 1))

  Qn <- sum(weights * residuals^2)

  # safeguard for optimizer
  if (!is.finite(Qn)) {

    if (!return_details) return(1e12)

    return(list(
      Qn         = 1e12,
      n_pairs    = nrow(raw_cov_df),
      rho        = rho_vec,
      sigma1     = sigma1,
      sigma2     = sigma2,
      fitted_cov = fitted_cov,
      residuals  = residuals,
      weights    = weights
    ))
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
