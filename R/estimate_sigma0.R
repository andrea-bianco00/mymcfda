# =============================================================================
# estimate_sigma0.R
# -----------------------------------------------------------------------------
# Runs the full SNPT pipeline (Block 1 + Block 2 + sigma0 block) on a single
# dataset and returns the estimated noise variance sigma^2_0_hat.
#
# Pipeline:
#   1. bandwidth_cv_mu()       -> h_mu
#   2. compute_mu_hat()        -> mu_hat on grid
#   3. compute_residuals()     -> Z_ij = {Y_ij - mu_hat(T_ij)}^2
#   4. bandwidth_cv_sigma()    -> h_sigma
#   5. compute_varsigma2_hat() -> varsigma2_hat on grid
#   6. compute_h0()            -> h0 (empirical rule)
#   7. compute_sigma0()        -> sigma^2_0_hat
#
# Dependencies:
#   - bandwidth_cv_mu()
#   - compute_mu_hat()
#   - compute_residuals()
#   - bandwidth_cv_sigma()
#   - compute_varsigma2_hat()
#   - compute_h0()
#   - compute_sigma0()
# =============================================================================


#' Run full SNPT pipeline and return sigma^2_0_hat
#'
#' @param T_list   List of length n. Observation times per subject.
#' @param Y_list   List of length n. Observed values per subject.
#' @param kappa    Integer. Number of CV folds. Default 5.
#' @param n_h      Integer. Number of bandwidth candidates. Default 20.
#' @param kernel   Character. Kernel type. Default "epanechnikov".
#' @param scheme   Character. Weighting scheme. Default "OBS".
#' @param seed_cv  Integer. Seed for CV fold assignment. Default 1.
#' @param n_grid   Integer. Number of grid points for mu_hat and varsigma2_hat.
#'                 Default 100.
#' @param verbose  Logical. If TRUE, print pipeline output. Default FALSE.
#'
#' @return A list with:
#'   $sigma0_hat   Numeric. Estimated noise variance.
#'   $h_mu         Numeric. Selected bandwidth for mu.
#'   $h_sigma      Numeric. Selected bandwidth for varsigma^2.
#'   $h0           Numeric. Selected bandwidth for sigma^2_0.
#'   $h0_adjusted  Logical. Whether h0 was adjusted by safety check.
#'   $res_mu       List. Output of compute_mu_hat().
#'   $Z_list       List. Squared residuals.
#'   $res_varsigma List. Output of compute_varsigma2_hat().
#'   $A0           Numeric. Estimated A0.
#'   $A2           Numeric. Estimated A2.
#'   $B            Numeric. Estimated B.
#'   $Delta        Numeric. Ridge term in sigma0 estimation.
#'   $error        Character or NULL. Error message if pipeline failed.

#' @export
estimate_sigma0 <- function(T_list, Y_list,
                            kappa   = 5,
                            n_h     = 20,
                            kernel  = "epanechnikov",
                            scheme  = "OBS",
                            seed_cv = 1,
                            n_grid  = 100,
                            verbose = FALSE) {
  
  maybe_silent <- function(expr) {
    if (verbose) {
      expr
    } else {
      suppressMessages(capture.output(expr))
    }
  }
  
  result <- tryCatch({
    
    # --- Step 1: bandwidth for mu ---------------------------------------------
    cv_mu <- NULL
    maybe_silent(
      cv_mu <- bandwidth_cv_mu(
        T_list = T_list,
        Y_list = Y_list,
        kappa  = kappa,
        n_h    = n_h,
        kernel = kernel,
        scheme = scheme,
        seed   = seed_cv
      )
    )
    h_mu <- cv_mu$h_opt
    
    # --- Step 2: mu estimation ------------------------------------------------
    res_mu <- NULL
    maybe_silent(
      res_mu <- compute_mu_hat(
        T_list = T_list,
        Y_list = Y_list,
        h_mu   = h_mu,
        kernel = kernel,
        scheme = scheme,
        n_t    = n_grid
      )
    )
    
    # --- Step 3: squared residuals --------------------------------------------
    Z_list <- compute_residuals(
      T_list = T_list,
      Y_list = Y_list,
      h_mu   = h_mu,
      kernel = kernel,
      scheme = scheme
    )
    
    # --- Step 4: bandwidth for varsigma^2 -------------------------------------
    cv_sigma <- NULL
    maybe_silent(
      cv_sigma <- bandwidth_cv_sigma(
        T_list = T_list,
        Z_list = Z_list,
        kappa  = kappa,
        n_h    = n_h,
        kernel = kernel,
        scheme = scheme,
        seed   = seed_cv
      )
    )
    h_sigma <- cv_sigma$h_opt
    
    # --- Step 5: varsigma^2 estimation ----------------------------------------
    res_varsigma <- NULL
    maybe_silent(
      res_varsigma <- compute_varsigma2_hat(
        T_list  = T_list,
        Z_list  = Z_list,
        h_sigma = h_sigma,
        kernel  = kernel,
        scheme  = scheme,
        n_t     = n_grid
      )
    )
    
    # --- Step 6: h0 via empirical rule ----------------------------------------
    h0_res <- NULL
    maybe_silent(
      h0_res <- compute_h0(
        T_list       = T_list,
        res_varsigma = res_varsigma
      )
    )
    
    # --- Step 7: sigma^2_0 estimation -----------------------------------------
    sigma0_res <- NULL
    maybe_silent(
      sigma0_res <- compute_sigma0(
        T_list = T_list,
        Y_list = Y_list,
        h0     = h0_res$h0
      )
    )
    
    # --- return results -------------------------------------------------------
    list(
      sigma0_hat   = sigma0_res$sigma0_hat,
      h_mu         = h_mu,
      h_sigma      = h_sigma,
      h0           = h0_res$h0,
      h0_adjusted  = h0_res$adjusted,
      res_mu       = res_mu,
      Z_list       = Z_list,
      res_varsigma = res_varsigma,
      A0           = sigma0_res$A0,
      A2           = sigma0_res$A2,
      B            = sigma0_res$B,
      Delta        = sigma0_res$Delta,
      error        = NULL
    )
    
  }, error = function(e) {
    list(
      sigma0_hat   = NA_real_,
      h_mu         = NA_real_,
      h_sigma      = NA_real_,
      h0           = NA_real_,
      h0_adjusted  = NA,
      res_mu       = NULL,
      Z_list       = NULL,
      res_varsigma = NULL,
      A0           = NA_real_,
      A2           = NA_real_,
      B            = NA_real_,
      Delta        = NA_real_,
      error        = conditionMessage(e)
    )
  })
  
  return(result)
}
