# =============================================================================
# kernel.R
# -----------------------------------------------------------------------------
# Pure kernel function K(u) for use in local linear smoothing.
#
# This file defines K(u) only — the mathematical kernel evaluated at a
# rescaled argument u. The scaled version K_h(u) = h^{-1} K(u/h) is defined
# separately and builds on this function.
#
# Condition (B1) of Lin & Wang (2022) requires K to be:
#   - symmetric:        K(u) = K(-u)
#   - density:          integral K(u) du = 1, K(u) >= 0
#   - compactly supported on [-1, 1]
#   - Lipschitz continuous
#
# Kernels marked with (*) satisfy condition (B1).
# Kernels marked with (!) do NOT satisfy (B1) — included for completeness.
# =============================================================================


#' Pure kernel function K(u)
#'
#' @param u       Numeric vector. Argument of the kernel.
#' @param kernel  Character string. One of:
#'                  "epanechnikov"  (*) optimal in MSE sense
#'                  "biweight"      (*) smoother than Epanechnikov
#'                  "triweight"     (*) even smoother
#'                  "tricube"       (*) used in LOESS
#'                  "triangular"    (*) tent-shaped
#'                  "cosine"        (*) cosine-based
#'                  "uniform"       (!) not Lipschitz at boundaries
#'                  "gaussian"      (!) infinite support
#'
#' @return Numeric vector of the same length as u with values K(u).
#'
#' @examples
#' u <- seq(-1.5, 1.5, by = 0.01)
#' plot(u, kernel_fun(u, "epanechnikov"), type = "l")
#' lines(u, kernel_fun(u, "biweight"),   col = "red")
#' lines(u, kernel_fun(u, "gaussian"),   col = "blue")

kernel_fun <- function(u, kernel = "epanechnikov") {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(u)) {
    stop("u must be a numeric vector.")
  }
  valid_kernels <- c("epanechnikov", "biweight", "triweight", "tricube",
                     "triangular", "cosine", "uniform", "gaussian")
  if (!kernel %in% valid_kernels) {
    stop(sprintf(
      "kernel must be one of: %s.", paste(valid_kernels, collapse = ", ")
    ))
  }
  
  # --- evaluate K(u) ----------------------------------------------------------
  k <- switch(kernel,
              
              # (*) Epanechnikov: K(u) = (3/4)(1 - u^2) * 1(|u| <= 1)
              # Optimal kernel in the MSE sense among compactly supported kernels.
              "epanechnikov" = ifelse(abs(u) <= 1, (3/4) * (1 - u^2), 0),
              
              # (*) Biweight (Quartic): K(u) = (15/16)(1 - u^2)^2 * 1(|u| <= 1)
              # Smoother than Epanechnikov, zero derivative at boundaries.
              "biweight" = ifelse(abs(u) <= 1, (15/16) * (1 - u^2)^2, 0),
              
              # (*) Triweight: K(u) = (35/32)(1 - u^2)^3 * 1(|u| <= 1)
              # Even smoother, faster decay away from 0.
              "triweight" = ifelse(abs(u) <= 1, (35/32) * (1 - u^2)^3, 0),
              
              # (*) Tricube: K(u) = (70/81)(1 - |u|^3)^3 * 1(|u| <= 1)
              # Used in LOESS smoothing.
              "tricube" = ifelse(abs(u) <= 1, (70/81) * (1 - abs(u)^3)^3, 0),
              
              # (*) Triangular: K(u) = (1 - |u|) * 1(|u| <= 1)
              # Simple tent-shaped kernel.
              "triangular" = ifelse(abs(u) <= 1, 1 - abs(u), 0),
              
              # (*) Cosine: K(u) = (pi/4) cos(pi/2 * u) * 1(|u| <= 1)
              "cosine" = ifelse(abs(u) <= 1, (pi/4) * cos((pi/2) * u), 0),
              
              # (!) Uniform: K(u) = (1/2) * 1(|u| <= 1)
              # Does NOT satisfy (B1): not Lipschitz continuous at u = +/-1.
              "uniform" = ifelse(abs(u) <= 1, 1/2, 0),
              
              # (!) Gaussian: K(u) = (1/sqrt(2*pi)) exp(-u^2/2)
              # Does NOT satisfy (B1): support is not compact.
              "gaussian" = dnorm(u)
  )
  
  return(k)
}
