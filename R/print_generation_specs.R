# =============================================================================
# print_generation_specs.R
# -----------------------------------------------------------------------------
# Print a human-readable summary of the full data-generation setup from the
# output of generate_full_snippet_data().
#
# This is a narrative companion to build_generation_specs_table().
# =============================================================================

#' @keywords internal
#' @noRd
print_generation_specs <- function(sim_data) {

  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(sim_data)) {
    stop("sim_data must be the output of generate_full_snippet_data().")
  }

  required_names <- c("config", "summary")
  if (!all(required_names %in% names(sim_data))) {
    stop("sim_data must contain at least: config and summary.")
  }

  cfg <- sim_data$config
  smy <- sim_data$summary

  # ---------------------------------------------------------------------------
  # HELPERS
  # ---------------------------------------------------------------------------
  fmt_value <- function(x) {
    if (is.null(x)) return("NULL")

    if (is.function(x)) {
      body_txt <- paste(deparse(body(x)), collapse = " ")
      return(paste0("function: ", body_txt))
    }

    if (is.atomic(x) && length(x) == 1) {
      return(as.character(x))
    }

    if (is.atomic(x) && length(x) > 1) {
      return(paste0("c(", paste(x, collapse = ", "), ")"))
    }

    if (is.list(x)) {
      if (length(x) == 0) return("list()")

      if (!is.null(names(x))) {
        parts <- mapply(
          function(nm, val) paste0(nm, "=", fmt_value(val)),
          names(x), x,
          SIMPLIFY = TRUE, USE.NAMES = FALSE
        )
        return(paste0("list(", paste(parts, collapse = ", "), ")"))
      } else {
        return(paste0(
          "list(",
          paste(vapply(x, fmt_value, character(1)), collapse = ", "),
          ")"
        ))
      }
    }

    paste(capture.output(str(x)), collapse = " ")
  }

  corr_formula_text <- function(corr_type, corr_params) {
    if (is.null(corr_type)) return("NULL")

    if (corr_type == "matern") {
      th1 <- if (!is.null(corr_params$theta1)) corr_params$theta1 else "?"
      th2 <- if (!is.null(corr_params$theta2)) corr_params$theta2 else "?"
      return(
        paste0(
          "rho(s,t) = Matern(|s-t|; theta1 = ", th1,
          ", theta2 = ", th2, ")"
        )
      )
    }

    if (corr_type == "power_exponential") {
      th1 <- if (!is.null(corr_params$theta1)) corr_params$theta1 else "?"
      th2 <- if (!is.null(corr_params$theta2)) corr_params$theta2 else "?"
      return(
        paste0(
          "rho(s,t) = exp( - ( |s-t| / ", th1, " )^", th2, " )"
        )
      )
    }

    if (corr_type == "rational_quadratic") {
      theta <- if (!is.null(corr_params$theta)) corr_params$theta else "?"
      alpha <- if (!is.null(corr_params$alpha)) corr_params$alpha else "?"
      return(
        paste0(
          "rho(s,t) = ( 1 + |s-t|^2 / (2 * ", alpha, " * ", theta,
          "^2) )^(-", alpha, ")"
        )
      )
    }

    "Unknown correlation formula"
  }

  noise_scale_text <- function(noise_structure, sigma0) {
    if (noise_structure == "homoscedastic") {
      return(paste0("constant SD = ", fmt_value(sigma0)))
    }
    if (noise_structure == "heteroscedastic") {
      return(paste0("time-varying SD = ", fmt_value(sigma0)))
    }
    fmt_value(sigma0)
  }

  seed_mode_text <- function(seed, seed_list) {
    if (!is.null(seed_list)) return("advanced block-specific seeding")
    if (!is.null(seed)) return("simple seed expanded into internal block seeds")
    "no seed provided"
  }

  cat_line <- function(title = "", width = 72, char = "=") {
    if (nchar(title) == 0) {
      cat(strrep(char, width), "\n", sep = "")
    } else {
      left <- paste0(" ", title, " ")
      pad <- max(0, width - nchar(left))
      cat(left, strrep(char, pad), "\n", sep = "")
    }
  }

  # ---------------------------------------------------------------------------
  # PRINT REPORT
  # ---------------------------------------------------------------------------
  cat("\n")
  cat_line("DATA GENERATION SPECIFICATIONS", char = "=")

  cat_line("GENERAL", char = "-")
  cat("n                  : ", fmt_value(cfg$n), "\n", sep = "")
  cat("m_mean             : ", fmt_value(cfg$m_mean), "\n", sep = "")
  cat("delta              : ", fmt_value(cfg$delta), "\n", sep = "")
  cat("road               : ", fmt_value(cfg$road), "\n", sep = "")

  cat_line("SEEDS", char = "-")
  cat("seed mode          : ", seed_mode_text(cfg$seed, cfg$seed_list), "\n", sep = "")
  cat("user seed          : ", fmt_value(cfg$seed), "\n", sep = "")
  cat("seed_list          : ", fmt_value(cfg$seed_list), "\n", sep = "")
  cat("resolved seeds     : ", fmt_value(cfg$resolved_seeds), "\n", sep = "")
  cat("design seed        : controls m_i, O_i, and T_ij generation\n", sep = "")
  cat("process seed       : controls centred process generation Z_i(T_i)\n", sep = "")
  cat("noise seed         : controls noise generation epsilon_ij\n", sep = "")

  cat_line("DESIGN", char = "-")
  cat("fO_type            : ", fmt_value(cfg$fO_type), "\n", sep = "")
  cat("fO_params          : ", fmt_value(cfg$fO_params), "\n", sep = "")
  cat("f0_type            : ", fmt_value(cfg$f0_type), "\n", sep = "")
  cat("f0_params          : ", fmt_value(cfg$f0_params), "\n", sep = "")

  cat_line("MEAN FUNCTION", char = "-")
  cat("mu(t)              : ", fmt_value(cfg$mu_fn), "\n", sep = "")

  if (identical(cfg$road, "kl")) {
    cat_line("PROCESS: KARHUNEN-LOEVE ROAD", char = "-")
    cat("representation     : Z_i(t) = sum_{k=1}^K xi_ik * phi_k(t)\n", sep = "")
    cat("score_type         : ", fmt_value(cfg$score_type), "\n", sep = "")
    cat("score_params       : ", fmt_value(cfg$score_params), "\n", sep = "")
    cat("eigenvalues        : ", fmt_value(cfg$eigenval_vec), "\n", sep = "")
    cat("n eigenfunctions   : ",
        if (is.null(cfg$eigenfn_list)) "NULL" else length(cfg$eigenfn_list),
        "\n", sep = "")
    cat("eigenfunctions     : ",
        if (is.null(cfg$eigenfn_list)) "NULL"
        else paste0("list of ", length(cfg$eigenfn_list), " functions"),
        "\n", sep = "")
  }

  if (identical(cfg$road, "decomposition")) {
    cat_line("PROCESS: DECOMPOSITION ROAD", char = "-")
    cat("representation     : C(s,t) = sigma_X(s) * rho(s,t) * sigma_X(t)\n", sep = "")
    cat("sigmaX_fn          : ", fmt_value(cfg$sigmaX_fn), "\n", sep = "")
    cat("corr_type          : ", fmt_value(cfg$corr_type), "\n", sep = "")
    cat("corr_params        : ", fmt_value(cfg$corr_params), "\n", sep = "")
    cat("corr_formula       : ", corr_formula_text(cfg$corr_type, cfg$corr_params), "\n", sep = "")
    cat("jitter             : ", fmt_value(cfg$jitter), "\n", sep = "")
  }

  cat_line("NOISE", char = "-")
  cat("noise_structure    : ", fmt_value(cfg$noise_structure), "\n", sep = "")
  cat("noise_scale        : ", noise_scale_text(cfg$noise_structure, cfg$sigma0), "\n", sep = "")
  cat("noise_mean         : ", fmt_value(cfg$noise_mean), "\n", sep = "")
  cat("noise_type         : ", fmt_value(cfg$noise_type), "\n", sep = "")
  cat("noise_params       : ", fmt_value(cfg$noise_params), "\n", sep = "")

  cat_line("SUMMARY", char = "-")
  cat("n_subjects         : ", fmt_value(smy$n_subjects), "\n", sep = "")
  cat("total observations : ", fmt_value(smy$total_observations), "\n", sep = "")
  cat("mean obs/subject   : ", fmt_value(smy$mean_obs_per_subject), "\n", sep = "")
  cat("road used          : ", fmt_value(smy$road), "\n", sep = "")
  cat("noise structure    : ", fmt_value(smy$noise_structure), "\n", sep = "")
  cat("has_noise          : ", fmt_value(smy$has_noise), "\n", sep = "")

  cat_line(char = "=")
  cat("\n")

  invisible(NULL)
}
