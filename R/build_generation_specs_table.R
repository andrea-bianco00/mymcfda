# =============================================================================
# build_generation_specs_table.R
# -----------------------------------------------------------------------------
# Build a detailed, human-readable table of the full data-generation setup
# from the output of generate_full_snippet_data().
#
# The goal is to provide a complete summary of how the simulated data were
# generated, including:
#   - general simulation setup
#   - seed logic and effective seeds used
#   - design specification
#   - mean function
#   - process-generation road and related parameters
#   - noise specification
#   - summary diagnostics
# =============================================================================


build_generation_specs_table <- function(sim_data) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(sim_data)) {
    stop("sim_data must be the output of generate_full_snippet_data().")
  }
  
  required_names <- c("config", "summary", "design_info", "process_info")
  if (!all(required_names %in% names(sim_data))) {
    stop("sim_data must contain: config, summary, design_info, process_info.")
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
        return(paste0("list(",
                      paste(vapply(x, fmt_value, character(1)), collapse = ", "),
                      ")"))
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
  
  seed_mode_text <- function(seed, seed_list, resolved_seeds) {
    if (!is.null(seed_list)) {
      return("Advanced block-specific seeding via seed_list")
    }
    if (!is.null(seed)) {
      return("Simple seed provided; internal block seeds derived automatically")
    }
    return("No seed provided")
  }
  
  noise_sigma_text <- function(noise_structure, sigma0) {
    if (noise_structure == "homoscedastic") {
      return(paste0("constant SD: ", fmt_value(sigma0)))
    }
    if (noise_structure == "heteroscedastic") {
      return(paste0("time-varying SD: ", fmt_value(sigma0)))
    }
    return(fmt_value(sigma0))
  }
  
  # ---------------------------------------------------------------------------
  # BUILD ROWS
  # ---------------------------------------------------------------------------
  rows <- list(
    c("General", "n", fmt_value(cfg$n)),
    c("General", "m_mean", fmt_value(cfg$m_mean)),
    c("General", "delta", fmt_value(cfg$delta)),
    c("General", "road", fmt_value(cfg$road)),
    
    c("Seeds", "seed_mode", seed_mode_text(cfg$seed, cfg$seed_list, cfg$resolved_seeds)),
    c("Seeds", "seed", fmt_value(cfg$seed)),
    c("Seeds", "seed_list", fmt_value(cfg$seed_list)),
    c("Seeds", "resolved_seeds", fmt_value(cfg$resolved_seeds)),
    c("Seeds", "design_seed_meaning", "controls m_i, O_i, T_ij generation"),
    c("Seeds", "process_seed_meaning", "controls centred process generation Z_i(T_i)"),
    c("Seeds", "noise_seed_meaning", "controls noise generation epsilon_ij"),
    
    c("Design", "fO_type", fmt_value(cfg$fO_type)),
    c("Design", "fO_params", fmt_value(cfg$fO_params)),
    c("Design", "f0_type", fmt_value(cfg$f0_type)),
    c("Design", "f0_params", fmt_value(cfg$f0_params)),
    
    c("Mean", "mu_fn", fmt_value(cfg$mu_fn)),
    
    c("Noise", "noise_structure", fmt_value(cfg$noise_structure)),
    c("Noise", "noise_scale", noise_sigma_text(cfg$noise_structure, cfg$sigma0)),
    c("Noise", "noise_mean", fmt_value(cfg$noise_mean)),
    c("Noise", "noise_type", fmt_value(cfg$noise_type)),
    c("Noise", "noise_params", fmt_value(cfg$noise_params)),
    
    c("Summary", "n_subjects", fmt_value(smy$n_subjects)),
    c("Summary", "total_observations", fmt_value(smy$total_observations)),
    c("Summary", "mean_obs_per_subject", fmt_value(smy$mean_obs_per_subject)),
    c("Summary", "road_used", fmt_value(smy$road)),
    c("Summary", "noise_structure_used", fmt_value(smy$noise_structure))
  )
  
  # ---------------------------------------------------------------------------
  # ROAD-SPECIFIC SECTION
  # ---------------------------------------------------------------------------
  if (identical(cfg$road, "kl")) {
    rows <- c(
      rows,
      list(
        c("Process (KL)", "representation",
          "Z_i(t) = sum_{k=1}^K xi_ik * phi_k(t)"),
        c("Process (KL)", "score_type", fmt_value(cfg$score_type)),
        c("Process (KL)", "score_params", fmt_value(cfg$score_params)),
        c("Process (KL)", "eigenval_vec", fmt_value(cfg$eigenval_vec)),
        c("Process (KL)", "number_of_eigenfunctions",
          if (is.null(cfg$eigenfn_list)) "NULL" else as.character(length(cfg$eigenfn_list))),
        c("Process (KL)", "eigenfn_list",
          if (is.null(cfg$eigenfn_list)) "NULL"
          else paste0("list of ", length(cfg$eigenfn_list), " functions"))
      )
    )
  }
  
  if (identical(cfg$road, "decomposition")) {
    rows <- c(
      rows,
      list(
        c("Process (Decomposition)", "representation",
          "C(s,t) = sigma_X(s) * rho(s,t) * sigma_X(t)"),
        c("Process (Decomposition)", "sigmaX_fn", fmt_value(cfg$sigmaX_fn)),
        c("Process (Decomposition)", "corr_type", fmt_value(cfg$corr_type)),
        c("Process (Decomposition)", "corr_params", fmt_value(cfg$corr_params)),
        c("Process (Decomposition)", "corr_formula",
          corr_formula_text(cfg$corr_type, cfg$corr_params)),
        c("Process (Decomposition)", "jitter", fmt_value(cfg$jitter))
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # BUILD TABLE
  # ---------------------------------------------------------------------------
  specs_table <- as.data.frame(
    do.call(rbind, rows),
    stringsAsFactors = FALSE
  )
  
  names(specs_table) <- c("section", "parameter", "value")
  rownames(specs_table) <- NULL
  
  return(specs_table)
}
