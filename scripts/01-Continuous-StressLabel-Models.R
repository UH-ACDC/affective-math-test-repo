#!/usr/bin/env Rscript
# ==============================================================================
# Affective Math Test — Physiological Signal Models (Continuous + Stress-Label)
#
# Author: Ioannis Pavlidis
# Affiliation: Affective and Data Computing Lab — University of Houston
#
# Repository script: scripts/01-Continuous-StressLabel-Models.R
#
# Roadmap (high-level modules)
#   0) Config + GitHub-safe paths
#   1) Packages
#   2) Utilities (formatting, safe fitting, AIC selection, R², coefficient extraction)
#   3) Data load + standardization (Q-level file → analysis frame D + channel map)
#   4) Model fitting
#        4a) Eq(10) continuous LMMs: FULL model + AIC reduction (Table 3)
#        4b) Eq(11) stress-label GLMMs: FULL model + AIC reduction (Table 4)
#   5) Table assembly + exports (CSV + TeX; plus AIC traces/summaries)
#   6) Figure 10 (emmeans by QTYPE; PDF + PNG)
#
# Outputs (GitHub naming convention)
#   • reports/Table3_Continuous_Models_N50.{csv,tex}
#   • reports/Table3_Continuous_Models_AICreduced_N50.{csv,tex}
#   • reports/Table3_Continuous_Models_AICtrace_N50.csv
#   • reports/Table3_Continuous_Models_KeptTerms_N50.csv
#   • reports/Table3_Continuous_Models_AICsummary_N50.{csv,tex}
#   • reports/Table4_StressLabel_Models_N50.{csv,tex}
#   • reports/Table4_StressLabel_Models_AICreduced_N50.{csv,tex}
#   • reports/Table4_StressLabel_Models_AICtrace_N50.csv
#   • reports/Table4_StressLabel_Models_KeptTerms_N50.csv
#   • reports/Table4_StressLabel_Models_AICsummary_N50.{csv,tex}
#   • figures/Figure10.pdf
#   • figures/Figure10.png
#
# Notes
#   • Paths resolve relative to the script location (RStudio, Rscript, CI).
#   • No absolute paths.
#   • Stress labels are assumed to be mean-threshold labels already present in the
#     processed Q-level dataset (e.g., Stress_fp, Stress_hraw, ...).
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# Module 0 — GitHub-safe paths
# ==============================================================================
get_script_dir <- function() {
  # 1) Rscript --file=...
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
  }

  # 2) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NA_character_)
    if (!is.na(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }

  # 3) Fallback: current working directory
  getwd()
}

make_paths <- function(project_root) {
  cand_data <- c(
    file.path(project_root, "data", "processed"),
    file.path(project_root, "data", "Processed"),
    file.path(project_root, "Data", "processed"),
    file.path(project_root, "Data", "Processed")
  )
  data_dir <- cand_data[dir.exists(cand_data)][1]
  if (is.na(data_dir)) stop("Cannot find processed data folder under {data,Data}/{processed,Processed}.", call. = FALSE)

  figures_dir <- if (dir.exists(file.path(project_root, "figures"))) file.path(project_root, "figures") else file.path(project_root, "Figures")

  list(
    root        = project_root,
    data_dir    = normalizePath(data_dir, winslash = "/", mustWork = FALSE),
    figures_dir = normalizePath(figures_dir, winslash = "/", mustWork = FALSE),
    reports_dir = normalizePath(file.path(project_root, "reports"), winslash = "/", mustWork = FALSE)
  )
}

ensure_dirs <- function(paths) {
  dir.create(paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$reports_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

# ==============================================================================
# Module 1 — Packages (only what we use)
# ==============================================================================
load_packages <- function() {
  suppressPackageStartupMessages({
    library(dplyr)     # data manipulation + bind_rows
    library(ggplot2)   # plotting + ggsave
    library(cowplot)   # figure layout
    library(lme4)      # mixed models (glmer, VarCorr)
  })

  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("Package 'lmerTest' is required for Eq(10) p-values. Install: install.packages('lmerTest')", call. = FALSE)
  }
  suppressPackageStartupMessages(library(lmerTest))  # Eq(10): p-values for fixed effects

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for Figure 10. Install: install.packages('emmeans')", call. = FALSE)
  }

  # Used for the right-side strip grobs (namespace calls below)
  if (!requireNamespace("gtable", quietly = TRUE)) {
    stop("Package 'gtable' is required for Figure 10 strip labels. Install: install.packages('gtable')", call. = FALSE)
  }
  invisible(TRUE)
}

# ==============================================================================
# Module 2 — Utilities
# ==============================================================================
pick_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  NA_character_
}

zscore <- function(x) as.numeric(scale(x))

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

coerce_stress01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(ifelse(is.na(x), NA_integer_, as.integer(x != 0)))
  xx <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(xx))
  out[xx %in% c("s","stress","stressed","1","true","t","rhs")] <- 1L
  out[xx %in% c("ns","no","nostress","0","false","f","rls")] <- 0L
  suppressWarnings({
    cand <- as.integer(xx)
    out[is.na(out) & !is.na(cand)] <- cand[is.na(out) & !is.na(cand)]
  })
  out
}

safe_AIC <- function(m) {
  if (is.null(m)) return(Inf)
  tryCatch(AIC(m), error = function(e) Inf)
}

safe_fit_last_error <- NULL
safe_fit <- function(expr) {
  safe_fit_last_error <<- NULL
  tryCatch(
    eval(expr, envir = parent.frame()),
    error = function(e) {
      safe_fit_last_error <<- conditionMessage(e)
      NULL
    }
  )
}

# Greedy backward AIC selection:
#   • start from FULL model (all candidate terms)
#   • try dropping each droppable term once
#   • keep the drop that reduces AIC the most
#   • repeat until no drop improves AIC
greedy_backward_AIC <- function(fit_fun, start_terms, protect_terms = character(0), AIC_EPS = 1e-8) {
  cur_terms <- start_terms
  cur_fit   <- fit_fun(cur_terms)
  cur_aic   <- safe_AIC(cur_fit)

  trace <- data.frame(
    step = 0, action = "start", dropped = "", AIC = cur_aic,
    terms = paste(cur_terms, collapse = " + "),
    stringsAsFactors = FALSE
  )

  step <- 0
  improved <- TRUE
  while (improved) {
    improved <- FALSE
    drop_candidates <- setdiff(cur_terms, protect_terms)
    if (length(drop_candidates) == 0) break

    best_aic   <- cur_aic
    best_fit   <- cur_fit
    best_terms <- cur_terms
    best_drop  <- NULL

    for (trm in drop_candidates) {
      trial_terms <- setdiff(cur_terms, trm)
      trial_fit   <- fit_fun(trial_terms)
      trial_aic   <- safe_AIC(trial_fit)

      if (trial_aic + AIC_EPS < best_aic) {
        best_aic   <- trial_aic
        best_fit   <- trial_fit
        best_terms <- trial_terms
        best_drop  <- trm
      }
    }

    if (!is.null(best_drop) && best_aic + AIC_EPS < cur_aic) {
      step <- step + 1
      cur_aic   <- best_aic
      cur_fit   <- best_fit
      cur_terms <- best_terms
      improved  <- TRUE
      trace <- rbind(trace, data.frame(
        step = step, action = "drop", dropped = best_drop,
        AIC = cur_aic, terms = paste(cur_terms, collapse = " + "),
        stringsAsFactors = FALSE
      ))
    }
  }

  list(model = cur_fit, terms = cur_terms, trace = trace)
}

# Manual Nakagawa-style R² for LMM (robust / avoids additional packages)
r2_lmm_nakagawa_manual <- function(m, grp_pid = "PID", grp_qid = "QID") {
  yhat_fixed <- as.numeric(predict(m, re.form = NA))
  var_fixed  <- stats::var(yhat_fixed, na.rm = TRUE)

  vc_df <- as.data.frame(lme4::VarCorr(m))
  var_part <- vc_df$vcov[vc_df$grp == grp_pid & vc_df$var1 == "(Intercept)"]
  var_item <- vc_df$vcov[vc_df$grp == grp_qid & vc_df$var1 == "(Intercept)"]
  var_part <- if (length(var_part) == 0) 0 else var_part[1]
  var_item <- if (length(var_item) == 0) 0 else var_item[1]

  var_resid <- vc_df$vcov[vc_df$grp == "Residual"]
  var_resid <- if (length(var_resid) == 0) stats::sigma(m)^2 else var_resid[1]

  denom <- var_fixed + var_part + var_item + var_resid
  r2m <- var_fixed / denom
  r2c <- (var_fixed + var_part + var_item) / denom
  c(R2_m = r2m, R2_c = r2c)
}

# Latent-scale R² for logit GLMM
r2_glmm_latent <- function(m) {
  var_res   <- (pi^2) / 3
  eta_fixed <- as.numeric(predict(m, re.form = NA, type = "link"))
  var_fixed <- var(eta_fixed, na.rm = TRUE)
  vc <- as.data.frame(VarCorr(m))
  var_u <- vc$vcov[vc$grp=="PID" & vc$var1=="(Intercept)"]; var_u <- ifelse(length(var_u)==0, 0, var_u[1])
  var_v <- vc$vcov[vc$grp=="QID" & vc$var1=="(Intercept)"]; var_v <- ifelse(length(var_v)==0, 0, var_v[1])
  denom <- var_fixed + var_u + var_v + var_res
  c(R2_m = var_fixed/denom, R2_c = (var_fixed+var_u+var_v)/denom)
}

extract_fix_lmm <- function(m) {
  if (is.null(m)) return(data.frame(term=character(0), b=numeric(0), se=numeric(0), stat=numeric(0), p=numeric(0)))
  s <- summary(m); co <- as.data.frame(s$coefficients)
  data.frame(
    term = rownames(co),
    b    = co[,"Estimate"],
    se   = co[,"Std. Error"],
    stat = co[,"t value"],
    p    = if ("Pr(>|t|)" %in% colnames(co)) co[,"Pr(>|t|)"] else NA_real_,
    row.names = NULL
  )
}

extract_fix_glmm <- function(m) {
  if (is.null(m)) return(data.frame(term=character(0), b=numeric(0), se=numeric(0), stat=numeric(0), p=numeric(0)))
  s <- summary(m); co <- as.data.frame(s$coefficients)
  data.frame(
    term = rownames(co),
    b    = co[,"Estimate"],
    se   = co[,"Std. Error"],
    stat = co[,"z value"],
    p    = co[,"Pr(>|z|)"],
    row.names = NULL
  )
}

# Safer PDF device: prefer cairo if available, otherwise pdf()
save_plot_dual <- function(p, filename_base, out_dir, width, height, dpi = 300, ...) {
  pdf_path <- file.path(out_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(out_dir, paste0(filename_base, ".png"))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(pdf_path, p, width = width, height = height, units = "in", device = pdf_device, bg = "white", ...)
  ggplot2::ggsave(png_path, p, width = width, height = height, units = "in", dpi = dpi, bg = "white", ...)
  invisible(list(pdf = pdf_path, png = png_path))
}

# ==============================================================================
# Module 3 — Data load + standardization
# ==============================================================================
load_and_prepare_data <- function(data_dir, N_TARGET) {
  # ---- Load Q-level data ----
  candidate_files <- c(
    file.path(data_dir, sprintf("Affective_Math_Qlevel_Data_N%s.csv", N_TARGET)),
    file.path(data_dir, "Affective_Math_Qlevel_Data_N50.csv")
  )
  in_file <- candidate_files[file.exists(candidate_files)][1]
  if (is.na(in_file) || !nzchar(in_file)) {
    stop(
      "Could not find Q-level CSV in: ", data_dir,
      "\nExpected something like: data/processed/Affective_Math_Qlevel_Data_N", N_TARGET, ".csv",
      "\nIf you just cloned the repo, run the preprocessing script first to generate processed data.",
      call. = FALSE
    )
  }
  Q <- read.csv(in_file, stringsAsFactors = FALSE)
  message("Using input: ", in_file)

  # ---- Column mapping (robust to name variants) ----
  pid_col   <- pick_col(Q, c("ParticipantID","Participant","PID"))
  qid_col   <- pick_col(Q, c("Question.Name","QuestionName","QName","QID","QuestionID"))
  qtype_col <- pick_col(Q, c("Question.Type","QuestionType","QType","QTYPE"))
  qtime_col <- pick_col(Q, c("QTime","QTIME","TimeOnQuestion","Time"))

  grade_col <- pick_col(Q, c("Grade","QGRADE","QGrade"))
  sex_col   <- pick_col(Q, c("Gender","Sex","SEX"))
  sai_col   <- pick_col(Q, c("SAI","zSAI"))

  # Continuous outcomes (normalized signals)
  y_fp   <- pick_col(Q, c("FPNorm","PPNorm"))
  y_hraw <- pick_col(Q, c("HR.AWNorm","HRAWNorm"))
  y_hre4 <- pick_col(Q, c("HR.E4Norm","HRE4Norm"))
  y_hrv  <- pick_col(Q, c("HRVNorm","NHRVNorm","NHRV"))

  # Stress-label outcomes (binary labels; mean-threshold per channel)
  s_fp   <- pick_col(Q, c("Stress.fp","Stress_FP","StressFp"))
  s_hraw <- pick_col(Q, c("Stress.hraw","Stress_HRAW","StressHraw"))
  s_hre4 <- pick_col(Q, c("Stress.hre4","Stress_HRE4","StressHre4"))
  s_hrv  <- pick_col(Q, c("Stress.nhrv","Stress_NHRV","StressNhrv"))

  if (any(is.na(c(pid_col, qid_col, qtype_col, qtime_col)))) stop("Missing required columns: PID/QID/QTYPE/QTIME", call. = FALSE)
  if (any(is.na(c(y_fp, y_hraw, y_hre4, y_hrv)))) stop("Missing one or more continuous outcome columns.", call. = FALSE)
  if (any(is.na(c(s_fp, s_hraw, s_hre4, s_hrv)))) stop("Missing one or more stress-label outcome columns.", call. = FALSE)

  # ---- Analysis frame D (standardize covariates) ----
  D <- Q %>%
    mutate(
      PID   = factor(.data[[pid_col]]),
      QID   = factor(.data[[qid_col]]),
      QTYPE = factor(.data[[qtype_col]], levels = c("V","A","W")),
      QTime = suppressWarnings(as.numeric(.data[[qtime_col]]))
    )

  D$QGRADE1 <- if (!is.na(grade_col)) as.integer(suppressWarnings(as.numeric(D[[grade_col]])) == 1) else NA_integer_

  if (!is.na(sex_col)) {
    sx <- toupper(trimws(as.character(D[[sex_col]])))
    D$SexF <- as.integer(sx %in% c("F","FEMALE","WOMAN","GIRL"))
  } else {
    D$SexF <- NA_integer_
  }

  # zSAI is always standardized here (even if the source is already zSAI)
  D$zSAI <- if (!is.na(sai_col)) zscore(suppressWarnings(as.numeric(D[[sai_col]]))) else NA_real_

  D$logQTime  <- ifelse(is.finite(D$QTime) & D$QTime > 0, log(D$QTime), NA_real_)
  D$zlogQTime <- as.numeric(scale(D$logQTime))

  channels <- list(
    NFP    = list(y = y_fp,   s = s_fp),
    NHR_AW = list(y = y_hraw, s = s_hraw),
    NHR_E4 = list(y = y_hre4, s = s_hre4),
    NHRV   = list(y = y_hrv,  s = s_hrv)
  )

  list(D = D, channels = channels)
}

# ==============================================================================
# Module 4 — Model fitting (FULL + AIC reduction)
# ==============================================================================
fit_eq10_one <- function(D, y_col, outcome_name, EQ10_REML = FALSE, AIC_EPS = 1e-8) {
  # Eq(10) uses a fixed-N (complete-case) frame so FULL vs AIC-reduced models are comparable by AIC.
  # FULL model happens here: `m_full <- fit_fun(full_terms)`
  # AIC reduction happens here: `sel <- greedy_backward_AIC(...)` and returns `m_red`.

  Dm0 <- D %>%
    mutate(Yraw = suppressWarnings(as.numeric(.data[[y_col]]))) %>%
    filter(
      is.finite(Yraw),
      !is.na(QTYPE), !is.na(PID), !is.na(QID),
      is.finite(zlogQTime),
      is.finite(QGRADE1),
      is.finite(SexF),
      is.finite(zSAI)
    ) %>%
    mutate(Yz = zscore(Yraw))

  base_terms <- c("QTYPE", "zlogQTime")          # forced in all models (paper behavior)
  opt_terms  <- c("QGRADE1", "SexF", "zSAI")     # AIC may drop
  full_terms <- c(base_terms, opt_terms)

  fit_fun <- function(terms_vec) {
    rhs <- paste(terms_vec, collapse = " + ")
    fml <- as.formula(paste0("Yz ~ ", rhs, " + (1|PID) + (1|QID)"))
    safe_fit(quote(
      lmerTest::lmer(
        fml, data = Dm0, REML = EQ10_REML,
        control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
    ))
  }

  # ---- FULL model (primary reporting; Table 3 + Figure 10a) ----
  m_full <- fit_fun(full_terms)
  err_full <- safe_fit_last_error
  aic_full <- safe_AIC(m_full); if (!is.finite(aic_full)) aic_full <- NA_real_

  # ---- AIC-reduced model (secondary reporting; "retained terms" files) ----
  sel <- greedy_backward_AIC(fit_fun, start_terms = full_terms, protect_terms = base_terms, AIC_EPS = AIC_EPS)
  m_red <- sel$model
  err_red <- if (is.null(m_red)) safe_fit_last_error else NA_character_
  aic_red <- safe_AIC(m_red); if (!is.finite(aic_red)) aic_red <- NA_real_
  delta_aic <- if (is.finite(aic_full) && is.finite(aic_red)) (aic_full - aic_red) else NA_real_

  list(
    model_full    = m_full,
    model_reduced = m_red,
    model_paper   = m_full,
    kept_terms    = sel$terms,
    trace         = sel$trace,
    frame_n       = nrow(Dm0),
    outcome       = outcome_name,
    aic_full      = aic_full,
    aic_reduced   = aic_red,
    delta_aic     = delta_aic,
    err_full      = err_full,
    err_reduced   = err_red
  )
}

fit_eq11_one <- function(D, s_col, outcome_name, AIC_EPS = 1e-8) {
  # Eq(11) uses a fixed-N (complete-case) frame so FULL vs AIC-reduced models are comparable by AIC.
  # FULL model happens here: `m_full <- fit_fun(full_terms)`
  # AIC reduction happens here: `sel <- greedy_backward_AIC(...)` and returns `m_red`.
  # If AIC selection fails, we fall back to a base model (QTYPE + zlogQTime).

  Dm0 <- D %>%
    mutate(S01 = coerce_stress01(.data[[s_col]])) %>%
    filter(
      !is.na(S01),
      !is.na(QTYPE), !is.na(PID), !is.na(QID),
      is.finite(zlogQTime),
      is.finite(QGRADE1),
      is.finite(SexF),
      is.finite(zSAI)
    )

  # Degenerate checks
  if (nrow(Dm0) == 0 || length(unique(Dm0$S01)) < 2 || nlevels(droplevels(Dm0$QTYPE)) < 2) {
    warning("Eq11 ", outcome_name, ": degenerate/empty frame after fixed-N filtering (n=", nrow(Dm0), ").")
    return(list(
      model_full=NULL, model_reduced=NULL, model_paper=NULL, kept_terms=character(0),
      trace=data.frame(step=0, action="start", dropped="", AIC=Inf, terms="", stringsAsFactors=FALSE),
      frame_n=nrow(Dm0), outcome=outcome_name,
      aic_full=NA_real_, aic_reduced=NA_real_, delta_aic=NA_real_,
      err_full="degenerate/empty frame", err_reduced=NA_character_
    ))
  }

  base_terms <- c("QTYPE", "zlogQTime")          # forced
  opt_terms  <- c("QGRADE1", "SexF", "zSAI")     # AIC may drop
  full_terms <- c(base_terms, opt_terms)

  fit_fun <- function(terms_vec) {
    rhs <- paste(terms_vec, collapse = " + ")
    fml <- as.formula(paste0("S01 ~ ", rhs, " + (1|PID) + (1|QID)"))
    safe_fit(quote(
      lme4::glmer(
        fml, data = Dm0,
        family = binomial(link="logit"),
        nAGQ = 1,
        control = lme4::glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
    ))
  }

  # ---- FULL model (primary reporting; Table 4 + Figure 10b) ----
  m_full <- fit_fun(full_terms)
  err_full <- safe_fit_last_error
  aic_full <- safe_AIC(m_full); if (!is.finite(aic_full)) aic_full <- NA_real_

  # ---- AIC-reduced model (secondary reporting) ----
  sel <- greedy_backward_AIC(fit_fun, start_terms = full_terms, protect_terms = base_terms, AIC_EPS = AIC_EPS)
  m_red <- sel$model
  err_red <- if (is.null(m_red)) safe_fit_last_error else NA_character_
  aic_red <- safe_AIC(m_red); if (!is.finite(aic_red)) aic_red <- NA_real_
  delta_aic <- if (is.finite(aic_full) && is.finite(aic_red)) (aic_full - aic_red) else NA_real_

  # Fallback if selection failed
  if (is.null(m_red)) {
    warning("Eq11 ", outcome_name, ": AIC selection produced NULL model; attempting fallback base GLMM...")
    f_base <- fit_fun(base_terms)
    if (!is.null(f_base)) {
      m_red <- f_base
      aic_red <- safe_AIC(m_red); if (!is.finite(aic_red)) aic_red <- NA_real_
      delta_aic <- if (is.finite(aic_full) && is.finite(aic_red)) (aic_full - aic_red) else NA_real_
      sel$terms <- base_terms
    }
  }

  list(
    model_full    = m_full,
    model_reduced = m_red,
    model_paper   = m_full,
    kept_terms    = sel$terms,
    trace         = sel$trace,
    frame_n       = nrow(Dm0),
    outcome       = outcome_name,
    aic_full      = aic_full,
    aic_reduced   = aic_red,
    delta_aic     = delta_aic,
    err_full      = err_full,
    err_reduced   = err_red
  )
}

# ==============================================================================
# Module 5 — Tables (assembly + CSV/TeX exports)
# ==============================================================================
make_wide_table3 <- function(tab) {
  terms <- c("(Intercept)","QTYPEA","QTYPEW","zlogQTime","QGRADE1","SexF","zSAI")
  labs  <- c("Intercept","QTYPE=A","QTYPE=W","zlog(QTIME)","QGRADE=1","SEX=F","zSAI")
  out <- data.frame(Predictor=labs, stringsAsFactors=FALSE)
  for (nm in names(tab)) {
    fx <- tab[[nm]]$fx
    idx <- match(terms, fx$term)
    out[[paste0(nm,"_b")]]  <- fx$b[idx]
    out[[paste0(nm,"_SE")]] <- fx$se[idx]
    out[[paste0(nm,"_t")]]  <- fx$stat[idx]
    out[[paste0(nm,"_p")]]  <- fx$p[idx]
  }
  footer <- data.frame(Predictor=c("Observations","R2_m","R2_c","sigma_u","sigma_v","sigma_eps"))
  for (nm in names(tab)) {
    footer[[paste0(nm,"_b")]] <- c(tab[[nm]]$n, tab[[nm]]$r2["R2_m"], tab[[nm]]$r2["R2_c"],
                                   tab[[nm]]$sd_u, tab[[nm]]$sd_v, tab[[nm]]$sd_e)
    footer[[paste0(nm,"_SE")]] <- NA
    footer[[paste0(nm,"_t")]]  <- NA
    footer[[paste0(nm,"_p")]]  <- NA
  }
  dplyr::bind_rows(out, footer)
}

make_wide_table4 <- function(tab) {
  terms <- c("(Intercept)","QTYPEA","QTYPEW","zlogQTime","QGRADE1","SexF","zSAI")
  labs  <- c("Intercept","QTYPE=A","QTYPE=W","zlog(QTIME)","QGRADE=1","SEX=F","zSAI")
  out <- data.frame(Predictor=labs, stringsAsFactors=FALSE)
  for (nm in names(tab)) {
    fx <- tab[[nm]]$fx
    idx <- match(terms, fx$term)
    out[[paste0(nm,"_b")]]  <- fx$b[idx]
    out[[paste0(nm,"_SE")]] <- fx$se[idx]
    out[[paste0(nm,"_z")]]  <- fx$stat[idx]
    out[[paste0(nm,"_p")]]  <- fx$p[idx]
  }
  footer <- data.frame(Predictor=c("Observations","R2_m","R2_c","sigma_u","sigma_v"))
  for (nm in names(tab)) {
    footer[[paste0(nm,"_b")]] <- c(tab[[nm]]$n, tab[[nm]]$r2["R2_m"], tab[[nm]]$r2["R2_c"],
                                   tab[[nm]]$sd_u, tab[[nm]]$sd_v)
    footer[[paste0(nm,"_SE")]] <- NA
    footer[[paste0(nm,"_z")]]  <- NA
    footer[[paste0(nm,"_p")]]  <- NA
  }
  dplyr::bind_rows(out, footer)
}

write_tex_table3 <- function(path, tab, delta_aic = NULL) {
  sigs <- names(tab)
  disp <- c("(Intercept)","QTYPEA","QTYPEW","zlogQTime","QGRADE1","SexF","zSAI")
  lab <- function(trm){
    if (trm=="(Intercept)") return("Intercept")
    if (trm=="QTYPEA")     return("$QTYPE=A$")
    if (trm=="QTYPEW")     return("$QTYPE=W$")
    if (trm=="zlogQTime")  return("$z\\log(QTIME)$")
    if (trm=="QGRADE1")    return("$QGRADE{=}1$")
    if (trm=="SexF")       return("$SEX{=}F$")
    if (trm=="zSAI")       return("$zSAI$")
    trm
  }

  con <- file(path, "wt"); on.exit(close(con), add=TRUE)
  cat("% --- AUTO-GENERATED (Table 3; Eq(10) continuous LMMs) ---\n", file=con)
  cat("\\begin{table*}[!t]\n\\centering\\scriptsize\n", file=con)
  cat("\\begin{tabular}{@{}l", file=con)
  for (k in sigs) cat(" r r r l", file=con)
  cat("@{}}\\toprule\n", file=con)

  cat("& ", file=con)
  cat(paste(sprintf("\\multicolumn{4}{c}{%s}", sigs), collapse=" & "), file=con)
  cat(" \\\\\n", file=con)
  cat("\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-13}\\cmidrule(lr){14-17}\n", file=con)

  cat("Predictor", file=con)
  for (k in sigs) cat(" & $b$ & SE & $t$ & $p$", file=con)
  cat(" \\\\\n\\midrule\n", file=con)

  for (trm in disp) {
    cat(lab(trm), file=con)
    for (k in sigs) {
      fx <- tab[[k]]$fx
      row <- fx[fx$term==trm, , drop=FALSE]
      if (nrow(row)==0) { cat(" &  &  &  & ", file=con); next }
      cat(sprintf(" & %.3f & %.3f & %.3f & %s",
                  row$b[1], row$se[1], row$stat[1], fmt_p(row$p[1])), file=con)
    }
    cat(" \\\\\n", file=con)
  }

  cat("\\midrule\n", file=con)

  footer <- c("Observations","$R^2_m$","$R^2_c$","$\\hat{\\sigma}_u$","$\\hat{\\sigma}_v$","$\\hat{\\sigma}_\\varepsilon$")
  for (nm in footer) {
    cat(nm, file=con)
    for (k in sigs) {
      if (nm=="Observations") cat(sprintf(" & %.0f &  &  & ", tab[[k]]$n), file=con)
      else if (nm=="$R^2_m$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$r2["R2_m"]), file=con)
      else if (nm=="$R^2_c$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$r2["R2_c"]), file=con)
      else if (nm=="$\\hat{\\sigma}_u$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$sd_u), file=con)
      else if (nm=="$\\hat{\\sigma}_v$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$sd_v), file=con)
      else cat(sprintf(" & %.3f &  &  & ", tab[[k]]$sd_e), file=con)
    }
    cat(" \\\\\n", file=con)
  }

  if (!is.null(delta_aic)) {
    cat("$\\Delta$AIC (full$-$reduced)", file=con)
    for (k in sigs) {
      da <- delta_aic[[k]]
      if (is.null(da) || is.na(da)) cat(" &  &  &  & ", file=con)
      else cat(sprintf(" & %.3f &  &  & ", da), file=con)
    }
    cat(" \\\\\n", file=con)
  }

  cat("\\bottomrule\n\\end{tabular}\n\\end{table*}\n", file=con)
}

write_tex_table4 <- function(path, tab, delta_aic = NULL) {
  sigs <- names(tab)
  disp <- c("(Intercept)","QTYPEA","QTYPEW","zlogQTime","QGRADE1","SexF","zSAI")
  lab <- function(trm){
    if (trm=="(Intercept)") return("Intercept")
    if (trm=="QTYPEA")     return("$QTYPE=A$")
    if (trm=="QTYPEW")     return("$QTYPE=W$")
    if (trm=="zlogQTime")  return("$z\\log(QTIME)$")
    if (trm=="QGRADE1")    return("$QGRADE{=}1$")
    if (trm=="SexF")       return("$SEX{=}F$")
    if (trm=="zSAI")       return("$zSAI$")
    trm
  }

  con <- file(path, "wt"); on.exit(close(con), add=TRUE)
  cat("% --- AUTO-GENERATED (Table 4; Eq(11) stress-label GLMMs) ---\n", file=con)
  cat("\\begin{table*}[!t]\n\\centering\\scriptsize\n", file=con)
  cat("\\begin{tabular}{@{}l", file=con)
  for (k in sigs) cat(" r r r l", file=con)
  cat("@{}}\\toprule\n", file=con)

  cat("& ", file=con)
  cat(paste(sprintf("\\multicolumn{4}{c}{%s}", sigs), collapse=" & "), file=con)
  cat(" \\\\\n", file=con)
  cat("\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-13}\\cmidrule(lr){14-17}\n", file=con)

  cat("Predictor", file=con)
  for (k in sigs) cat(" & $b$ & SE & $z$ & $p$", file=con)
  cat(" \\\\\n\\midrule\n", file=con)

  for (trm in disp) {
    cat(lab(trm), file=con)
    for (k in sigs) {
      fx <- tab[[k]]$fx
      row <- fx[fx$term==trm, , drop=FALSE]
      if (nrow(row)==0) { cat(" &  &  &  & ", file=con); next }
      cat(sprintf(" & %.3f & %.3f & %.3f & %s",
                  row$b[1], row$se[1], row$stat[1], fmt_p(row$p[1])), file=con)
    }
    cat(" \\\\\n", file=con)
  }

  cat("\\midrule\n", file=con)

  footer <- c("Observations","$R^2_m$","$R^2_c$","$\\hat{\\sigma}_u$","$\\hat{\\sigma}_v$","$\\sigma_\\varepsilon$")
  for (nm in footer) {
    cat(nm, file=con)
    for (k in sigs) {
      if (nm=="Observations") cat(sprintf(" & %.0f &  &  & ", tab[[k]]$n), file=con)
      else if (nm=="$R^2_m$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$r2["R2_m"]), file=con)
      else if (nm=="$R^2_c$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$r2["R2_c"]), file=con)
      else if (nm=="$\\hat{\\sigma}_u$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$sd_u), file=con)
      else if (nm=="$\\hat{\\sigma}_v$") cat(sprintf(" & %.3f &  &  & ", tab[[k]]$sd_v), file=con)
      else cat(" & $\\pi/\\sqrt{3}$ &  &  & ", file=con)
    }
    cat(" \\\\\n", file=con)
  }

  if (!is.null(delta_aic)) {
    cat("$\\Delta$AIC (full$-$reduced)", file=con)
    for (k in sigs) {
      da <- delta_aic[[k]]
      if (is.null(da) || is.na(da)) cat(" &  &  &  & ", file=con)
      else cat(sprintf(" & %.3f &  &  & ", da), file=con)
    }
    cat(" \\\\\n", file=con)
  }

  cat("\\bottomrule\n\\end{tabular}\n\\end{table*}\n", file=con)
}

# AIC-retained summaries (compact)
label_term <- function(trm) {
  if (trm=="QTYPEA") return("QTYPE=A")
  if (trm=="QTYPEW") return("QTYPE=W")
  if (trm=="zlogQTime") return("zlog(QTIME)")
  if (trm=="QGRADE1") return("QGRADE=1")
  if (trm=="SexF") return("SEX=F")
  if (trm=="zSAI") return("zSAI")
  trm
}

make_aic_summary <- function(fits_list, kind = c("eq10","eq11")) {
  kind <- match.arg(kind)
  rows <- list()
  cand <- c("QTYPEA","QTYPEW","zlogQTime","QGRADE1","SexF","zSAI")

  for (nm in names(fits_list)) {
    m <- fits_list[[nm]]$model_reduced
    if (is.null(m)) next
    fx <- if (kind == "eq10") extract_fix_lmm(m) else extract_fix_glmm(m)
    fx <- fx[fx$term %in% cand, , drop=FALSE]
    if (nrow(fx)==0) next
    fx$Outcome   <- nm
    fx$Predictor <- vapply(fx$term, label_term, FUN.VALUE=character(1))
    fx$DeltaAIC  <- fits_list[[nm]]$delta_aic
    rows[[nm]] <- fx[, c("Outcome","Predictor","b","se","stat","p","DeltaAIC")]
  }
  if (length(rows)==0) return(data.frame())
  do.call(rbind, rows)
}

write_tex_aic_summary <- function(path, df, stat_name=c("t","z")) {
  stat_name <- match.arg(stat_name)
  con <- file(path, "wt"); on.exit(close(con), add=TRUE)
  cat("% --- AUTO-GENERATED: AIC-retained terms only ---\n", file=con)
  cat("\\begin{table*}[!t]\n\\centering\\scriptsize\n", file=con)
  cat("\\begin{tabular}{@{}l l r r r l r@{}}\\toprule\n", file=con)
  cat("Outcome & Predictor & $b$ & SE & $", stat_name, "$ & $p$ & $\\Delta$AIC \\\\\n\\midrule\n", sep="", file=con)

  if (nrow(df)==0) {
    cat("\\multicolumn{7}{c}{(no models fit)}\\\\\n", file=con)
  } else {
    df <- df[order(df$Outcome, df$Predictor), ]
    outcomes <- unique(df$Outcome)
    for (o in outcomes) {
      sub <- df[df$Outcome==o, ]
      for (i in seq_len(nrow(sub))) {
        da <- if (i==1) sub$DeltaAIC[i] else NA_real_
        da_txt <- ifelse(is.na(da), "", sprintf("%.3f", da))
        cat(sprintf("%s & %s & %.3f & %.3f & %.3f & %s & %s \\\\\n",
                    sub$Outcome[i], sub$Predictor[i],
                    sub$b[i], sub$se[i], sub$stat[i], fmt_p(sub$p[i]), da_txt),
            file=con)
      }
      cat("\\addlinespace\n", file=con)
    }
  }

  cat("\\bottomrule\\end{tabular}\n\\end{table*}\n", file=con)
}

# ==============================================================================
# Module 6 — Figure 10 (emmeans panels)
# ==============================================================================
make_figure10 <- function(fits10, fits11) {

  sig_levels <- c("NS","*","**","***")
  sig_colors <- c("NS"="#000000", "*"="#00BFC4", "**"="#E69F00", "***"="#FF0000")
  VREF_COLOR <- "#DADADA"

  p_to_star <- function(p) {
    if (is.na(p)) return("NS")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    "NS"
  }

  qtype_emm_df_cont <- function(fit) {
    if (is.null(fit)) return(NULL)
    emm <- emmeans::emmeans(fit, ~ QTYPE)
    sm  <- as.data.frame(summary(emm, infer = c(TRUE, TRUE)))
    sm$pred <- sm$emmean
    if ("lower.CL" %in% names(sm) && "upper.CL" %in% names(sm)) {
      sm$lo <- sm$lower.CL; sm$hi <- sm$upper.CL
    } else if ("asymp.LCL" %in% names(sm) && "asymp.UCL" %in% names(sm)) {
      sm$lo <- sm$asymp.LCL; sm$hi <- sm$asymp.UCL
    } else {
      stop("emmeans (continuous) did not contain CI bounds.", call. = FALSE)
    }
    sm$QTYPE <- factor(sm$QTYPE, levels = c("V","A","W"))
    sm[order(sm$QTYPE), ]
  }

  qtype_emm_df_stress <- function(fit) {
    if (is.null(fit)) return(NULL)
    emm <- emmeans::emmeans(fit, ~ QTYPE, type = "response")
    sm  <- as.data.frame(summary(emm, infer = c(TRUE, TRUE), type = "response"))
    if ("prob" %in% names(sm)) sm$pred <- sm$prob
    else if ("response" %in% names(sm)) sm$pred <- sm$response
    else if ("emmean" %in% names(sm)) sm$pred <- sm$emmean
    else stop("emmeans (stress) did not contain a probability column.", call. = FALSE)

    if ("lower.CL" %in% names(sm) && "upper.CL" %in% names(sm)) {
      sm$lo <- sm$lower.CL; sm$hi <- sm$upper.CL
    } else if ("asymp.LCL" %in% names(sm) && "asymp.UCL" %in% names(sm)) {
      sm$lo <- sm$asymp.LCL; sm$hi <- sm$asymp.UCL
    } else {
      stop("emmeans (stress) did not contain CI bounds.", call. = FALSE)
    }
    sm$QTYPE <- factor(sm$QTYPE, levels = c("V","A","W"))
    sm[order(sm$QTYPE), ]
  }

  qtype_pvals_vs_V <- function(fit, kind = c("cont","stress")) {
    kind <- match.arg(kind)
    emm <- if (kind == "stress") emmeans::emmeans(fit, ~ QTYPE, type = "response") else emmeans::emmeans(fit, ~ QTYPE)
    ctr <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "V")
    sm  <- as.data.frame(summary(ctr, infer = c(TRUE, TRUE), adjust = "none"))

    pA <- NA_real_; pW <- NA_real_
    if (nrow(sm) > 0) {
      if (any(grepl("A", sm$contrast) & grepl("V", sm$contrast))) pA <- sm$p.value[which(grepl("A", sm$contrast) & grepl("V", sm$contrast))[1]]
      if (any(grepl("W", sm$contrast) & grepl("V", sm$contrast))) pW <- sm$p.value[which(grepl("W", sm$contrast) & grepl("V", sm$contrast))[1]]
    }
    list(pA = pA, pW = pW)
  }

  ylims_cont   <- c(-0.5, 0.5)
  ylims_stress <- c(0.1, 0.7)
  YLABEL_STRESS <- expression(italic(P)*"("*italic(S)==plain(RHS)*")")

  panel_theme <- ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 11, face = "plain"),
      axis.text    = ggplot2::element_text(size = 11),
      plot.margin  = ggplot2::margin(0, 1, 0, 0),
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = ggplot2::element_blank()
    )

  make_qtype_cont_plot <- function(fit, ylab, show_x_text = TRUE) {
    if (is.null(fit)) return(ggplot2::ggplot() + ggplot2::theme_void())
    df <- qtype_emm_df_cont(fit)
    pv <- qtype_pvals_vs_V(fit, kind = "cont")

    df$sig <- "NS"
    df$sig[df$QTYPE == "A"] <- p_to_star(pv$pA)
    df$sig[df$QTYPE == "W"] <- p_to_star(pv$pW)
    df$sig <- factor(df$sig, levels = sig_levels)
    df$col_key <- ifelse(df$QTYPE == "V", "Vref", as.character(df$sig))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = QTYPE, y = pred, color = col_key)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), linewidth = 1.25) +
      ggplot2::geom_point(shape = 16, stroke = 0, size = 3.2) +
      ggplot2::scale_color_manual(
        values = c(sig_colors, Vref = VREF_COLOR),
        limits = c(sig_levels, "Vref"),
        breaks = sig_levels,
        drop = FALSE
      ) +
      ggplot2::coord_cartesian(ylim = ylims_cont) +
      ggplot2::labs(y = ylab) +
      panel_theme +
      ggplot2::scale_x_discrete(
        labels = function(x) parse(text = paste0("italic(", x, ")")),
        expand = ggplot2::expansion(mult = 0.18)
      ) +
      ggplot2::scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
      ggplot2::theme(axis.ticks.x = ggplot2::element_line(linewidth = 0.7))

    if (!show_x_text) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    p
  }

  make_qtype_stress_plot <- function(fit, show_x_text = TRUE) {
    if (is.null(fit)) return(ggplot2::ggplot() + ggplot2::theme_void())
    df <- qtype_emm_df_stress(fit)
    pv <- qtype_pvals_vs_V(fit, kind = "stress")

    df$sig <- "NS"
    df$sig[df$QTYPE == "A"] <- p_to_star(pv$pA)
    df$sig[df$QTYPE == "W"] <- p_to_star(pv$pW)
    df$sig <- factor(df$sig, levels = sig_levels)
    df$col_key <- ifelse(df$QTYPE == "V", "Vref", as.character(df$sig))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = QTYPE, y = pred, color = col_key)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), linewidth = 1.25) +
      ggplot2::geom_point(shape = 16, stroke = 0, size = 3.2) +
      ggplot2::scale_color_manual(
        values = c(sig_colors, Vref = VREF_COLOR),
        limits = c(sig_levels, "Vref"),
        breaks = sig_levels,
        drop = FALSE
      ) +
      ggplot2::coord_cartesian(ylim = ylims_stress) +
      ggplot2::labs(y = YLABEL_STRESS) +
      panel_theme +
      ggplot2::scale_x_discrete(
        labels = function(x) parse(text = paste0("italic(", x, ")")),
        expand = ggplot2::expansion(mult = 0.18)
      ) +
      ggplot2::theme(axis.ticks.x = ggplot2::element_line(linewidth = 0.7))

    if (!show_x_text) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    p
  }

  acol <- function(lab) {
    cowplot::ggdraw() +
      cowplot::draw_label(lab, fontface = "bold", size = 13, x = 0.98, y = 0.98, hjust = 1, vjust = 1) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  }

  # Right strip helper uses gtable + grid; keep calls explicit.
  right_strip_from_plot <- function(p, label_txt, fill = "grey90") {
    g <- ggplot2::ggplotGrob(p)
    idx <- which(g$layout$name == "panel")
    if (length(idx) == 0) {
      return(cowplot::ggdraw() +
               cowplot::draw_grob(grid::rectGrob(gp = grid::gpar(fill = fill, col = NA)), x = 0, y = 0, width = 1, height = 1) +
               cowplot::draw_label(label_txt, angle = -90, fontface = "bold", size = 12, x = 0.5, y = 0.5) +
               ggplot2::theme_void())
    }
    top <- min(g$layout$t[idx]); bot <- max(g$layout$b[idx])

    strip <- gtable::gtable(widths = grid::unit(1, "npc"), heights = g$heights)
    rect  <- grid::rectGrob(gp = grid::gpar(fill = fill, col = NA))
    strip <- gtable::gtable_add_grob(strip, rect, t = top, b = bot, l = 1, r = 1, clip = "on")

    cowplot::ggdraw() +
      cowplot::draw_grob(strip, x = 0, y = 0, width = 1, height = 1) +
      cowplot::draw_label(label_txt, angle = -90, fontface = "bold", size = 12, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5) +
      ggplot2::theme_void()
  }

  legend_items <- data.frame(sig = factor(sig_levels, levels = sig_levels), x_dot = c(1.6, 3.6, 5.5, 7.4), y = 0)
  legend_manual <- ggplot2::ggplot() +
    ggplot2::geom_point(data = legend_items, ggplot2::aes(x = x_dot, y = y, color = sig), size = 8.8) +
    ggplot2::geom_text(data = legend_items, ggplot2::aes(x = x_dot + 0.48, y = y, label = as.character(sig)),
                       hjust = 0, vjust = 0.35, fontface = "bold", size = 6.5) +
    ggplot2::scale_color_manual(values = sig_colors, limits = sig_levels, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0.5, 8.7), ylim = c(-0.6, 0.6), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))

  legend <- cowplot::as_grob(legend_manual)
  header <- cowplot::ggdraw() + cowplot::draw_label(expression(bolditalic(QTYPE)), size = 15) + ggplot2::theme_void()

  # a1-a4 (continuous)
  pa1 <- make_qtype_cont_plot(fits10$NFP$model_full,    ylab = expression(bar(italic(NFP))),            show_x_text = FALSE)
  pa2 <- make_qtype_cont_plot(fits10$NHR_AW$model_full, ylab = expression(bar(italic(NHR))[plain(AW)]), show_x_text = FALSE)
  pa3 <- make_qtype_cont_plot(fits10$NHR_E4$model_full, ylab = expression(bar(italic(NHR))[plain(E4)]), show_x_text = FALSE)
  pa4 <- make_qtype_cont_plot(fits10$NHRV$model_full,   ylab = expression(bar(italic(NHRV))),           show_x_text = TRUE)

  # b1-b4 (stress labels)
  pb1 <- make_qtype_stress_plot(fits11$NFP$model_full,    show_x_text = FALSE)
  pb2 <- make_qtype_stress_plot(fits11$NHR_AW$model_full, show_x_text = FALSE)
  pb3 <- make_qtype_stress_plot(fits11$NHR_E4$model_full, show_x_text = FALSE)
  pb4 <- make_qtype_stress_plot(fits11$NHRV$model_full,   show_x_text = TRUE)

  spacer <- cowplot::ggdraw() + cowplot::draw_grob(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))) + ggplot2::theme_void()
  ROW_GAP <- 0.10
  row_gap <- ggplot2::ggplot() + ggplot2::theme_void()
  REL_W <- c(0.06, 0.94)

  rowa1 <- cowplot::plot_grid(acol("a1"), pa1, ncol = 2, rel_widths = REL_W)
  rowa2 <- cowplot::plot_grid(acol("a2"), pa2, ncol = 2, rel_widths = REL_W)
  rowa3 <- cowplot::plot_grid(acol("a3"), pa3, ncol = 2, rel_widths = REL_W)
  rowa4 <- cowplot::plot_grid(acol("a4"), pa4, ncol = 2, rel_widths = REL_W)

  rowb1 <- cowplot::plot_grid(acol("b1"), pb1, ncol = 2, rel_widths = REL_W)
  rowb2 <- cowplot::plot_grid(acol("b2"), pb2, ncol = 2, rel_widths = REL_W)
  rowb3 <- cowplot::plot_grid(acol("b3"), pb3, ncol = 2, rel_widths = REL_W)
  rowb4 <- cowplot::plot_grid(acol("b4"), pb4, ncol = 2, rel_widths = REL_W)

  left_col <- cowplot::plot_grid(rowa1, row_gap, rowa2, row_gap, rowa3, row_gap, rowa4,
                                 ncol = 1, rel_heights = c(1, ROW_GAP, 1, ROW_GAP, 1, ROW_GAP, 1),
                                 align = "v", axis = "lr")
  right_col <- cowplot::plot_grid(rowb1, row_gap, rowb2, row_gap, rowb3, row_gap, rowb4,
                                  ncol = 1, rel_heights = c(1, ROW_GAP, 1, ROW_GAP, 1, ROW_GAP, 1),
                                  align = "v", axis = "lr")

  strip_col <- cowplot::plot_grid(
    right_strip_from_plot(pb1, "FP"),    row_gap,
    right_strip_from_plot(pb2, "HR AW"), row_gap,
    right_strip_from_plot(pb3, "HR E4"), row_gap,
    right_strip_from_plot(pb4, "HRV"),
    ncol = 1,
    rel_heights = c(1, ROW_GAP, 1, ROW_GAP, 1, ROW_GAP, 1),
    align = "v", axis = "lr"
  )

  middle_lr <- cowplot::plot_grid(left_col, spacer, right_col, strip_col, nrow = 1,
                                  rel_widths = c(1, 0.14, 1, 0.10), align = "h", axis = "tb")

  cowplot::plot_grid(header, middle_lr, legend, ncol = 1, rel_heights = c(0.10, 1, 0.14))
}

# ==============================================================================
# Main runner
# ==============================================================================
main <- function(N_TARGET = 50, AIC_EPS = 1e-8) {
  N_TARGET <- as.integer(N_TARGET)
  AIC_EPS  <- as.numeric(AIC_EPS)

  # Paper behavior switches (kept explicit here for visibility)
  FORCE_TIME_EQ10 <- TRUE   # Eq(10) keeps zlogQTime in all channels (protect term)
  FORCE_TIME_EQ11 <- FALSE  # Eq(11) time can drop (not used; retained for paper parity)
  EQ10_REML <- FALSE        # ML fits for Eq(10), matching paper

  # ---- Paths ----
  script_dir   <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
  paths        <- make_paths(project_root)
  ensure_dirs(paths)

  message("Script directory:  ", script_dir)
  message("Project root:      ", paths$root)
  message("Data directory:    ", paths$data_dir)
  message("Figures directory: ", paths$figures_dir)
  message("Reports directory: ", paths$reports_dir, "\n")

  # ---- Packages ----
  load_packages()

  # ---- Data load + standardization ----
  data_obj <- load_and_prepare_data(paths$data_dir, N_TARGET)
  D <- data_obj$D
  channels <- data_obj$channels

  # ---- Fit models ----
  message("\nFitting Eq(10) continuous models (FULL + AIC reduction)...")
  fits10 <- lapply(names(channels), function(nm) fit_eq10_one(D, channels[[nm]]$y, nm, EQ10_REML = EQ10_REML, AIC_EPS = AIC_EPS))
  names(fits10) <- names(channels)

  message("Fitting Eq(11) stress-label models (FULL + AIC reduction)...")
  fits11 <- lapply(names(channels), function(nm) fit_eq11_one(D, channels[[nm]]$s, nm, AIC_EPS = AIC_EPS))
  names(fits11) <- names(channels)

  # ---- Build Table blocks (FULL and AIC-reduced) ----
  tab3_full <- list(); tab3_red <- list()
  for (nm in names(fits10)) {
    m_full <- fits10[[nm]]$model_full
    m_red  <- fits10[[nm]]$model_reduced

    if (is.null(m_full)) {
      tab3_full[[nm]] <- list(fx=data.frame(), n=NA_real_, r2=c(R2_m=NA_real_, R2_c=NA_real_), sd_u=NA_real_, sd_v=NA_real_, sd_e=NA_real_)
    } else {
      fx <- extract_fix_lmm(m_full)
      r2 <- r2_lmm_nakagawa_manual(m_full, grp_pid="PID", grp_qid="QID")
      vc <- as.data.frame(VarCorr(m_full))
      sd_u <- sqrt(vc$vcov[vc$grp=="PID" & vc$var1=="(Intercept)"][1])
      var_v <- vc$vcov[vc$grp=="QID" & vc$var1=="(Intercept)"]
      sd_v <- ifelse(length(var_v)==0, 0, sqrt(var_v[1]))
      sd_e <- sigma(m_full)
      tab3_full[[nm]] <- list(fx=fx, n=nobs(m_full), r2=r2, sd_u=sd_u, sd_v=sd_v, sd_e=sd_e)
    }

    if (is.null(m_red)) {
      tab3_red[[nm]] <- list(fx=data.frame(), n=NA_real_, r2=c(R2_m=NA_real_, R2_c=NA_real_), sd_u=NA_real_, sd_v=NA_real_, sd_e=NA_real_)
    } else {
      fx <- extract_fix_lmm(m_red)
      r2 <- r2_lmm_nakagawa_manual(m_red, grp_pid="PID", grp_qid="QID")
      vc <- as.data.frame(VarCorr(m_red))
      sd_u <- sqrt(vc$vcov[vc$grp=="PID" & vc$var1=="(Intercept)"][1])
      var_v <- vc$vcov[vc$grp=="QID" & vc$var1=="(Intercept)"]
      sd_v <- ifelse(length(var_v)==0, 0, sqrt(var_v[1]))
      sd_e <- sigma(m_red)
      tab3_red[[nm]] <- list(fx=fx, n=nobs(m_red), r2=r2, sd_u=sd_u, sd_v=sd_v, sd_e=sd_e)
    }
  }

  tab4_full <- list(); tab4_red <- list()
  for (nm in names(fits11)) {
    m_full <- fits11[[nm]]$model_full
    m_red  <- fits11[[nm]]$model_reduced

    if (is.null(m_full)) {
      tab4_full[[nm]] <- list(fx=data.frame(), n=NA_real_, r2=c(R2_m=NA_real_, R2_c=NA_real_), sd_u=NA_real_, sd_v=NA_real_)
    } else {
      fx <- extract_fix_glmm(m_full)
      nobs_m <- tryCatch(nobs(m_full), error=function(e) NA_real_)
      r2 <- tryCatch(r2_glmm_latent(m_full), error=function(e) c(R2_m=NA_real_, R2_c=NA_real_))
      sd_u <- NA; sd_v <- NA
      if (inherits(m_full, "merMod")) {
        vc <- as.data.frame(VarCorr(m_full))
        if (any(vc$grp=="PID")) sd_u <- sqrt(vc$vcov[vc$grp=="PID" & vc$var1=="(Intercept)"][1])
        if (any(vc$grp=="QID")) sd_v <- sqrt(vc$vcov[vc$grp=="QID" & vc$var1=="(Intercept)"][1])
      }
      tab4_full[[nm]] <- list(fx=fx, n=nobs_m, r2=r2, sd_u=sd_u, sd_v=sd_v)
    }

    if (is.null(m_red)) {
      tab4_red[[nm]] <- list(fx=data.frame(), n=NA_real_, r2=c(R2_m=NA_real_, R2_c=NA_real_), sd_u=NA_real_, sd_v=NA_real_)
    } else {
      fx <- extract_fix_glmm(m_red)
      nobs_m <- tryCatch(nobs(m_red), error=function(e) NA_real_)
      r2 <- tryCatch(r2_glmm_latent(m_red), error=function(e) c(R2_m=NA_real_, R2_c=NA_real_))
      sd_u <- NA; sd_v <- NA
      if (inherits(m_red, "merMod")) {
        vc <- as.data.frame(VarCorr(m_red))
        if (any(vc$grp=="PID")) sd_u <- sqrt(vc$vcov[vc$grp=="PID" & vc$var1=="(Intercept)"][1])
        if (any(vc$grp=="QID")) sd_v <- sqrt(vc$vcov[vc$grp=="QID" & vc$var1=="(Intercept)"][1])
      }
      tab4_red[[nm]] <- list(fx=fx, n=nobs_m, r2=r2, sd_u=sd_u, sd_v=sd_v)
    }
  }

  # ---- AIC traces + kept-terms exports ----
  trace3 <- dplyr::bind_rows(lapply(names(fits10), function(nm) cbind(Outcome=nm, fits10[[nm]]$trace)))
  trace4 <- dplyr::bind_rows(lapply(names(fits11), function(nm) cbind(Outcome=nm, fits11[[nm]]$trace)))

  kept3 <- data.frame(
    Outcome   = names(fits10),
    FrameN    = sapply(fits10, `[[`, "frame_n"),
    AIC_full  = sapply(fits10, `[[`, "aic_full"),
    AIC_red   = sapply(fits10, `[[`, "aic_reduced"),
    DeltaAIC  = sapply(fits10, `[[`, "delta_aic"),
    KeptTerms = sapply(fits10, function(x) paste(x$kept_terms, collapse=" + ")),
    ErrFull   = sapply(fits10, function(x) ifelse(is.null(x$err_full), NA_character_, x$err_full)),
    ErrRed    = sapply(fits10, function(x) ifelse(is.null(x$err_reduced), NA_character_, x$err_reduced)),
    stringsAsFactors = FALSE
  )

  kept4 <- data.frame(
    Outcome   = names(fits11),
    FrameN    = sapply(fits11, `[[`, "frame_n"),
    AIC_full  = sapply(fits11, `[[`, "aic_full"),
    AIC_red   = sapply(fits11, `[[`, "aic_reduced"),
    DeltaAIC  = sapply(fits11, `[[`, "delta_aic"),
    KeptTerms = sapply(fits11, function(x) paste(x$kept_terms, collapse=" + ")),
    ErrFull   = sapply(fits11, function(x) ifelse(is.null(x$err_full), NA_character_, x$err_full)),
    ErrRed    = sapply(fits11, function(x) ifelse(is.null(x$err_reduced), NA_character_, x$err_reduced)),
    stringsAsFactors = FALSE
  )

  write.csv(trace3, file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_AICtrace_N%s.csv", N_TARGET)), row.names=FALSE)
  write.csv(trace4, file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_AICtrace_N%s.csv", N_TARGET)), row.names=FALSE)
  write.csv(kept3,  file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_KeptTerms_N%s.csv", N_TARGET)), row.names=FALSE)
  write.csv(kept4,  file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_KeptTerms_N%s.csv", N_TARGET)), row.names=FALSE)

  # ---- Wide CSV exports (FULL + AIC reduced) ----
  t3_csv <- file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_N%s.csv", N_TARGET))
  t4_csv <- file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_N%s.csv", N_TARGET))
  write.csv(make_wide_table3(tab3_full), t3_csv, row.names=FALSE)
  write.csv(make_wide_table4(tab4_full), t4_csv, row.names=FALSE)

  t3_csv_red <- file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_AICreduced_N%s.csv", N_TARGET))
  t4_csv_red <- file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_AICreduced_N%s.csv", N_TARGET))

  w3_red <- make_wide_table3(tab3_red)
  w4_red <- make_wide_table4(tab4_red)

  da_row3 <- data.frame(Predictor="DeltaAIC_full_minus_reduced", stringsAsFactors=FALSE)
  for (nm in names(tab3_red)) {
    da_row3[[paste0(nm,"_b")]]  <- fits10[[nm]]$delta_aic
    da_row3[[paste0(nm,"_SE")]] <- NA
    da_row3[[paste0(nm,"_t")]]  <- NA
    da_row3[[paste0(nm,"_p")]]  <- NA
  }
  w3_red <- dplyr::bind_rows(w3_red, da_row3)

  da_row4 <- data.frame(Predictor="DeltaAIC_full_minus_reduced", stringsAsFactors=FALSE)
  for (nm in names(tab4_red)) {
    da_row4[[paste0(nm,"_b")]]  <- fits11[[nm]]$delta_aic
    da_row4[[paste0(nm,"_SE")]] <- NA
    da_row4[[paste0(nm,"_z")]]  <- NA
    da_row4[[paste0(nm,"_p")]]  <- NA
  }
  w4_red <- dplyr::bind_rows(w4_red, da_row4)

  write.csv(w3_red, t3_csv_red, row.names=FALSE)
  write.csv(w4_red, t4_csv_red, row.names=FALSE)

  # ---- TeX exports ----
  t3_tex <- file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_N%s.tex", N_TARGET))
  t4_tex <- file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_N%s.tex", N_TARGET))
  write_tex_table3(t3_tex, tab3_full)
  write_tex_table4(t4_tex, tab4_full)

  t3_tex_red <- file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_AICreduced_N%s.tex", N_TARGET))
  t4_tex_red <- file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_AICreduced_N%s.tex", N_TARGET))
  delta3 <- lapply(names(fits10), function(nm) fits10[[nm]]$delta_aic); names(delta3) <- names(fits10)
  delta4 <- lapply(names(fits11), function(nm) fits11[[nm]]$delta_aic); names(delta4) <- names(fits11)
  write_tex_table3(t3_tex_red, tab3_red, delta_aic = delta3)
  write_tex_table4(t4_tex_red, tab4_red, delta_aic = delta4)

  # ---- AIC summary exports ----
  aic3_sum <- make_aic_summary(fits10, kind="eq10")
  aic4_sum <- make_aic_summary(fits11, kind="eq11")

  t3_aic_sum_csv <- file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_AICsummary_N%s.csv", N_TARGET))
  t4_aic_sum_csv <- file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_AICsummary_N%s.csv", N_TARGET))
  write.csv(aic3_sum, t3_aic_sum_csv, row.names=FALSE)
  write.csv(aic4_sum, t4_aic_sum_csv, row.names=FALSE)

  t3_aic_sum_tex <- file.path(paths$reports_dir, sprintf("Table3_Continuous_Models_AICsummary_N%s.tex", N_TARGET))
  t4_aic_sum_tex <- file.path(paths$reports_dir, sprintf("Table4_StressLabel_Models_AICsummary_N%s.tex", N_TARGET))
  write_tex_aic_summary(t3_aic_sum_tex, aic3_sum, stat_name="t")
  write_tex_aic_summary(t4_aic_sum_tex, aic4_sum, stat_name="z")

  # ---- Figure 10 ----
  fig10 <- make_figure10(fits10, fits11)
  saved10 <- save_plot_dual(fig10, "Figure10", paths$figures_dir, width = 7.16, height = 7.40, dpi = 300, limitsize = FALSE)

  message("\nWrote: ", t3_csv)
  message("Wrote: ", t3_tex)
  message("Wrote: ", t4_csv)
  message("Wrote: ", t4_tex)
  message("Wrote: ", saved10$pdf)
  message("Wrote: ", saved10$png)
  message("\nDONE.")
  invisible(list(paths = paths, fits10 = fits10, fits11 = fits11))
}

if (identical(environment(), globalenv())) {
  main()
}
