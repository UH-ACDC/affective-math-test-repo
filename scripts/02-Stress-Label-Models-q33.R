#!/usr/bin/env Rscript
# ==============================================================================
# Affective Math Test — Stress-Label GLMMs (q33 threshold)
#
# Author: Ioannis Pavlidis
# Affiliation: Affective and Data Computing Lab — University of Houston
#
# Repository script: scripts/02-Stress-Label-Models-q33_GitHub.R
#
# What this script does
#   • Fits stress-label mixed-effects logistic regression models (Eq. 11 form)
#     for four physiological channels, where stress labels are defined using a
#     lower/upper tercile threshold (q33; “high-stress” vs “low-stress”).
#   • Reports FULL models (paper reporting) and an AIC-reduced variant
#     (secondary; backward AIC with protected terms).
#   • Exports tables (CSV + LaTeX) to /reports.
#
# Outputs (GitHub naming convention; words in Initial Caps)
#   • reports/STable1_StressLabel_Models_q33_N50.{csv,tex}
#   • reports/STable1_StressLabel_Models_q33_AICreduced_N50.{csv,tex}
#   • reports/STable1_StressLabel_Models_q33_AICtrace_N50.csv
#   • reports/STable1_StressLabel_Models_q33_KeptTerms_N50.csv
#
# Notes
#   • Paths resolve relative to the script location (RStudio, Rscript, CI).
#   • No absolute paths.
#   • Do NOT clear the workspace inside `main()`; it can remove function
#     arguments in some source() workflows.
# ==============================================================================

options(stringsAsFactors = FALSE)

main <- function(N_TARGET = 50, AIC_EPS = 1e-8) {

  # ---------------------------------------------------------------------------
  # 0) Config (types + modeling behavior)
  # ---------------------------------------------------------------------------
  N_TARGET <- as.integer(N_TARGET)
  AIC_EPS  <- as.numeric(AIC_EPS)

  # Protected (always-kept) terms for AIC selection
  PROTECT_TERMS <- c("QTYPE", "zlogQTime")

  # ---------------------------------------------------------------------------
  # 1) Paths (GitHub-safe)
  # ---------------------------------------------------------------------------
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
    if (is.na(data_dir)) {
      stop("Cannot find processed data folder under {data,Data}/{processed,Processed}.", call. = FALSE)
    }

    list(
      root        = project_root,
      data_dir    = normalizePath(data_dir, winslash = "/", mustWork = FALSE),
      reports_dir = normalizePath(file.path(project_root, "reports"), winslash = "/", mustWork = FALSE)
    )
  }

  ensure_dirs <- function(paths) {
    dir.create(paths$reports_dir, recursive = TRUE, showWarnings = FALSE)
    invisible(TRUE)
  }

  script_dir   <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
  paths        <- make_paths(project_root)
  ensure_dirs(paths)

  data_dir    <- paths$data_dir
  reports_dir <- paths$reports_dir

  message("Script directory:  ", script_dir)
  message("Project root:      ", paths$root)
  message("Data directory:    ", data_dir)
  message("Reports directory: ", reports_dir, "\n")

  # ---------------------------------------------------------------------------
  # 2) Libraries
  # ---------------------------------------------------------------------------
  suppressPackageStartupMessages({
    library(dplyr)
    library(lme4)
  })

  # ---------------------------------------------------------------------------
  # 3) Utility helpers
  # ---------------------------------------------------------------------------
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

  # Convert a stress label to 0/1 (robust to factors/strings)
  coerce_stress01 <- function(x) {
    if (is.logical(x)) return(as.integer(x))
    if (is.numeric(x)) return(ifelse(is.na(x), NA_integer_, as.integer(x != 0)))
    xx <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(xx))
    out[xx %in% c("s","stress","stressed","1","true","t","rhs")] <- 1L
    out[xx %in% c("ns","no","nostress","0","false","f","rls")]   <- 0L
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

  # Greedy backward AIC (robust if some trial fits fail)
  greedy_backward_AIC <- function(fit_fun, start_terms, protect_terms = character(0)) {
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

  # Latent-scale R2 for logit GLMM (Nakagawa-style; simple and stable)
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

  extract_fix_glmm <- function(m) {
    if (is.null(m)) {
      return(data.frame(term=character(0), b=numeric(0), se=numeric(0), stat=numeric(0), p=numeric(0)))
    }
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

  # ---------------------------------------------------------------------------
  # 4) Data load + standardization
  # ---------------------------------------------------------------------------
  candidate_files <- c(
    file.path(data_dir, sprintf("Affective_Math_Qlevel_Data_N%s.csv", N_TARGET)),
    file.path(data_dir, "Affective_Math_Qlevel_Data_N50.csv")
  )
  in_file <- candidate_files[file.exists(candidate_files)][1]
  if (is.na(in_file) || !nzchar(in_file)) {
    stop("Could not find Q-level CSV in: ", data_dir,
         "\nExpected something like: data/processed/Affective_Math_Qlevel_Data_N", N_TARGET, ".csv",
         "\nIf you just cloned the repo, run the preprocessing script first to generate processed data.",
         call. = FALSE)
  }
  Q <- read.csv(in_file, stringsAsFactors = FALSE)
  message("Using input: ", in_file, "\n")

  # Required identifiers/covariates
  pid_col   <- pick_col(Q, c("ParticipantID","Participant","PID"))
  qid_col   <- pick_col(Q, c("Question.Name","QuestionName","QName","QID","QuestionID"))
  qtype_col <- pick_col(Q, c("Question.Type","QuestionType","QType","QTYPE"))
  qtime_col <- pick_col(Q, c("QTime","QTIME","TimeOnQuestion","Time"))

  grade_col <- pick_col(Q, c("Grade","QGRADE","QGrade"))
  sex_col   <- pick_col(Q, c("Gender","Sex","SEX"))
  sai_col   <- pick_col(Q, c("SAI","zSAI"))

  if (any(is.na(c(pid_col, qid_col, qtype_col, qtime_col)))) {
    stop("Missing required columns: PID/QID/QTYPE/QTIME", call. = FALSE)
  }

  # Stress-label columns (q33 variant)
  s_fp   <- pick_col(Q, c("Stress_NFP_q33","Stress.fp_q33","Stress_FP_q33","StressFp_q33"))
  s_hraw <- pick_col(Q, c("Stress_HR.AW_q33","Stress.hraw_q33","Stress_HRAW_q33","StressHraw_q33","Stress_NHR_AW_q33"))
  s_hre4 <- pick_col(Q, c("Stress_HR.E4_q33","Stress.hre4_q33","Stress_HRE4_q33","StressHre4_q33","Stress_NHR_E4_q33"))
  s_hrv  <- pick_col(Q, c("Stress_NHRV_q33","Stress.nhrv_q33","Stress_NHRV_q33","StressNhrv_q33"))

  if (any(is.na(c(s_fp, s_hraw, s_hre4, s_hrv)))) {
    stop("Missing one or more q33 stress-label columns. Found: ",
         "\n  FP:   ", s_fp,
         "\n  HRAW: ", s_hraw,
         "\n  HRE4: ", s_hre4,
         "\n  HRV:  ", s_hrv,
         call. = FALSE)
  }

  # Working frame with normalized covariates used by Eq(11)
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

  D$zSAI <- if (!is.na(sai_col)) zscore(suppressWarnings(as.numeric(D[[sai_col]]))) else NA_real_

  D$logQTime  <- ifelse(is.finite(D$QTime) & D$QTime > 0, log(D$QTime), NA_real_)
  D$zlogQTime <- as.numeric(scale(D$logQTime))

  channels <- list(
    NFP    = list(s = s_fp),
    NHR_AW = list(s = s_hraw),
    NHR_E4 = list(s = s_hre4),
    NHRV   = list(s = s_hrv)
  )

  # ---------------------------------------------------------------------------
  # 5) Model fitting (FULL + AIC-reduced)
  # ---------------------------------------------------------------------------
  fit_eq11_one <- function(s_col, outcome_name) {

    # Fixed-N (complete-case) frame for AIC comparisons:
    # full vs reduced must use identical rows
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

    # Degenerate checks (avoid glmer failures)
    if (nrow(Dm0) == 0) {
      warning("Eq11 ", outcome_name, ": empty frame after fixed-N filtering.")
      return(list(
        model_full=NULL, model_reduced=NULL, kept_terms=character(0),
        trace=data.frame(step=0, action="start", dropped="", AIC=Inf, terms="", stringsAsFactors=FALSE),
        frame_n=0, outcome=outcome_name,
        aic_full=NA_real_, aic_reduced=NA_real_, delta_aic=NA_real_,
        err_full="empty frame", err_reduced=NA_character_
      ))
    }
    if (length(unique(Dm0$S01)) < 2) {
      warning("Eq11 ", outcome_name, ": degenerate stress label (one class). n=", nrow(Dm0))
      return(list(
        model_full=NULL, model_reduced=NULL, kept_terms=character(0),
        trace=data.frame(step=0, action="start", dropped="", AIC=Inf, terms="", stringsAsFactors=FALSE),
        frame_n=nrow(Dm0), outcome=outcome_name,
        aic_full=NA_real_, aic_reduced=NA_real_, delta_aic=NA_real_,
        err_full="degenerate S01", err_reduced=NA_character_
      ))
    }

    base_terms <- PROTECT_TERMS
    opt_terms  <- c("QGRADE1", "SexF", "zSAI")
    full_terms <- c(base_terms, opt_terms)

    fit_fun <- function(terms_vec) {
      rhs <- paste(terms_vec, collapse = " + ")
      fml <- as.formula(paste0("S01 ~ ", rhs, " + (1|PID) + (1|QID)"))
      safe_fit(quote(
        glmer(
          fml, data = Dm0,
          family = binomial(link="logit"),
          nAGQ = 1,
          control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
        )
      ))
    }

    # FULL model (primary reporting)
    m_full  <- fit_fun(full_terms)
    err_full <- safe_fit_last_error
    aic_full <- safe_AIC(m_full); if (!is.finite(aic_full)) aic_full <- NA_real_

    # AIC reduction (secondary)
    sel    <- greedy_backward_AIC(fit_fun, start_terms = full_terms, protect_terms = base_terms)
    m_red  <- sel$model
    err_red <- if (is.null(m_red)) safe_fit_last_error else NA_character_
    aic_red <- safe_AIC(m_red); if (!is.finite(aic_red)) aic_red <- NA_real_
    delta_aic <- if (is.finite(aic_full) && is.finite(aic_red)) (aic_full - aic_red) else NA_real_

    # Fallback: if selection returns NULL, attempt base model (protected terms only)
    if (is.null(m_red)) {
      warning("Eq11 ", outcome_name, ": AIC selection produced NULL; attempting fallback base GLMM...")
      m_base <- fit_fun(base_terms)
      if (!is.null(m_base)) {
        m_red <- m_base
        aic_red <- safe_AIC(m_red); if (!is.finite(aic_red)) aic_red <- NA_real_
        delta_aic <- if (is.finite(aic_full) && is.finite(aic_red)) (aic_full - aic_red) else NA_real_
        sel$terms <- base_terms
      }
    }

    list(
      model_full    = m_full,
      model_reduced = m_red,
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

  cat("Fitting stress-label GLMMs (q33)...\n")
  fits <- lapply(names(channels), function(nm) fit_eq11_one(channels[[nm]]$s, nm))
  names(fits) <- names(channels)

  # ---------------------------------------------------------------------------
  # 6) Build Table blocks (FULL and AIC-reduced)
  # ---------------------------------------------------------------------------
  build_table_block <- function(m) {
    if (is.null(m)) {
      return(list(
        fx=data.frame(term=character(0), b=numeric(0), se=numeric(0), stat=numeric(0), p=numeric(0)),
        n=NA_real_, r2=c(R2_m=NA_real_, R2_c=NA_real_), sd_u=NA_real_, sd_v=NA_real_
      ))
    }
    fx <- extract_fix_glmm(m)
    nobs_m <- tryCatch(nobs(m), error=function(e) NA_real_)
    r2 <- tryCatch(r2_glmm_latent(m), error=function(e) c(R2_m=NA_real_, R2_c=NA_real_))
    sd_u <- NA_real_; sd_v <- NA_real_
    if (inherits(m, "merMod")) {
      vc <- as.data.frame(VarCorr(m))
      if (any(vc$grp=="PID")) sd_u <- sqrt(vc$vcov[vc$grp=="PID" & vc$var1=="(Intercept)"][1])
      if (any(vc$grp=="QID")) sd_v <- sqrt(vc$vcov[vc$grp=="QID" & vc$var1=="(Intercept)"][1])
    }
    list(fx=fx, n=nobs_m, r2=r2, sd_u=sd_u, sd_v=sd_v)
  }

  tab_full <- lapply(names(fits), function(nm) build_table_block(fits[[nm]]$model_full))
  names(tab_full) <- names(fits)

  tab_red <- lapply(names(fits), function(nm) build_table_block(fits[[nm]]$model_reduced))
  names(tab_red) <- names(fits)

  # ---------------------------------------------------------------------------
  # 7) Exports (CSV + LaTeX)
  # ---------------------------------------------------------------------------
  # ---- AIC traces + kept terms (debug/inspection) ----
  trace <- bind_rows(lapply(names(fits), function(nm) cbind(Outcome=nm, fits[[nm]]$trace)))
  kept <- data.frame(
    Outcome   = names(fits),
    FrameN    = sapply(fits, `[[`, "frame_n"),
    AIC_full  = sapply(fits, `[[`, "aic_full"),
    AIC_red   = sapply(fits, `[[`, "aic_reduced"),
    DeltaAIC  = sapply(fits, `[[`, "delta_aic"),
    KeptTerms = sapply(fits, function(x) paste(x$kept_terms, collapse=" + ")),
    ErrFull   = sapply(fits, function(x) ifelse(is.null(x$err_full), NA_character_, x$err_full)),
    ErrRed    = sapply(fits, function(x) ifelse(is.null(x$err_reduced), NA_character_, x$err_reduced)),
    stringsAsFactors=FALSE
  )

  trace_csv <- file.path(reports_dir, sprintf("STable1_StressLabel_Models_q33_AICtrace_N%s.csv", N_TARGET))
  kept_csv  <- file.path(reports_dir, sprintf("STable1_StressLabel_Models_q33_KeptTerms_N%s.csv", N_TARGET))
  write.csv(trace, trace_csv, row.names=FALSE)
  write.csv(kept,  kept_csv,  row.names=FALSE)

  # ---- Wide CSV tables (one block per channel; stable shape) ----
  make_wide_table <- function(tab) {
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
    footer <- data.frame(Predictor=c("Observations","R2_m","R2_c","sigma_u","sigma_v","sigma_eps"))
    for (nm in names(tab)) {
      footer[[paste0(nm,"_b")]] <- c(tab[[nm]]$n, tab[[nm]]$r2["R2_m"], tab[[nm]]$r2["R2_c"],
                                     tab[[nm]]$sd_u, tab[[nm]]$sd_v, pi/sqrt(3))
      footer[[paste0(nm,"_SE")]] <- NA
      footer[[paste0(nm,"_z")]]  <- NA
      footer[[paste0(nm,"_p")]]  <- NA
    }
    bind_rows(out, footer)
  }

  t_csv_full <- file.path(reports_dir, sprintf("STable1_StressLabel_Models_q33_N%s.csv", N_TARGET))
  t_csv_red  <- file.path(reports_dir, sprintf("STable1_StressLabel_Models_q33_AICreduced_N%s.csv", N_TARGET))
  write.csv(make_wide_table(tab_full), t_csv_full, row.names=FALSE)

  w_red <- make_wide_table(tab_red)
  da_row <- data.frame(Predictor="DeltaAIC_full_minus_reduced", stringsAsFactors=FALSE)
  for (nm in names(tab_red)) {
    da_row[[paste0(nm,"_b")]]  <- fits[[nm]]$delta_aic
    da_row[[paste0(nm,"_SE")]] <- NA
    da_row[[paste0(nm,"_z")]]  <- NA
    da_row[[paste0(nm,"_p")]]  <- NA
  }
  w_red <- bind_rows(w_red, da_row)
  write.csv(w_red, t_csv_red, row.names=FALSE)

  # ---- LaTeX writer (booktabs only; style in main LaTeX) ----
  write_tex_table <- function(path, tab, delta_aic = NULL) {
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
    cat("% --- AUTO-GENERATED ---\n", file=con)
    cat("\\begin{table*}[!t]\n", file=con)
    cat("\\centering\\scriptsize\n", file=con)
    cat("\\begin{tabular}{@{}l", file=con)
    for (k in sigs) cat(" r r r l", file=con)
    cat("@{}}\\toprule\n", file=con)

    # Header
    cat("& ", file=con)
    cat(paste(sprintf("\\multicolumn{4}{c}{%s}", sigs), collapse=" & "), file=con)
    cat(" \\\\\n", file=con)
    cat("\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-13}\\cmidrule(lr){14-17}\n", file=con)

    # Column labels
    cat("Predictor", file=con)
    for (k in sigs) cat(" & $b$ & SE & $z$ & $p$", file=con)
    cat(" \\\\\n", file=con)
    cat("\\midrule\n", file=con)

    # Fixed effects
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

    # Footer
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

    # ΔAIC row (optional; values under b-column)
    if (!is.null(delta_aic)) {
      cat("$\\Delta$AIC (full$-$reduced)", file=con)
      for (k in sigs) {
        da <- delta_aic[[k]]
        if (is.null(da) || is.na(da)) cat(" &  &  &  & ", file=con)
        else cat(sprintf(" & %.3f &  &  & ", da), file=con)
      }
      cat(" \\\\\n", file=con)
    }

    cat("\\bottomrule\n", file=con)
    cat("\\end{tabular}\n", file=con)
    cat("\\end{table*}\n", file=con)
  }

  t_tex_full <- file.path(reports_dir, sprintf("STable1_StressLabel_Models_q33_N%s.tex", N_TARGET))
  t_tex_red  <- file.path(reports_dir, sprintf("STable1_StressLabel_Models_q33_AICreduced_N%s.tex", N_TARGET))

  delta <- lapply(names(fits), function(nm) fits[[nm]]$delta_aic)
  names(delta) <- names(fits)

  write_tex_table(t_tex_full, tab_full)
  write_tex_table(t_tex_red,  tab_red, delta_aic = delta)

  # ---------------------------------------------------------------------------
  # 8) Console summary
  # ---------------------------------------------------------------------------
  cat("Wrote: ", t_csv_full, "\n", sep="")
  cat("Wrote: ", t_tex_full, "\n", sep="")
  cat("Wrote: ", t_csv_red,  "\n", sep="")
  cat("Wrote: ", t_tex_red,  "\n", sep="")
  cat("Wrote: ", trace_csv,  "\n", sep="")
  cat("Wrote: ", kept_csv,   "\n", sep="")
  cat("DONE.\n")
}

# Run when sourced/executed at top-level
if (identical(environment(), globalenv())) {
  main()
}
