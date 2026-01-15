#!/usr/bin/env Rscript
# ==============================================================================
# Affective Math Test — Valence + Performance Mixed-Effects Models
#
# Author: Ioannis Pavlidis
# Affiliation: Affective and Data Computing Lab — University of Houston
#
# Repository script: scripts/03-Valence-Performance-Models.R
#
# What this script does
#   • Fits valence mixed-effects models (Table 5; Figure 11, left column)
#   • Fits performance mixed-effects logistic models for correctness (Table 6; Figure 11, right column)
#   • Exports tables to /reports and Figure 11 to /figures (or /Figures fallback)
# Outputs (GitHub naming convention)
#   • reports/Table5_Valence_Models_N50.{csv,tex}
#   • reports/Table6_Performance_Model_N50.{csv,tex}
#   • reports/Table6_Performance_Model_AICtrace_N50.csv
#   • figures/Figure11.{pdf,png}
#
# Notes
#   • Paths resolve relative to the script location (RStudio, Rscript, CI).
#   • No absolute paths.
# ==============================================================================
options(stringsAsFactors = FALSE)

# -----------------------------
# Utility helpers (GitHub-safe paths)
# -----------------------------
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

ensure_dirs <- function(...) {
  dirs <- list(...)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

# Safer PDF device: prefer cairo if available, otherwise fall back to default pdf()
save_plot_dual <- function(p, filename_base, out_dir, width, height, dpi = 300, ...) {
  pdf_path <- file.path(out_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(out_dir, paste0(filename_base, ".png"))

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggplot2::ggsave(pdf_path, p, width = width, height = height, units = "in", device = pdf_device, bg = "white", ...)
  ggplot2::ggsave(png_path, p, width = width, height = height, units = "in", dpi = dpi, bg = "white", ...)
  invisible(list(pdf = pdf_path, png = png_path))
}

# ==============================================================================
# Main
# ==============================================================================
main <- function(N_TARGET = 50, AIC_EPS = 1e-8) {

  # Do NOT clear the workspace inside `main()`. That can remove function arguments
  # (e.g., `N_TARGET`) from the function environment in some user workflows.
  options(stringsAsFactors = FALSE)

  # =============================================================================
  # Module 0 — Config
  #   * N_TARGET: selects the processed dataset file (default N=50)
  #   * AIC_EPS : retained for a consistent CLI signature across scripts
  # =============================================================================
  N_TARGET <- as.integer(N_TARGET)
  AIC_EPS  <- as.numeric(AIC_EPS)

  # =============================================================================
  # Module 1 — Project paths (GitHub-safe)
  #   Resolves repo root relative to this script location.
  # =============================================================================
  script_dir <- get_script_dir()
  setwd(script_dir)
  repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

  DATA_DIR    <- file.path(repo_root, "data", "processed")
  REPORTS_DIR <- file.path(repo_root, "reports")
  FIG_DIR1    <- file.path(repo_root, "Figures")
  FIG_DIR2    <- file.path(repo_root, "figures")   # fallback in case repo uses lowercase
  FIGURES_DIR <- if (dir.exists(FIG_DIR1) || !dir.exists(FIG_DIR2)) FIG_DIR1 else FIG_DIR2

  dir.create(REPORTS_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("Repo root:\n  ", repo_root, "\n")
  cat("Input dir:\n  ", DATA_DIR, "\n")
  cat("Outputs:\n  Figures -> ", FIGURES_DIR, "\n  Tables  -> ", REPORTS_DIR, "\n\n", sep="")


  # ==============================================================================
  # Module 2 — Libraries
  #   - Load required packages (quietly)
  #   - Fail fast with clear install messages
  # ==============================================================================

  # -----------------------------
  # LIBS
  # -----------------------------
  pkgs <- c("dplyr","readr","stringr","forcats","ggplot2","cowplot",
            "lme4","lmerTest","ggeffects","performance","knitr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse=", "))
  suppressPackageStartupMessages({
    library(dplyr); library(readr); library(stringr); library(forcats)
    library(ggplot2); library(cowplot)
    library(lme4); library(lmerTest)
    library(ggeffects); library(performance); library(knitr)
  })

  # -----------------------------
  # INPUT
  # -----------------------------
  candidate_files <- c(
    file.path(DATA_DIR, sprintf("Affective_Math_Qlevel_Data_N%s.csv", N_TARGET)),
    file.path(DATA_DIR, sprintf("Affective_Math_Qlevel_Data_N%d.csv", N_TARGET)),
    file.path(DATA_DIR, "Affective_Math_Qlevel_Data_N50.csv"),
    file.path(DATA_DIR, "Affective_Math_Qlevel_Data.csv")
  )
  INPUT_CSV <- candidate_files[file.exists(candidate_files)][1]
  if (is.na(INPUT_CSV) || !nzchar(INPUT_CSV)) {
    stop("Could not find Qlevel CSV in ", DATA_DIR, "\nTried:\n", paste(candidate_files, collapse="\n"))
  }
  cat("Using input:\n  ", INPUT_CSV, "\n\n", sep="")
  Qlevel <- readr::read_csv(INPUT_CSV, show_col_types = FALSE)

  # -----------------------------
  # Helpers
  # -----------------------------
  coerce01 <- function(x) {
    if (is.logical(x)) return(as.integer(x))
    if (is.numeric(x)) return(as.integer(ifelse(is.na(x), NA, x != 0)))
    x <- as.character(x)
    x <- trimws(tolower(x))
    out <- rep(NA_integer_, length(x))
    out[x %in% c("1","true","t","yes","y","correct","right")] <- 1L
    out[x %in% c("0","false","f","no","n","incorrect","wrong")] <- 0L
    suppressWarnings({
      num <- as.numeric(x)
      out[is.na(out) & !is.na(num)] <- as.integer(num != 0)
    })
    out
  }

  zlog_time <- function(qt) {
    qt <- as.numeric(qt)
    if (any(qt <= 0, na.rm = TRUE)) return(as.numeric(scale(log1p(qt))))
    as.numeric(scale(log(qt)))
  }

  # p formatting (paper style)
  pfmt <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("<0.001")
    # fixed decimals, trim trailing zeros
    s <- sprintf("%.3f", p)
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
    s
  }

  # ==============================================================================
  # Module 3 — Helper functions (AIC, extraction, plotting helpers)
  #   - Formatting helpers (p-values, labels)
  #   - Safe model fitting + greedy backward AIC
  #   - Effect extraction + plotting helpers used later
  # ==============================================================================

  # significance bins (matches your Figure 11 legend)
  p_to_sig <- function(p) {
    if (is.na(p)) return("NS")
    if (p <= 0.001) return("***")
    if (p <= 0.01)  return("**")
    if (p <= 0.05)  return("*")
    "NS"
  }

  sig_levels <- c("REF","NS","*","**","***")
  sig_cols <- c("REF"="gray60", "NS"="black", "*"="#00D5D5", "**"="orange", "***"="red")

  # robust R2 for glmer(logit) if performance fails (Nakagawa latent-scale)
  r2_glmm_latent_robust <- function(m) {
    out <- c(R2_m = NA_real_, R2_c = NA_real_)
    if (is.null(m) || !inherits(m, "merMod")) return(out)
    r2 <- tryCatch(performance::r2_nakagawa(m), error=function(e) NULL)
    if (!is.null(r2)) {
      Rm <- suppressWarnings(as.numeric(r2$R2_marginal))
      Rc <- suppressWarnings(as.numeric(r2$R2_conditional))
      if (is.finite(Rm)) out["R2_m"] <- Rm
      if (is.finite(Rc)) out["R2_c"] <- Rc
      # if Rc missing but random var ~0, collapse; otherwise compute below
      if (is.finite(out["R2_m"]) && is.finite(out["R2_c"])) return(out)
    }
    # manual fallback
    mu <- tryCatch(predict(m, type="link"), error=function(e) NULL)
    if (is.null(mu)) return(out)
    varF <- stats::var(as.numeric(mu), na.rm=TRUE)
    vc <- tryCatch(lme4::VarCorr(m), error=function(e) NULL)
    varRE <- 0
    if (!is.null(vc)) {
      for (g in names(vc)) {
        v <- as.numeric(vc[[g]])
        v <- v[is.finite(v)]
        if (length(v) > 0) varRE <- varRE + sum(v)
      }
    }
    varE <- (pi^2)/3
    denom_m <- varF + varE
    denom_c <- varF + varRE + varE
    if (is.finite(varF) && denom_m > 0) out["R2_m"] <- varF / denom_m
    if (is.finite(varF) && denom_c > 0) out["R2_c"] <- (varF + varRE) / denom_c
    out
  }

  # alias used elsewhere in older scripts
  r2_latent_glmm <- r2_glmm_latent_robust

  # backward AIC for glmer fixed effects only; retain keep terms
  backward_AIC_glmer <- function(mod, keep = character(0)) {
    aic0  <- AIC(mod)
    steps <- data.frame(step = 0, dropped = NA_character_, AIC = aic0)
    repeat {
      tt <- attr(terms(mod), "term.labels")
      fixed_terms <- tt[!grepl("\\|", tt)]
      cand <- setdiff(fixed_terms, keep)
      if (length(cand) == 0) break
      aics <- sapply(cand, function(term) {
        m2 <- try(suppressWarnings(update(mod, as.formula(paste(". ~ . -", term)))), silent = TRUE)
        if (inherits(m2, "try-error")) return(Inf)
        AIC(m2)
      })
      best_term <- cand[which.min(aics)]
      best_aic  <- min(aics)
      if (!is.finite(best_aic) || best_aic >= (AIC(mod) - AIC_EPS)) break
      mod <- suppressWarnings(update(mod, as.formula(paste(". ~ . -", best_term))))
      steps <- rbind(steps, data.frame(step = nrow(steps), dropped = best_term, AIC = AIC(mod)))
    }
    list(model = mod, steps = steps)
  }

  # extract fixed effect table for lmer (t + p) or glmer (z + p)
  coef_table <- function(fit, model_name, stat = c("t","z")) {
    stat <- match.arg(stat)
    sm <- summary(fit)
    cm <- sm$coefficients
    if (stat == "t") {
      out <- data.frame(
        Model = model_name,
        term  = rownames(cm),
        b     = cm[, "Estimate"],
        SE    = cm[, "Std. Error"],
        stat  = cm[, "t value"],
        p     = cm[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    } else {
      out <- data.frame(
        Model = model_name,
        term  = rownames(cm),
        b     = cm[, "Estimate"],
        SE    = cm[, "Std. Error"],
        stat  = cm[, "z value"],
        p     = cm[, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
    names(out)[names(out)=="stat"] <- stat
    out
  }

  # ---------------------------------------------------------------------------
  # Table 6 row ordering (paper layout)
  #   In the performance tables, place the SEX row near the bottom of the table,
  #   immediately before zSAI.
  # ---------------------------------------------------------------------------
  reorder_perf_terms <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    ord <- c(
      "(Intercept)",
      "Question.TypeA", "Question.TypeW",
      "zlogQTime", "zNFP",
      "SexF", "zSAI"
    )
    df$..ord <- match(df$term, ord)
    df$..ord[is.na(df$..ord)] <- Inf
    df <- df[order(df$..ord, df$term), , drop = FALSE]
    df$..ord <- NULL
    df
  }

  # =============================================================================
  # Module 4 — Data load + standardization (Q-level)
  #   - Read processed Q-level file (N_TARGET)
  #   - Harmonize column names (PID/QID/QTYPE/QTIME/etc.)
  #   - Create analysis variables (zlogQTime, zSAI, factor levels)
  #   - Drop incomplete rows *within each model frame* to ensure valid AIC comparisons
  #   We construct two analysis frames from the processed question-level file:
  #   (i)  Valence frame (Table 5): HAPPY, AFRAID
  #   (ii) Performance frame (Table 6): correctness (QGRADE=1)
  #
  #   This is the single place where we:
  #     • select/rename columns
  #     • coerce data types
  #     • create standardized predictors (zlogQTime, zNFP, zSAI)
  #     • apply complete-case filtering (fixed-N per model family)
  # =============================================================================

  # -----------------------------
  # 4A) Valence analysis frame
  #   - Question.Type in {V,A,W}
  #   - logHappy = log(HAPPY) (positive values only)
  # -----------------------------
  Val <- Qlevel %>%
    dplyr::select(ParticipantID, Question.Name, Question.Type, Happy, Afraid) %>%
    mutate(
      Question.Type = factor(Question.Type, levels = c("V","A","W")),
      Happy = as.numeric(Happy),
      Afraid = as.numeric(Afraid),
      logHappy = ifelse(is.finite(Happy) & Happy > 0, log(Happy), NA_real_)
    ) %>%
    filter(is.finite(logHappy), is.finite(Afraid)) %>%
    droplevels()

  # =============================================================================
  # Module 5 — Model fits
  #   - FULL model = the paper model we report (Tables)
  #   - AIC-reduced model = secondary diagnostic (same fixed-N frame)
  #   - QTYPE + zlog(QTIME) are protected during AIC selection
  #   5A) Valence LMMs (Table 5)
  #       Eq(12)-style: outcome ~ QTYPE + (1|PID) + (1|QID)
  #       These are the FULL models reported in the paper.
  # =============================================================================

  m_happy <- lmer(logHappy ~ Question.Type + (1|ParticipantID) + (1|Question.Name), data = Val)
  m_afraid<- lmer(Afraid   ~ Question.Type + (1|ParticipantID) + (1|Question.Name), data = Val)

  # -----------------------------
  # 5B) Performance analysis frame
  #   - Grade01: correctness indicator (QGRADE=1)
  #   - zlogQTime: standardized log(time-on-question)
  #   - zNFP: standardized facial perspiration (FPNorm)
  #   - SexF: SEX=F indicator
  #   - zSAI: standardized state anxiety
  # -----------------------------
  Perf <- Qlevel %>%
    dplyr::select(ParticipantID, Question.Name, Question.Type, Gender, QTime, Grade, FPNorm, SAI) %>%
    mutate(
      Question.Type = factor(Question.Type, levels = c("V","A","W")),
      QTime  = as.numeric(QTime),
      Grade01 = coerce01(Grade),
      zlogQTime = zlog_time(QTime),
      zNFP = as.numeric(scale(as.numeric(FPNorm))),
      # Keep both raw and standardized SAI.
      # Models use zSAI, but Figure 11b2 shows the raw SAI scale on the TOP axis.
      SAI_raw = as.numeric(SAI),
      zSAI = as.numeric(scale(SAI_raw)),
      SexF   = as.integer(as.character(Gender) %in% c('F','Female','female','FEMALE')),
      Question.Name = as.factor(Question.Name),
      ParticipantID = as.factor(ParticipantID)
    ) %>%
    filter(is.finite(zlogQTime), is.finite(zNFP), is.finite(zSAI), !is.na(Grade01)) %>%
    droplevels()

  # Raw SAI mapping for the secondary (top) axis in panel b2.
  sai_mu <- mean(Perf$SAI_raw, na.rm = TRUE)
  sai_sd <- stats::sd(Perf$SAI_raw, na.rm = TRUE)

  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

  # =============================================================================
  # 5C) Performance GLMMs (Table 6)
  #     FULL model: includes SEX, QTYPE contrasts, zlogQTime, zNFP, zSAI
  #     AIC reduction: backward AIC starting from the full model, while *protecting*
  #     QTYPE and zlogQTime (they are always retained).
  # =============================================================================

  m_perf_full <- glmer(
    Grade01 ~ 1 + SexF + Question.Type + zlogQTime + zNFP + zSAI +
      (1|ParticipantID) + (1|Question.Name),
    data = Perf, family = binomial, control = ctrl
  )

  sel <- backward_AIC_glmer(m_perf_full, keep = c("Question.Type","zlogQTime"))
  m_perf_reduced <- sel$model
  aic_steps <- sel$steps
  # (m_perf_full is the FULL model; m_perf_reduced is the AIC-reduced model)
  m_perf <- m_perf_reduced  # backward-compat alias

  # -----------------------------
  # TABLE 5 (Valence): CSV + TeX (paper style)
  # -----------------------------
  tab5_h <- coef_table(m_happy,  "log(Pbar_HAPPY)", stat="t") %>% mutate(Outcome="log(Pbar_HAPPY)")
  tab5_a <- coef_table(m_afraid, "Pbar_AFRAID",     stat="t") %>% mutate(Outcome="Pbar_AFRAID")

  # Keep only intercept + QTYPE contrasts
  keep_terms_val <- c("(Intercept)","Question.TypeA","Question.TypeW")
  tab5_h <- tab5_h %>% filter(term %in% keep_terms_val)
  tab5_a <- tab5_a %>% filter(term %in% keep_terms_val)

  # R2 + variance components (SDs)
  r2_h <- performance::r2_nakagawa(m_happy)
  r2_a <- performance::r2_nakagawa(m_afraid)

  sd_h <- as.data.frame(VarCorr(m_happy)) %>%
    mutate(sd = sqrt(vcov)) %>% select(grp, sd)
  sd_a <- as.data.frame(VarCorr(m_afraid)) %>%
    mutate(sd = sqrt(vcov)) %>% select(grp, sd)

  get_sd <- function(sd_df, grp_name) {
    x <- sd_df$sd[sd_df$grp == grp_name]
    if (length(x)==0) return(NA_real_)
    x[1]
  }

  summ5 <- function(m, r2, sd_df) {
    data.frame(
      Observations = nobs(m),
      R2_m = as.numeric(r2$R2_marginal),
      R2_c = as.numeric(r2$R2_conditional),
      sigma_u = get_sd(sd_df, "ParticipantID"),
      sigma_v = get_sd(sd_df, "Question.Name"),
      sigma_e = sigma(m),
      stringsAsFactors = FALSE
    )
  }

  s5_h <- summ5(m_happy, r2_h, sd_h)
  s5_a <- summ5(m_afraid, r2_a, sd_a)

  # Write CSV in a simple long format + a compact "paper wide" export
  tab5_csv <- bind_rows(
    tab5_h %>% transmute(Outcome, Predictor=term, b, SE, t, p),
    tab5_a %>% transmute(Outcome, Predictor=term, b, SE, t, p)
  )

  # ==============================================================================
  # Module 6 — Export tables (CSV + TeX)
  #   - Table 5: Valence models
  #   - Table 6: Performance models
  #   - File names use InitialCaps words (paper convention)
  # ==============================================================================

  tab5_csv_path <- file.path(REPORTS_DIR, sprintf("Table5_Valence_Models_N%s.csv", N_TARGET))
  write.csv(tab5_csv, tab5_csv_path, row.names = FALSE)

  # LaTeX writer (Table 5 as two side-by-side blocks)
  write_table5_tex <- function(path) {

    fmt_row <- function(term, b, se, t, p) {
      # predictor label like paper
      pred <- term
      pred <- sub("^\\(Intercept\\)$", "(Intercept)", pred)
      pred <- sub("^Question\\.TypeA$", "QTYPE[A]", pred)
      pred <- sub("^Question\\.TypeW$", "QTYPE[W]", pred)
      paste0(pred, " & ",
             sprintf("%.3f", b), " & ",
             sprintf("%.3f", se), " & ",
             sprintf("%.3f", t), " & ",
             pfmt(p), " \\\\")
    }

    lines <- c()
    lines <- c(lines, "\\begin{table}[!t]")
    lines <- c(lines, "\\centering")
    lines <- c(lines, "\\caption{Valence mixed-effects models by question format (reference: Video). Outcomes are $\\log(\\overline{P}_{\\mathrm{HAPPY}})$ and $\\overline{P}_{\\mathrm{AFRAID}}$. Tables report fixed effects (estimate $b$, SE, $t$, $p$), the number of observations, Nakagawa--Schielzeth $R^2$ (marginal $R^2_m$ / conditional $R^2_c$), and variance components as standard deviations $\\hat\\sigma_u$, $\\hat\\sigma_v$, and $\\hat\\sigma_\\epsilon$.}")
    lines <- c(lines, "\\label{tab:Table5}")
    lines <- c(lines, "\\vspace{0.5em}")
    lines <- c(lines, "\\begin{tabular}{lrrrr \\hspace{1.2em} lrrrr}")
    lines <- c(lines, "\\toprule")
    lines <- c(lines, "\\multicolumn{5}{c}{$\\log(\\overline{P}_{\\mathrm{HAPPY}})$} & \\multicolumn{5}{c}{$\\overline{P}_{\\mathrm{AFRAID}}$}\\\\")
    lines <- c(lines, "\\cmidrule(lr){1-5} \\cmidrule(lr){6-10}")
    lines <- c(lines, "Predictor & $b$ & SE & $t$ & $p$ & Predictor & $b$ & SE & $t$ & $p$\\\\")
    lines <- c(lines, "\\midrule")

    # fixed rows
    for (k in keep_terms_val) {
      rh <- tab5_h[tab5_h$term==k,]
      ra <- tab5_a[tab5_a$term==k,]
      lines <- c(lines,
                 paste0(
                   fmt_row(rh$term, rh$b, rh$SE, rh$t, rh$p),
                   " & ",
                   sub(" \\\\","", fmt_row(ra$term, ra$b, ra$SE, ra$t, ra$p))
                 )
      )
    }

    # footer / summaries
    footer_pair <- function(lbl, v1, v2, fmt="%.3f") {
      if (lbl=="Observations") {
        a <- as.integer(v1); b <- as.integer(v2)
        return(sprintf("%s & \\multicolumn{4}{r}{%d} & %s & \\multicolumn{4}{r}{%d}\\\\", lbl, a, lbl, b))
      }
      sprintf("%s & \\multicolumn{4}{r}{%s} & %s & \\multicolumn{4}{r}{%s}\\\\",
              lbl, sprintf(fmt, v1), lbl, sprintf(fmt, v2))
    }

    lines <- c(lines, "\\midrule")
    lines <- c(lines, footer_pair("Observations", s5_h$Observations, s5_a$Observations, fmt="%.0f"))
    lines <- c(lines, footer_pair("$R^2_m$", s5_h$R2_m, s5_a$R2_m))
    lines <- c(lines, footer_pair("$R^2_c$", s5_h$R2_c, s5_a$R2_c))
    lines <- c(lines, footer_pair("$\\hat\\sigma_u$", s5_h$sigma_u, s5_a$sigma_u))
    lines <- c(lines, footer_pair("$\\hat\\sigma_v$", s5_h$sigma_v, s5_a$sigma_v))
    lines <- c(lines, footer_pair("$\\hat\\sigma_\\epsilon$", s5_h$sigma_e, s5_a$sigma_e))
    lines <- c(lines, "\\bottomrule")
    lines <- c(lines, "\\end{tabular}")
    lines <- c(lines, "\\end{table}")

    writeLines(lines, con = path)
  }

  tab5_tex_path <- file.path(REPORTS_DIR, sprintf("Table5_Valence_Models_N%s.tex", N_TARGET))
  write_table5_tex(tab5_tex_path)

  # -----------------------------
  # TABLE 6 (Performance): FULL + AIC-reduced (paper workflow)
  #   - Report FULL model coefficients in Table 6 outputs
  #   - Then run AIC reduction while ALWAYS retaining QTYPE and zlogQTime
  # -----------------------------

  # FULL table
  tab6_full <- coef_table(m_perf_full, "Performance: logit(P(QGRADE=1))", stat="z")

  # AIC-reduced table (for transparency / checks; not necessarily used in main paper)
  tab6_red  <- coef_table(m_perf_reduced, "Performance: logit(P(QGRADE=1)) [AIC-reduced]", stat="z")

  # Paper table order: place SEX just before zSAI
  tab6_full <- reorder_perf_terms(tab6_full)
  tab6_red  <- reorder_perf_terms(tab6_red)

  # AIC trace
  aic_trace_path <- file.path(REPORTS_DIR, sprintf("Table6_Performance_Model_AICtrace_N%s.csv", N_TARGET))
  write.csv(aic_steps, aic_trace_path, row.names = FALSE)

  # CSV writers
  tab6_full_csv <- tab6_full %>% transmute(Predictor=term, b, SE, z, p)
  tab6_full_csv_path <- file.path(REPORTS_DIR, sprintf("Table6_Performance_Model_N%s.csv", N_TARGET))
  write.csv(tab6_full_csv, tab6_full_csv_path, row.names = FALSE)
  tab6_csv_path <- tab6_full_csv_path  # alias used in end-of-script prints

  tab6_red_csv <- tab6_red %>% transmute(Predictor=term, b, SE, z, p)
  tab6_red_csv_path <- file.path(REPORTS_DIR, sprintf("Table6_Performance_Model_AICreduced_N%s.csv", N_TARGET))
  write.csv(tab6_red_csv, tab6_red_csv_path, row.names = FALSE)

  # Stats rows (use FULL model; matches the "report FULL" rule)
  sd6 <- VarCorr(m_perf_full)
  get_sd <- function(vc, grp) {
    if (is.null(vc[[grp]])) return(0)
    s <- attr(vc[[grp]], "stddev")
    if (length(s) == 0 || !is.finite(s[1])) return(0)
    as.numeric(s[1])
  }
  s6 <- data.frame(
    N       = nobs(m_perf_full),
    R2_m    = r2_latent_glmm(m_perf_full)["R2_m"],
    R2_c    = r2_latent_glmm(m_perf_full)["R2_c"],
    sigma_u = get_sd(sd6, "ParticipantID"),
    sigma_v = get_sd(sd6, "Question.Name"),
    sigma_e = pi/sqrt(3),
    stringsAsFactors = FALSE
  )

  # TeX writer (Table-4-like simplicity; no sci notation)
  write_table6_tex <- function(path, tab, caption, label) {

    pfmt <- function(p) {
      if (!is.finite(p)) return("")
      if (p < 0.001) return("$<0.001$")
      sprintf("%.3f", p)
    }

    fmt_pred <- function(term) {
      pred <- term
      pred <- sub("^\\(Intercept\\)$", "(Intercept)", pred)
      pred <- sub("^Question\\.TypeA$", "QTYPE[A]", pred)
      pred <- sub("^Question\\.TypeW$", "QTYPE[W]", pred)
      pred <- sub("^zlogQTime$", "$z\\log(QTIME)$", pred)
      pred <- sub("^zNFP$", "$zNFP$", pred)
      pred <- sub("^zSAI$", "$zSAI$", pred)
      pred <- sub("^SexF$", "Sex[F]", pred)
      pred
    }

    fmt_row <- function(term, b, se, z, p) {
      paste0(fmt_pred(term), " & ",
             sprintf("%.3f", b), " & ",
             sprintf("%.3f", se), " & ",
             sprintf("%.3f", z), " & ",
             pfmt(p), " \\\\")
    }

    lines <- c()
    lines <- c(lines, "\\begin{table}[!t]")
    lines <- c(lines, "\\centering")
    lines <- c(lines, paste0("\\caption{", caption,
                             " Random intercepts are included for participant ($\\hat\\sigma_u$) and question ($\\hat\\sigma_v$), with $\\sigma_\\epsilon=\\pi/\\sqrt{3}$.}"))
    lines <- c(lines, paste0("\\label{", label, "}"))
    lines <- c(lines, "\\vspace{0.5ex}")
    lines <- c(lines, "\\begin{tabular}{lrrrr}")
    lines <- c(lines, "\\toprule")
    lines <- c(lines, "Predictor & $b$ & $SE$ & $z$ & $p$ \\\\")
    lines <- c(lines, "\\midrule")

    for (i in seq_len(nrow(tab))) {
      lines <- c(lines, fmt_row(tab$term[i], tab$b[i], tab$SE[i], tab$z[i], tab$p[i]))
    }

    lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf("$N$ & \\multicolumn{4}{r}{%d}\\\\", as.integer(s6$N)))
    lines <- c(lines, sprintf("$R^2_m$ & \\multicolumn{4}{r}{%.3f}\\\\", s6$R2_m))
    lines <- c(lines, sprintf("$R^2_c$ & \\multicolumn{4}{r}{%.3f}\\\\", s6$R2_c))
    lines <- c(lines, sprintf("$\\hat\\sigma_u$ & \\multicolumn{4}{r}{%.3f}\\\\", s6$sigma_u))
    lines <- c(lines, sprintf("$\\hat\\sigma_v$ & \\multicolumn{4}{r}{%.3f}\\\\", s6$sigma_v))
    lines <- c(lines, "$\\sigma_\\epsilon$ & \\multicolumn{4}{r}{$\\pi/\\sqrt{3}$}\\\\")
    lines <- c(lines, "\\bottomrule")
    lines <- c(lines, "\\end{tabular}")
    lines <- c(lines, "\\end{table}")

    writeLines(lines, con = path)
  }

  # FULL TeX
  tab6_tex_path <- file.path(REPORTS_DIR, sprintf("Table6_Performance_Model_N%s.tex", N_TARGET))
  write_table6_tex(tab6_tex_path, tab6_full, "Performance mixed-effects logistic regression (FULL model).", "tab:Table6")

  # AIC-reduced TeX (separate label)
  tab6_red_tex_path <- file.path(REPORTS_DIR, sprintf("Table6_Performance_Model_AICreduced_N%s.tex", N_TARGET))
  write_table6_tex(tab6_red_tex_path, tab6_red, "Performance mixed-effects logistic regression (AIC-reduced; QTYPE and $z\\log(QTIME)$ retained).", "tab:Table6AICreduced")

  # -----------------------------
  # FIGURE 11

  # -----------------------------
  theme_set(theme_classic())
  theme_update(plot.title = element_text(hjust = 0.5))

  # -------------------------
  # Panel look (paper style): rectangular border + bold/italic tick labels
  # -------------------------
  panel_theme <- theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line    = element_blank(),  # prevents left/bottom double-thick edges
    axis.ticks   = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks.length = grid::unit(2.5, "mm"),
    axis.text.x  = element_text(size = 16, face = "bold.italic"),
    axis.text.y  = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold.italic"),
    legend.position = "none",
    plot.margin = margin(8, 10, 8, 10)
  )

  # Create QTYPE marginal-mean panel with per-level coloring (REF vs A/W sig)
  make_qtype_panel <- function(mod, outcome, y_limits, y_breaks, ylab_expr, title_expr) {

    dfp <- ggeffects::ggpredict(mod, terms = "Question.Type", bias_correction = TRUE) %>% as.data.frame()
    dfp$x <- factor(dfp$x, levels = c("V","A","W"))

    # determine per-level sig for contrasts vs V using model coef p-values
    coefmat <- summary(mod)$coefficients
    pA <- if ("Question.TypeA" %in% rownames(coefmat)) coefmat["Question.TypeA", ncol(coefmat)] else NA_real_
    pW <- if ("Question.TypeW" %in% rownames(coefmat)) coefmat["Question.TypeW", ncol(coefmat)] else NA_real_

    lev_sig <- c(V="REF", A=p_to_sig(pA), W=p_to_sig(pW))
    dfp$sig <- factor(lev_sig[as.character(dfp$x)], levels = sig_levels)

    ggplot(dfp, aes(x = x, y = predicted)) +
      geom_linerange(aes(ymin = conf.low, ymax = conf.high, color = sig), linewidth = 1.2) +
      geom_point(aes(color = sig), size = 3.5) +
      scale_color_manual(values = sig_cols, breaks = sig_levels, drop = FALSE) +
      scale_x_discrete(limits = c("V","A","W")) +
      scale_y_continuous(breaks = y_breaks) +
      coord_cartesian(ylim = y_limits) +
      labs(title = title_expr, x = "", y = ylab_expr) +
      theme_classic() +
      panel_theme
  }

  # a1, a2
  p_a1 <- make_qtype_panel(
    m_happy, "logHappy",
    y_limits = c(-5.2, -3.9),
    y_breaks = seq(-5.0, -4.0, by = 0.5),
    ylab_expr = expression(log(bar(italic(P))[HAPPY])),
    title_expr = expression(italic(QTYPE))
  )

  p_a2 <- make_qtype_panel(
    m_afraid, "Afraid",
    y_limits = c(0.155, 0.245),
    y_breaks = seq(0.16, 0.24, by = 0.02),
    ylab_expr = expression(bar(italic(P))[AFRAID]),
    title_expr = expression(italic(QTYPE))
  )

  # b1: performance by QTYPE
  p_b1 <- make_qtype_panel(
    m_perf, "Grade01",
    y_limits = c(0.3, 0.7),
    y_breaks = seq(0.3, 0.7, by = 0.1),
    ylab_expr = expression(italic(P)(italic(QGRADE) == 1)),
    title_expr = expression(italic(QTYPE))
  )

  # b2: SAI curve (if present)
  terms_perf <- attr(terms(m_perf), "term.labels")
  has_sai <- any(grepl("\\bzSAI\\b", terms_perf)) || any(grepl("\\bSAI\\b", terms_perf))

  blank_panel <- function(msg) {
    cowplot::ggdraw() + cowplot::draw_label(msg, fontface="italic", size=12) + theme_void()
  }

  if (has_sai) {
    coefmat <- summary(m_perf)$coefficients
    sai_term <- rownames(coefmat)[grepl("SAI", rownames(coefmat)) & !grepl("Intercept", rownames(coefmat))]
    sai_term <- if (length(sai_term)>0) sai_term[1] else NA_character_
    p_sai <- if (!is.na(sai_term)) coefmat[sai_term, "Pr(>|z|)"] else NA_real_
    sig_sai <- p_to_sig(p_sai)
    col_sai <- sig_cols[[sig_sai]]
    if (is.null(col_sai) || is.na(col_sai)) col_sai <- sig_cols[["NS"]]

    sai_pred_term <- if (!is.na(sai_term)) paste0(sai_term, " [all]") else NA_character_
    sai_df <- if (!is.na(sai_pred_term)) ggeffects::ggpredict(m_perf, terms = sai_pred_term, bias_correction = TRUE) %>% as.data.frame() else data.frame()

    if (nrow(sai_df) == 0) {
      p_b2 <- blank_panel("SAI not in model")
    } else {
      p_b2 <- ggplot(sai_df, aes(x = x, y = predicted)) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
        geom_line(color = col_sai, linewidth = 2) +
        geom_vline(xintercept = mean(Perf[[sai_term]], na.rm=TRUE), linetype="dashed", color="gray", linewidth=1) +
        scale_y_continuous(breaks = seq(0.3, 0.7, by = 0.1)) +
        coord_cartesian(ylim = c(0.3, 0.7)) +
        labs(title = expression(italic(SAI)), x = "", y = "") +
        theme_classic() +
        panel_theme +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")

      # Restore the *raw* SAI scale on the top x-axis (paper style).
      # The model is fit on zSAI, but we display raw SAI units as a secondary axis.
      if (identical(sai_term, "zSAI") && "SAI_raw" %in% names(Perf)) {
        sai_mu <- mean(Perf$SAI_raw, na.rm = TRUE)
        sai_sd <- stats::sd(Perf$SAI_raw, na.rm = TRUE)

        # Use stable breaks from the predicted x-range (zSAI scale)
        br_z <- pretty(range(sai_df$x, na.rm = TRUE), n = 5)
        p_b2 <- p_b2 +
          scale_x_continuous(
            breaks = br_z,
            sec.axis = sec_axis(
              trans  = ~ . * sai_sd + sai_mu,
              name   = NULL,
              labels = function(raw) sprintf("%.0f", raw)
            )
          ) +
          theme(
            axis.title.x.top = element_blank(),
            axis.text.x.top  = element_text(size = 16, face = "bold")
            ,axis.ticks.x.top = element_line(linewidth = 0.7)
          )
      }
    }
  } else {
    p_b2 <- blank_panel("SAI not retained")
  }

  # legend strip (REF/NS/*/**/***)
  legend_items <- data.frame(
    sig = factor(sig_levels, levels = sig_levels),
    lab = c("REF","NS","*","**","***"),
    x   = c(2.0, 3.5, 5.0, 6.0, 7.5)
  )

  legend_strip <- ggplot(legend_items, aes(x = x, y = 0, color = sig)) +
    geom_point(size = 5, show.legend = FALSE) +
    geom_text(
      aes(x = x + 0.35, label = lab),
      color = "black", fontface = "bold", size = 7,
      hjust = 0, vjust = 0.35, show.legend = FALSE
    ) +
    scale_color_manual(values = sig_cols, drop = FALSE) +
    coord_cartesian(xlim = c(-0.3, 9.6), clip = "off") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
    guides(color = "none")

  # layout with spacer columns
  spacer_plot <- cowplot::ggdraw() + theme_void()

  row_a <- cowplot::plot_grid(p_a1, spacer_plot, p_a2, nrow=1, rel_widths=c(1,0.12,1), align="hv")
  row_b <- cowplot::plot_grid(p_b1, spacer_plot, p_b2, nrow=1, rel_widths=c(1,0.12,1), align="hv")

  valence_labeled <- cowplot::ggdraw(row_a) +
    cowplot::draw_plot_label(label=c("a1","a2"),
                             x=c(0.01, 0.53), y=c(0.99, 0.99),
                             hjust=0, vjust=1, fontface="bold", size=16)

  perf_labeled <- cowplot::ggdraw(row_b) +
    cowplot::draw_plot_label(label=c("b1","b2"),
                             x=c(0.01, 0.53), y=c(0.99, 0.99),
                             hjust=0, vjust=1, fontface="bold", size=16)

  # ==============================================================================
  # Module 7 — Figure 11 (paper-style panels)
  #   - Built from FULL models so it matches Tables 5–6
  # ==============================================================================

  fig11 <- cowplot::plot_grid(
    valence_labeled,
    perf_labeled,
    legend_strip,
    nrow = 3,
    rel_heights = c(1, 1, 0.18)
  )

  saved11 <- save_plot_dual(
    p             = fig11,
    filename_base = "Figure11",
    out_dir       = FIGURES_DIR,
    width         = 7.2,
    height        = 5.6,
    dpi           = 300,
    limitsize     = FALSE
  )

  message("Wrote: ", saved11$pdf)
  message("Wrote: ", saved11$png)

  # -----------------------------
  # Also write AIC trace (helpful; not used in paper tables)
  # -----------------------------
  aic_path <- file.path(REPORTS_DIR, sprintf("Table6_Performance_Model_AICtrace_N%s.csv", N_TARGET))
  write.csv(aic_steps, aic_path, row.names = FALSE)

cat("DONE.\n")
cat("Figure:\n  ", saved11$pdf, "\n  ", saved11$png, "\n", sep = "")
  cat("Tables:\n  ", tab5_csv_path, "\n  ", tab5_tex_path, "\n  ", tab6_csv_path, "\n  ", tab6_tex_path, "\n", sep="")
  cat("AIC trace:\n  ", aic_path, "\n\n", sep="")

}

# ==============================================================================
# Module 8 — Script entry point
#   - Running via source() in RStudio or Rscript executes main()
# ==============================================================================

if (identical(environment(), globalenv())) {
  main()
}