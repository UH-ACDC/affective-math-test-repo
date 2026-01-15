#!/usr/bin/env Rscript
# ==============================================================================
# Affective Math Test — Video-story Calmness Scoring (Table 7)
#
# Author: Ioannis Pavlidis
# Affiliation: Affective and Data Computing Lab — University of Houston
#
# Repository script: scripts/04-VideoStory_Calmness_Scoring.R
#
# What this script does
#   • Fits four per-channel linear mixed models (LMMs) on standardized signals:
#       - NFP, NHR_AW, NHR_E4, NHRV
#   • Computes participant-specific QTYPE means (V/A/W) at zlog(QTIME)=0 using
#     fixed effects + participant random effects
#   • Defines per-channel "calmness" as the (signed) contrast between Video and
#     the mean of (Abstract, Word)
#   • Averages the four per-channel calmness scores into an overall score and
#     ranks participants
#   • Exports Table 7 (CSV + LaTeX) to /reports
#
# Output files
#   • reports/Table7_Overall_Video_Calmness_Ranking.csv
#   • reports/Table7_Overall_Video_Calmness_Ranking.tex
#
# IMPORTANT
#   This GitHub-ready script intentionally reproduces the *exact computation*
#   used in the original lab script (04-VideoStory_Calmness_Scoring.R):
#     - outcome z-scoring per channel
#     - model: Yz ~ QTYPE + zlogQTime + (1 + QTYPE | PID) + (1 | QID)
#     - calmness direction: arousal channels use -Delta; HRV uses +Delta
# ==============================================================================

options(stringsAsFactors = FALSE)

main <- function(N_TARGET = 50) {

  # ---------------------------------------------------------------------------
  # Module 0 — Project paths (GitHub-safe)
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

    # 3) fallback
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
      root        = normalizePath(project_root, winslash = "/", mustWork = FALSE),
      data_dir    = normalizePath(data_dir, winslash = "/", mustWork = FALSE),
      reports_dir = normalizePath(file.path(project_root, "reports"), winslash = "/", mustWork = FALSE)
    )
  }

  ensure_dirs <- function(...) {
    dirs <- list(...)
    for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    invisible(TRUE)
  }

  script_dir   <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
  paths        <- make_paths(project_root)
  ensure_dirs(paths$reports_dir)

  cat("Repo root:\n  ", paths$root, "\n", sep = "")
  cat("Input dir:\n  ", paths$data_dir, "\n", sep = "")
  cat("Outputs:\n  Tables -> ", paths$reports_dir, "\n\n", sep = "")

  # ---------------------------------------------------------------------------
  # Module 1 — Libraries + small helpers
  # ---------------------------------------------------------------------------
  suppressPackageStartupMessages({
    library(dplyr)
    library(lme4)
  })

  pick_col <- function(df, candidates) {
    for (nm in candidates) if (nm %in% names(df)) return(nm)
    NA_character_
  }

  zscore <- function(x) as.numeric(scale(x))

  # ---------------------------------------------------------------------------
  # Module 2 — Load + standardize data (Q-level)
  # ---------------------------------------------------------------------------
  candidate_files <- c(
    file.path(paths$data_dir, sprintf("Affective_Math_Qlevel_Data_N%s.csv", as.integer(N_TARGET))),
    file.path(paths$data_dir, "Affective_Math_Qlevel_Data_N50.csv")
  )
  in_file <- candidate_files[file.exists(candidate_files)][1]
  if (is.na(in_file) || !nzchar(in_file)) {
    stop(
      "Could not find Q-level CSV in: ", paths$data_dir,
      "\nExpected something like: data/processed/Affective_Math_Qlevel_Data_N", as.integer(N_TARGET), ".csv",
      call. = FALSE
    )
  }

  cat("Using input:\n  ", in_file, "\n\n", sep = "")
  Q <- read.csv(in_file, stringsAsFactors = FALSE)

  # Required IDs
  pid_col   <- pick_col(Q, c("ParticipantID", "Participant", "PID"))
  qid_col   <- pick_col(Q, c("Question.Name", "QuestionName", "QName", "QID", "QuestionID"))
  qtype_col <- pick_col(Q, c("Question.Type", "QuestionType", "QType", "QTYPE"))
  qtime_col <- pick_col(Q, c("QTime", "QTIME", "TimeOnQuestion", "Time"))

  if (any(is.na(c(pid_col, qid_col, qtype_col, qtime_col)))) {
    stop("Missing required columns: PID/QID/QTYPE/QTIME", call. = FALSE)
  }

  # Channel columns (normalized continuous signals)
  y_fp   <- pick_col(Q, c("FPNorm", "PPNorm"))
  y_hraw <- pick_col(Q, c("HR.AWNorm", "HRAWNorm"))
  y_hre4 <- pick_col(Q, c("HR.E4Norm", "HRE4Norm"))
  y_hrv  <- pick_col(Q, c("HRVNorm", "NHRVNorm", "NHRV"))

  if (any(is.na(c(y_fp, y_hraw, y_hre4, y_hrv)))) {
    stop("Missing one or more continuous channel columns (FP/HR_AW/HR_E4/HRV).", call. = FALSE)
  }

  D <- Q %>%
    mutate(
      PID   = factor(.data[[pid_col]]),
      QID   = factor(.data[[qid_col]]),
      QTYPE = factor(.data[[qtype_col]], levels = c("V", "A", "W")),
      QTime = suppressWarnings(as.numeric(.data[[qtime_col]]))
    )

  D$logQTime  <- ifelse(is.finite(D$QTime) & D$QTime > 0, log(D$QTime), NA_real_)
  D$zlogQTime <- as.numeric(scale(D$logQTime))

  # ---------------------------------------------------------------------------
  # Module 3 — Per-channel LMMs (exactly as original)
  #   Model: Yz ~ QTYPE + zlogQTime + (1 + QTYPE | PID) + (1 | QID)
  # ---------------------------------------------------------------------------
  fit_channel_model <- function(df, y_col, channel_name) {
    dm <- df %>%
      mutate(Yraw = suppressWarnings(as.numeric(.data[[y_col]]))) %>%
      filter(
        is.finite(Yraw),
        !is.na(PID), !is.na(QID), !is.na(QTYPE),
        is.finite(zlogQTime)
      ) %>%
      mutate(Yz = zscore(Yraw))

    if (nrow(dm) < 10) {
      warning("Too few rows for ", channel_name, " after filtering (n=", nrow(dm), ").")
    }

    fml <- Yz ~ QTYPE + zlogQTime + (1 + QTYPE | PID) + (1 | QID)

    m <- lme4::lmer(
      fml,
      data = dm,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )

    list(model = m, data = dm)
  }

  cat("Fitting channel models and computing calmness scores...\n")
  m_fp   <- fit_channel_model(D, y_fp,   "NFP")
  m_hraw <- fit_channel_model(D, y_hraw, "NHR_AW")
  m_hre4 <- fit_channel_model(D, y_hre4, "NHR_E4")
  m_hrv  <- fit_channel_model(D, y_hrv,  "NHRV")

  # ---------------------------------------------------------------------------
  # Module 4 — Participant-specific QTYPE means (at zlogQTime = 0)
  # ---------------------------------------------------------------------------
  participant_format_means <- function(m, participant_var = "PID") {
    fix <- lme4::fixef(m)
    re  <- lme4::ranef(m)[[participant_var]]

    # Random-effect columns are typically: (Intercept), QTYPEA, QTYPEW
    colA <- if ("QTYPEA" %in% colnames(re)) "QTYPEA" else NA_character_
    colW <- if ("QTYPEW" %in% colnames(re)) "QTYPEW" else NA_character_

    ids <- rownames(re)

    out <- lapply(ids, function(pid) {
      b0 <- fix["(Intercept)"] + re[pid, "(Intercept)"]
      bA <- if (!is.na(colA)) (fix["QTYPEA"] + re[pid, colA]) else fix["QTYPEA"]
      bW <- if (!is.na(colW)) (fix["QTYPEW"] + re[pid, colW]) else fix["QTYPEW"]

      muV <- b0
      muA <- b0 + bA
      muW <- b0 + bW

      data.frame(PID = pid, meanV = muV, meanA = muA, meanW = muW, stringsAsFactors = FALSE)
    })

    dplyr::bind_rows(out)
  }

  # ---------------------------------------------------------------------------
  # Module 5 — Calmness scoring + ranking
  #   Delta = mu_V - mean(mu_A, mu_W)
  #   Calmness = -Delta for arousal channels; +Delta for HRV
  # ---------------------------------------------------------------------------
  calmness_from_means <- function(df_means, direction = c("arousal", "hrv")) {
    direction <- match.arg(direction)

    df_means %>%
      mutate(
        Delta = meanV - (meanA + meanW) / 2,
        Calmness = if (direction == "arousal") -Delta else Delta
      ) %>%
      select(PID, Calmness)
  }

  means_fp   <- participant_format_means(m_fp$model)
  means_hraw <- participant_format_means(m_hraw$model)
  means_hre4 <- participant_format_means(m_hre4$model)
  means_hrv  <- participant_format_means(m_hrv$model)

  calm_fp   <- calmness_from_means(means_fp,   direction = "arousal")
  calm_hraw <- calmness_from_means(means_hraw, direction = "arousal")
  calm_hre4 <- calmness_from_means(means_hre4, direction = "arousal")
  calm_hrv  <- calmness_from_means(means_hrv,  direction = "hrv")

  overall <- calm_fp %>%
    rename(Calm_FP = Calmness) %>%
    inner_join(rename(calm_hraw, Calm_HR_AW = Calmness), by = "PID") %>%
    inner_join(rename(calm_hre4, Calm_HR_E4 = Calmness), by = "PID") %>%
    inner_join(rename(calm_hrv,  Calm_HRV   = Calmness), by = "PID") %>%
    mutate(
      CalmnessScore = rowMeans(across(starts_with("Calm_")), na.rm = FALSE)
    ) %>%
    arrange(desc(CalmnessScore)) %>%
    mutate(Rank = row_number()) %>%
    select(Participant = PID, CalmnessScore, Rank)

  # ---------------------------------------------------------------------------
  # Module 6 — Export (CSV + LaTeX)
  # ---------------------------------------------------------------------------
  out_csv <- file.path(paths$reports_dir, "Table7_Overall_Video_Calmness_Ranking_N50.csv")
  out_tex <- file.path(paths$reports_dir, "Table7_Overall_Video_Calmness_Ranking_N50.tex")

  write.csv(overall, out_csv, row.names = FALSE)

  latex_table_5cols <- function(df_rank, tex_file, caption, label, digits = 3) {
    stopifnot(all(c("Participant", "CalmnessScore") %in% names(df_rank)))

    fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)

    # Make 10x5 (50 entries) by filling down columns (as in the paper)
    n <- nrow(df_rank)
    ncols <- 5
    nrows <- ceiling(n / ncols)

    # Build a 10x5 matrix of "Sxx (score)" strings
    cells <- rep("", nrows * ncols)
    for (i in seq_len(n)) {
      r <- ((i - 1) %% nrows) + 1
      c <- ((i - 1) %/% nrows) + 1
      cells[(c - 1) * nrows + r] <- paste0(df_rank$Participant[i], " (", fmt(df_rank$CalmnessScore[i]), ")")
    }
    mat <- matrix(cells, nrow = nrows, ncol = ncols, byrow = FALSE)

    con <- file(tex_file, open = "wt")
    on.exit(close(con), add = TRUE)

    cat("% Auto-generated from 04-VideoStory_Calmness_Scoring_GitHub.R\n", file = con)
    cat("\\begin{table}[!thb]\n", file = con)
    cat("  \\caption{\\textbf{", caption, "}}\n", sep = "", file = con)
    cat("  \\label{", label, "}\n", sep = "", file = con)
    cat("  \\centering\n", file = con)
    cat("  \\footnotesize\n", file = con)
    cat("  \\setlength{\\tabcolsep}{8pt}\n", file = con)
    cat("  \\renewcommand{\\arraystretch}{1.2}\n", file = con)
    cat("  \\begin{tabular}{lllll}\n", file = con)
    cat("    \\toprule\n", file = con)

    for (r in seq_len(nrows)) {
      row_cells <- mat[r, ]
      cat("    ", paste(row_cells, collapse = " & "), " \\\\\n", sep = "", file = con)
    }

    cat("    \\bottomrule\n", file = con)
    cat("  \\end{tabular}\n", file = con)
    cat("\\end{table}\n", file = con)
  }

  # Helper for a safe basename in the LaTeX comment (no need for rlang)
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

  latex_table_5cols(
    df_rank  = overall,
    tex_file = out_tex,
    caption  = "Overall calmness score ranking for all 50 participants.",
    label    = "tab:Table7",
    digits   = 3
  )

  cat("Wrote: ", out_csv, "\n", sep = "")
  cat("Wrote: ", out_tex, "\n\n", sep = "")
  cat("DONE.\n")
}

if (identical(environment(), globalenv())) {
  main()
}
