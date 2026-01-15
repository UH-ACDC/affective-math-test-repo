# ==============================================================================
# Affective Math Test — Exploratory Data Analysis, Preprocessing & Figure Generation
#
# Author: Ioannis Pavlidis
# Affiliation: Affective and Data Computing Lab — University of Houston
#
#
# Repository script: scripts/00-Exploratory-Data-Analysis.R
#
# What this script does
#   1) Loads the processed baseline + test datasets (CSV) from /data/processed
#   2) Performs exploratory summaries and transformations used in the paper
#   3) Generates the paper figures (Figures 2b–8) in /figures
#   4) Exports intermediate tables in /reports
#
# Notes on terminology
#   • The dataset filenames use “_Test.csv”. Some legacy variable names still use “Exam”
#     (e.g., Df_Exam) from early prototypes. We keep those variable names to minimize
#     diff-noise, but user-facing labels and comments use “test” throughout.
#
# Reproducibility
#   • The script resolves paths relative to the script location, so it can be run from
#     anywhere (RStudio, Rscript, CI).
#   • No setwd() calls are required.
# ==============================================================================

# ---- Packages (keep lean; attach only what is used) ----
suppressPackageStartupMessages({
  library(tidyverse)     # dplyr/tidyr/readr/ggplot2/etc.
  library(cowplot)       # multi-panel layouts + shared legends
  library(ggrastr)       # rasterize dense points in vector outputs
  library(latex2exp)     # TeX() in labels
  library(scales)        # label formatting helpers
  library(rstatix)       # summary stats helpers used in several modules
  library(bestNormalize) # transformation helpers used in EDA (e.g., normalizing signals)
})

# ---- Utility helpers ----
get_script_dir <- function() {
  # Works for: RStudio, Rscript, and most CI runners.
  # 1) Rscript --file=... provides a --file= argument in commandArgs()
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")) )
  }

  # 2) RStudio: rstudioapi if available
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NA_character_)
    if (!is.na(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }

  # 3) Fallback: current working directory
  return(getwd())
}

assert_exists <- function(path, what = "file") {
  if (!file.exists(path)) stop(sprintf("Missing %s: %s", what, path), call. = FALSE)
  invisible(TRUE)
}

make_paths <- function(project_root) {
  list(
    root        = project_root,
    scripts_dir = file.path(project_root, "scripts"),
    data_dir    = file.path(project_root, "data", "processed"),
    figures_dir = file.path(project_root, "figures"),
    reports_dir = file.path(project_root, "reports")
  )
}

ensure_dirs <- function(paths) {
  dir.create(paths$figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$reports_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

# Safer PDF device: prefer cairo if available, otherwise fall back to default pdf()
save_plot_dual <- function(p, filename_base, out_dir, width, height, dpi = 300, ...) {
  pdf_path <- file.path(out_dir, paste0(filename_base, ".pdf"))
  png_path <- file.path(out_dir, paste0(filename_base, ".png"))

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

  ggsave(pdf_path, p, width = width, height = height, units = "in", device = pdf_device, bg = "white", ...)
  ggsave(png_path, p, width = width, height = height, units = "in", dpi = dpi, bg = "white", ...)
  invisible(list(pdf = pdf_path, png = png_path))
}

# ---- Main entrypoint ----
main <- function(N = 50, is_save = TRUE) {
  script_dir   <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
  paths        <- make_paths(project_root)

  message("Script directory:  ", script_dir)
  message("Project root:      ", paths$root)
  message("Data directory:    ", paths$data_dir)
  message("Figures directory: ", paths$figures_dir)
  message("Reports directory: ", paths$reports_dir)

  # Create output folders if missing
  ensure_dirs(paths)

  # ==============================================================================
  # MODULE 1: Configuration, paths, themes, and figure controls
  # ==============================================================================

  # Convenience aliases expected by downstream (legacy) code blocks
  data_dir    <- paste0(paths$data_dir, .Platform$file.sep)
  plot_dir    <- paste0(paths$figures_dir, .Platform$file.sep)
  reports_dir <- paste0(paths$reports_dir, .Platform$file.sep)

  # Script-level config used throughout
  is_save <- is_save
  N       <- N

  # External plotting helpers (defines my_theme, etc.)
  assert_exists(file.path(script_dir, "fm_plot_func.R"), what = "helper script")
  source(file.path(script_dir, "fm_plot_func.R"))

  theme_set(theme_classic())
  theme_update(plot.title = element_text(hjust = 0.5))

  my_theme_bold.italic <- my_theme + theme(text = element_text(size = 14, face = "bold"))


  # -----------------------------
  # Figure 2b (SUS) toggle
  # -----------------------------
  do_fig2b <- TRUE
  # ==============================================================================
  # MODULE 2: Load baseline and test datasets
  # ==============================================================================


  file_path_BL   <- sprintf("%sAffective_Math_Dataset_N%s_BL.csv",   data_dir, N)
  file_path_Exam <- sprintf("%sAffective_Math_Dataset_N%s_Test.csv", data_dir, N)

  # NOTE: The input file uses the term *Test* (Affective_Math_Dataset_*_Test.csv).
  # Legacy variable names from the original pipeline use 'Exam' (Df_Exam) — we keep them
  # for backward compatibility, but the paper refers to the *test* throughout.

  Df_BL   <- read.csv(file_path_BL,   header = TRUE, sep = ",")
  Df_Exam <- read.csv(file_path_Exam, header = TRUE, sep = ",", stringsAsFactors = TRUE)

  # Friendly aliases (optional; improves readability in this script)
  Df_Baseline <- Df_BL
  Df_TestRaw  <- Df_Exam


  names(Df_Exam)
  colSums(!is.na(Df_Exam))


  Df_BL.SubSet <- subset(Df_BL,
                         select = c(ParticipantID, Time, Timestamp, Perspiration, HR.E4, HR.AW, HRV.IBI)
  )

  idx <- is.na(Df_BL.SubSet$HR.E4) | is.na(Df_BL.SubSet$HR.AW)
  Df_BL.SubSet$HR.E4[idx] <- NA
  Df_BL.SubSet$HR.AW[idx] <- NA

  Df_BL.mean <- Df_BL.SubSet %>%
    group_by(ParticipantID) %>%
    summarise(
      Perspiration = round(mean(Perspiration, na.rm = TRUE), 4),
      HR.E4 = mean(HR.E4, na.rm = TRUE),
      HR.AW = mean(HR.AW, na.rm = TRUE),
      HRV.IBI = mean(HRV.IBI, na.rm = TRUE),
      .groups = "drop"
    )

  # ==============================================================================
  # MODULE 3: Figure 2b — SUS score distribution
  # ==============================================================================
  if (isTRUE(do_fig2b)) {

    # Participant-level SUS (SUS.Score repeats across rows; take first per participant)
    df_sus <- Df_Exam %>%
      dplyr::filter(!is.na(SUS.Score)) %>%
      dplyr::group_by(ParticipantID) %>%
      dplyr::summarise(
        SUS = as.numeric(as.character(first(SUS.Score))),
        .groups = "drop"
      )

    # Define the function
    make_fig2b_plot <- function(df_sus) {
      stopifnot("SUS" %in% names(df_sus))

      n   <- sum(is.finite(df_sus$SUS))
      mu  <- mean(df_sus$SUS, na.rm = TRUE)
      sdv <- sd(df_sus$SUS, na.rm = TRUE)

      # Use plotmath via atop() to get true multi-line math text (no parse errors)
      ann_txt <- sprintf(
        "atop(italic(N)==%d, atop(Mean==%.2f, SD==%.2f))",
        n, mu, sdv
      )

      ggplot(df_sus, aes(x = 1, y = SUS)) +
        geom_boxplot(width = 0.28, outlier.shape = NA) +

        # lighter, less obtrusive points
        geom_jitter(
          width = 0.06, height = 0,
          size  = 1.0,
          alpha = 0.35,
          color = "gray55"
        )+

        # threshold + mean marker (contrast)
        geom_hline(yintercept = 68, linetype = "dashed", linewidth = 0.8, color = "red") +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 3.2, color = "red") +

        # IMPORTANT: ggpp masks annotate(); call ggplot2::annotate explicitly
        ggplot2::annotate(
          "label",
          x = 1.62, y = 8,
          label = ann_txt,
          parse = TRUE,
          hjust = 1, vjust = 0,
          size = 4.2,          # smaller text
          label.size = 0.2     # thinner border
        ) +

        coord_cartesian(xlim = c(0.65, 1.75), ylim = c(0, 100), clip = "off") +
        scale_x_continuous(breaks = NULL) +
        labs(x = NULL, y = "SUS score") +
        theme_classic() +
        theme(
          axis.ticks.x = element_blank(),
          axis.text.x  = element_blank(),
          plot.margin  = margin(5.5, 18, 5.5, 5.5),
          text         = element_text(size = 14, face = "bold"),
          axis.text    = element_text(size = 12),
          axis.title   = element_text(size = 14, face = "bold")
        )
    }

    # Create plot
    p2b <- make_fig2b_plot(df_sus)

    # Save (use your robust saver so you don’t get blank PDFs)
    if (isTRUE(is_save)) {
      save_plot_dual(
        p2b,
        filename_base = "Figure2b",
        out_dir = paths$figures_dir,
        width = 3.2,
        height = 2.6,
        dpi = 300
      )
    } else {
      print(p2b)
    }
  }

  # ==============================================================================
  # MODULE 4: Test data cleaning + sensor alignment (AW/E4 matching)
  # ==============================================================================
  Df_Exam <- Df_Exam[!Df_Exam$Question.Type == "Example", ]
  Df_Exam <- Df_Exam[Df_Exam$Attempt == 1, ]
  Df_Exam <- Df_Exam[!is.na(Df_Exam$Question.Name), ]

  Df_Exam.SubSet <- subset(Df_Exam, select = -c(Timestamp))

  colSums(!is.na(Df_Exam.SubSet))

  Df_Exam.SubSet$HRV.IBINorm <- NA
  Df_Exam.SubSet$HRV.IBINorm <- round(Df_Exam.SubSet$HRV.IBI / Df_Exam.SubSet$HR.E4, 5)

  Df_Exam.SubSet["HRE4Perf"] <- Df_Exam.SubSet$HR.E4
  Df_Exam.SubSet[is.na(Df_Exam.SubSet$HR.AW), ]$HRE4Perf <- NA

  print("Before matching")
  colSums(!is.na(Df_Exam.SubSet %>% dplyr::select(
    ParticipantID, Time, Question.Name, Question.Type, Perspiration, HR.E4, HR.AW, HRV.IBI, HRV.IBINorm
  )))

  idx <- is.na(Df_Exam.SubSet$HR.E4) | is.na(Df_Exam.SubSet$HR.AW)
  Df_Exam.SubSet$HR.E4[idx] <- NA
  Df_Exam.SubSet$HR.AW[idx] <- NA

  rownames(Df_Exam.SubSet) <- 1:nrow(Df_Exam.SubSet)
  Df_Exam.SubSet <- droplevels(Df_Exam.SubSet)

  print("After matching")
  colSums(!is.na(Df_Exam.SubSet %>% dplyr::select(
    ParticipantID, Time, Question.Name, Question.Type, Perspiration, HR.E4, HR.AW, HRV.IBI, HRV.IBINorm
  )))
  # Keep a copy of the AW–E4 matched (pre-Cook) data for Table 1 counts
  Df_Exam.SubSet.matched <- Df_Exam.SubSet



  # ==============================================================================
  # MODULE 5: Outlier handling — Cook's distance (Figure 4) on HR.E4 vs HR.AW
  # ==============================================================================
  signal.lm.HRs <- lm(HR.E4 ~ HR.AW, data = Df_Exam.SubSet)

  cooksD    <- cooks.distance(signal.lm.HRs)
  cooksD.95 <- quantile(cooksD, prob = .95)

  influential          <- cooksD[cooksD > cooksD.95]
  names_of_influential <- names(influential)

  df_outlier <- Df_Exam.SubSet[names_of_influential, ]

  Df_Exam.SubSet.afterCook <- Df_Exam.SubSet
  Df_Exam.SubSet.afterCook[names_of_influential, c("HR.E4","HR.AW")] <- NA

  rownames(Df_Exam.SubSet.afterCook) <- seq_len(nrow(Df_Exam.SubSet.afterCook))

  cat(sprintf("Before Cook's Distance:"))
  colSums(!is.na(Df_Exam.SubSet))

  cat(sprintf("After Cook's Distance:"))
  colSums(!is.na(Df_Exam.SubSet.afterCook))


  hr_limits <- c(50, 140)

  # --- Correlation labels (lowercase p; subscripts upright like Fig. 6) ---
  df_bfr_cor <- Df_Exam.SubSet %>% dplyr::select(HR.E4, HR.AW) %>% tidyr::drop_na()
  ct_bfr     <- suppressWarnings(cor.test(df_bfr_cor$HR.E4, df_bfr_cor$HR.AW, method = "pearson"))
  p_bfr_txt  <- ifelse(is.finite(ct_bfr$p.value) && ct_bfr$p.value < 0.001,
                      "<0.001",
                      sprintf("==%.3f", ct_bfr$p.value))
  lab_bfr    <- sprintf("italic(n)==%d~~italic(p)%s~~italic(R)==%.2f",
                      nrow(df_bfr_cor),
                      p_bfr_txt,
                      unname(ct_bfr$estimate))

  df_aftr_cor <- Df_Exam.SubSet.afterCook %>% dplyr::select(HR.E4, HR.AW) %>% tidyr::drop_na()
  ct_aftr     <- suppressWarnings(cor.test(df_aftr_cor$HR.E4, df_aftr_cor$HR.AW, method = "pearson"))
  p_aftr_txt  <- ifelse(is.finite(ct_aftr$p.value) && ct_aftr$p.value < 0.001,
                       "<0.001",
                       sprintf("==%.3f", ct_aftr$p.value))
  lab_aftr    <- sprintf("italic(n)==%d~~italic(p)%s~~italic(R)==%.2f",
                       nrow(df_aftr_cor),
                       p_aftr_txt,
                       unname(ct_aftr$estimate))

  bfrCD.plot <- ggplot(Df_Exam.SubSet, aes(x = HR.E4, y = HR.AW)) +
    ggrastr::geom_point_rast(alpha = 0.25, size = 1.2, raster.dpi = 300) +
    geom_abline(colour = "red") +
    # Upright subscripts via plotmath expressions (match Fig. 6)
    ylab(expression(italic(HR)[AW]~"[bpm]")) +
    xlab(expression(italic(HR)[E4]~"[bpm]")) +
    theme(
      text = element_text(size = 18, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 16)
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = lab_bfr,
      parse = TRUE,
      hjust = -0.05, vjust = 1.25,
      size = 6,
      fontface = "bold"
    ) +
    coord_cartesian(xlim = hr_limits, ylim = hr_limits)

  aftrCD.plot <- ggplot(Df_Exam.SubSet.afterCook, aes(x = HR.E4, y = HR.AW)) +
    ggrastr::geom_point_rast(alpha = 0.25, size = 1.2, raster.dpi = 300) +
    geom_abline(colour = "red") +
    ylab("") + xlab(expression(italic(HR)[E4]~"[bpm]")) +
    theme(
      text = element_text(size = 18, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = lab_aftr,
      parse = TRUE,
      hjust = -0.05, vjust = 1.25,
      size = 6,
      fontface = "bold"
    ) +
    coord_cartesian(xlim = hr_limits, ylim = hr_limits)

  bfr.aftr <- cowplot::plot_grid(
    bfrCD.plot, aftrCD.plot,
    nrow = 1,
    scale = c(.97, .97),
    labels = c("a", "b"),
    label_size = 20
  )

  # Draw/save Figure 4 safely (cowplot can occasionally error under some ggplot2 versions)
  tryCatch({
    print(bfr.aftr)
    if (is_save == TRUE) {
      # Save as both PDF and PNG
      ggsave(sprintf("%sFigure4.pdf", plot_dir), bfr.aftr, width = 10, height = 5, units = "in", device = cairo_pdf,
        bg = "white"
      )
      ggsave(sprintf("%sFigure4.png", plot_dir), bfr.aftr, width = 10, height = 5, units = "in", dpi = 300)
    }
  }, error = function(e) {
    warning("Skipping Figure4 (Cook's distance) due to plotting error: ", conditionMessage(e))
  })


  # ==============================================================================
  # MODULE 6: Baseline normalization + bestNormalize transforms (Figure 5)
  # ==============================================================================
  Df_Exam.SubSet <- Df_Exam.SubSet.afterCook

  Df_Exam.SubSet$FPNorm    <- NA
  Df_Exam.SubSet$HR.E4Norm <- NA
  Df_Exam.SubSet$HR.AWNorm <- NA

  for (p in unique(Df_Exam.SubSet$ParticipantID)) {

    tmpExam_pp    <- Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, ]$Perspiration
    tmpBL_pp.mean <- Df_BL.mean[Df_BL.mean$ParticipantID == p, ]$Perspiration

    tmpExam_HR.E4     <- Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, ]$HR.E4
    tmpBL_HR.E4_Mean  <- Df_BL.mean[Df_BL.mean$ParticipantID == p, ]$HR.E4

    tmpExam_HR.AW     <- Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, ]$HR.AW
    tmpBL_HR.AW_Mean  <- Df_BL.mean[Df_BL.mean$ParticipantID == p, ]$HR.AW

    Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, ]$FPNorm    <- tmpExam_pp    - tmpBL_pp.mean
    Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, ]$HR.E4Norm <- tmpExam_HR.E4 - tmpBL_HR.E4_Mean
    Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, ]$HR.AWNorm <- tmpExam_HR.AW - tmpBL_HR.AW_Mean
  }


  cat("-- FPNorm Best Normalization( Correction) \n")
  FP_bst_method <- bestNormalize(Df_Exam.SubSet$FPNorm, na.rm = TRUE)
  print(FP_bst_method)
  Df_Exam.SubSet$FPNorm.Corrected <- FP_bst_method$x.t

  cat("-- HR.E4Norm Best Normalization( Correction) \n")
  HR_E4_bst_method <- bestNormalize(Df_Exam.SubSet$HR.E4Norm, na.rm = TRUE)
  print(HR_E4_bst_method)
  Df_Exam.SubSet$HR.E4Norm.Corrected <- HR_E4_bst_method$x.t

  cat("-- HR.AWNorm Best Normalization( Correction) \n")
  HR_AW_bst_method <- bestNormalize(Df_Exam.SubSet$HR.AWNorm, na.rm = TRUE)
  print(HR_AW_bst_method)
  Df_Exam.SubSet$HR.AWNorm.Corrected <- HR_AW_bst_method$x.t

  HRV_bst_method <- bestNormalize(Df_Exam.SubSet$HRV.IBINorm, na.rm = TRUE)
  cat("-- HRVNorm Best Normalization( Correction) \n")
  print(HRV_bst_method)
  Df_Exam.SubSet$NHRV.Corrected <- HRV_bst_method$x.t



  vars <- c(
    "FPNorm",
    "FPNorm.Corrected",
    "HR.AWNorm",
    "HR.AWNorm.Corrected",
    "HR.E4Norm",
    "HR.E4Norm.Corrected",
    "HRV.IBINorm",
    "NHRV.Corrected"
  )

  titles <- list(
    expression(italic(CFP)),
    expression(italic(FPN)),
    expression(italic(CHR)[AW]),
    expression(italic(NHR)[AW]),
    expression(italic(CHR)[E4]),
    expression(italic(NHR)[E4]),
    expression(italic(CHRV)),
    expression(italic(NHRV))
  )

  ncol <- 2
  nrow <- 4

  qq_plots <- lapply(seq_along(vars), function(i) {
    var_i   <- vars[i]
    title_i <- titles[[i]]

    row_i <- ceiling(i / ncol)
    col_i <- ifelse(i %% ncol == 1, 1, 2)

    p <- ggplot(Df_Exam.SubSet, aes_string(sample = var_i)) +
      ggrastr::rasterise(geom_qq(color = "black", size = 0.8), dpi = 300) +
      ggplot2::stat_qq_line(color = "red", linewidth = 0.8) +
      my_theme1 +
      ggtitle(title_i) +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

    n <- sum(!is.na(Df_Exam.SubSet[[var_i]]))
    p <- p + annotate("text", x = 0, y = Inf, hjust = 0.5, vjust = 1.5, parse = TRUE, size = 5,
                      label = paste0("~italic(n)", " == ", n))

    if (col_i > 1) p <- p + theme(axis.title.y = element_blank())
    if (row_i < nrow) p <- p + theme(axis.title.x = element_blank())

    # Bring the bottom-row x-axis title ("Theoretical Quantiles") a tad lower,
    # i.e., add more separation from the x tick labels.
    if (row_i == nrow) {
      p <- p + theme(axis.title.x = element_text(margin = margin(t = 10)))
    }

    p
  })

  final <- cowplot::plot_grid(plotlist = qq_plots, ncol = ncol, nrow = nrow, align = "hv")
  print(final)

  if (is_save == TRUE) {
    # PNG (raster)
    ggsave(
      filename = file.path(plot_dir, "Figure5.png"),
      plot     = final,
      width    = 7,
      height   = 9,
      dpi      = 300
    )

    # PDF: render via base pdf() + grid.draw() for maximum robustness.
    # (Some systems/viewers can intermittently write a blank page when saving
    # a cowplot::plot_grid object via ggsave() with certain PDF devices.)
    grDevices::pdf(
      file   = file.path(plot_dir, "Figure5.pdf"),
      width  = 7,
      height = 9,
      useDingbats = FALSE
    )
    grid::grid.newpage()
    grid::grid.draw(ggplotGrob(final))
    grDevices::dev.off()
  }



  summ_stat_Participant <- Df_Exam.SubSet %>%
    group_by(Question.Type) %>%
    get_summary_stats(
      Time, Perspiration, HR.AW, HR.E4, HRV.IBI,
      show = c("mean", "sd", "median", "min", "max" ,"n")
    )

  summ_stat_Participant


  # ==============================================================================
  # MODULE 7: Question-level aggregation (Qlevel) + sanity checks + Table 1 export
  # ==============================================================================
  Qlevel <- Df_Exam.SubSet %>%
    group_by(ParticipantID, Question.Name, Question.Type) %>%
    summarise(
      Gender       = first(Gender),
      Perspiration = mean(Perspiration, na.rm = TRUE),
      FPNorm       = mean(FPNorm.Corrected, na.rm = TRUE),
      HR.AW        = mean(HR.AW, na.rm = TRUE),
      HR.AWNorm    = mean(HR.AWNorm.Corrected, na.rm = TRUE),
      HR.E4        = mean(HR.E4, na.rm = TRUE),
      HR.E4Norm    = mean(HR.E4Norm.Corrected, na.rm = TRUE),
      HRV.IBI      = mean(HRV.IBI, na.rm = TRUE),
      HRVNorm      = mean(NHRV.Corrected, na.rm = TRUE),
      SAI          = first(SAI.Score),
      QType        = first(Question.Type),
      QOrder       = first(Question.Order),
      QNumber      = first(Question.Number),
      QTime        = n(),
      Grade        = first(as.character(Attempt.Correctness)),
      SUS          = first(SUS.Score),
      .groups = "drop"
    )

  Qlevel$Gender <- as.character(Qlevel$Gender)
  Qlevel$Gender <- relevel(factor(Qlevel$Gender), ref = "M")
  Qlevel$Grade  <- relevel(factor(Qlevel$Grade),  ref = "0")
  Qlevel$Question.Type <- factor(Qlevel$Question.Type, levels = c("V", "A", "W"))

  print("Number of participants in Qlevel data frame")
  colSums(!is.na(Qlevel))


  all_QNames <- unique(Df_Exam.SubSet$Question.Name)
  all_PIDs   <- unique(Df_Exam.SubSet$ParticipantID)

  missing_QNames <- data.frame(ParticipantID = character(), Missing_QName = character(), stringsAsFactors = FALSE)

  for (p in all_PIDs) {
    QNames <- unique(Df_Exam.SubSet[Df_Exam.SubSet$ParticipantID == p, "Question.Name"])
    missing_names <- setdiff(all_QNames, QNames)
    if (length(missing_names) > 0) {
      missing_QNames <- rbind(
        missing_QNames,
        data.frame(ParticipantID = p, Missing_QName = paste(missing_names, collapse = ", "), stringsAsFactors = FALSE)
      )
    }
  }

  missing_QNames



  # Helper for safe quantiles
  qfun <- function(x, p) as.numeric(quantile(x, probs = p, na.rm = TRUE, type = 7))

  # 1) thresholds computed on NORMALIZED Qlevel signals
  participant_stats <- Qlevel %>%
    group_by(ParticipantID) %>%
    summarise(
      # mean thresholds
      FPNorm_mean    = mean(FPNorm,    na.rm = TRUE),
      HR.E4Norm_mean = mean(HR.E4Norm, na.rm = TRUE),
      HR.AWNorm_mean = mean(HR.AWNorm, na.rm = TRUE),
      HRVNorm_mean   = mean(HRVNorm,   na.rm = TRUE),

      # "33% stress" thresholds: top-33% => q67 for PP/HR; bottom-33% => q33 for HRV
      FPNorm_q67     = qfun(FPNorm,    0.67),
      HR.E4Norm_q67  = qfun(HR.E4Norm, 0.67),
      HR.AWNorm_q67  = qfun(HR.AWNorm, 0.67),
      HRVNorm_q33    = qfun(HRVNorm,   0.33),

      # Ns (optional)
      FPNorm_n       = sum(!is.na(FPNorm)),
      HR.E4Norm_n    = sum(!is.na(HR.E4Norm)),
      HR.AWNorm_n    = sum(!is.na(HR.AWNorm)),
      HRVNorm_n      = sum(!is.na(HRVNorm)),
      .groups = "drop"
    )

  # 2) create BOTH label sets
  Qlevel_Simpl_Stress <- Qlevel %>%
    left_join(participant_stats, by = "ParticipantID") %>%
    mutate(
      # Mean-split (keep original names for backward compatibility)
      Stress.fp   = factor(if_else(FPNorm     > FPNorm_mean,     "S", "NS"), levels = c("NS","S")),
      Stress.hre4 = factor(if_else(HR.E4Norm  > HR.E4Norm_mean,  "S", "NS"), levels = c("NS","S")),
      Stress.hraw = factor(if_else(HR.AWNorm  > HR.AWNorm_mean,  "S", "NS"), levels = c("NS","S")),
      Stress.nhrv = factor(if_else(HRVNorm    < HRVNorm_mean,    "S", "NS"), levels = c("NS","S")),

      # Quantile "33% stress" labels
      Stress.fp_q33   = factor(if_else(FPNorm     > FPNorm_q67,    "S", "NS"), levels = c("NS","S")),
      Stress.hre4_q33 = factor(if_else(HR.E4Norm  > HR.E4Norm_q67, "S", "NS"), levels = c("NS","S")),
      Stress.hraw_q33 = factor(if_else(HR.AWNorm  > HR.AWNorm_q67, "S", "NS"), levels = c("NS","S")),
      Stress.nhrv_q33 = factor(if_else(HRVNorm    < HRVNorm_q33,   "S", "NS"), levels = c("NS","S")),

      Gender = factor(Gender, levels = c("F","M")),
      QType  = factor(QType,  levels = c("V","A","W"))
    )

  names(Qlevel_Simpl_Stress)


  colNamesFACS <- c("F_Angry", "F_Disgusted", "F_Afraid", "F_Happy", "F_Sad", "F_Surprised", "F_Neutral")

  Df_Exam.FACS <- Df_Exam.SubSet %>%
    dplyr::select(ParticipantID, Question.Name, Question.Type, dplyr::all_of(colNamesFACS))

  Df_Exam.FACS <- Df_Exam.FACS %>%
    mutate(
      F_MaxEmotion = colNamesFACS[
        max.col(dplyr::select(., dplyr::all_of(colNamesFACS)), ties.method = "first")
      ]
    )

  Df_Exam.SubSet$F_MaxEmotion <- as.factor(Df_Exam.FACS$F_MaxEmotion)

  # Keep only rows with complete FACS vector
  Df_Exam.FACS <- tidyr::drop_na(Df_Exam.FACS)



  FACS.Qlevel <- Df_Exam.FACS %>%
    group_by(ParticipantID, Question.Name, Question.Type) %>%
    summarise(
      Angry      = mean(F_Angry,     na.rm = TRUE),
      Disgusted  = mean(F_Disgusted, na.rm = TRUE),
      Afraid     = mean(F_Afraid,    na.rm = TRUE),
      Happy      = mean(F_Happy,     na.rm = TRUE),
      Sad        = mean(F_Sad,       na.rm = TRUE),
      Surprised  = mean(F_Surprised, na.rm = TRUE),
      Neutral    = mean(F_Neutral,   na.rm = TRUE),
      .groups = "drop"
    )

  FACS.Qlevel$Question.Type <- factor(FACS.Qlevel$Question.Type, levels = c("V", "A", "W"))

  print("Number of participants in Qlevel FACS data frame")
  table(FACS.Qlevel$ParticipantID)


  Qlevel_Stress_FACS <- Qlevel_Simpl_Stress %>%
    left_join(FACS.Qlevel, by = c("ParticipantID", "Question.Name", "Question.Type"))

  names(Qlevel_Stress_FACS)

  # Original output used by modeling scripts
  write.csv(Qlevel_Stress_FACS, file = sprintf("%sAffective_Math_Qlevel_Data_N%s.csv", data_dir, N), row.names = FALSE)

  # ============================================================
  # Table 1: Channel-wise sample counts across the processing pipeline
  # (Exports CSV + TeX under /reports; `reports_dir` is defined in the main config above.)
  # ============================================================

  nn <- function(x) sum(!is.na(x))

  # Sec-Level counts (raw exam file, before cross-channel matching / Cook's filter)
  sec_fp   <- nn(Df_Exam$Perspiration)
  sec_hre4 <- nn(Df_Exam$HR.E4)
  sec_hraw <- nn(Df_Exam$HR.AW)
  sec_hrv  <- nn(Df_Exam$HRV.IBI)
  sec_e    <- nrow(Df_Exam)  # exam timeline rows (E-channel in Table 1)

  # Cross-Channel counts (AW–E4 matched)
  cross_fp   <- sec_fp
  cross_hre4 <- nn(Df_Exam.SubSet.matched$HR.E4)
  cross_hraw <- nn(Df_Exam.SubSet.matched$HR.AW)
  cross_hrv  <- sec_hrv
  cross_e    <- sec_e

  # Cook's Filter (AW–E4 only)
  cooks_fp   <- cross_fp
  cooks_hre4 <- nn(Df_Exam.SubSet.afterCook$HR.E4)
  cooks_hraw <- nn(Df_Exam.SubSet.afterCook$HR.AW)
  cooks_hrv  <- cross_hrv
  cooks_e    <- cross_e

  # Q-Level counts (question-aggregated)
  q_fp   <- nn(Qlevel$FPNorm)
  q_hre4 <- nn(Qlevel$HR.E4Norm)
  q_hraw <- nn(Qlevel$HR.AWNorm)
  q_hrv  <- nn(Qlevel$HRVNorm)

  # Facial-expression availability at Q-Level (joined into Qlevel_Stress_FACS)
  fac_cols <- c("Angry","Disgusted","Afraid","Happy","Sad","Surprised","Neutral")
  q_e <- if (all(fac_cols %in% names(Qlevel_Stress_FACS))) nn(Qlevel_Stress_FACS$Neutral) else NA_integer_

  table1_counts <- tibble::tibble(
    Channel       = c("FP","HRE4","HRAW","HRV","E"),
    Sec_Level     = c(sec_fp, sec_hre4, sec_hraw, sec_hrv, sec_e),
    Cross_Channel = c(cross_fp, cross_hre4, cross_hraw, cross_hrv, cross_e),
    Cooks_Filter  = c(cooks_fp, cooks_hre4, cooks_hraw, cooks_hrv, cooks_e),
    Q_Level       = c(q_fp, q_hre4, q_hraw, q_hrv, q_e)
  )

  out_table1 <- file.path(reports_dir, sprintf("Table1_Channel_Wise_Sample_Counts_N%s.csv", N))
  readr::write_csv(table1_counts, out_table1)
  message("Wrote Table 1 counts CSV: ", out_table1)

  # ---- TeX (booktabs; no caption) ----
  out_table1_tex <- file.path(reports_dir, sprintf("Table1_Channel_Wise_Sample_Counts_N%s.tex", N))

  tex1 <- c(
    "\\begin{table}[t]", "\\centering",
    "\\small",
    "\\setlength{\\tabcolsep}{6pt}",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Channel & Sec. Level & Cross-Channel & Cook's Filter & Q-Level \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(table1_counts))) {
    r <- table1_counts[i, ]
    tex1 <- c(tex1,
              sprintf("%s & %d & %d & %d & %s \\\\",
                      r$Channel,
                      r$Sec_Level,
                      r$Cross_Channel,
                      r$Cooks_Filter,
                      ifelse(is.na(r$Q_Level), "--", as.character(r$Q_Level))
              )
    )
  }

  tex1 <- c(tex1,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}"
  )

  writeLines(tex1, out_table1_tex)
  message("Wrote Table 1 TeX: ", out_table1_tex)


  # Put the y-axis label as a top-centered label (just outside the panel)
  top_panel_label <- function(p, lab, size = 12, top_margin = 18) {
    p +
      labs(y = NULL, title = lab) +
      theme(
        axis.title.y = element_blank(),
        plot.title = element_text(
          hjust = 0.5, vjust = 0, face = "bold", size = size,
          margin = margin(b = 3)
        ),
        plot.title.position = "panel",   # <<<<<< key change
        plot.margin = margin(t = top_margin, r = 6, b = 2, l = 6)
      )
  }

  # Panel tag (a/b) in the top-left OUTSIDE the panel rectangle (in plot margin area)
  add_panel_tag <- function(p, tag, size = 12, x = 0.01, y = 0.995) {
    cowplot::ggdraw(p) +
      cowplot::draw_label(
        tag,
        x = x, y = y,
        hjust = 0, vjust = 1,
        fontface = "bold",
        size = size
      )
  }

  # ==============================================================================
  # MODULE 8: Paper figures (Figure 6–7)
  # ==============================================================================

  # =========================
  # Figure 6: Question-level physiological time series (example participant)
  # =========================
  pid_fig6 <- "S006"

  df_fig6 <- Qlevel %>%
    dplyr::filter(ParticipantID == pid_fig6) %>%
    dplyr::arrange(QNumber) %>%
    dplyr::mutate(
      QNumber = as.numeric(QNumber),
      QType   = as.character(QType)
    )

  if (nrow(df_fig6) > 0) {

    # Baseline means from baseline segment
    base_fp   <- Df_BL.mean$Perspiration[match(pid_fig6, Df_BL.mean$ParticipantID)]
    base_hraw <- Df_BL.mean$HR.AW[match(pid_fig6, Df_BL.mean$ParticipantID)]
    base_hre4 <- Df_BL.mean$HR.E4[match(pid_fig6, Df_BL.mean$ParticipantID)]

    # Pastel tones (Word/Video/Abstract)
    qtype_cols_leg <- c(
      "W" = "#9A7BC0",
      "V" = "#F08DB4",
      "A" = "#A7DB96"
    )

    # Legend (dots) at bottom
    legend_df <- data.frame(QType = factor(c("W","V","A"), levels = c("W","V","A")))

    p_leg <- ggplot(legend_df, aes(x = QType, y = 1, color = QType)) +
      geom_point(size = 3) +
      scale_color_manual(
        values = qtype_cols_leg,
        breaks = c("W","V","A"),
        labels = c("Word","Video","Abstract")
      ) +
      guides(color = guide_legend(title = NULL, nrow = 1, byrow = TRUE,
                                  override.aes = list(size = 4))) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.box = "horizontal"
      )

    leg_grob <- cowplot::get_legend(p_leg)

    # Identify contiguous format periods
    seg_df <- df_fig6 %>%
      dplyr::mutate(seg = cumsum(.data$QType != dplyr::lag(.data$QType, default = first(.data$QType)))) %>%
      dplyr::group_by(seg, QType) %>%
      dplyr::summarise(
        x_min = min(.data$QNumber, na.rm = TRUE) - 0.5,
        x_max = max(.data$QNumber, na.rm = TRUE) + 0.5,
        .groups = "drop"
      )

    stopifnot(is.numeric(seg_df$x_min), is.numeric(seg_df$x_max))

    # Helper: add thin format rectangles at bottom
    add_format_rects <- function(p, y_vals, top_pad = 0.10, bottom_pad = 0.18) {
      y_rng  <- range(y_vals, na.rm = TRUE)
      y_span <- y_rng[2] - y_rng[1]

      y_bar_min <- y_rng[1] - 0.16 * y_span
      y_bar_max <- y_rng[1] - 0.12 * y_span

      p +
        geom_rect(
          data = seg_df,
          aes(xmin = x_min, xmax = x_max, ymin = y_bar_min, ymax = y_bar_max, fill = QType),
          inherit.aes = FALSE,
          color = NA
        ) +
        scale_fill_manual(values = qtype_cols_leg, drop = FALSE, guide = "none") +
        coord_cartesian(
          ylim = c(y_rng[1] - bottom_pad * y_span, y_rng[2] + top_pad * y_span),
          clip = "off"
        ) +
        theme(legend.position = "none")
    }

    # -------------------------
    # Panel (a): FP
    # -------------------------
    n_fp <- sum(is.finite(df_fig6$Perspiration))   # <-- add this

    p_fp_ts <- ggplot(df_fig6, aes(x = QNumber, y = Perspiration)) +
      geom_line(color = "gray20", linewidth = 0.5) +
      geom_hline(yintercept = base_fp, linetype = "dashed", color = "dodgerblue3", linewidth = 0.8) +
      labs(x = "Question Number", y = NULL) +
      my_theme1

    p_fp_ts <- add_format_rects(p_fp_ts, df_fig6$Perspiration, top_pad = 0.10, bottom_pad = 0.18)

    p_fp_ts <- p_fp_ts +
      scale_y_continuous(
        breaks = c(0.006, 0.008, 0.010, 0.012),
        labels = scales::number_format(accuracy = 0.001)
      )

    # <-- add this annotation (top-left INSIDE the panel)
    p_fp_ts <- p_fp_ts +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = sprintf("italic(n)[FP]==%d", n_fp),
        parse = TRUE,
        hjust = 1.08,   # pull left from right border (increase if needed)
        vjust = 1.45,   # pull down from top border (increase if needed)
        size = 4.3,
        fontface = "bold"
      )

    p_fp_ts <- top_panel_label(p_fp_ts, expression(italic(FP) * " [" * degree * "C"^2 * "]"), top_margin = 18)
    p_fp_ts <- add_panel_tag(p_fp_ts, "a", size = 14, x = 0.03, y = 0.98)

    # -------------------------
    # Panel (b): HR_E4 line + HR_AW points
    # -------------------------
    n_aw <- sum(is.finite(df_fig6$HR.AW))
    n_e4 <- sum(is.finite(df_fig6$HR.E4))

    yr_hr <- range(c(df_fig6$HR.E4, df_fig6$HR.AW), na.rm = TRUE)
    hr_breaks <- seq(floor(yr_hr[1] / 2) * 2, ceiling(yr_hr[2] / 2) * 2, by = 2)

    x_mid <- mean(range(df_fig6$QNumber, na.rm = TRUE))

    p_hr_ts <- ggplot(df_fig6, aes(x = QNumber)) +
      geom_line(aes(y = HR.E4), color = "gray20", linewidth = 0.5) +
      geom_point(aes(y = HR.AW), color = "red3", size = 1.6, alpha = 0.95) +
      geom_hline(yintercept = base_hre4, linetype = "dashed", color = "dodgerblue3", linewidth = 0.8) +
      geom_hline(yintercept = base_hraw, linetype = "solid",  color = "red3",  linewidth = 0.9) +
      annotate(
        "text",
        x = x_mid, y = Inf,                    # <- anchor to top
        label = sprintf("italic(n)[AW]==%d~~italic(n)[E4]==%d", n_aw, n_e4),
        parse = TRUE,
        size = 4.3, fontface = "bold",
        hjust = 0.5,
        vjust = 1.35                           # <- pushes DOWN so it won't touch the border
      ) +
      scale_y_continuous(breaks = hr_breaks) +
      labs(x = "Question Number", y = NULL) +
      my_theme1

    p_hr_ts <- add_format_rects(p_hr_ts, c(df_fig6$HR.E4, df_fig6$HR.AW), top_pad = 0.26, bottom_pad = 0.18)
    p_hr_ts <- top_panel_label(p_hr_ts, expression(italic(HR)[AW] * " & " * italic(HR)[E4] * " [bpm]"), top_margin = 18)
    p_hr_ts <- add_panel_tag(p_hr_ts, "b", size = 14, x = 0.03, y = 0.98)

    # Combine + legend
    gap <- ggplot() + theme_void()

    fig6_main <- cowplot::plot_grid(
      p_fp_ts,
      gap,          # <- spacer
      p_hr_ts,
      nrow = 1,
      rel_widths = c(1, 0.06, 1),   # increase 0.06 -> 0.08/0.10 for more gap
      align = "hv"
    )

    fig6 <- cowplot::plot_grid(
      fig6_main,
      leg_grob,
      ncol = 1,
      rel_heights = c(1, 0.18)
    )

    print(fig6)

    if (is_save == TRUE) {
      ggsave(file.path(plot_dir, "Figure6.pdf"), fig6, width = 7.0, height = 3.0, units = "in")
      ggsave(file.path(plot_dir, "Figure6.png"), fig6, width = 7.0, height = 3.0, units = "in", dpi = 300)
    }
  }

  # =========================
  # Figure 7: Channel-specific relative stress labeling (example participant)
  # =========================
  pid_fig7 <- "S002"

  stress_cols <- c(
    "RLS" = "#B2DF8A",  # light green
    "RHS" = "red3"      # red
  )

  df_fig7 <- Qlevel_Stress_FACS %>%
    filter(ParticipantID == pid_fig7) %>%
    arrange(QNumber)

  if (nrow(df_fig7) > 0) {

    make_hist <- function(values, reverse = FALSE, ylab = NULL, show_legend = FALSE, bins = 41) {

      vals <- values[is.finite(values)]
      thr  <- mean(vals, na.rm = TRUE)

      h <- hist(vals, breaks = bins, plot = FALSE)
      binwidth <- diff(h$breaks)[1]

      dd <- data.frame(
        xmin  = h$breaks[-length(h$breaks)],
        xmax  = h$breaks[-1],
        xmid  = (h$breaks[-length(h$breaks)] + h$breaks[-1]) / 2,
        count = h$counts
      )

      dd$lab <- if (reverse) {
        ifelse(dd$xmid < thr, "RHS", "RLS")
      } else {
        ifelse(dd$xmid > thr, "RHS", "RLS")
      }
      dd$lab <- factor(dd$lab, levels = c("RLS", "RHS"))

      ggplot(dd, aes(x = xmid, y = count, fill = lab)) +
        geom_col(
          width = 0.65 * binwidth,
          color = "black", linewidth = 0.15
        ) +
        geom_vline(xintercept = thr, linetype = "dotted", color = "dodgerblue3", linewidth = 0.9) +
        scale_fill_manual(values = stress_cols, drop = FALSE) +
        labs(x = NULL, y = ylab, fill = pid_fig7) +
        guides(fill = guide_legend(
          title.position = "left",
          direction = "horizontal",
          nrow = 1,
          byrow = TRUE
        )) +
        theme_classic(base_size = 10) +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          axis.line = element_blank(),
          legend.position = "none",

          # --- spacing controls you asked for ---
          legend.title = element_text(size = 10, margin = margin(r = 10)),  # S002 -> green key gap
          legend.spacing.x = grid::unit(24, "pt"),                               # RLS -> next (red) key gap

          legend.text  = element_text(size = 10),
          plot.margin  = margin(2, 2, 2, 2)
        )
    }

    make_qq <- function(values, conf = 0.95) {
      vals <- values[is.finite(values)]
      n <- length(vals)
      if (n < 3) stop("Not enough finite values for QQ plot.")

      # Order statistics
      samp <- sort(vals)
      i    <- seq_len(n)
      p    <- (i - 0.5) / n
      theo <- qnorm(p)

      # qqline-style fit (through Q1/Q3), like base R's qqline()
      qy <- quantile(samp, c(0.25, 0.75), names = FALSE)
      qx <- qnorm(c(0.25, 0.75))
      b  <- (qy[2] - qy[1]) / (qx[2] - qx[1])
      a  <- qy[1] - b * qx[1]

      # Approx. pointwise envelope for normal order stats (in theo-space),
      # then mapped through the same straight line y = a + b*theo
      zcrit <- qnorm((1 + conf) / 2)

      # Var approximation for normal quantiles of order stats:
      # Var(Z_p) ≈ p(1-p) / (n * φ(z)^2)
      se_theo <- sqrt(p * (1 - p) / (n * dnorm(theo)^2))

      y_line <- a + b * theo
      y_low  <- a + b * (theo - zcrit * se_theo)
      y_high <- a + b * (theo + zcrit * se_theo)

      dd <- data.frame(theo = theo, samp = samp, y_line = y_line, y_low = y_low, y_high = y_high)

      x_mid <- mean(range(theo, na.rm = TRUE))

      ggplot(dd, aes(x = theo, y = samp)) +
        geom_ribbon(aes(ymin = y_low, ymax = y_high),
                    fill = "gray85", color = NA) +
        geom_line(aes(y = y_line), color = "black", linewidth = 0.7) +
        ggrastr::geom_point_rast(size = 1, raster.dpi = 300) +
        annotate(
          "text",
          x = x_mid, y = Inf,
          label = sprintf("italic(n)==%d", n),
          parse = TRUE,
          hjust = 0.5,
          vjust = 1.65,   # pushes it a bit down from the top border
          size = 3.3
        ) +
        labs(x = NULL, y = NULL) +
        scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
        scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
        theme_classic(base_size = 10) +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          axis.line = element_blank(),
          plot.margin = margin(2, 2, 2, 2)
        )
    }

    # Build panels (left: histogram with mean threshold; right: QQ plot)
    p_fp_hist   <- make_hist(df_fig7$FPNorm,    reverse = FALSE,
                             ylab = expression(bar(italic(NFP))), show_legend = TRUE)
    p_fp_qq     <- make_qq(df_fig7$FPNorm)

    p_hraw_hist <- make_hist(df_fig7$HR.AWNorm, reverse = FALSE,
                             ylab = expression(bar(italic(NHR))[AW]))
    p_hraw_qq   <- make_qq(df_fig7$HR.AWNorm)

    p_hre4_hist <- make_hist(df_fig7$HR.E4Norm, reverse = FALSE,
                             ylab = expression(bar(italic(NHR))[E4]))
    p_hre4_qq   <- make_qq(df_fig7$HR.E4Norm)

    p_hrv_hist  <- make_hist(df_fig7$HRVNorm,   reverse = TRUE,
                             ylab = expression(bar(italic(NHRV))))
    p_hrv_qq    <- make_qq(df_fig7$HRVNorm)

    # --- Force 1 decimal on NHR x-axes (e.g., 1.0 instead of 1) ---
    fmt_1dp <- function(x) sprintf("%.1f", x)

    p_hraw_hist <- p_hraw_hist + scale_x_continuous(labels = fmt_1dp)

    p_hre4_hist <- p_hre4_hist + scale_x_continuous(labels = fmt_1dp)

    # --- Custom legend for Figure 7 (precise spacing control) ---
    legend_fig7 <- {
      # spacing knobs (in "data" x-units)
      x_title <- 2.5
      x_key1  <- 3.5   # S002 -> green key distance (increase to add more space)
      x_key2  <- 4.5   # RLS  -> red key distance (increase to add more space)
      dx_lab  <- 0.25  # gap between key and its label

      leg_df <- data.frame(
        key = factor(c("RLS","RHS"), levels = c("RLS","RHS")),
        xk  = c(x_key1, x_key2),
        xl  = c(x_key1 + dx_lab, x_key2 + dx_lab),
        lab = c("RLS","RHS")
      )

      ggplot() +
        annotate("text", x = x_title, y = 1, label = pid_fig7,
                 hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold") +
        geom_tile(
          data = leg_df,
          aes(x = xk, y = 1, fill = key),
          width = 0.42, height = 0.42,
          color = "black", linewidth = 0.2
        ) +
        geom_text(
          data = leg_df,
          aes(x = xl, y = 1, label = lab),
          hjust = 0, vjust = 0.5, size = 3.6, fontface = "bold"
        ) +
        scale_fill_manual(values = stress_cols, drop = FALSE, guide = "none") +
        coord_cartesian(xlim = c(0, 8), ylim = c(0.5, 1.5), clip = "off") +
        theme_void() +
        theme(plot.margin = margin(0, 2, 0, 2))
    }

    row1 <- cowplot::plot_grid(p_fp_hist,   p_fp_qq,   nrow = 1, rel_widths = c(1, 1))
    row2 <- cowplot::plot_grid(p_hraw_hist, p_hraw_qq, nrow = 1, rel_widths = c(1, 1))
    row3 <- cowplot::plot_grid(p_hre4_hist, p_hre4_qq, nrow = 1, rel_widths = c(1, 1))
    row4 <- cowplot::plot_grid(p_hrv_hist,  p_hrv_qq,  nrow = 1, rel_widths = c(1, 1))

    body_fig7 <- cowplot::plot_grid(row1, row2, row3, row4, ncol = 1, align = "hv")

    fig7 <- cowplot::plot_grid(legend_fig7, body_fig7, ncol = 1, rel_heights = c(0.12, 1))
    print(fig7)

    if (is_save == TRUE) {
      ggsave(file.path(plot_dir, "Figure7.pdf"), fig7, width = 7.0, height = 6.6, units = "in")
      ggsave(file.path(plot_dir, "Figure7.png"), fig7, width = 7.0, height = 6.6, units = "in", dpi = 300)
    }
  }


  # ==============================================================================
  # MODULE 9: Test-level dual-label export + Table 2 export + Figure 8
  # ==============================================================================
  Df_Exam_Dual <- Df_Exam.SubSet %>%
    left_join(
      Qlevel_Stress_FACS %>%
        dplyr::select(
          ParticipantID, Question.Name, Question.Type,
          Stress.fp, Stress.hre4, Stress.hraw, Stress.nhrv,
          Stress.fp_q33, Stress.hre4_q33, Stress.hraw_q33, Stress.nhrv_q33
        ),
      by = c("ParticipantID", "Question.Name", "Question.Type")
    )

  out_exam_dual <- sprintf("%sAffective_Math_Dataset_N%s_Test-Dual-Labels.csv", data_dir, N)
  write.csv(Df_Exam_Dual, file = out_exam_dual, row.names = FALSE)
  message("Wrote dual-label Test file: ", out_exam_dual)

  # ============================================================
  # Figure 8: (a) QTIME by matched k, (b) QGRADE by matched k, (c) SAI
  # Stacked vertically, single-column width (two-column paper)
  # FIXES:
  #  - removes duplicate pipelines
  #  - prevents clipping of k labels + a/b/c labels (adds headroom + clip="off")
  #  - keeps rectangular frames
  #  - enforces exactly 12 ks and exactly 3 highlighted ks
  # ============================================================


  # ---- Shared: math-italic tick labels for V/A/W ----
  x_labs_math <- c(expression(italic(V)),
                   expression(italic(A)),
                   expression(italic(W)))

  # ---- Theme: rectangular frame, no axis lines, extra margins for stacked plots ----
  theme_col_box <- theme_classic(base_size = 9) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      plot.margin = margin(6, 6, 6, 6)   # more breathing room so borders/labels don’t collide
    )

  # ============================================================
  # Panel (b): QGRADE by matched index k=1..12 (worst-3 Video ks highlighted)
  # ============================================================

  perf_q <- Qlevel %>%
    mutate(
      correct = as.integer(as.character(Grade)),
      QType   = factor(Question.Type, levels = c("V","A","W")),
      k       = as.integer(trimws(as.character(QOrder)))  # <-- k is the matched index (1..12)
    ) %>%
    filter(!is.na(k), !is.na(QType))

  perf_by_k <- perf_q %>%
    group_by(k, QType) %>%
    summarise(mean_acc = mean(correct, na.rm = TRUE), .groups = "drop") %>%
    mutate(x = as.numeric(QType))  # V=1, A=2, W=3

  # HARD sanity checks (avoid “mystery” ks)
  stopifnot(n_distinct(perf_by_k$k) == 12)
  stopifnot(nrow(perf_by_k) == 36)

  delta_tbl <- perf_by_k %>%
    select(k, QType, mean_acc) %>%
    pivot_wider(names_from = QType, values_from = mean_acc) %>%
    mutate(
      mean_AW = (A + W) / 2,
      delta_V_vs_AW = V - mean_AW
    )

  stopifnot(nrow(delta_tbl) == 12)
  stopifnot(all(is.finite(delta_tbl$delta_V_vs_AW)))

  # ============================================================
  # Table 2: Performance within matched triplets (ranked by Δ_V(k))
  #   Acc_F(k) = mean(QGRADE == 1) × 100
  #   Δ_V(k)   = Acc_V(k) - (Acc_A(k) + Acc_W(k))/2
  # Exports CSV + TeX under ../reports/ (no caption).
  # ============================================================

  table2_items <- Qlevel %>%
    transmute(
      QType = factor(Question.Type, levels = c("V","A","W")),
      k     = as.integer(trimws(as.character(QOrder))),
      item  = as.character(Question.Name)
    ) %>%
    filter(!is.na(k), !is.na(QType), !is.na(item)) %>%
    group_by(k, QType) %>%
    summarise(item = dplyr::first(unique(item)), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = QType, values_from = item) %>%
    rename(Vitem = V, Aitem = A, Witem = W)

  table2_n <- perf_q %>%
    group_by(k, QType) %>%
    summarise(n = sum(!is.na(correct)), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = QType, values_from = n) %>%
    # Paper column order: A, W, V
    mutate(n_A_W_V = paste0(A, "/", W, "/", V)) %>%
    select(k, n_A_W_V)

  table2_df <- delta_tbl %>%
    transmute(
      k       = k,
      Acc_A   = 100 * A,
      Acc_W   = 100 * W,
      Acc_V   = 100 * V,
      Delta_V = 100 * delta_V_vs_AW
    ) %>%
    mutate(
      Overall = (Acc_A + Acc_V + Acc_W) / 3
    ) %>%
    left_join(table2_items, by = "k") %>%
    left_join(table2_n,     by = "k") %>%
    arrange(Delta_V, k) %>%
    mutate(Rank = dplyr::row_number()) %>%
    # Paper column order: A, W, V (items, Acc, n)
    relocate(Rank, k, Aitem, Witem, Vitem, Acc_A, Acc_W, Acc_V, Delta_V, Overall, n_A_W_V) %>%
    mutate(
      across(c(Acc_A, Acc_V, Acc_W, Delta_V, Overall), ~ round(., 1))
    )

  stopifnot(nrow(table2_df) == 12)

  # CSV: use paper-like column headers
  table2_csv_df <- table2_df %>%
    rename(
      `A item`     = Aitem,
      `W item`     = Witem,
      `V item`     = Vitem,
      `Acc(A)`     = Acc_A,
      `Acc(W)`     = Acc_W,
      `Acc(V)`     = Acc_V,
      `Delta_V`    = Delta_V,
      `n(A/W/V)`   = n_A_W_V
    )

  out_table2_csv <- file.path(reports_dir, sprintf("Table2_Performance_Matched_Triplets_N%s.csv", N))
  readr::write_csv(table2_csv_df, out_table2_csv)
  message("Wrote Table 2 CSV: ", out_table2_csv)

  # ---- TeX (booktabs; no caption) ----
  out_table2_tex <- file.path(reports_dir, sprintf("Table2_Performance_Matched_Triplets_N%s.tex", N))

  fmt1 <- function(x) sprintf("%.1f", x)

  tex_lines <- c(
    "\\begin{table}[t]",  "\\centering",
    "\\small",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{rrlllrrrrrl}",
    "\\toprule",
    "Rank & $k$ & A item & W item & V item & $\\mathrm{Acc}(A)$ & $\\mathrm{Acc}(W)$ & $\\mathrm{Acc}(V)$ & $\\Delta_V$ & Overall & $n(A/W/V)$ \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(table2_df))) {
    r <- table2_df[i, ]
    tex_lines <- c(tex_lines,
      sprintf("%d & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        r$Rank, r$k, r$Aitem, r$Witem, r$Vitem,
        fmt1(r$Acc_A), fmt1(r$Acc_W), fmt1(r$Acc_V),
        fmt1(r$Delta_V), fmt1(r$Overall), r$n_A_W_V
      )
    )
  }

  tex_lines <- c(tex_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(tex_lines, out_table2_tex)
  message("Wrote Table 2 TeX: ", out_table2_tex)

  bad_k <- delta_tbl %>%
    slice_min(order_by = delta_V_vs_AW, n = 3, with_ties = FALSE) %>%  # EXACTLY 3
    arrange(delta_V_vs_AW, k) %>%
    pull(k)

  stopifnot(length(bad_k) == 3)

  cat("Figure8 highlight ks (worst by Δk = V - mean(A,W)):", paste(bad_k, collapse = ", "), "\n")

  perf_by_k <- perf_by_k %>%
    mutate(
      off    = scales::rescale(k, to = c(-0.12, 0.12)),
      x_off  = x + off,
      is_bad = k %in% bad_k
    )

  mean_by_format <- perf_by_k %>%
    group_by(QType) %>%
    summarise(mean_acc = mean(mean_acc, na.rm = TRUE), .groups = "drop") %>%
    mutate(x = as.numeric(QType))

  # y headroom so k labels aren’t clipped
  p_grade_by_k <- ggplot(perf_by_k, aes(x = x_off, y = mean_acc, group = k)) +
    geom_line(data = subset(perf_by_k, !is_bad), alpha = 0.25, linewidth = 0.5, color = "gray40") +
    geom_point(data = subset(perf_by_k, !is_bad), alpha = 0.45, size = 1.6, color = "gray20") +
    geom_line(data = subset(perf_by_k,  is_bad), linewidth = 0.9, color = "black") +
    geom_point(data = subset(perf_by_k, is_bad), size = 2.0, color = "black") +
    geom_text(
      data = subset(perf_by_k, QType == "W"),
      aes(label = k),
      size = 2.6,
      vjust = 0,                       # anchor at text baseline
      nudge_y = 0.035                  # nudge upward (in data units below via coord expansion)
    ) +
    geom_line(
      data = mean_by_format,
      aes(x = x, y = mean_acc, group = 1),
      inherit.aes = FALSE,
      linewidth = 0.8,
      linetype = "dashed",
      color = "gray30"
    ) +
    geom_point(
      data = mean_by_format,
      aes(x = x, y = mean_acc),
      inherit.aes = FALSE,
      size = 2.0,
      shape = 21,
      stroke = 0.6,
      color = "gray30",
      fill = "white"
    ) +
    scale_x_continuous(breaks = 1:3, labels = x_labs_math) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0.02, 0.14))  # TOP headroom for k labels
    ) +
    labs(x = NULL, y = expression(italic(QGRADE) == 1)) +
    coord_cartesian(clip = "off") +
    theme_col_box

  # ============================================================
  # Panel (a): QTIME by matched index k=1..12 (median)
  # ============================================================

  qtime_by_k <- Qlevel %>%
    mutate(
      k     = as.integer(trimws(as.character(QOrder))),
      QType = factor(Question.Type, levels = c("V","A","W")),
      QTime = as.numeric(QTime)
    ) %>%
    filter(!is.na(k), !is.na(QType), is.finite(QTime)) %>%
    group_by(k, QType) %>%
    summarise(
      med_time = median(QTime, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      x     = as.numeric(QType),
      off   = scales::rescale(k, to = c(-0.12, 0.12)),
      x_off = x + off
    )

  stopifnot(n_distinct(qtime_by_k$k) == 12)
  stopifnot(nrow(qtime_by_k) == 36)

  med_by_format <- qtime_by_k %>%
    group_by(QType) %>%
    summarise(med_time = median(med_time, na.rm = TRUE), .groups = "drop") %>%
    mutate(x = as.numeric(QType))

  p_time_by_k <- ggplot(qtime_by_k, aes(x = x_off, y = med_time, group = k)) +
    geom_line(alpha = 0.25, linewidth = 0.5, color = "gray40") +
    geom_point(alpha = 0.45, size = 1.6, color = "gray20") +
    geom_text(
      data = subset(qtime_by_k, QType == "W"),
      aes(label = k),
      size = 2.6,
      vjust = 0,
      nudge_y = 0.035 * diff(range(qtime_by_k$med_time, na.rm = TRUE))  # proportional nudge
    ) +
    geom_line(
      data = med_by_format,
      aes(x = x, y = med_time, group = 1),
      inherit.aes = FALSE,
      linewidth = 0.8,
      linetype = "dashed",
      color = "gray30"
    ) +
    geom_point(
      data = med_by_format,
      aes(x = x, y = med_time),
      inherit.aes = FALSE,
      size = 2.0,
      shape = 21,
      stroke = 0.6,
      color = "gray30",
      fill = "white"
    ) +
    scale_x_continuous(breaks = 1:3, labels = x_labs_math) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +  # TOP headroom for k labels
    labs(x = NULL, y = expression(italic(QTIME)~"[s]")) +
    coord_cartesian(clip = "off") +
    theme_col_box

  # ============================================================
  # Panel (c): SAI distribution (participant-level)
  # ============================================================

  sai_df <- Qlevel %>%
    group_by(ParticipantID) %>%
    summarise(SAI = first(SAI), .groups = "drop") %>%
    filter(is.finite(SAI))

  n_sai <- nrow(sai_df)

  p_sai <- ggplot(sai_df, aes(x = "", y = SAI)) +
    geom_violin(fill = "gray80", color = "black", linewidth = 0.4, trim = FALSE) +
    labs(x = NULL, y = expression(italic(SAI))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    coord_cartesian(clip = "off") +
    theme_col_box +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    annotate(
      "text",
      x = 1, y = Inf,
      label = sprintf("italic(N)==%d", n_sai),
      parse = TRUE,
      vjust = 1.6, hjust = 1.1, size = 3
    )

  # ============================================================
  # Arrange (a,b,c) horizontally for full-page width (two-column figure)
  # (labels placed safely inside so they don’t get clipped)
  # ============================================================

  fig_abc <- cowplot::plot_grid(
    p_time_by_k, p_grade_by_k, p_sai,
    nrow = 1,
    labels = c("a", "b", "c"),
    label_size = 11,
    label_x = 0.02,   # inside-left
    label_y = 0.98,   # inside-top
    hjust = 0, vjust = 1,
    align = "h",
    axis = "tb",
    rel_widths = c(1.05, 1.05, 1.0)
  )

  print(fig_abc)

  if (is_save == TRUE) {
    out_pdf <- file.path(plot_dir, "Figure8.pdf")
    out_png <- file.path(plot_dir, "Figure8.png")

    # IEEE-like full-page width for two-column figure ~7.16 in
    ggsave(out_pdf, fig_abc, width = 7.16, height = 3.2, units = "in", device   = cairo_pdf,
           bg       = "white"
    )
    ggsave(out_png, fig_abc, width = 7.16, height = 3.2, units = "in", dpi = 300)
  }
  # =============================================================================
  # Figure 9a — Facially displayed emotions (stacked proportions by participant × format)
  # =============================================================================
  # Paper: Fig. 9a (stacked proportions of the seven CNN emotions displayed on faces)
  # Uses question-level emotion probabilities (FACS.Qlevel) computed above.

  if (exists("FACS.Qlevel")) {

    fac_cols <- c("Angry","Disgusted","Afraid","Happy","Sad","Surprised","Neutral")
    needed   <- c("ParticipantID","Question.Type", fac_cols)

    if (!all(needed %in% names(FACS.Qlevel))) {
      warning("Figure 9a not generated: FACS.Qlevel is missing required columns: ",
              paste(setdiff(needed, names(FACS.Qlevel)), collapse = ", "))
    } else {

      # Long-format data for stacked bars (matches the original Figure 9a layout)
      long_data <- FACS.Qlevel %>%
        dplyr::select(dplyr::all_of(needed)) %>%
        tidyr::pivot_longer(
          cols      = dplyr::all_of(fac_cols),
          names_to  = "Emotion",
          values_to = "Probability"
        ) %>%
        dplyr::mutate(
          Question.Type  = factor(Question.Type, levels = c("V", "A", "W")),
          Emotion        = factor(Emotion, levels = fac_cols),
          ParticipantID  = factor(ParticipantID)
        )

      strip_labels <- c(
        "V" = "VIDEO",
        "A" = "ABSTRACT",
        "W" = "WORD"
      )

      ## --- Build Figure 9a (paper-matching aesthetics) ---

      # Enforce canonical Participant order S001..S050 (prevents accidental reorder)
      pid_levels <- sprintf("S%03d", 1:200)
      pid_levels <- pid_levels[pid_levels %in% unique(as.character(long_data$ParticipantID))]
      long_data <- long_data %>%
        dplyr::mutate(ParticipantID = factor(as.character(ParticipantID), levels = pid_levels))

      # Compact, paper-like theme for this figure only (do not rely on global theme_set)
      theme_fig9a <- theme_classic(base_size = 9) +
        theme(
          legend.position      = "top",
          legend.direction     = "horizontal",
          legend.title         = element_blank(),
          legend.text          = element_text(size = 8, face = "bold"),
          legend.key.size      = unit(0.32, "cm"),
          legend.margin        = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin    = margin(0, 0, 2, 0, unit = "pt"),
          axis.title           = element_blank(),
          axis.text.x          = element_text(size = 8, angle = 75, hjust = 1, vjust = 1),
          axis.text.y          = element_text(size = 8),
          axis.ticks.length    = unit(1.5, "pt"),
          panel.grid           = element_blank(),
          panel.spacing        = unit(0.10, "lines"),
          strip.placement      = "outside",
          strip.background     = element_rect(fill = "white", color = "black", linewidth = 0.35),
          strip.text.y.right   = element_text(size = 9, face = "bold"),
          plot.margin          = margin(2, 2, 2, 2, unit = "pt")
        )

      Fig9a <- ggplot(long_data, aes(x = ParticipantID, y = Probability, fill = Emotion)) +
        geom_col(position = "fill", width = 0.95) +
        facet_wrap(
          ~ Question.Type,
          nrow = 3,
          strip.position = "right",
          labeller = as_labeller(strip_labels)
        ) +
        scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00), labels = scales::label_number(accuracy = 0.01)) +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        scale_fill_manual(
          values = c(
            "Angry"     = "red",
            "Disgusted" = "green",
            "Afraid"    = "orange",
            "Happy"     = "yellow",
            "Sad"       = "blue",
            "Surprised" = "purple",
            "Neutral"   = "grey"
          ),
          breaks = fac_cols,
          labels = c("ANGRY", "DISGUSTED", "AFRAID", "HAPPY", "SAD", "SURPRISED", "NEUTRAL"),
          guide  = guide_legend(nrow = 1)
        ) +
        theme_fig9a

      print(Fig9a)

      if (isTRUE(is_save)) {
        # Match paper naming convention
        ggsave(
          filename = file.path(plot_dir, "Figure9a.pdf"),
          plot     = Fig9a,
          width    = 7.2,
          height   = 4.2,
          units    = "in",
          device   = cairo_pdf,
          bg       = "white"
        )
        ggsave(
          filename = file.path(plot_dir, "Figure9a.png"),
          plot     = Fig9a,
          width    = 7.2,
          height   = 4.2,
          units    = "in",
          dpi      = 300,
          bg       = "white"
        )
        message("Saved Figure9a to: ", normalizePath(file.path(plot_dir, "Figure9a.pdf"), winslash = "/", mustWork = FALSE))
      }
    }
  } else {
    warning("Figure 9a not generated: FACS.Qlevel object not found in this script run.")
  }

}

# Run when executed as a script (not when the file is imported by another script)
if (identical(environment(), globalenv())) {
  main()
}