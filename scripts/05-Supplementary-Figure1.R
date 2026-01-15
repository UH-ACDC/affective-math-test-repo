#!/usr/bin/env Rscript
# =============================================================================
# Affective Math Test — Supplementary Figure 1 (All Participants)
#
# Author: Ioannis Pavlidis
# Affiliation: Affective and Data Computing Lab — University of Houston
#
# Repository script: scripts/05-Supplementary-Figure1R
#
# Creates a 50-page PDF (1 page per participant) that matches the look/logic of
# Figure 6 in 00-Exploratory-Data-Analysis_Fig9a_v5.R, but iterates over *all*
# ParticipantID values in the Q-level dataset.
#
# Expected repo structure (run from ./scripts):
#   scripts/02-Supplementary-Figure1_AllParticipants_Fig6Style.R   (this file)
#   scripts/fm_plot_func.R
#   data/processed/Affective_Math_Qlevel_Data_N50.csv
#   data/processed/Affective_Math_Dataset_N50_BL.csv
#
# Output:
#   figures/Supplementary_Figure1.pdf
# =============================================================================

# ---- Robust script dir (GitHub-safe) ----
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(p))
  }
  # Fallback for Rscript
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path_idx <- grep(file_arg, cmd_args)
  if (length(path_idx) > 0) {
    p <- sub(file_arg, "", cmd_args[path_idx[1]])
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  return(getwd())
}

script_dir <- get_script_dir()
setwd(script_dir)
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

# ---- Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

# ---- Theme helpers ----
# Prefer the project's shared plotting helper if available
if (file.exists(file.path(script_dir, "fm_plot_func.R"))) {
  source(file.path(script_dir, "fm_plot_func.R"))
} else {
  theme_set(theme_classic())
  my_theme1 <- theme_classic() +
    theme(
      text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11)
    )
}

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
      plot.title.position = "panel",
      plot.margin = margin(t = top_margin, r = 6, b = 2, l = 6)
    )
}

# Panel tag (a/b) in the top-left (in plot margin area)
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

# ---- Config ----
N <- 50
DATA_DIR <- file.path(repo_root, "data", "processed")
FIG_DIR  <- file.path(repo_root, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

in_qlevel <- file.path(DATA_DIR, sprintf("Affective_Math_Qlevel_Data_N%s.csv", N))
in_bl     <- file.path(DATA_DIR, sprintf("Affective_Math_Dataset_N%s_BL.csv", N))

stopifnot(file.exists(in_qlevel))
stopifnot(file.exists(in_bl))

# ---- Load data ----
Qlevel <- read.csv(in_qlevel, stringsAsFactors = FALSE)
Df_BL  <- read.csv(in_bl, stringsAsFactors = FALSE)

# Baseline means (same logic as main script)
Df_BL.SubSet <- Df_BL %>%
  dplyr::select(ParticipantID, Time, Timestamp, Perspiration, HR.E4, HR.AW, HRV.IBI)

idx <- is.na(Df_BL.SubSet$HR.E4) | is.na(Df_BL.SubSet$HR.AW)
Df_BL.SubSet$HR.E4[idx] <- NA
Df_BL.SubSet$HR.AW[idx] <- NA

Df_BL.mean <- Df_BL.SubSet %>%
  dplyr::group_by(ParticipantID) %>%
  dplyr::summarise(
    Perspiration = mean(Perspiration, na.rm = TRUE),
    HR.E4        = mean(HR.E4,        na.rm = TRUE),
    HR.AW        = mean(HR.AW,        na.rm = TRUE),
    HRV.IBI      = mean(HRV.IBI,      na.rm = TRUE),
    .groups = "drop"
  )

# ---- Aesthetics (match Supplementary Figure 1) ----
# Word/Video/Abstract colors
qtype_cols <- c(
  "W" = "#9A7BC0",  # Word (purple)
  "V" = "#F08DB4",  # Video (pink)
  "A" = "#A7DB96"   # Abstract (green)
)

# Legend (rectangles) at bottom
make_legend_row <- function(pid) {
  legend_df <- data.frame(QType = factor(c("W","V","A"), levels = c("W","V","A")))

  p_leg <- ggplot(legend_df, aes(x = QType, y = 1, fill = QType)) +
    geom_tile(width = 0.9, height = 0.9, color = NA) +
    scale_fill_manual(
      values = qtype_cols,
      breaks = c("W","V","A"),
      labels = c("WORD","VIDEO","ABSTRACT")
    ) +
    guides(fill = guide_legend(title = NULL, nrow = 1, byrow = TRUE,
                              override.aes = list(alpha = 1))) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(face = "bold", size = 10)
    )

  leg_grob <- cowplot::get_legend(p_leg)

  pid_grob <- cowplot::ggdraw() +
    cowplot::draw_label(pid, fontface = "bold", size = 11, x = 2.50, hjust = 1)
  
  cowplot::plot_grid(pid_grob, leg_grob, nrow = 1, rel_widths = c(0.13, 0.87))
}

# Background rectangles: 1 per question (full height)
make_rect_df <- function(df_pid) {
  df_pid %>%
    dplyr::mutate(
      QNumber = as.numeric(QNumber),
      QType   = as.character(QType)
    ) %>%
    dplyr::filter(is.finite(QNumber)) %>%
    dplyr::transmute(
      xmin  = QNumber - 0.5,
      xmax  = QNumber + 0.5,
      QType = factor(QType, levels = c("W","V","A"))
    )
}

# Build one page (Figure-6-like) for a given participant
make_page <- function(pid) {
  df_pid <- Qlevel %>%
    dplyr::filter(ParticipantID == pid) %>%
    dplyr::arrange(as.numeric(QNumber)) %>%
    dplyr::mutate(
      QNumber = as.numeric(QNumber),
      QType   = as.character(QType)
    )

  if (nrow(df_pid) == 0) return(NULL)

  rect_df <- make_rect_df(df_pid)

  # Baseline means from baseline segment
  base_fp   <- Df_BL.mean$Perspiration[match(pid, Df_BL.mean$ParticipantID)]
  base_hraw <- Df_BL.mean$HR.AW[match(pid, Df_BL.mean$ParticipantID)]
  base_hre4 <- Df_BL.mean$HR.E4[match(pid, Df_BL.mean$ParticipantID)]

  # Panel (a): FP
  n_fp <- sum(is.finite(df_pid$Perspiration))

  p_fp <- ggplot(df_pid, aes(x = QNumber, y = Perspiration)) +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = QType),
      inherit.aes = FALSE,
      alpha = 0.35,
      color = NA
    ) +
    geom_vline(
      xintercept = seq(min(df_pid$QNumber, na.rm = TRUE) - 0.5,
                       max(df_pid$QNumber, na.rm = TRUE) + 0.5,
                       by = 1),
      linewidth = 0.25,
      alpha = 0.20
    ) +
    geom_line(color = "gray20", linewidth = 0.55) +
    geom_hline(yintercept = base_fp, linetype = "dashed", color = "dodgerblue3", linewidth = 0.8) +
    scale_fill_manual(values = qtype_cols, drop = FALSE, guide = "none") +
    labs(x = "QNUMBER", y = NULL) +
    my_theme1 +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.margin  = margin(t = 18, r = 6, b = 2, l = 6)
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = sprintf("italic(n)[FP]==%d", n_fp),
      parse = TRUE,
      hjust = 1.08,
      vjust = 1.40,
      size = 4.2,
      fontface = "bold"
    )

  p_fp <- top_panel_label(p_fp, expression(italic(FP) * " [" * degree * "C"^2 * "]"), top_margin = 18)
  p_fp <- add_panel_tag(p_fp, "a", size = 12, x = 0.03, y = 0.98)

  # Panel (b): HR
  n_aw <- sum(is.finite(df_pid$HR.AW))
  n_e4 <- sum(is.finite(df_pid$HR.E4))

  yr_hr <- range(c(df_pid$HR.E4, df_pid$HR.AW), na.rm = TRUE)
  hr_breaks <- seq(floor(yr_hr[1] / 5) * 5, ceiling(yr_hr[2] / 5) * 5, by = 5)
  x_mid <- mean(range(df_pid$QNumber, na.rm = TRUE))
  
  # Build a compact label when counts match
  hr_n_label <- if (!is.na(n_aw) && !is.na(n_e4) && n_aw == n_e4) {
    sprintf("italic(n)[AW]==italic(n)[E4]~~'='~~%d", n_aw)   # n_AW = n_E4 = 36
  } else {
    sprintf("italic(n)[AW]==%d~~italic(n)[E4]==%d", n_aw, n_e4)
  }

  p_hr <- ggplot(df_pid, aes(x = QNumber)) +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = QType),
      inherit.aes = FALSE,
      alpha = 0.35,
      color = NA
    ) +
    geom_vline(
      xintercept = seq(min(df_pid$QNumber, na.rm = TRUE) - 0.5,
                       max(df_pid$QNumber, na.rm = TRUE) + 0.5,
                       by = 1),
      linewidth = 0.25,
      alpha = 0.20
    ) +
    geom_line(aes(y = HR.E4), color = "gray20", linewidth = 0.55) +
    geom_point(aes(y = HR.AW), color = "red3", size = 1.5, alpha = 0.9) +
    geom_hline(yintercept = base_hre4, linetype = "dashed", color = "dodgerblue3", linewidth = 0.8) +
    geom_hline(yintercept = base_hraw, linetype = "solid",  color = "red3", linewidth = 0.9) +
    scale_fill_manual(values = qtype_cols, drop = FALSE, guide = "none") +
    scale_y_continuous(breaks = hr_breaks) +
    labs(x = "QNUMBER", y = NULL) +
    my_theme1 +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      axis.title.x = element_text(margin = margin(t = 6)),
      plot.margin  = margin(t = 18, r = 6, b = 2, l = 6)
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = hr_n_label,
      parse = TRUE,
      hjust = 1.08,
      vjust = 1.40,
      size = 4.2,
      fontface = "bold"
    )

  p_hr <- top_panel_label(p_hr, expression(italic(HR)[AW] * " & " * italic(HR)[E4] * " [bpm]"), top_margin = 18)
  p_hr <- add_panel_tag(p_hr, "b", size = 12, x = 0.03, y = 0.98)

  # Combine with center gap + legend
  gap <- ggplot() + theme_void()

  main_row <- cowplot::plot_grid(
    p_fp,
    gap,
    p_hr,
    nrow = 1,
    rel_widths = c(1, 0.06, 1),
    align = "hv"
  )

  cowplot::plot_grid(
    main_row,
    make_legend_row(pid),
    ncol = 1,
    rel_heights = c(1, 0.18)
  )
}

# ---- Render multi-page PDF ----
all_pids <- sort(unique(Qlevel$ParticipantID))

out_pdf <- file.path(FIG_DIR, "Supplementary_Figure1.pdf")

grDevices::pdf(out_pdf, width = 7.0, height = 3.0, useDingbats = FALSE)
for (pid in all_pids) {
  p <- make_page(pid)
  if (!is.null(p)) print(p)
}
# NOTE: dev.off() is provided by grDevices (not graphics)
grDevices::dev.off()

message("Wrote: ", out_pdf)
