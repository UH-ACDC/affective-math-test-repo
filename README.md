# affective-math-test-repo

Companion repository for our paper on story-based math problems (Video / Abstract / Word formats) and affective responses during test-taking.  
It contains reproducible analysis scripts and the processed dataset used to generate the paper’s results and supplementary material.

## Repository structure

- `scripts/` — runnable analysis scripts (paper + supplementary)
- `data/processed/` — processed input data used by the scripts
- `figures/` — generated figures (created by scripts if missing)
- `reports/` — generated tables / model outputs (created by scripts if missing)

## Requirements

- R (recent R 4.x recommended)
- Typical packages used across scripts include: `tidyverse`, `lme4`, `ggplot2`, `cowplot`, `emmeans` (and related dependencies).

## How to run

From the repository root:

```bash
Rscript scripts/00-Exploratory-Data-Analysis.R
Rscript scripts/01-Continuous-StressLabel-Models.R
Rscript scripts/02-Stress-Label-Models-q33.R
Rscript scripts/03-Valence-Performance-Models.R
Rscript scripts/04-VideoStory-Calmness-Scoring.R
Rscript scripts/05-Supplementary-Figure1.R

## Outputs

Outputs are written to:

- `figures/` (plots and figure PDFs)
- `reports/` (tables, model summaries, intermediate artifacts)

## Recommended run order

Most users can run scripts independently; however, for a “full reproduction” workflow we recommend:

1. `00-Exploratory-Data-Analysis.R` (descriptives / core figures)
2. `01-Continuous-StressLabel-Models.R` (continuous + stress-label modeling outputs)
3. `02-Stress-Label-Models-q33.R` (q33/top-33% stress-label variants)
4. `03-Valence-Performance-Models.R` (valence + performance models)
5. `04-VideoStory-Calmness-Scoring.R` (calmness scoring analyses)
6. `05-Supplementary-Figure1.R` (supplementary figure reproduction)

## Script index (what each script does)

- `scripts/00-Exploratory-Data-Analysis.R`  
  **Purpose:** Loads the processed datasets, performs descriptive summaries, and generates core exploratory figures.  
  **Inputs:** `data/processed/*.csv` (see `data/processed/`).  
  **Outputs:** Figures to `figures/` and supporting summary artifacts to `reports/`.

- `scripts/01-Continuous-StressLabel-Models.R`  
  **Purpose:** Fits mixed-effects models for continuous physiology-derived stress signals and participant-relative stress labels.  
  **Inputs:** Processed datasets in `data/processed/`.  
  **Outputs:** Model tables and diagnostics to `reports/`; associated plots/figures to `figures/`.

- `scripts/02-Stress-Label-Models-q33.R`  
  **Purpose:** Stress-label mixed-effects models using the top-33% (q33) criterion (variant analysis).  
  **Inputs:** Processed datasets in `data/processed/`.  
  **Outputs:** q33-specific model tables/exports to `reports/` (and figures to `figures/` if enabled).

- `scripts/03-Valence-Performance-Models.R`  
  **Purpose:** Mixed-effects models linking affective/behavioral variables (including valence and timing measures) to performance outcomes.  
  **Inputs:** Processed datasets in `data/processed/`.  
  **Outputs:** Model tables to `reports/` and figures to `figures/`.

- `scripts/04-VideoStory-Calmness-Scoring.R`  
  **Purpose:** Computes and analyzes calmness scoring for the video/story stimuli and produces the associated outputs.  
  **Inputs:** Processed datasets in `data/processed/`.  
  **Outputs:** Exports to `reports/` (and figures to `figures/` if enabled).

- `scripts/05-Supplementary-Figure1.R`  
  **Purpose:** Reproduces Supplementary Figure 1 from the paper.  
  **Inputs:** Processed datasets in `data/processed/`.  
  **Outputs:** `figures/Supplementary_Figure1.pdf` (or similarly named figure file produced by the script).

- `scripts/fm_plot_func.R`  
  **Purpose:** Helper plotting utilities sourced by one or more scripts.

## Notes

- The scripts create `figures/` and `reports/` automatically if they do not exist.
- If you encounter a missing-package error, install the required package(s) in R and rerun the script.

## License

- **Code:** MIT License (see `LICENSE`)
- **Data:** CC BY 4.0 (attribution required) for contents of `data/`
