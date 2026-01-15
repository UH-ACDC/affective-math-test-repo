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
