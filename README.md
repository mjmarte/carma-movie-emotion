# Altered emotional reactivity in aphasia during naturalistic movie-viewing

Analysis code for the paper.

People with post-stroke aphasia and healthy controls watched eight short movie clips, rated emotional valence continuously,
rated discrete emotions after each clip (PANAS-SF), and answered comprehension and antonym generation items. The code computes
the emotion and language features, fits the group models, runs the principal component analysis, and evaluates the
classification (SVM) and severity prediction (SVR) models with the reference-dependent features recomputed within each fold.

## System requirements

- R 4.4.2 with the packages in `r-packages.txt` (the versions tested are listed beside each package).
- Tested on macOS 26.6 (Apple silicon). No non-standard hardware is required.

## Installation

```r
install.packages(sub("\\s*#.*", "", grep("^[^#]", readLines("r-packages.txt"), value = TRUE)))
```

To install the exact versions tested, use `remotes::install_version()` with the versions in `r-packages.txt`.

## Usage

1. Copy `config.example.json` to `config.json` and set the path of each input and of the output folder. `config.json` is not
   tracked, and every script refers to its inputs only by key (`R/config.R`).
2. From the repository root, run `Rscript run_all.R`, which runs the scripts in `src/` in numbered order (about 6 minutes;
   the SVM takes most of it). Each script writes its results to the output folder; the figure scripts also write vector PDFs (`figure1.pdf` to
   `figure7.pdf`) next to the PNGs.

## Reproducing the paper

| Script | In the paper |
|---|---|
| `01_realtime_ratings.R` | Analytic sample (57 PWA, 43 HC); real-time valence ratings |
| `02_holistic_msz.R` | Holistic MSZ (Methods 4.4) |
| `03_isc.R` | Time-resolved ISC (Methods 4.4) |
| `04_sample_entropy.R` | Sample entropy (Methods 4.4) |
| `05_task_accuracy.R` | Comprehension and antonym accuracy and their models (Results 2.3) |
| `06_realtime_msz.R` | Real-time MSZ (Methods 4.4), and the version referenced to all participants used in Results 2.3 and 2.4 |
| `07_feature_correlations.R` | Feature selection (Supplementary text); Supplementary Figure S2 |
| `08_pca.R` | Principal component analysis (Results 2.4); Figure 4 |
| `09_pca_validation.R` | Bootstrap stability of the components (Results 2.4) |
| `10_group_models.R` | Group effects on the emotion features (Results 2.3) and on the components (Results 2.4) |
| `11_demographics.R` | Table 1; Results 2.1 |
| `12_figure5.R` | PC1 correlations with AQ, PHQ-8 and GAD-7 (Results 2.5); Figure 5 |
| `13_figures_1_2.R` | Figures 1 and 2 |
| `14_table2.R` | Table 2; Results 2.2 |
| `16_supplementary.R` | Holistic MSZ with an HC-only reference (Supplementary text); inputs to Supplementary Table S2; education, sex, race and mood by group (Results 2.1, Table 1); joint models of PC1 on severity and mood |
| `17_svm.R` | Classification: Table 3, Results 2.6 |
| `18_svr.R` | Severity prediction: Results 2.6 |
| `19_figures_6_7.R` | Figures 6 and 7 |
| `20_sensitivity_analyses.R` | Comprehension sensitivity analyses (Results 2.7); Supplementary Table S1 |
| `21_familiarity.R` | Prior familiarity with the clips; Supplementary Table S3 |
| `22_figure_s1.R` | Supplementary Figure S1 |
| `23_table_s2.R` | Supplementary Table S2 |
| `24_input_method.R` | Response-method check (Supplementary text) |
| `25_outlier_checks.R` | Influential-observation checks (Methods 4.6) |
| `26_figure3.R` | Figure 3 |
| `27_reporting_statistics.R` | Degrees of freedom, 95% CIs, effect sizes and exact p values reported in the paper and Supplementary Tables |

`R/foldwise_features.R` rebuilds the reference-dependent features (MSZ, ISC) from the training participants of each fold for the
SVM and SVR; `R/theme_marte.R` is the plotting theme.

## Data

Participant data are not included. Requests go to the corresponding author.

## Citation

Marte, M. J., Gillis, B., Wei, X., Galvin, C., Rigolo, L., Tie, Y., Kiran, S., & Liebenthal, E. Altered emotional reactivity in
aphasia during naturalistic movie-viewing.

Code archived at Zenodo: https://doi.org/10.5281/zenodo.23325284
