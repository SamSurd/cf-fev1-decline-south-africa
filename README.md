# Analysis code: Determinants of FEV₁ decline in cystic fibrosis in South Africa

R code for every analysis in the revised manuscript and its supplementary material:

> Surdut SP, Verstraete J, Zampoli M, Calligaro G. Determinants of FEV₁ decline in cystic fibrosis in South Africa: a longitudinal cohort. *Journal of Cystic Fibrosis* (manuscript JCF-D-26-00377, revised November 2026).

The scripts read the South African Cystic Fibrosis Registry (SACFR) export for 2018–2024, build the analytic cohort (459 patients, 2,292 annual FEV₁ percent predicted values), fit the linear mixed-effects and event-anchored models, and write every table, figure and estimate reported.

## Data

**No data are included.** The registry data are held by the SACFR and are available on reasonable request to the registry steering committee, subject to ethics approval. To run the code, place the registry export at `data_raw/CF_registry_master.xlsx` (sheet `Sheet1`: one row per patient, with yearly columns such as `fev1_percent_pred_local_ref_2018`). See `data_raw/README.md`.

The scripts write patient-level derived data to `outputs/data_derived/`. Treat that folder like the raw data and do not share it.

## Requirements

- R 4.5.2. Package versions used for the published results are in `R_session_info.txt`: tidyverse 2.0.0 (dplyr 1.2.0, ggplot2 4.0.2, tidyr 1.3.2), lme4 1.1.38, lmerTest 3.2.0, performance 0.16.0, broom.mixed 0.2.9.6, arrow 22.0.0.1, readxl 1.4.5, janitor 2.2.1, gt 1.3.0, ragg 1.5.0.
- pandoc (3.6.3 was used), needed by gt to write the `.docx` tables. If pandoc is not on the PATH, set `RSTUDIO_PANDOC` to the folder that contains it; `run_all.sh` does this automatically when RStudio is installed on macOS.

Several scripts install missing CRAN packages on first run.

## How to run

From this folder:

```
./run_all.sh
```

This runs the scripts in the order below (about a minute on a laptop), writes a log for each to `outputs/logs/`, and stops at the first failure. Each script can also be run on its own with `Rscript <script>` once the scripts before it have run. Random seeds are set in each script, so the same data and package versions reproduce the reported estimates exactly.

## Scripts

| Order | Script | What it does | Used for |
|---|---|---|---|
| 1 | `A01.R` | Reads the registry export, builds the long analytic dataset and derived variables, applies the inclusion criteria and writes `outputs/data_derived/analysis_long.parquet` | All analyses |
| 2 | `A02.R` | Baseline-time models M0 (unadjusted), M1 (main effects) and M2 (with slope modifiers); random intercepts and slopes, with uncorrelated and then intercept-only fallbacks | Results 3.2–3.3; Tables S2–S3 |
| 3 | `Supplementary Scripts/Build_Exac_Anchored_Tidies.R` | Exacerbation event models as specified in the original submission (unadjusted) | Table S6 |
| 4 | `A03.R` | Model tables: fixed and random effects, R², diagnostics | Model output |
| 5 | `A03.1.R` | Tables 2A and 2B in the submitted layout (inputs to Figure 1) | Figure 1 |
| 6 | `A04.R` | Baseline characteristics by modulator exposure | Results 3.1 |
| 7 | `A05.2.R` | Table 1: baseline characteristics by healthcare sector | Table 1 |
| 8 | `A08_revision_fixes.R` | CFTR modulator models C0–C6 (C2 primary), population-averaged annual change, within-model post-event slopes, chronic infection models E0–E3 (E1 primary) and event classes | Results 3.2–3.4; Figure 1; Tables S3–S5 |
| 9 | `A09_revision_outputs.R` | Infection models censored at modulator initiation, follow-up density, cohort flow, respiratory samples, 3120+1G>A and eligibility, Benjamini–Hochberg p-values | Results 3.1; Tables S2–S4, S9; Figure S1 |
| 10 | `A10_fev1_reference_check.R` | Source of each FEV₁ value (GLI or local reference) and models with a source indicator | Table S7 |
| 11 | `A12_revision_sensitivity.R` | Patients with denser follow-up, R², CFTR confounding analyses, co-infection, age at acquisition | Tables S4, S5, S7, S8 |
| 12 | `A13_revision_figures.R` | Figure 2 and the cohort flowchart | Figure 2; Figure S1 |
| 13 | `A14_sampling_sensitivity.R` | Infection models in well-sampled patients | Table S4 |
| 14 | `A15_exacerbations.R` | Exacerbation coding, previous-year models and adjusted event models | Results 3.2–3.3; Table S6; Figure 1 |
| 15 | `A16_figure1.R` | Figure 1 | Figure 1 |
| 16 | `A18_chronic_status_unknown.R` | Chronic infection status recorded per registry year | Results 3.1; Table S9 |
| 17 | `A21_revision_consistency.R` | R² and intraclass correlation, GLI vs local reference, modulator exposure by commencement date | Table 1; Tables S1, S7 |
| 18 | `A22_figureS2_infections.R` | FEV₁ around onset of incident chronic infection | Figure S2 |
| 19 | `A23_table1_extra.R` | Chronic *P. aeruginosa* at baseline, by sector | Table 1 |
| 20 | `A19_column_audit.R` | Quality check: audit of every year-specific registry column read by `A01.R` | Data checks |
| 21 | `A24_age_check.R` | Quality check: definitions of age at baseline | Data checks |

Outputs are written to `outputs/tables/`, `outputs/models/` and `outputs/revision/` (with figures in `outputs/revision/figures/`).

## Corrections made during peer review

Three errors in how `A01.R` read yearly registry fields were corrected on 5 November 2026, and the corrected lines are marked `revision (5 Nov 2026)`: the 2023 local-reference FEV₁ column, the 2019 *P. aeruginosa* and MSSA status columns and the 2018–2019 MRSA, *B. cepacia* complex and diabetes columns, and carrying CFRD and pancreatic insufficiency forward once documented. `A08_revision_fixes.R` documents the re-specified models in its header. All reported results come from the corrected code.

## Licence

MIT (see `LICENSE`).

## Contact

Samuel P. Surdut, Department of Medicine, University of Cape Town (srdsam001@myuct.ac.za, medicine@surdut.com).
