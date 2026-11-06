#!/usr/bin/env Rscript
# ==============================================================================
# A24 - WHICH AGE IS "AGE AT BASELINE"? (aggregate check, JCF-D-26-00377)
# ------------------------------------------------------------------------------
# Table 1 (A05.2) takes age_years from each patient's first record; Table S1
# (A21) and the models take baseline_age_years. Their means differ slightly
# (18.6 vs about 18.4). This compares the two on the same baseline row.
# Prints aggregates only. Writes: outputs/revision/row24_age_check.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow) })
dat <- arrow::read_parquet("outputs/data_derived/analysis_long.parquet") %>% as_tibble() %>% ungroup()
bl <- as.data.frame(dat, stringsAsFactors = FALSE)
bl <- bl[order(bl$patient_id, bl$visit_date, na.last = TRUE), , drop = FALSE]
bl <- bl[!duplicated(bl$patient_id), , drop = FALSE] %>% as_tibble()
d <- bl %>% transmute(age_years = as.numeric(age_years), baseline_age_years = as.numeric(baseline_age_years),
                      diff = age_years - baseline_age_years,
                      visit_year = as.integer(format(as.Date(visit_date), "%Y")),
                      same_date = as.Date(visit_date) == as.Date(baseline_date))
out <- tibble(
  patients = nrow(d),
  mean_age_years = mean(d$age_years, na.rm = TRUE), sd_age_years = sd(d$age_years, na.rm = TRUE),
  mean_baseline_age_years = mean(d$baseline_age_years, na.rm = TRUE), sd_baseline_age_years = sd(d$baseline_age_years, na.rm = TRUE),
  n_missing_age_years = sum(is.na(d$age_years)), n_missing_baseline_age = sum(is.na(d$baseline_age_years)),
  mean_diff = mean(d$diff, na.rm = TRUE), median_diff = median(d$diff, na.rm = TRUE),
  n_absdiff_gt_0.1 = sum(abs(d$diff) > 0.1, na.rm = TRUE), n_absdiff_gt_0.5 = sum(abs(d$diff) > 0.5, na.rm = TRUE),
  n_absdiff_gt_1 = sum(abs(d$diff) > 1, na.rm = TRUE), max_absdiff = max(abs(d$diff), na.rm = TRUE),
  first_row_is_baseline_date = sum(d$same_date, na.rm = TRUE))
print(as.data.frame(t(out)))
by_year <- d %>% group_by(visit_year) %>% summarise(n = n(), mean_diff = mean(diff, na.rm = TRUE), n_gt_0.5 = sum(abs(diff) > 0.5, na.rm = TRUE), .groups = "drop")
print(as.data.frame(by_year))
write_csv(out, "outputs/revision/row24_age_check.csv")
