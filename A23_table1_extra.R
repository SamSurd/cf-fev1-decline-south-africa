#!/usr/bin/env Rscript
# ==============================================================================
# A23 - EXTRA TABLE 1 ROW: CHRONIC P. AERUGINOSA AT BASELINE, BY SECTOR
# ------------------------------------------------------------------------------
# Table 1 (A05.2) shows chronic infections as ever recorded during follow-up;
# Results 3.1 quotes chronic P. aeruginosa at baseline (150/459, 32.7%, from
# A05 Table1_Baseline_Characteristics). This adds that row by sector, using the
# same baseline row (first record per patient by visit date), sector grouping
# and tests as A05.2 (chi-square without continuity correction, Fisher's exact
# if any expected count < 5). Aggregate output only.
# Writes: outputs/revision/row23_table1_baseline_pa.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow) })
dat <- arrow::read_parquet("outputs/data_derived/analysis_long.parquet") %>% as_tibble() %>% ungroup()
bl <- as.data.frame(dat, stringsAsFactors = FALSE)
bl <- bl[order(bl$patient_id, bl$visit_date, na.last = TRUE), , drop = FALSE]
bl <- bl[!duplicated(bl$patient_id), , drop = FALSE] %>% as_tibble()
sec <- tolower(trimws(as.character(bl$sector)))
bl$sector_group <- case_when(is.na(sec) ~ NA_character_,
                             str_detect(sec, "public") ~ "Public Only",
                             str_detect(sec, "private|combined") ~ "Private/Combined",
                             TRUE ~ NA_character_)
bl <- bl %>% filter(!is.na(sector_group)) %>%
  mutate(pa_cat = as.character(baseline_pa_cat), pa_cat = if_else(is.na(pa_cat), "Negative", pa_cat),
         pa_bl = pa_cat == "Chronic")
n_pct <- function(k, n) sprintf("%d (%.1f%%)", as.integer(k), 100 * k / n)
fmt_p <- function(p) if (!is.finite(p)) "-" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
p_cat <- function(x, g) {
  tab <- table(x, g)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (any(chi$expected < 5)) fisher.test(tab)$p.value else chi$p.value
}
pu <- bl$sector_group == "Public Only"; pr <- bl$sector_group == "Private/Combined"
out <- tibble(domain = "Chronic infections", variable = "P. aeruginosa at baseline", level = "Chronic",
              Public_Only = n_pct(sum(bl$pa_bl[pu]), sum(pu)),
              Private_Combined = n_pct(sum(bl$pa_bl[pr]), sum(pr)),
              Overall = n_pct(sum(bl$pa_bl), nrow(bl)),
              p_value = fmt_p(p_cat(bl$pa_bl, bl$sector_group)),
              p_value_3cat = fmt_p(p_cat(bl$pa_cat, bl$sector_group)))
lv <- bl %>% count(sector_group, pa_cat) %>% group_by(sector_group) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
print(as.data.frame(out)); print(as.data.frame(lv))
write_csv(out, "outputs/revision/row23_table1_baseline_pa.csv")
