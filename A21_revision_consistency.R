#!/usr/bin/env Rscript
# ==============================================================================
# A21 - CONSISTENCY ITEMS FOR THE RESUBMISSION (corrected data, after fix 2)
# ------------------------------------------------------------------------------
# 1. M2 (REML): marginal and conditional R2 and the adjusted ICC from the same
#    fit (the text previously took the ICC from the ML fit in A02).
# 2. GLI 2012 vs local-reference % predicted for the same test, in modelled
#    patient-years where the registry holds both values.
# 3. CFTR modulator exposure by commencement date (the definition used by the
#    modulator models and Supplementary Methods S1), replacing the A01 flag used
#    by A05.2 (Table 1) and A09 (row 10): Table 1 row by sector and the
#    Supplementary Table S1 comparison of treated and untreated patients.
# Run after A08 (uses its helper functions). Aggregate outputs only.
# Writes: outputs/revision/A21_summary.txt, row22_r2_icc.csv, row22_gli_vs_local.csv,
#         row22_table1_modulator.csv, row22_baseline_by_modulator.csv
# ==============================================================================
set.seed(20251111)
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(lme4); library(lmerTest); library(readxl); library(janitor); library(performance)
})
out_dir <- "outputs/revision"
summary_file <- file.path(out_dir, "A21_summary.txt")
fit_log <- list()
for (e in parse("A08_revision_fixes.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, envir = globalenv())
}
file.create(summary_file)
say <- function(...) { msg <- paste0(...); message(msg); con <- file(summary_file, open = "a"); writeLines(msg, con); close(con) }
show_tbl <- function(x) { for (l in capture.output(print(as.data.frame(x), row.names = FALSE))) say(l) }
n_pct <- function(k, n) sprintf("%d (%.1f%%)", as.integer(k), 100 * k / n)

# Data preparation identical to A12 ---------------------------------------------
dat <- arrow::read_parquet("outputs/data_derived/analysis_long.parquet") %>% as_tibble() %>%
  mutate(
    sex           = factor(sex, levels = c("Female", "Male")),
    nutrition_cat = factor(nutrition_cat, levels = c("Normal", "Underweight", "Overweight", "Unknown")),
    exac_cat      = factor(exac_cat, levels = c("0", "1", "2+", "Unknown")),
    sector        = factor(sector, levels = c("Public", "Private", "Combined")),
    ethnicity_cat = factor(ethnicity_cat, levels = c("Black African", "Mixed Race", "Indian", "Asian", "Caucasian", "Other", "Unknown")),
    f508_genotype_group = factor(f508_genotype_group, levels = c("Other", "Heterozygous", "Homozygous")),
    ethnicity_caucasian = dplyr::if_else(ethnicity_cat == "Caucasian", 1L, 0L, missing = 0L),
    sector_public_only  = dplyr::if_else(sector == "Public", 1L, 0L, missing = 0L),
    f508_any = dplyr::if_else(f508_genotype_group %in% c("Heterozygous", "Homozygous"), 1L, 0L, missing = 0L))
dat_model <- dat %>% filter(nutrition_cat != "Unknown" | is.na(nutrition_cat)) %>% droplevels() %>%
  group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup()
pancreas_start <- dat_model %>% filter(pancreatic_insufficiency == 1L) %>% group_by(patient_id) %>%
  summarise(pancreas_start_date = suppressWarnings(min(visit_date, na.rm = TRUE)), .groups = "drop")
cfrd_start <- dat_model %>% filter(cfrd_indicator == 1L) %>% group_by(patient_id) %>%
  summarise(cfrd_start_date = suppressWarnings(min(visit_date, na.rm = TRUE)), .groups = "drop")
dat_model <- dat_model %>% left_join(pancreas_start, by = "patient_id") %>% left_join(cfrd_start, by = "patient_id")
flags_main <- compute_use_flags(dat_model)
say("A21 consistency items, run ", format(Sys.time(), "%Y-%m-%d %H:%M"), "; analytic cohort ", n_distinct(dat$patient_id),
    " patients; model data ", nrow(dat_model), " rows, ", n_distinct(dat_model$patient_id), " patients")

# 1. R2 and ICC from the same REML fit --------------------------------------------
m2_all <- fit_triple_fallback(build_rhs(flags_main, "m2"), dat_model, TRUE, "M2 all")
r2 <- performance::r2_nakagawa(m2_all)
ic <- tryCatch(performance::icc(m2_all), error = function(e) NULL)
icc_adj <- if (!is.null(ic)) ic$ICC_adjusted else NA_real_
icc_from_r2 <- (r2$R2_conditional - r2$R2_marginal) / (1 - r2$R2_marginal)
row22_r2 <- tibble(model = "M2 (REML)", marginal_R2 = r2$R2_marginal, conditional_R2 = r2$R2_conditional,
                   icc_adjusted = icc_adj, icc_from_r2 = icc_from_r2)
say(""); say("=== 1. M2 (REML) R2 and adjusted ICC ===")
show_tbl(row22_r2 %>% mutate(across(where(is.numeric), ~round(.x, 4))))
write_csv(row22_r2, file.path(out_dir, "row22_r2_icc.csv"))

# 2. GLI vs local-reference % predicted for the same test ------------------------
s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_df <- s1[grep("^patcode", names(s1), value = TRUE)]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
num <- function(v) if (v %in% names(s1)) suppressWarnings(as.numeric(s1[[v]])) else rep(NA_real_, nrow(s1))
src <- bind_rows(lapply(2018:2024, function(y) {
  gli <- dplyr::coalesce(num(paste0("gli_fev1_pp_", y)), num(paste0("gli_percpredfev", y)))
  loc <- dplyr::coalesce(num(paste0("fev1_percent_pred_local_ref_", y)), num(paste0("fev1_pp_", y)),
                         num(paste0("fev1_local_", y)), num(paste0("x", y, "_fev1_percent_pred_local_ref")))
  tibble(patient_id = pid, visit_year = y, gli = gli, local = loc)
})) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, visit_year, .keep_all = TRUE)
both <- dat_model %>% mutate(visit_year = as.integer(visit_year)) %>% select(patient_id, visit_year) %>%
  left_join(src, by = c("patient_id", "visit_year")) %>% filter(!is.na(gli), !is.na(local)) %>% mutate(d = gli - local)
gl <- both %>% group_by(visit_year) %>%
  summarise(n_both = n(), mean_gli_minus_local = mean(d), sd = sd(d), median = median(d),
            q1 = quantile(d, 0.25), q3 = quantile(d, 0.75), .groups = "drop") %>%
  bind_rows(both %>% summarise(visit_year = NA_integer_, n_both = n(), mean_gli_minus_local = mean(d), sd = sd(d),
                               median = median(d), q1 = quantile(d, 0.25), q3 = quantile(d, 0.75)))
say(""); say("=== 2. GLI minus local-reference % predicted, same test (modelled patient-years with both; last row = all years) ===")
show_tbl(gl %>% mutate(across(where(is.double), ~round(.x, 2))))
write_csv(gl, file.path(out_dir, "row22_gli_vs_local.csv"))

# 3. Modulator exposure by commencement date (rule as A08/A12/A17) ----------------
to_date2 <- function(x) { if (inherits(x, c("POSIXct", "POSIXt", "Date"))) return(as.Date(x)); suppressWarnings(lubridate::dmy(as.character(x), quiet = TRUE)) }
dcols <- intersect(c("date_cftrm2018", "date_cftrm2019", "date_cftrm2020", "date_cftrm2021", "date_cftrm2022", "date_cftrm2023",
                     "date_cftrm_commenced_2024", "treatment_start_date_cftr", "treatment_start_date_cftr_new"), names(s1))
dmat <- sapply(dcols, function(v) { d <- to_date2(s1[[v]]); d[!is.na(d) & (d < as.Date("2010-01-01") | d > as.Date("2024-12-31"))] <- NA; as.numeric(d) })
en <- suppressWarnings(apply(dmat, 1, min, na.rm = TRUE)); en[!is.finite(en)] <- NA
earliest <- as.Date(en, origin = "1970-01-01")
fcols <- grep("^cftr_modifier_?20[0-9]{2}$", names(s1), value = TRUE); fy <- as.integer(stringr::str_extract(fcols, "20[0-9]{2}"))
Yf <- sapply(fcols, function(v) as.character(s1[[v]]) %in% "Yes")
fyes <- apply(Yf, 1, function(r) if (any(r)) min(fy[r]) else NA_integer_)
imp <- !is.na(fyes) & (is.na(earliest) | as.integer(format(earliest, "%Y")) > fyes)
srev <- earliest; srev[imp] <- as.Date(paste0(fyes[imp], "-07-01"))
starts <- tibble(patient_id = pid, cftr_start_rev = srev) %>% filter(!is.na(patient_id)) %>%
  group_by(patient_id) %>% summarise(cftr_start_rev = suppressWarnings(min(cftr_start_rev, na.rm = TRUE)), .groups = "drop") %>%
  mutate(cftr_start_rev = dplyr::if_else(is.finite(as.numeric(cftr_start_rev)), cftr_start_rev, as.Date(NA)))
chr <- function(v) { x <- if (v %in% names(s1)) as.character(s1[[v]]) else rep(NA_character_, nrow(s1)); x[is.na(x)] <- ""; x }
m1 <- paste(chr("mutation1"), chr("othermutation1")); m2 <- paste(chr("mutation2"), chr("othermutation2"))
gen <- tibble(patient_id = pid, carrier_3120 = grepl("3120", m1) | grepl("3120", m2),
              eligible = dplyr::na_if(trimws(chr("cftrm_eligible")), "")) %>%
  filter(!is.na(patient_id)) %>% distinct(patient_id, .keep_all = TRUE)
bl <- dat %>% arrange(patient_id, visit_date) %>% group_by(patient_id) %>% slice_head(n = 1) %>% ungroup() %>%
  left_join(gen, by = "patient_id") %>% left_join(starts, by = "patient_id") %>%
  mutate(group = if_else(!is.na(cftr_start_rev), "Ever", "Never"),
         public_only = sector == "Public",
         chronic_pa_bl = baseline_pa_cat %in% "Chronic",
         chronic_mssa_bl = sa_status %in% "Chronic",
         underweight_bl = baseline_nutrition_cat %in% "Underweight",
         eligible_yes = eligible %in% "Yes")
say(""); say("=== 3. Modulator exposure: commencement date vs the A01 flag (analytic cohort, ", nrow(bl), " patients) ===")
show_tbl(bl %>% count(A01_flag = ever_cftr_yes, start_date = group))
p_cat <- function(x, g) {
  tab <- table(x, g); if (nrow(tab) < 2) return(NA_real_)
  ex <- suppressWarnings(chisq.test(tab)$expected)
  if (any(ex < 5)) fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value else chisq.test(tab)$p.value
}
p_cont <- function(x, g) tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
sec2 <- if_else(bl$public_only %in% TRUE, "Public_Only", "Private_Combined")
t1 <- tibble(domain = "CFTR modulator therapy", variable = "Ever on modulator (commencement date)", level = "Yes",
             Public_Only = n_pct(sum(bl$group == "Ever" & sec2 == "Public_Only"), sum(sec2 == "Public_Only")),
             Private_Combined = n_pct(sum(bl$group == "Ever" & sec2 == "Private_Combined"), sum(sec2 == "Private_Combined")),
             Overall = n_pct(sum(bl$group == "Ever"), nrow(bl)),
             p_value = signif(p_cat(bl$group == "Ever", sec2), 3))
say(""); say("Table 1 row (replaces the A01-flag row):"); show_tbl(t1)
say("Public-sector-only share of treated patients: ", n_pct(sum(bl$group == "Ever" & bl$public_only %in% TRUE), sum(bl$group == "Ever")))
write_csv(t1, file.path(out_dir, "row22_table1_modulator.csv"))
row_cont <- function(v, label) { x <- bl[[v]]
  tibble(characteristic = label,
         Never = sprintf("%.1f (%.1f)", mean(x[bl$group == "Never"], na.rm = TRUE), sd(x[bl$group == "Never"], na.rm = TRUE)),
         Ever  = sprintf("%.1f (%.1f)", mean(x[bl$group == "Ever"], na.rm = TRUE),  sd(x[bl$group == "Ever"], na.rm = TRUE)),
         p = p_cont(x, bl$group)) }
row_bin <- function(v, label) { x <- as.logical(bl[[v]])
  tibble(characteristic = label,
         Never = n_pct(sum(x[bl$group == "Never"], na.rm = TRUE), sum(!is.na(x[bl$group == "Never"]))),
         Ever  = n_pct(sum(x[bl$group == "Ever"], na.rm = TRUE),  sum(!is.na(x[bl$group == "Ever"]))),
         p = p_cat(x, bl$group)) }
row_fac <- function(v, label) { x <- droplevels(factor(bl[[v]])); p <- p_cat(x, bl$group)
  bind_rows(lapply(levels(x), function(l) tibble(characteristic = paste0(label, ": ", l),
    Never = n_pct(sum(x[bl$group == "Never"] == l, na.rm = TRUE), sum(!is.na(x[bl$group == "Never"]))),
    Ever  = n_pct(sum(x[bl$group == "Ever"] == l, na.rm = TRUE),  sum(!is.na(x[bl$group == "Ever"]))),
    p = if (l == levels(x)[1]) p else NA_real_))) }
row22 <- bind_rows(
  tibble(characteristic = "Patients", Never = as.character(sum(bl$group == "Never")), Ever = as.character(sum(bl$group == "Ever")), p = NA_real_),
  row_cont("age_years", "Age at baseline, years, mean (SD)"),  # revision: age at the first eligible FEV1, as in Table 1 (A05.2); baseline_age_years is the first registry visit
  row_fac("sex", "Sex"),
  row_bin("public_only", "Public sector only"),
  row_fac("ethnicity_cat", "Ethnicity"),
  row_fac("f508_genotype_group", "F508del genotype"),
  row_bin("eligible_yes", "Registry: modulator-eligible"),
  row_bin("carrier_3120", "3120+1G>A carrier"),
  row_cont("baseline_fev1_pp", "Baseline FEV1pp, mean (SD)"),
  row_bin("underweight_bl", "Underweight at baseline"),
  row_bin("baseline_cfrd_indicator", "CFRD at baseline"),
  row_bin("baseline_pancreatic_insuff", "Pancreatic insufficiency at baseline"),
  row_bin("chronic_pa_bl", "Chronic P. aeruginosa at baseline"),
  row_bin("chronic_mssa_bl", "Chronic MSSA at baseline")
) %>% mutate(p = signif(p, 3))
say(""); say("Supplementary Table S1 (treated = any recorded commencement date, 2010-2024):"); show_tbl(row22)
write_csv(row22, file.path(out_dir, "row22_baseline_by_modulator.csv"))
say(""); say("Saved row22_r2_icc.csv, row22_gli_vs_local.csv, row22_table1_modulator.csv, row22_baseline_by_modulator.csv")
