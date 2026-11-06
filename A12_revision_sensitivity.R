#!/usr/bin/env Rscript
# ==============================================================================
# A12 - REMAINING SENSITIVITY ANALYSES (JCF-D-26-00377 revision, corrected data)
# ------------------------------------------------------------------------------
# Row 8   Main models restricted to patients with >=4 FEV1 values over >=3 years
# Row 9   R2 (Nakagawa) for the reporting models; M2 without baseline FEV1
#         category as a level covariate (kept as a slope modifier)
# Row 11  CFTR model: eligible private/combined patients only; calendar-year
#         and pandemic-period adjustment
# Row 13  Co-infection state (neither / MSSA only / P. aeruginosa only / both)
# Row 14  Age at incident chronic infection: summary and stratified models
# Run after A08 (uses its helper functions). Writes outputs/revision/ only.
# ==============================================================================

set.seed(20251111)
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(lme4); library(lmerTest); library(readxl); library(janitor); library(performance)
})
out_dir <- "outputs/revision"
summary_file <- file.path(out_dir, "A12_summary.txt")
fit_log <- list()
for (e in parse("A08_revision_fixes.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, envir = globalenv())
}
file.create(summary_file)

# Data preparation identical to A02/A08 --------------------------------------
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
dat_model <- dat_model %>% left_join(pancreas_start, by = "patient_id") %>% left_join(cfrd_start, by = "patient_id") %>%
  mutate(covid_period = as.integer(visit_year %in% c(2020L, 2021L)))
flags_main <- compute_use_flags(dat_model)
say("A12 sensitivity analyses, run ", format(Sys.time(), "%Y-%m-%d %H:%M"), "; model data ", nrow(dat_model), " rows, ",
    n_distinct(dat_model$patient_id), " patients")

# Revised CFTR start dates and registry eligibility (same rule as A08) ------
s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_df <- s1[grep("^patcode", names(s1), value = TRUE)]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
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
chr <- function(v) { x <- if (v %in% names(s1)) as.character(s1[[v]]) else rep(NA_character_, nrow(s1)); x[is.na(x)] <- ""; x }
anyf508 <- grepl("F508|Phe508del", paste(chr("mutation1"), chr("othermutation1"), chr("mutation2"), chr("othermutation2")))
pt_s1 <- tibble(patient_id = pid, cftr_start_rev = srev, eligible_reg = trimws(chr("cftrm_eligible")), any_f508_reg = anyf508) %>%
  filter(!is.na(patient_id)) %>% distinct(patient_id, .keep_all = TRUE) %>%
  mutate(cftr_start_rev = dplyr::if_else(is.finite(as.numeric(cftr_start_rev)), cftr_start_rev, as.Date(NA)),
         eligible_any = eligible_reg %in% "Yes" | any_f508_reg)

add_cftr <- function(df) {
  df %>% left_join(pt_s1, by = "patient_id") %>%
    mutate(cftr_post = as.integer(!is.na(cftr_start_rev) & visit_date >= cftr_start_rev),
           cftr_tsi = dplyr::if_else(cftr_post == 1L, as.numeric(visit_date - cftr_start_rev) / 365.25, 0))
}
drop_prevalent_users <- function(df) {
  prev <- df %>% group_by(patient_id) %>%
    summarise(pre = sum(!is.na(cftr_start_rev) & cftr_post == 0L), post = sum(cftr_post == 1L), .groups = "drop") %>%
    filter(post > 0 & pre == 0) %>% pull(patient_id)
  df %>% filter(!patient_id %in% prev)
}
rhs_c2 <- function(flags) paste(c("time_cs", "cftr_post", "cftr_tsi", event_level_covariates(flags)), collapse = " + ")
build_e1 <- function(df, spec) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index")
  df %>% mutate(event_start = .data[[spec$start_col]], anchor_date = dplyr::coalesce(event_start, baseline_date),
                !!tv := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
                !!pv := as.integer(!is.na(event_start) & visit_date >= event_start))
}
specs <- list(
  list(name = "pa",   start_col = "pa_chronic_start_date",   status_col = "pa_status", own = "chronic_paeruginosa", label = "P. aeruginosa"),
  list(name = "mssa", start_col = "mssa_chronic_start_date", status_col = "sa_status", own = "chronic_mssa",        label = "MSSA"))
rhs_e1 <- function(flags, spec) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index")
  rhs_drop(build_rhs_event(flags, tv, pv), spec$own)
}
key_terms <- c(sector_public_only = "Public sector only", nutrition_catUnderweight = "Underweight", cfrd_indicator = "CFRD",
               `exac_cat2+` = ">=2 exacerbations/year", cftr_modulator_indicator = "On CFTR modulator (level)",
               baseline_age_years_c = "Baseline age (per year)", sexMale = "Male sex", ethnicity_caucasian = "Caucasian ethnicity")
level_rows <- function(mod, label) {
  fe <- names(lme4::fixef(mod))
  bind_rows(lapply(names(key_terms)[names(key_terms) %in% fe], function(t) lc(mod, setNames(1, t), key_terms[[t]]))) %>%
    mutate(model = label, .before = 1)
}
avg_slope_row <- function(mod, label) {
  sd <- slope_design(mod)
  lc(mod, colMeans(sd$D), "Population-averaged slope (pp/yr)") %>% mutate(model = label, .before = 1)
}

# ------------------------------------------------------------------------------
# ROW 8: patients with >=4 FEV1 values over >=3 years
# ------------------------------------------------------------------------------
say("")
say("=== ROW 8: patients with >=4 FEV1 values spanning >=3 years ===")
dense_ids <- dat_model %>% filter(!is.na(fev1_pp)) %>% group_by(patient_id) %>%
  summarise(n = n(), span = as.numeric(max(visit_date) - min(visit_date)) / 365.25, .groups = "drop") %>%
  filter(n >= 4, span >= 3) %>% pull(patient_id)
dm8 <- dat_model %>% filter(patient_id %in% dense_ids) %>% droplevels()
f8 <- compute_use_flags(dm8)
say("Dense subset: ", length(dense_ids), " patients, ", nrow(dm8), " FEV1 values")
m2_all <- fit_triple_fallback(build_rhs(flags_main, "m2"), dat_model, TRUE, "M2 all")
m2_8 <- fit_triple_fallback(build_rhs(f8, "m2"), dm8, TRUE, "M2 dense")
dc_all <- drop_prevalent_users(add_cftr(dat_model))
dc8 <- dc_all %>% filter(patient_id %in% dense_ids)
c2_all <- fit_triple_fallback(rhs_c2(flags_main), dc_all, TRUE, "C2 all")
c2_8 <- fit_triple_fallback(rhs_c2(f8), dc8, TRUE, "C2 dense")
row8 <- bind_rows(avg_slope_row(m2_all, "All patients"), avg_slope_row(m2_8, "Dense subset"),
                  level_rows(m2_all, "All patients") %>% filter(contrast %in% c("Public sector only", "Underweight", "CFRD")),
                  level_rows(m2_8, "Dense subset") %>% filter(contrast %in% c("Public sector only", "Underweight", "CFRD")),
                  seg_rows(c2_all, "All patients", "time_cs", "cftr_tsi", "cftr_post") %>% mutate(contrast = paste("CFTR:", contrast)),
                  seg_rows(c2_8, "Dense subset", "time_cs", "cftr_tsi", "cftr_post") %>% mutate(contrast = paste("CFTR:", contrast)))
for (spec in specs) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index")
  ea <- fit_triple_fallback(rhs_e1(flags_main, spec), build_e1(dat_model, spec), TRUE, paste(spec$label, "E1 all"))
  e8 <- fit_triple_fallback(rhs_e1(f8, spec), build_e1(dm8, spec), TRUE, paste(spec$label, "E1 dense"))
  row8 <- bind_rows(row8,
    seg_rows(ea, "All patients", tv, paste0(tv, ":", pv), pv) %>% mutate(contrast = paste0(spec$label, ": ", contrast)),
    seg_rows(e8, "Dense subset", tv, paste0(tv, ":", pv), pv) %>% mutate(contrast = paste0(spec$label, ": ", contrast)))
}
row8 <- row8 %>% filter(!grepl("Difference from pre-event trend", contrast))
show_tbl(row8 %>% select(model, contrast, estimate, conf.low, conf.high, p.value) %>% arrange(contrast, model))
write_csv(row8, file.path(out_dir, "row8_dense_subset.csv"))

# ------------------------------------------------------------------------------
# ROW 9: R2 and M2 without baseline FEV1 category as a level covariate
# ------------------------------------------------------------------------------
say("")
say("=== ROW 9: R2 (Nakagawa) and M2 without the baseline FEV1 level term ===")
m2_nobl <- fit_triple_fallback(rhs_drop(build_rhs(flags_main, "m2"), "baseline_fev1_cat"), dat_model, TRUE, "M2 without baseline FEV1 level term")
e1_models <- lapply(specs, function(spec) fit_triple_fallback(rhs_e1(flags_main, spec), build_e1(dat_model, spec), TRUE, paste(spec$label, "E1")))
r2row <- function(mod, label) {
  r <- tryCatch(performance::r2_nakagawa(mod), error = function(e) NULL)
  tibble(model = label, marginal_R2 = if (is.null(r)) NA_real_ else round(r$R2_marginal, 3),
         conditional_R2 = if (is.null(r)) NA_real_ else round(r$R2_conditional, 3))
}
row9_r2 <- bind_rows(r2row(m2_all, "M2 (REML)"), r2row(m2_nobl, "M2 without baseline FEV1 level term"),
                     r2row(c2_all, "C2 CFTR segmented"), r2row(e1_models[[1]], "E1 P. aeruginosa"), r2row(e1_models[[2]], "E1 MSSA"))
show_tbl(row9_r2)
row9_est <- bind_rows(level_rows(m2_all, "M2"), level_rows(m2_nobl, "M2 without baseline FEV1 level term"),
                      avg_slope_row(m2_all, "M2"), avg_slope_row(m2_nobl, "M2 without baseline FEV1 level term"))
show_tbl(row9_est %>% select(model, contrast, estimate, conf.low, conf.high, p.value) %>% arrange(contrast, model))
write_csv(row9_r2, file.path(out_dir, "row9_r2.csv"))
write_csv(row9_est, file.path(out_dir, "row9_without_baseline_fev1.csv"))

# ------------------------------------------------------------------------------
# ROW 11: CFTR confounding checks
# ------------------------------------------------------------------------------
say("")
say("=== ROW 11: CFTR model, eligible private/combined only; calendar year; pandemic period ===")
dc_elig <- dc_all %>% filter(eligible_any, sector != "Public") %>% droplevels()
f11 <- compute_use_flags(dc_elig)
cnt <- function(df) {
  g <- df %>% group_by(patient_id) %>% summarise(tr = any(cftr_post == 1L), .groups = "drop")
  sprintf("%d patients (%d with on-treatment FEV1)", nrow(g), sum(g$tr))
}
say("Primary C2 sample: ", cnt(dc_all), "; eligible private/combined: ", cnt(dc_elig))
c2_elig <- fit_triple_fallback(rhs_c2(f11), dc_elig, TRUE, "C2 eligible private/combined")
c2_cal  <- fit_triple_fallback(paste(rhs_c2(flags_main), "+ calendar_year"), dc_all, TRUE, "C2 + calendar year")
c2_cov  <- fit_triple_fallback(paste(rhs_c2(flags_main), "+ covid_period"), dc_all, TRUE, "C2 + pandemic period 2020-21")
row11 <- bind_rows(seg_rows(c2_all, "C2 primary", "time_cs", "cftr_tsi", "cftr_post"),
                   seg_rows(c2_elig, "Eligible private/combined only", "time_cs", "cftr_tsi", "cftr_post"),
                   seg_rows(c2_cal, "+ calendar year", "time_cs", "cftr_tsi", "cftr_post"),
                   seg_rows(c2_cov, "+ pandemic period 2020-21", "time_cs", "cftr_tsi", "cftr_post")) %>%
  filter(!grepl("Difference from pre-event trend", contrast))
show_tbl(row11 %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
write_csv(row11, file.path(out_dir, "row11_cftr_confounding.csv"))

# ------------------------------------------------------------------------------
# ROW 13: co-infection state
# ------------------------------------------------------------------------------
say("")
say("=== ROW 13: co-infection state (time-varying) ===")
dm13 <- dat_model %>%
  mutate(infection_state = factor(dplyr::case_when(
    chronic_paeruginosa == 1L & chronic_mssa == 1L ~ "Both",
    chronic_paeruginosa == 1L ~ "PA_only",
    chronic_mssa == 1L ~ "MSSA_only",
    TRUE ~ "Neither"), levels = c("Neither", "MSSA_only", "PA_only", "Both")))
say("Visits by state:")
show_tbl(dm13 %>% count(infection_state, name = "visits") %>%
           left_join(dm13 %>% distinct(patient_id, infection_state) %>% count(infection_state, name = "patients_ever"), by = "infection_state"))
cov13 <- setdiff(event_level_covariates(flags_main), c("chronic_paeruginosa", "chronic_mssa"))
rhs13 <- paste(c("time_cs", "infection_state", cov13, "cftr_modulator_indicator", "time_cs:infection_state"), collapse = " + ")
m13 <- fit_triple_fallback(rhs13, dm13, TRUE, "Co-infection model")
row13 <- bind_rows(lapply(c("MSSA_only", "PA_only", "Both"), function(s) bind_rows(
  lc(m13, setNames(1, paste0("infection_state", s)), paste0(s, ": level vs neither (pp)")),
  lc(m13, setNames(1, paste0("time_cs:infection_state", s)), paste0(s, ": slope vs neither (pp/yr)")),
  lc(m13, setNames(c(1, 1), c("time_cs", paste0("time_cs:infection_state", s))), paste0(s, ": slope (pp/yr)")))),
  lc(m13, c(time_cs = 1), "Neither: slope (pp/yr)")) %>% mutate(model = "Co-infection model", .before = 1)
show_tbl(row13 %>% select(contrast, estimate, conf.low, conf.high, p.value))
write_csv(row13, file.path(out_dir, "row13_coinfection.csv"))

# ------------------------------------------------------------------------------
# ROW 14: age at incident chronic infection
# ------------------------------------------------------------------------------
say("")
say("=== ROW 14: age at incident chronic infection ===")
row14 <- list(); row14_age <- list()
for (spec in specs) {
  ev <- classify_events(dat_model, spec$start_col, spec$status_col)
  inc <- ev %>% filter(event_class == "incident")
  age_on <- dat_model %>% inner_join(inc %>% select(patient_id, onset), by = "patient_id") %>%
    filter(visit_date >= onset) %>% group_by(patient_id) %>% slice_min(visit_date, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(patient_id, age_onset = age_years)
  med <- median(age_on$age_onset, na.rm = TRUE)
  row14_age[[spec$label]] <- tibble(organism = spec$label, incident = nrow(inc), age_at_onset = med_iqr(age_on$age_onset, 1),
                                    range = sprintf("%.1f-%.1f", min(age_on$age_onset, na.rm = TRUE), max(age_on$age_onset, na.rm = TRUE)),
                                    split_at = round(med, 1))
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index")
  d1 <- build_e1(dat_model, spec)
  never <- ev %>% filter(event_class == "never") %>% pull(patient_id)
  young <- age_on %>% filter(age_onset < med) %>% pull(patient_id)
  old   <- age_on %>% filter(age_onset >= med) %>% pull(patient_id)
  my <- fit_triple_fallback(rhs_e1(flags_main, spec), d1 %>% filter(patient_id %in% c(never, young)), TRUE, paste(spec$label, "incident, younger"))
  mo <- fit_triple_fallback(rhs_e1(flags_main, spec), d1 %>% filter(patient_id %in% c(never, old)),   TRUE, paste(spec$label, "incident, older"))
  row14[[spec$label]] <- bind_rows(
    seg_rows(my, paste0(spec$label, ": onset age < ", round(med, 1), " y (n=", length(young), ")"), tv, paste0(tv, ":", pv), pv),
    seg_rows(mo, paste0(spec$label, ": onset age >= ", round(med, 1), " y (n=", length(old), ")"), tv, paste0(tv, ":", pv), pv)) %>%
    filter(contrast %in% c("Level change at event (pp)", "Slope change after event (pp/yr)", "Post-event slope, net (pp/yr)"))
}
show_tbl(bind_rows(row14_age))
show_tbl(bind_rows(row14) %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
write_csv(bind_rows(row14_age), file.path(out_dir, "row14_age_at_onset.csv"))
write_csv(bind_rows(row14), file.path(out_dir, "row14_age_stratified_models.csv"))

say("")
say("=== Fit log ===")
show_tbl(bind_rows(fit_log) %>% select(model, re_structure, n_patients, n_obs, n_warnings))
say("Done.")
