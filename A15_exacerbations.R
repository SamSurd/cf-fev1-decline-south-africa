#!/usr/bin/env Rscript
# ==============================================================================
# A15 - EXACERBATIONS (checklist row 18, corrected data)
# ------------------------------------------------------------------------------
# 1. Coding check: exac_cat is the SAME-YEAR count of IV antibiotic courses
#    (registry ivi_number<year>; a missing count is coded 0, as in A01).
# 2. Lagged sensitivity: M2 refitted with the PREVIOUS registry year's
#    exacerbation category (prior-year counts from Sheet1 through A01's own
#    builder), alone and together with the same-year category, on visits that
#    have a prior-year registry record. L0 = M2 on the same visits.
# 3. Exacerbation event models: the submitted unadjusted models (tidies written
#    by Supplementary Scripts/Build_Exac_Anchored_Tidies.R) and adjusted
#    versions in the E1 framework (A08 covariates, triple fallback), anchored at
#    the first FEV1 year with >=1 and with >=2 exacerbations.
# 4. Which model covariates vary within patients (for the Figure 1 grouping).
# Run after A08. Writes outputs/revision/ only.
# ==============================================================================
set.seed(20251111)
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(lme4); library(lmerTest); library(readxl); library(janitor) })
out_dir <- "outputs/revision"; summary_file <- file.path(out_dir, "A15_summary.txt"); fit_log <- list()
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
  group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup() %>%
  mutate(visit_year = as.integer(visit_year))
flags_main <- compute_use_flags(dat_model)
dm <- dat_model %>% filter(!is.na(fev1_pp))
say("A15 exacerbation analyses, run ", format(Sys.time(), "%Y-%m-%d %H:%M"), "; model data ", nrow(dat_model), " rows (",
    nrow(dm), " with FEV1), ", n_distinct(dm$patient_id), " patients")

# ------------------------------------------------------------------------------
# 1. Coding check
# ------------------------------------------------------------------------------
say(""); say("=== 1. Exacerbation coding (modelled visits with FEV1) ===")
show_tbl(dm %>% count(exac_cat) %>% mutate(pct = round(100 * n / sum(n), 1)))
if ("exacerbation_count" %in% names(dm)) {
  say("Visits with no recorded IV-course count (coded 0): ", sum(is.na(dm$exacerbation_count)), " of ", nrow(dm),
      " (", round(100 * mean(is.na(dm$exacerbation_count)), 1), "%)")
}
say("exac_cat is taken from the same registry year as the FEV1 value (A01: ivi_number<year> -> exacerbations_per_year -> exac_cat).")

# 1b. The registry also has a yes/no field, ivi_antibiotics<year>, which A01 does not use.
s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_df <- s1[grep("^patcode", names(s1), value = TRUE)]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
yn_cols <- grep("^ivi_antibiotics_?20[0-9]{2}$", names(s1), value = TRUE)
nm_cols <- grep("^ivi_number_?20[0-9]{2}$", names(s1), value = TRUE)
ivi <- bind_rows(lapply(2018:2024, function(y) {
  yc <- yn_cols[grepl(as.character(y), yn_cols)]; nc <- nm_cols[grepl(as.character(y), nm_cols)]
  yn <- if (length(yc)) tolower(trimws(as.character(s1[[yc[1]]]))) else rep(NA_character_, nrow(s1))
  nn <- if (length(nc)) suppressWarnings(as.numeric(as.character(s1[[nc[1]]]))) else rep(NA_real_, nrow(s1))
  tibble(patient_id = pid, visit_year = as.integer(y), ivi_yn = yn, ivi_n = nn)
})) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, visit_year, .keep_all = TRUE)
dmi <- dat_model %>% left_join(ivi, by = c("patient_id", "visit_year")) %>%
  mutate(yn = dplyr::case_when(ivi_yn %in% c("yes", "y", "1") ~ "Yes", ivi_yn %in% c("no", "n", "0") ~ "No",
                               is.na(ivi_yn) | ivi_yn == "" ~ "Blank", TRUE ~ "Other"),
         count_recorded = !is.na(ivi_n),
         exac_rec = factor(dplyr::case_when(!is.na(ivi_n) & ivi_n >= 2 ~ "2+", !is.na(ivi_n) & ivi_n == 1 ~ "1",
                                            !is.na(ivi_n) & ivi_n == 0 ~ "0", yn == "No" ~ "0",
                                            yn == "Yes" ~ "yes_nocount", TRUE ~ "unknown"),
                           levels = c("0", "1", "2+", "yes_nocount", "unknown")))
say("Sheet1 yes/no columns: ", paste(yn_cols, collapse = ", "))
if ("exacerbation_count" %in% names(dmi)) say("Sheet1 ivi_number agrees with A01 exacerbation_count on ",
  sum(dplyr::coalesce(dmi$ivi_n, -1) == dplyr::coalesce(dmi$exacerbation_count, -1)), " of ", nrow(dmi), " modelled visits")
show_tbl(dmi %>% count(ivi_antibiotics = yn, count_recorded, exac_cat))
show_tbl(dmi %>% count(exac_rec) %>% mutate(pct = round(100 * n / sum(n), 1)))

# ------------------------------------------------------------------------------
# 4 (run early). Within-patient variation of the model covariates
# ------------------------------------------------------------------------------
say(""); say("=== 4. Patients whose covariate value changes during follow-up ===")
cov_vars <- intersect(c("baseline_age_years_c", "sex", "ethnicity_caucasian", "sector_public_only", "f508_any", "baseline_pa_any",
                        "baseline_fev1_cat", "nutrition_cat", "cfrd_indicator", "pancreatic_insufficiency", "chronic_paeruginosa",
                        "chronic_mssa", "chronic_mrsa", "chronic_bcepacia", "exac_cat", "cftr_modulator_indicator"), names(dm))
variation <- map_dfr(cov_vars, function(v) {
  k <- dm %>% group_by(patient_id) %>% summarise(k = dplyr::n_distinct(.data[[v]], na.rm = TRUE), .groups = "drop")
  tibble(variable = v, patients_varying = sum(k$k > 1), patients = nrow(k))
})
show_tbl(variation)
write_csv(variation, file.path(out_dir, "row18_covariate_variation.csv"))

# ------------------------------------------------------------------------------
# 2. Lagged exacerbations
# ------------------------------------------------------------------------------
say(""); say("=== 2. Prior-year exacerbations ===")
a01 <- new.env()
for (e in parse("A01.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, envir = a01)
}
raw <- suppressMessages(a01$make_cf_longitudinal_from_sheet1("data_raw/CF_registry_master.xlsx")) %>% janitor::clean_names()
stopifnot(all(c("patcode", "year", "exacerbations_per_year") %in% names(raw)))
say("Registry patient-years from A01's builder: ", nrow(raw), "; IV-course count recorded in ", sum(!is.na(raw$exacerbations_per_year)))
prior <- raw %>%
  transmute(patient_id = as.character(patcode), visit_year = as.integer(year) + 1L,
            n_iv = suppressWarnings(as.numeric(exacerbations_per_year))) %>%
  filter(!is.na(patient_id)) %>%
  group_by(patient_id, visit_year) %>%
  summarise(prev_count = if (all(is.na(n_iv))) 0 else max(n_iv, na.rm = TRUE), .groups = "drop")  # record exists; missing count = 0 (A01 rule)
d_lag <- dat_model %>% left_join(prior, by = c("patient_id", "visit_year")) %>%
  filter(!is.na(prev_count)) %>%
  mutate(exac_prev_cat = factor(dplyr::case_when(prev_count == 0 ~ "0", prev_count == 1 ~ "1", TRUE ~ "2+"), levels = c("0", "1", "2+"))) %>%
  group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup()
say("Visits with a prior-year registry record: ", sum(!is.na(d_lag$fev1_pp)), " of ", nrow(dm), " FEV1 visits; ",
    n_distinct(d_lag$patient_id), " patients")
show_tbl(d_lag %>% filter(!is.na(fev1_pp)) %>% count(prior_year = exac_prev_cat, same_year = exac_cat) %>%
           pivot_wider(names_from = same_year, values_from = n, values_fill = 0))

exac_terms <- c("exac_cat1", "exac_cat2+", "exac_prev_cat1", "exac_prev_cat2+")
term_rows <- function(mod, label) {
  if (is.null(mod)) return(tibble())
  tt <- intersect(exac_terms, names(lme4::fixef(mod)))
  bind_rows(lapply(tt, function(t) lc(mod, setNames(1, t), t))) %>%
    mutate(model = label, .before = 1) %>%
    mutate(n_patients = unname(lme4::ngrps(mod)[["patient_id"]]), n_obs = nobs(mod))
}
rhs_full <- build_rhs(flags_main, "m2")
rhs_lag  <- build_rhs(compute_use_flags(d_lag), "m2")
m2_full <- fit_triple_fallback(rhs_full, dat_model, TRUE, "M2 full cohort (same-year, as reported)")
l0 <- fit_triple_fallback(rhs_lag, d_lag, TRUE, "L0 M2 same-year, lag subset")
l1 <- fit_triple_fallback(paste(rhs_drop(rhs_lag, "exac_cat"), "+ exac_prev_cat"), d_lag, TRUE, "L1 M2 prior-year only, lag subset")
l2 <- fit_triple_fallback(paste(rhs_lag, "+ exac_prev_cat"), d_lag, TRUE, "L2 M2 same-year + prior-year, lag subset")
lag_tab <- bind_rows(term_rows(m2_full, "M2 full cohort (same-year, as reported)"),
                     term_rows(l0, "L0 M2 same-year, lag subset"),
                     term_rows(l1, "L1 M2 prior-year only, lag subset"),
                     term_rows(l2, "L2 M2 same-year + prior-year, lag subset"))
show_tbl(lag_tab %>% select(model, contrast, estimate, conf.low, conf.high, p.value, n_patients, n_obs))
pa_slope <- function(mod, label) { sd <- slope_design(mod); lc(mod, colMeans(sd$D), "Population-averaged slope (pp/yr)") %>% mutate(model = label, .before = 1) }
show_tbl(bind_rows(pa_slope(m2_full, "M2 full"), pa_slope(l0, "L0"), pa_slope(l2, "L2")) %>% select(model, estimate, conf.low, conf.high))

say(""); say("=== 2b. Missing IV-course counts: true zeros from the yes/no field ===")
rec_terms <- c("exac_rec1", "exac_rec2+", "exac_recyes_nocount", "exac_recunknown")
rec_rows <- function(mod, label) {
  if (is.null(mod)) return(tibble())
  tt <- intersect(rec_terms, names(lme4::fixef(mod)))
  bind_rows(lapply(tt, function(t) lc(mod, setNames(1, t), t))) %>% mutate(model = label, .before = 1) %>%
    mutate(n_patients = unname(lme4::ngrps(mod)[["patient_id"]]), n_obs = nobs(mod))
}
rhs_rec <- paste(rhs_drop(build_rhs(flags_main, "m2"), "exac_cat"), "+ exac_rec")
s_all <- fit_triple_fallback(rhs_rec, dmi, TRUE, "S1 M2, exacerbations recoded (unknown as own category)")
d_known <- dmi %>% filter(exac_rec %in% c("0", "1", "2+")) %>% droplevels() %>% group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup()
s_known <- fit_triple_fallback(paste(rhs_drop(build_rhs(compute_use_flags(d_known), "m2"), "exac_cat"), "+ exac_rec"), d_known, TRUE,
                               "S2 M2, visits with known exacerbation status only")
rec_tab <- bind_rows(rec_rows(s_all, "S1 M2, exacerbations recoded (unknown as own category)"),
                     rec_rows(s_known, "S2 M2, visits with known exacerbation status only"))
show_tbl(rec_tab %>% select(model, contrast, estimate, conf.low, conf.high, p.value, n_patients, n_obs))

# ------------------------------------------------------------------------------
# 3. Exacerbation event models
# ------------------------------------------------------------------------------
say(""); say("=== 3a. Submitted unadjusted exacerbation event models (re-run tidies, corrected data) ===")
crude <- bind_rows(lapply(c(`First year with exactly 1 (submitted, unadjusted)` = ".model_exac1_anchored_tidy.csv",
                            `First year with >=2 (submitted, unadjusted)` = ".model_exac2plus_anchored_tidy.csv"),
                          function(f) read_csv(file.path("outputs/tables", f), show_col_types = FALSE)), .id = "model") %>%
  filter(term != "(Intercept)")
show_tbl(crude %>% select(model, term, estimate, conf.low, conf.high, p.value))

say(""); say("=== 3b. Adjusted exacerbation event models (E1 framework: A08 covariates without exac_cat) ===")
exac_specs <- list(list(name = "exac1", cats = c("1", "2+"), label = "First year with >=1 exacerbation"),
                   list(name = "exac2", cats = "2+",        label = "First year with >=2 exacerbations"))
ev_rows <- list(); ev_counts <- list()
for (sp in exac_specs) {
  tv <- paste0("time_", sp$name, "_index_cs"); pv <- paste0(sp$name, "_post_index"); iv <- paste0(tv, ":", pv)
  first_ev <- dm %>% filter(exac_cat %in% sp$cats) %>% group_by(patient_id) %>%
    summarise(event_start = min(visit_date), .groups = "drop")
  d1 <- dat_model %>% left_join(first_ev, by = "patient_id") %>%
    mutate(anchor_date = dplyr::coalesce(event_start, baseline_date),
           !!tv := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
           !!pv := as.integer(!is.na(event_start) & visit_date >= event_start))
  cls <- d1 %>% filter(!is.na(fev1_pp)) %>% group_by(patient_id) %>%
    summarise(event = any(!is.na(event_start)), first_visit = min(visit_date), ev = dplyr::first(event_start),
              n_pre = sum(.data[[pv]] == 0L), n_post = sum(.data[[pv]] == 1L), .groups = "drop") %>%
    mutate(class = dplyr::case_when(!event ~ "never", ev <= first_visit ~ "at first FEV1 year", TRUE ~ "later"))
  ev_counts[[sp$name]] <- cls %>% group_by(event = sp$label, class) %>%
    summarise(patients = dplyr::n(), median_pre = median(n_pre), median_post = median(n_post), .groups = "drop")
  m <- fit_triple_fallback(rhs_drop(build_rhs_event(compute_use_flags(d1), tv, pv), "exac_cat"), d1, TRUE, paste(sp$label, "adjusted"))
  ev_rows[[sp$name]] <- seg_rows(m, paste(sp$label, "(adjusted)"), tv, iv, pv) %>%
    mutate(n_patients = unname(lme4::ngrps(m)[["patient_id"]]), n_obs = nobs(m))
  # As E2 for infections: drop patients whose first such year is their first FEV1 year
  keep_inc <- cls %>% filter(class != "at first FEV1 year") %>% pull(patient_id)
  d_inc <- d1 %>% filter(patient_id %in% keep_inc)
  m_inc <- fit_triple_fallback(rhs_drop(build_rhs_event(compute_use_flags(d_inc), tv, pv), "exac_cat"), d_inc, TRUE,
                               paste(sp$label, "adjusted, later events only"))
  ev_rows[[paste0(sp$name, "_inc")]] <- seg_rows(m_inc, paste(sp$label, "(adjusted, later events only)"), tv, iv, pv) %>%
    mutate(n_patients = unname(lme4::ngrps(m_inc)[["patient_id"]]), n_obs = nobs(m_inc))
}
show_tbl(bind_rows(ev_counts))
ev_tab <- bind_rows(ev_rows)
show_tbl(ev_tab %>% filter(contrast %in% c("Pre-event slope (pp/yr)", "Level change at event (pp)", "Slope change after event (pp/yr)",
                                           "Post-event slope, net (pp/yr)")) %>%
           select(model, contrast, estimate, conf.low, conf.high, p.value))

say(""); say("=== Fits ===")
show_tbl(bind_rows(fit_log) %>% select(model, re_structure, n_patients, n_obs, n_warnings))

out <- bind_rows(lag_tab %>% mutate(part = "lagged"),
                 rec_tab %>% mutate(part = "missing counts recoded"),
                 ev_tab %>% mutate(part = "event, adjusted"),
                 crude %>% transmute(part = "event, submitted unadjusted", model, contrast = term, estimate, se = std.error,
                                     df, conf.low, conf.high, p.value))
write_csv(out, file.path(out_dir, "row18_exacerbations.csv"))
write_csv(bind_rows(ev_counts), file.path(out_dir, "row18_exacerbation_event_counts.csv"))
say(""); say("Wrote row18_exacerbations.csv, row18_exacerbation_event_counts.csv, row18_covariate_variation.csv")
