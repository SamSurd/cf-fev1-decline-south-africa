#!/usr/bin/env Rscript
# ==============================================================================
# A09 - REVISION OUTPUTS, BATCH 2 (JCF-D-26-00377)
# ------------------------------------------------------------------------------
# Row 5   Infection event models (all-patients specification, as A08 "E1")
#         censored at CFTR modulator start
# Row 6   Follow-up density: FEV1 values per patient and follow-up span
# Row 7   Sequential cohort flow for the STROBE flowchart
# Row 10  Baseline characteristics by modulator exposure
# Row 12  Respiratory samples per patient-year
# Row 16  3120+1G>A carriage and modulator eligibility
# Row 17  Benjamini-Hochberg adjusted p-values for the reporting models
#
# Run from the project folder after A08:  Rscript A09_revision_outputs.R
# Writes: outputs/revision/ only. Existing outputs are not touched.
# ==============================================================================

set.seed(20251111)

suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(lme4); library(lmerTest)
  library(readxl); library(janitor); library(broom.mixed)
})

data_file <- "outputs/data_derived/analysis_long.parquet"
raw_file  <- "data_raw/CF_registry_master.xlsx"
out_dir   <- "outputs/revision"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path(out_dir, "A09_summary.txt")
file.create(summary_file)
say <- function(...) {
  msg <- paste0(...)
  message(msg)
  con <- file(summary_file, open = "a"); writeLines(msg, con); close(con)
}
show_tbl <- function(x, digits = 3) {
  txt <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  for (l in txt) say(l)
}
med_iqr <- function(x, d = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  q <- quantile(x, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf(paste0("%.", d, "f (IQR %.", d, "f-%.", d, "f)"), q[1], q[2], q[3])
}
n_pct <- function(k, n) sprintf("%d (%.1f%%)", as.integer(k), 100 * k / n)

# ------------------------------------------------------------------------------
# 1. Data preparation: identical to A02.R / A08
# ------------------------------------------------------------------------------
dat <- arrow::read_parquet(data_file) %>% as_tibble()
dat <- dat %>%
  mutate(
    sex           = factor(sex, levels = c("Female", "Male")),
    nutrition_cat = factor(nutrition_cat, levels = c("Normal", "Underweight", "Overweight", "Unknown")),
    exac_cat      = factor(exac_cat, levels = c("0", "1", "2+", "Unknown")),
    sector        = factor(sector, levels = c("Public", "Private", "Combined")),
    ethnicity_cat = factor(ethnicity_cat, levels = c("Black African", "Mixed Race", "Indian", "Asian", "Caucasian", "Other", "Unknown")),
    f508_genotype_group = factor(f508_genotype_group, levels = c("Other", "Heterozygous", "Homozygous")),
    ethnicity_caucasian = dplyr::if_else(ethnicity_cat == "Caucasian", 1L, 0L, missing = 0L),
    sector_public_only  = dplyr::if_else(sector == "Public", 1L, 0L, missing = 0L),
    f508_any = dplyr::if_else(f508_genotype_group %in% c("Heterozygous", "Homozygous"), 1L, 0L, missing = 0L)
  )
dat_model <- dat %>%
  filter(nutrition_cat != "Unknown" | is.na(nutrition_cat)) %>%
  droplevels() %>%
  group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup()

say("A09 revision outputs, run ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Parquet: ", nrow(dat), " rows, ", n_distinct(dat$patient_id), " patients, ", sum(!is.na(dat$fev1_pp)), " rows with FEV1pp")
say("Model data: ", nrow(dat_model), " rows, ", n_distinct(dat_model$patient_id), " patients")

# ------------------------------------------------------------------------------
# 2. Helpers copied from A08 (fitting, contrasts, event models, start dates)
# ------------------------------------------------------------------------------
compute_use_flags <- function(df) {
  nlev <- function(x) length(unique(na.omit(x)))
  list(use_blf_cat = nlev(df$baseline_fev1_cat) >= 2, use_sex = nlev(df$sex) >= 2,
       use_nut = nlev(df$nutrition_cat) >= 2, use_exac = nlev(df$exac_cat) >= 2,
       use_sector = nlev(df$sector_public_only) >= 2, use_eth = nlev(df$ethnicity_caucasian) >= 2,
       use_geno = nlev(df$f508_any) >= 2)
}
fit_log <- list()
fit_triple_fallback <- function(rhs_formula, data, use_reml = FALSE, label) {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  forms <- list(
    correlated   = as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 + time_cs | patient_id)")),
    uncorrelated = as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 + time_cs || patient_id)")),
    intercept    = as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 | patient_id)")))
  for (nm in names(forms)) {
    w <- character(0)
    mod <- tryCatch(
      withCallingHandlers(lmerTest::lmer(forms[[nm]], data = data, REML = use_reml, control = ctrl),
                          warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }),
      error = function(e) { say("    ", label, ": ", nm, " fit failed: ", conditionMessage(e)); NULL })
    if (is.null(mod)) next
    sing <- lme4::isSingular(mod, tol = 1e-4)
    if (!sing || nm == "intercept") {
      fit_log[[label]] <<- tibble(model = label, re_structure = paste0(nm, if (sing) " (singular)" else ""),
                                  n_patients = unname(lme4::ngrps(mod)[["patient_id"]]), n_obs = nobs(mod),
                                  n_warnings = length(w), first_warning = if (length(w)) w[1] else "")
      return(mod)
    }
  }
  NULL
}
lc <- function(mod, w, label) {
  fe <- lme4::fixef(mod)
  if (!all(names(w) %in% names(fe))) stop("Unknown term(s) in contrast: ", paste(setdiff(names(w), names(fe)), collapse = ", "))
  L <- setNames(numeric(length(fe)), names(fe)); L[names(w)] <- w
  r <- tryCatch(as.data.frame(lmerTest::contest1D(mod, unname(L), confint = TRUE)), error = function(e) NULL)
  need <- c("Estimate", "Std. Error", "df", "lower", "upper", "Pr(>|t|)")
  if (!is.null(r) && all(need %in% names(r))) {
    tibble(contrast = label, estimate = r[["Estimate"]], se = r[["Std. Error"]], df = r[["df"]],
           conf.low = r[["lower"]], conf.high = r[["upper"]], p.value = r[["Pr(>|t|)"]])
  } else {
    V <- as.matrix(vcov(mod)); est <- sum(L * fe); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    tibble(contrast = label, estimate = est, se = se, df = NA_real_, conf.low = est - 1.96 * se,
           conf.high = est + 1.96 * se, p.value = 2 * pnorm(-abs(est / se)))
  }
}
seg_rows <- function(mod, label, t_term, slope_term, post_term) {
  bind_rows(
    lc(mod, setNames(1, t_term), "Pre-event slope (pp/yr)"),
    lc(mod, setNames(1, post_term), "Level change at event (pp)"),
    lc(mod, setNames(1, slope_term), "Slope change after event (pp/yr)"),
    lc(mod, setNames(c(1, 1), c(t_term, slope_term)), "Post-event slope, net (pp/yr)"),
    lc(mod, setNames(c(1, 2), c(post_term, slope_term)), "Difference from pre-event trend at +2 yr (pp)")
  ) %>% mutate(model = label, .before = 1)
}
event_level_covariates <- function(flags) {
  terms <- c("pancreatic_insufficiency", "cfrd_indicator", "baseline_fev1_pp_c", "baseline_age_years_c")
  if (flags$use_nut)    terms <- c(terms, "nutrition_cat")
  if (flags$use_exac)   terms <- c(terms, "exac_cat")
  if (flags$use_sector) terms <- c(terms, "sector")
  if (flags$use_eth)    terms <- c(terms, "ethnicity_cat")
  if (flags$use_geno)   terms <- c(terms, "f508_genotype_group")
  c(terms, "chronic_paeruginosa", "chronic_bcepacia", "chronic_mssa", "chronic_mrsa")
}
build_rhs_event <- function(flags, time_var_cs, post_var) {
  paste(c(time_var_cs, post_var, event_level_covariates(flags), "cftr_modulator_indicator",
          paste0(time_var_cs, ":", post_var)), collapse = " + ")
}
rhs_drop <- function(rhs, drop) paste(setdiff(strsplit(rhs, " + ", fixed = TRUE)[[1]], drop), collapse = " + ")
flags_main <- compute_use_flags(dat_model)

s1 <- suppressMessages(readxl::read_excel(raw_file, sheet = "Sheet1")) %>% janitor::clean_names()
id_cols <- grep("^patcode", names(s1), value = TRUE)
id_df <- s1[id_cols]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)

# Revised CFTR start dates (same rule as A08)
to_date <- function(x) {
  if (inherits(x, c("POSIXct", "POSIXt", "Date"))) return(as.Date(x))
  suppressWarnings(lubridate::dmy(as.character(x), quiet = TRUE))
}
date_cols <- intersect(c("date_cftrm2018", "date_cftrm2019", "date_cftrm2020", "date_cftrm2021", "date_cftrm2022",
                         "date_cftrm2023", "date_cftrm_commenced_2024", "treatment_start_date_cftr",
                         "treatment_start_date_cftr_new"), names(s1))
lo <- as.Date("2010-01-01"); hi <- as.Date("2024-12-31")
dm <- sapply(date_cols, function(v) { d <- to_date(s1[[v]]); d[!is.na(d) & (d < lo | d > hi)] <- NA; as.numeric(d) })
earliest_num <- suppressWarnings(apply(dm, 1, min, na.rm = TRUE)); earliest_num[!is.finite(earliest_num)] <- NA
earliest <- as.Date(earliest_num, origin = "1970-01-01")
flag_cols <- grep("^cftr_modifier_?20[0-9]{2}$", names(s1), value = TRUE)
flag_year <- as.integer(stringr::str_extract(flag_cols, "20[0-9]{2}"))
Y <- sapply(flag_cols, function(v) as.character(s1[[v]]) %in% "Yes")
first_yes <- apply(Y, 1, function(r) if (any(r)) min(flag_year[r]) else NA_integer_)
impute <- !is.na(first_yes) & (is.na(earliest) | as.integer(format(earliest, "%Y")) > first_yes)
start_rev <- earliest; start_rev[impute] <- as.Date(paste0(first_yes[impute], "-07-01"))
starts <- tibble(patient_id = pid, cftr_start_rev = start_rev) %>%
  filter(!is.na(patient_id)) %>% group_by(patient_id) %>%
  summarise(cftr_start_rev = suppressWarnings(min(cftr_start_rev, na.rm = TRUE)), .groups = "drop") %>%
  mutate(cftr_start_rev = dplyr::if_else(is.finite(as.numeric(cftr_start_rev)), cftr_start_rev, as.Date(NA)))

# ------------------------------------------------------------------------------
# 3. ROW 5: infection models censored at modulator start
# ------------------------------------------------------------------------------
say("")
say("=== ROW 5: infection models (all patients) censored at CFTR modulator start ===")
event_specs <- list(
  list(name = "pa",   start_col = "pa_chronic_start_date",   own = "chronic_paeruginosa", label = "P. aeruginosa"),
  list(name = "mssa", start_col = "mssa_chronic_start_date", own = "chronic_mssa",        label = "MSSA"))
row5 <- list(); e1_models <- list(); cens_counts <- list()
for (spec in event_specs) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index"); iv <- paste0(tv, ":", pv)
  d1 <- dat_model %>%
    mutate(event_start = .data[[spec$start_col]],
           anchor_date = dplyr::coalesce(event_start, baseline_date),
           !!tv := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
           !!pv := as.integer(!is.na(event_start) & visit_date >= event_start))
  rhs1 <- rhs_drop(build_rhs_event(flags_main, tv, pv), spec$own)
  m_e1 <- fit_triple_fallback(rhs1, d1, TRUE, paste(spec$label, "E1 all patients"))
  e1_models[[spec$label]] <- m_e1

  d1c <- d1 %>% left_join(starts, by = "patient_id") %>%
    filter(is.na(cftr_start_rev) | visit_date < cftr_start_rev) %>%
    group_by(patient_id) %>% filter(sum(!is.na(fev1_pp)) >= 2) %>% ungroup()
  m_e1c <- fit_triple_fallback(rhs_drop(rhs1, "cftr_modulator_indicator"), d1c, TRUE,
                               paste(spec$label, "E1 censored at modulator start"))
  cens_counts[[spec$label]] <- tibble(
    organism = spec$label,
    patients_after_censoring = n_distinct(d1c$patient_id), fev1_after_censoring = sum(!is.na(d1c$fev1_pp)),
    patients_with_post_onset_fev1 = d1c %>% filter(.data[[pv]] == 1L, !is.na(fev1_pp)) %>% distinct(patient_id) %>% nrow())
  row5[[spec$label]] <- bind_rows(
    seg_rows(m_e1,  paste(spec$label, "E1 all patients"), tv, iv, pv),
    seg_rows(m_e1c, paste(spec$label, "E1 censored at modulator start"), tv, iv, pv))
}
row5_tab <- bind_rows(row5)
show_tbl(bind_rows(cens_counts))
show_tbl(row5_tab %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
write_csv(row5_tab, file.path(out_dir, "row5_censored_event_models.csv"))

# ------------------------------------------------------------------------------
# 4. ROW 6: follow-up density
# ------------------------------------------------------------------------------
say("")
say("=== ROW 6: follow-up density ===")
fu <- function(df) df %>% filter(!is.na(fev1_pp)) %>% group_by(patient_id) %>%
  summarise(n_fev1 = n(), span = as.numeric(max(visit_date) - min(visit_date)) / 365.25, .groups = "drop")
fu_summary <- function(f, label) tibble(
  cohort = label, patients = nrow(f), fev1_values = sum(f$n_fev1),
  fev1_per_patient = med_iqr(f$n_fev1, 0),
  pct_ge3 = round(100 * mean(f$n_fev1 >= 3), 1), pct_ge4 = round(100 * mean(f$n_fev1 >= 4), 1),
  pct_ge5 = round(100 * mean(f$n_fev1 >= 5), 1),
  followup_years = med_iqr(f$span, 1), pct_span_ge3y = round(100 * mean(f$span >= 3), 1))
f439 <- fu(dat); f438 <- fu(dat_model)
row6 <- bind_rows(fu_summary(f439, "Analytic cohort"), fu_summary(f438, "Modelled"))
show_tbl(row6)
say("FEV1 values per patient (analytic cohort):")
show_tbl(as.data.frame(table(n_fev1 = f439$n_fev1)))
write_csv(row6, file.path(out_dir, "row6_followup_density.csv"))

# ------------------------------------------------------------------------------
# 5. ROW 7: sequential cohort flow (re-uses A01's own functions)
# ------------------------------------------------------------------------------
say("")
say("=== ROW 7: sequential cohort flow ===")
flow_ok <- tryCatch({
a01_exprs <- parse("A01.R")
for (e in a01_exprs) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = globalenv())
  }
}
raw_long <- suppressMessages(make_cf_longitudinal_from_sheet1(raw_file)) %>% janitor::clean_names()
lg <- raw_long %>% transmute(patient_id = as.character(patcode), visit_year = year,
                             visit_date = as.Date(visit_date),
                             age_years = suppressWarnings(as.numeric(age_at_visit)),
                             fev1_pp = suppressWarnings(as.numeric(fev1pp)))
lung_cols <- grep("^date_lung_tx", names(s1), value = TRUE)
tx_vec <- purrr::reduce(lapply(lung_cols, function(nm) tryCatch(as.Date(s1[[nm]]), error = function(e) rep(as.Date(NA), nrow(s1)))), dplyr::coalesce)
tx <- tibble(patient_id = pid, tx_date = tx_vec) %>% filter(!is.na(patient_id), !is.na(tx_date)) %>%
  group_by(patient_id) %>% summarise(tx_date = min(tx_date), .groups = "drop")
lg2 <- lg %>% left_join(tx, by = "patient_id") %>% filter(is.na(tx_date) | visit_date < tx_date)
lg3 <- lg2 %>% filter(!is.na(age_years) & age_years >= 6)
ids4 <- lg3 %>% filter(!is.na(fev1_pp)) %>% distinct(patient_id) %>% pull(patient_id)
ids5 <- lg3 %>% filter(!is.na(fev1_pp) & !is.na(visit_year)) %>% group_by(patient_id) %>%
  summarise(n_years = n_distinct(visit_year), .groups = "drop") %>% filter(n_years >= 2) %>% pull(patient_id)
flow <- tibble(
  step = c("Registered in SACFR extract (Sheet1)", "With at least one annual registry record",
           "After removing records on or after lung transplantation", "With a record at age >=6 years",
           "With an FEV1pp value at age >=6 years", "With FEV1pp in >=2 calendar years at age >=6 (analytic cohort)",
           "Modelled (records with known nutritional status; >=2 records)"),
  patients = c(n_distinct(pid[!is.na(pid)]), n_distinct(lg$patient_id), n_distinct(lg2$patient_id),
               n_distinct(lg3$patient_id), length(ids4), length(ids5), n_distinct(dat_model$patient_id)),
  fev1_values = c(NA, sum(!is.na(lg$fev1_pp)), sum(!is.na(lg2$fev1_pp)), sum(!is.na(lg3$fev1_pp)),
                  sum(!is.na(lg3$fev1_pp)), sum(!is.na(lg3$fev1_pp[lg3$patient_id %in% ids5])), sum(!is.na(dat_model$fev1_pp)))
) %>% mutate(excluded_at_step = dplyr::lag(patients) - patients)
show_tbl(flow)
say("Check: analytic cohort IDs equal parquet IDs: ", setequal(ids5, unique(dat$patient_id)))
write_csv(flow, file.path(out_dir, "row7_cohort_flow.csv"))
TRUE
}, error = function(e) { say("ROW 7 failed: ", conditionMessage(e)); FALSE })

# ------------------------------------------------------------------------------
# 6. ROW 16 inputs: genotype detail and registry eligibility (patient level)
# ------------------------------------------------------------------------------
chr <- function(v) { x <- if (v %in% names(s1)) as.character(s1[[v]]) else rep(NA_character_, nrow(s1)); x[is.na(x)] <- ""; x }
m1 <- paste(chr("mutation1"), chr("othermutation1")); m2 <- paste(chr("mutation2"), chr("othermutation2"))
gen <- tibble(
  patient_id = pid,
  c3120_1 = grepl("3120", m1), c3120_2 = grepl("3120", m2),
  f508_1 = grepl("F508|Phe508del", m1), f508_2 = grepl("F508|Phe508del", m2),
  eligible = dplyr::na_if(trimws(chr("cftrm_eligible")), "")
) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, .keep_all = TRUE) %>%
  mutate(carrier_3120 = c3120_1 | c3120_2, homo_3120 = c3120_1 & c3120_2,
         any_f508 = f508_1 | f508_2, carrier_3120_with_f508 = carrier_3120 & any_f508)

# ------------------------------------------------------------------------------
# 7. ROW 10: baseline characteristics by modulator exposure (analytic cohort)
# ------------------------------------------------------------------------------
say("")
say("=== ROW 10: baseline characteristics by modulator exposure (Ever = A01 flag, as Table 1) ===")
bl <- dat %>% arrange(patient_id, visit_date) %>% group_by(patient_id) %>% slice_head(n = 1) %>% ungroup() %>%
  left_join(gen, by = "patient_id") %>%
  mutate(group = if_else(ever_cftr_yes == 1L, "Ever", "Never"),
         public_only = sector == "Public",
         chronic_pa_bl = baseline_pa_cat %in% "Chronic",
         chronic_mssa_bl = sa_status %in% "Chronic",
         underweight_bl = baseline_nutrition_cat %in% "Underweight",
         eligible_yes = eligible %in% "Yes")
p_cont <- function(x, g) tryCatch(t.test(x ~ g)$p.value, error = function(e) NA_real_)
p_cat <- function(x, g) {
  tab <- table(x, g)
  if (nrow(tab) < 2) return(NA_real_)
  ex <- suppressWarnings(chisq.test(tab)$expected)
  if (any(ex < 5)) fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value else chisq.test(tab)$p.value
}
row_cont <- function(v, label) {
  x <- bl[[v]]
  tibble(characteristic = label,
         Never = sprintf("%.1f (%.1f)", mean(x[bl$group == "Never"], na.rm = TRUE), sd(x[bl$group == "Never"], na.rm = TRUE)),
         Ever  = sprintf("%.1f (%.1f)", mean(x[bl$group == "Ever"], na.rm = TRUE),  sd(x[bl$group == "Ever"], na.rm = TRUE)),
         p = p_cont(x, bl$group))
}
row_bin <- function(v, label) {
  x <- as.logical(bl[[v]])
  tibble(characteristic = label,
         Never = n_pct(sum(x[bl$group == "Never"], na.rm = TRUE), sum(!is.na(x[bl$group == "Never"]))),
         Ever  = n_pct(sum(x[bl$group == "Ever"], na.rm = TRUE),  sum(!is.na(x[bl$group == "Ever"]))),
         p = p_cat(x, bl$group))
}
row_fac <- function(v, label) {
  x <- droplevels(factor(bl[[v]]))
  p <- p_cat(x, bl$group)
  bind_rows(lapply(levels(x), function(l) tibble(
    characteristic = paste0(label, ": ", l),
    Never = n_pct(sum(x[bl$group == "Never"] == l, na.rm = TRUE), sum(!is.na(x[bl$group == "Never"]))),
    Ever  = n_pct(sum(x[bl$group == "Ever"] == l, na.rm = TRUE),  sum(!is.na(x[bl$group == "Ever"]))),
    p = if (l == levels(x)[1]) p else NA_real_)))
}
row10 <- bind_rows(
  tibble(characteristic = "Patients", Never = as.character(sum(bl$group == "Never")), Ever = as.character(sum(bl$group == "Ever")), p = NA_real_),
  row_cont("baseline_age_years", "Age at baseline, years, mean (SD)"),
  row_bin("public_only", "Public sector only"),
  row_fac("ethnicity_cat", "Ethnicity"),
  row_fac("f508_genotype_group", "F508del genotype"),
  row_bin("eligible_yes", "Registry: modulator-eligible"),
  row_bin("carrier_3120", "3120+1G>A carrier"),
  row_fac("sex", "Sex"),
  row_cont("baseline_fev1_pp", "Baseline FEV1pp, mean (SD)"),
  row_bin("underweight_bl", "Underweight at baseline"),
  row_bin("baseline_cfrd_indicator", "CFRD at baseline"),
  row_bin("baseline_pancreatic_insuff", "Pancreatic insufficiency at baseline"),
  row_bin("chronic_pa_bl", "Chronic P. aeruginosa at baseline"),
  row_bin("chronic_mssa_bl", "Chronic MSSA at baseline")
) %>% mutate(p = signif(p, 3))
show_tbl(row10)
write_csv(row10, file.path(out_dir, "row10_baseline_by_modulator.csv"))

# ------------------------------------------------------------------------------
# 8. ROW 12: respiratory samples per patient-year
# ------------------------------------------------------------------------------
say("")
say("=== ROW 12: respiratory samples per patient-year (registry category 5 = five or more) ===")
num_cols <- grep("^resp_sample_number_?20[0-9]{2}$", names(s1), value = TRUE)
any_cols <- grep("^resp_samples_?20[0-9]{2}$", names(s1), value = TRUE)
samp <- bind_rows(lapply(2018:2024, function(y) {
  nc <- num_cols[grepl(as.character(y), num_cols)]; ac <- any_cols[grepl(as.character(y), any_cols)]
  nraw <- if (length(nc)) tolower(trimws(as.character(s1[[nc[1]]]))) else rep(NA_character_, nrow(s1))
  araw <- if (length(ac)) tolower(trimws(as.character(s1[[ac[1]]]))) else rep(NA_character_, nrow(s1))
  num <- suppressWarnings(as.numeric(nraw)); num[nraw %in% "5 or more"] <- 5
  anyv <- dplyr::case_when(araw %in% c("yes", "1") ~ TRUE, araw %in% c("no", "0") ~ FALSE, TRUE ~ NA)
  tibble(patient_id = pid, visit_year = y,
         samples = dplyr::case_when(anyv %in% FALSE ~ 0, !is.na(num) ~ pmin(num, 5), TRUE ~ NA_real_))
})) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, visit_year, .keep_all = TRUE)
py <- dat %>% filter(!is.na(fev1_pp)) %>% distinct(patient_id, visit_year, sector) %>%
  left_join(samp, by = c("patient_id", "visit_year")) %>%
  mutate(sector_group = if_else(sector == "Public", "Public only", "Private/Combined"))
samp_summary <- function(d, label) tibble(
  group = label, patient_years = nrow(d), with_sample_data = n_pct(sum(!is.na(d$samples)), nrow(d)),
  samples_median_iqr = med_iqr(d$samples, 0),
  zero_samples = n_pct(sum(d$samples == 0, na.rm = TRUE), sum(!is.na(d$samples))),
  ge4_samples = n_pct(sum(d$samples >= 4, na.rm = TRUE), sum(!is.na(d$samples))))
row12 <- bind_rows(samp_summary(py, "All patient-years"),
                   samp_summary(filter(py, sector_group == "Public only"), "Public only"),
                   samp_summary(filter(py, sector_group == "Private/Combined"), "Private/Combined"))
show_tbl(row12)
say("Distribution of samples per patient-year (5 = five or more):")
show_tbl(as.data.frame(table(samples = py$samples, useNA = "ifany")))
say("Share with >=4 samples by registry year:")
show_tbl(py %>% group_by(visit_year) %>% summarise(patient_years = n(), with_data = sum(!is.na(samples)),
                                                   pct_ge4 = round(100 * mean(samples >= 4, na.rm = TRUE), 1), .groups = "drop"))
write_csv(row12, file.path(out_dir, "row12_resp_samples.csv"))

# ------------------------------------------------------------------------------
# 9. ROW 16: 3120+1G>A and modulator eligibility (analytic cohort)
# ------------------------------------------------------------------------------
say("")
say("=== ROW 16: 3120+1G>A carriage and registry modulator eligibility (analytic cohort) ===")
pt16 <- bl %>% select(patient_id, ethnicity_cat, sector, carrier_3120, homo_3120, carrier_3120_with_f508, any_f508, eligible)
g16 <- function(d, label) tibble(
  group = label, patients = nrow(d),
  carrier_3120 = n_pct(sum(d$carrier_3120, na.rm = TRUE), nrow(d)),
  homozygous_3120 = n_pct(sum(d$homo_3120, na.rm = TRUE), nrow(d)),
  carrier_with_F508del = n_pct(sum(d$carrier_3120_with_f508, na.rm = TRUE), nrow(d)),
  no_F508del_allele = n_pct(sum(!d$any_f508, na.rm = TRUE), nrow(d)),
  registry_ineligible = n_pct(sum(d$eligible %in% "No"), sum(!is.na(d$eligible))),
  ineligible_carrying_3120 = n_pct(sum(d$eligible %in% "No" & d$carrier_3120), max(1, sum(d$eligible %in% "No"))))
row16 <- bind_rows(
  g16(pt16, "All"),
  bind_rows(lapply(levels(droplevels(pt16$ethnicity_cat)), function(l) g16(filter(pt16, ethnicity_cat == l), paste("Ethnicity:", l)))),
  g16(filter(pt16, sector == "Public"), "Sector: Public only"),
  g16(filter(pt16, sector != "Public"), "Sector: Private/Combined"))
show_tbl(row16)
write_csv(row16, file.path(out_dir, "row16_3120_eligibility.csv"))

# ------------------------------------------------------------------------------
# 10. ROW 17: Benjamini-Hochberg adjusted p-values (within each model)
# ------------------------------------------------------------------------------
say("")
say("=== ROW 17: Benjamini-Hochberg adjusted p-values (within model, as A02) ===")
bh <- function(mod, label) broom.mixed::tidy(mod, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(model = label, p.value.BH = p.adjust(p.value, method = "BH")) %>%
  select(model, term, estimate, p.value, p.value.BH)
m2_tidy <- read_csv("outputs/tables/.model_m2_reml_tidy.csv", show_col_types = FALSE) %>%
  filter(term != "(Intercept)") %>%
  transmute(model = "M2 (submitted, REML)", term, estimate, p.value, p.value.BH = p.adjust(p.value, method = "BH"))
c2 <- readRDS(file.path(out_dir, "models", "C2_cftr_segmented.rds"))
row17 <- bind_rows(m2_tidy, bh(c2, "C2 CFTR segmented"),
                   bh(e1_models[["P. aeruginosa"]], "E1 P. aeruginosa"), bh(e1_models[["MSSA"]], "E1 MSSA"))
say("Terms with p < 0.05 and their BH-adjusted values:")
show_tbl(row17 %>% filter(p.value < 0.05) %>% mutate(across(c(estimate, p.value, p.value.BH), ~ signif(.x, 3))))
write_csv(row17, file.path(out_dir, "row17_fdr_pvalues.csv"))

say("")
say("=== Fit log ===")
show_tbl(bind_rows(fit_log) %>% select(model, re_structure, n_patients, n_obs, n_warnings))
say("Done. Outputs in ", normalizePath(out_dir))
