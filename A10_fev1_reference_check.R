#!/usr/bin/env Rscript
# ==============================================================================
# A10 - FEV1 REFERENCE-EQUATION CHECK (JCF-D-26-00377 revision)
# ------------------------------------------------------------------------------
# A01 builds fev1_pp as coalesce(GLI % predicted, local-reference % predicted).
# GLI is available for ~25% of patient-years in 2018-2023 and ~all of 2024, and
# the 2023 local-reference column (x2023_fev1_percent_pred_local_ref) is not
# matched by A01's name list. This script quantifies:
#   1. which source each analysed FEV1 value came from, by year;
#   2. how the main slope estimates move when a source indicator is added;
#   3. whether raw spirometry + height exist to recompute one reference throughout.
# Writes: outputs/revision/A10_* only.
# ==============================================================================

set.seed(20251111)
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(lme4); library(lmerTest); library(readxl); library(janitor)
})
out_dir <- "outputs/revision"
summary_file <- file.path(out_dir, "A10_summary.txt")
fit_log <- list()

# Re-use the helper functions defined in A08 (function definitions only)
for (e in parse("A08_revision_fixes.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = globalenv())
  }
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
dat_model <- dat_model %>% left_join(pancreas_start, by = "patient_id") %>% left_join(cfrd_start, by = "patient_id")
flags_main <- compute_use_flags(dat_model)

# 1. Source of each FEV1 value (A01 rule: GLI if present, else local) ---------
s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_cols <- grep("^patcode", names(s1), value = TRUE)
id_df <- s1[id_cols]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
num <- function(v) if (v %in% names(s1)) suppressWarnings(as.numeric(s1[[v]])) else rep(NA_real_, nrow(s1))
src <- bind_rows(lapply(2018:2024, function(y) {
  gli <- dplyr::coalesce(num(paste0("gli_fev1_pp_", y)), num(paste0("gli_percpredfev", y)))
  loc <- dplyr::coalesce(num(paste0("fev1_percent_pred_local_ref_", y)), num(paste0("fev1_pp_", y)),
                         num(paste0("fev1_local_", y)), num(paste0("x", y, "_fev1_percent_pred_local_ref")))
  tibble(patient_id = pid, visit_year = y, gli = gli, local = loc)
})) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, visit_year, .keep_all = TRUE)

dm <- dat_model %>% left_join(src, by = c("patient_id", "visit_year")) %>%
  mutate(fev1_src_gli = as.integer(!is.na(gli) & abs(fev1_pp - gli) < 0.01))
say("Analysed FEV1 values by source and year (A01 rule):")
show_tbl(dm %>% group_by(visit_year) %>%
           summarise(n = n(), from_gli = sum(fev1_src_gli), pct_gli = round(100 * mean(fev1_src_gli), 1),
                     local_available = sum(!is.na(local)), .groups = "drop"))
say("Values matching neither source: ", sum(!(abs(dm$fev1_pp - dm$gli) < 0.01 | abs(dm$fev1_pp - dm$local) < 0.01), na.rm = TRUE))
miss23 <- src %>% filter(visit_year == 2023, !is.na(local), is.na(gli)) %>% semi_join(distinct(dat, patient_id), by = "patient_id") %>% nrow()
say("2023 patient-years in the analytic cohort with a local-reference value but no GLI value (dropped by A01): ", miss23)

# 2. Do the key estimates move when the source is adjusted for? ---------------
say("")
say("=== Key estimates without and with a GLI-source indicator ===")
res <- list()
m1a <- fit_triple_fallback(build_rhs(flags_main, "m1"), dm, FALSE, "M1")
m1b <- fit_triple_fallback(paste(build_rhs(flags_main, "m1"), "+ fev1_src_gli"), dm, FALSE, "M1 + source")
res$m1 <- bind_rows(lc(m1a, c(time_cs = 1), "Adjusted annual slope (M1)") %>% mutate(model = "M1", .before = 1),
                    lc(m1b, c(time_cs = 1), "Adjusted annual slope (M1)") %>% mutate(model = "M1 + source", .before = 1),
                    lc(m1b, c(fev1_src_gli = 1), "GLI-source offset (pp)") %>% mutate(model = "M1 + source", .before = 1))

# CFTR segmented model (A08 C2 specification)
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
starts <- tibble(patient_id = pid, cftr_start_rev = srev) %>% filter(!is.na(patient_id)) %>% group_by(patient_id) %>%
  summarise(cftr_start_rev = suppressWarnings(min(cftr_start_rev, na.rm = TRUE)), .groups = "drop") %>%
  mutate(cftr_start_rev = dplyr::if_else(is.finite(as.numeric(cftr_start_rev)), cftr_start_rev, as.Date(NA)))
dc <- dm %>% left_join(starts, by = "patient_id") %>%
  mutate(cftr_post = as.integer(!is.na(cftr_start_rev) & visit_date >= cftr_start_rev),
         cftr_tsi = dplyr::if_else(cftr_post == 1L, as.numeric(visit_date - cftr_start_rev) / 365.25, 0))
prev_users <- dc %>% group_by(patient_id) %>%
  summarise(pre = sum(!is.na(cftr_start_rev) & cftr_post == 0L), post = sum(cftr_post == 1L), .groups = "drop") %>%
  filter(post > 0 & pre == 0) %>% pull(patient_id)
dc <- dc %>% filter(!patient_id %in% prev_users)
rhs_c2 <- paste(c("time_cs", "cftr_post", "cftr_tsi", event_level_covariates(flags_main)), collapse = " + ")
c2a <- fit_triple_fallback(rhs_c2, dc, TRUE, "C2")
c2b <- fit_triple_fallback(paste(rhs_c2, "+ fev1_src_gli"), dc, TRUE, "C2 + source")
res$c2 <- bind_rows(seg_rows(c2a, "C2", "time_cs", "cftr_tsi", "cftr_post"),
                    seg_rows(c2b, "C2 + source", "time_cs", "cftr_tsi", "cftr_post")) %>%
  filter(contrast %in% c("Pre-event slope (pp/yr)", "Level change at event (pp)", "Slope change after event (pp/yr)", "Post-event slope, net (pp/yr)"))

# Infection models, all patients (A08 E1 specification)
for (spec in list(list(name = "pa", start_col = "pa_chronic_start_date", own = "chronic_paeruginosa", label = "P. aeruginosa"),
                  list(name = "mssa", start_col = "mssa_chronic_start_date", own = "chronic_mssa", label = "MSSA"))) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index"); iv <- paste0(tv, ":", pv)
  d1 <- dm %>% mutate(event_start = .data[[spec$start_col]], anchor_date = dplyr::coalesce(event_start, baseline_date),
                      !!tv := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
                      !!pv := as.integer(!is.na(event_start) & visit_date >= event_start))
  rhs1 <- rhs_drop(build_rhs_event(flags_main, tv, pv), spec$own)
  ea <- fit_triple_fallback(rhs1, d1, TRUE, paste(spec$label, "E1"))
  eb <- fit_triple_fallback(paste(rhs1, "+ fev1_src_gli"), d1, TRUE, paste(spec$label, "E1 + source"))
  res[[spec$name]] <- bind_rows(seg_rows(ea, paste(spec$label, "E1"), tv, iv, pv),
                                seg_rows(eb, paste(spec$label, "E1 + source"), tv, iv, pv)) %>%
    filter(contrast %in% c("Pre-event slope (pp/yr)", "Slope change after event (pp/yr)", "Post-event slope, net (pp/yr)"))
}
res_tab <- bind_rows(res)
show_tbl(res_tab %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
write_csv(res_tab, file.path(out_dir, "A10_source_adjusted_estimates.csv"))

# 3. Can one reference be recomputed for every year? -------------------------
say("")
say("=== Raw inputs for recomputing GLI (non-missing per year) ===")
inputs <- bind_rows(lapply(2018:2024, function(y) {
  pick <- function(pat) { v <- grep(pat, names(s1), value = TRUE); v <- v[grepl(as.character(y), v)]; if (length(v)) v[1] else NA_character_ }
  cols <- c(fev1_litres = pick("fev1.*lit|lit.*fev1|^valfev1"), height = pick("height|^ht_|_ht"), weight = pick("weight|^wt_|_wt"))
  tibble(year = y, fev1_litres_col = cols[["fev1_litres"]], n_fev1_litres = if (is.na(cols[["fev1_litres"]])) 0L else sum(!is.na(num(cols[["fev1_litres"]]))),
         height_col = cols[["height"]], n_height = if (is.na(cols[["height"]])) 0L else sum(!is.na(num(cols[["height"]]))))
}))
show_tbl(inputs)
say("Height-like columns found: ", paste(head(grep("height|^ht|_ht", names(s1), value = TRUE), 20), collapse = ", "))
say("Done.")
