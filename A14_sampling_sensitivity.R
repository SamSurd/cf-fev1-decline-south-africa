#!/usr/bin/env Rscript
# ==============================================================================
# A14 - INFECTION MODELS IN WELL-SAMPLED PATIENTS (R1.6 sensitivity, corrected data)
# Well-sampled = >=4 respiratory samples in at least half of the patient's
# FEV1 years (registry category 5 = five or more). Writes outputs/revision/ only.
# ==============================================================================
set.seed(20251111)
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(lme4); library(lmerTest); library(readxl); library(janitor) })
out_dir <- "outputs/revision"; summary_file <- file.path(out_dir, "A14_summary.txt"); fit_log <- list()
for (e in parse("A08_revision_fixes.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, envir = globalenv())
}
file.create(summary_file)
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
flags_main <- compute_use_flags(dat_model)

s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_df <- s1[grep("^patcode", names(s1), value = TRUE)]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
num_cols <- grep("^resp_sample_number_?20[0-9]{2}$", names(s1), value = TRUE)
any_cols <- grep("^resp_samples_?20[0-9]{2}$", names(s1), value = TRUE)
samp <- bind_rows(lapply(2018:2024, function(y) {
  nc <- num_cols[grepl(as.character(y), num_cols)]; ac <- any_cols[grepl(as.character(y), any_cols)]
  nraw <- if (length(nc)) tolower(trimws(as.character(s1[[nc[1]]]))) else rep(NA_character_, nrow(s1))
  araw <- if (length(ac)) tolower(trimws(as.character(s1[[ac[1]]]))) else rep(NA_character_, nrow(s1))
  numv <- suppressWarnings(as.numeric(nraw)); numv[nraw %in% "5 or more"] <- 5
  anyv <- dplyr::case_when(araw %in% c("yes", "1") ~ TRUE, araw %in% c("no", "0") ~ FALSE, TRUE ~ NA)
  tibble(patient_id = pid, visit_year = y, samples = dplyr::case_when(anyv %in% FALSE ~ 0, !is.na(numv) ~ pmin(numv, 5), TRUE ~ NA_real_))
})) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, visit_year, .keep_all = TRUE)
well <- dat_model %>% distinct(patient_id, visit_year) %>% left_join(samp, by = c("patient_id", "visit_year")) %>%
  group_by(patient_id) %>% summarise(share_ge4 = mean(samples >= 4, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(share_ge4), share_ge4 >= 0.5) %>% pull(patient_id)
say("Well-sampled patients: ", length(well), " of ", n_distinct(dat_model$patient_id))

specs <- list(list(name = "pa", start_col = "pa_chronic_start_date", own = "chronic_paeruginosa", label = "P. aeruginosa"),
              list(name = "mssa", start_col = "mssa_chronic_start_date", own = "chronic_mssa", label = "MSSA"))
rows <- list()
for (spec in specs) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index"); iv <- paste0(tv, ":", pv)
  d1 <- dat_model %>% mutate(event_start = .data[[spec$start_col]], anchor_date = dplyr::coalesce(event_start, baseline_date),
                             !!tv := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
                             !!pv := as.integer(!is.na(event_start) & visit_date >= event_start))
  d1w <- d1 %>% filter(patient_id %in% well) %>% droplevels()
  rhs <- rhs_drop(build_rhs_event(compute_use_flags(d1w), tv, pv), spec$own)
  ma <- fit_triple_fallback(rhs_drop(build_rhs_event(flags_main, tv, pv), spec$own), d1, TRUE, paste(spec$label, "all patients"))
  mw <- fit_triple_fallback(rhs, d1w, TRUE, paste(spec$label, "well-sampled"))
  rows[[spec$label]] <- bind_rows(seg_rows(ma, paste(spec$label, "all patients"), tv, iv, pv),
                                  seg_rows(mw, paste(spec$label, "well-sampled"), tv, iv, pv)) %>%
    filter(contrast %in% c("Level change at event (pp)", "Slope change after event (pp/yr)", "Post-event slope, net (pp/yr)"))
}
res <- bind_rows(rows)
show_tbl(res %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
show_tbl(bind_rows(fit_log) %>% select(model, re_structure, n_patients, n_obs))
write_csv(res, file.path(out_dir, "row12_infection_models_well_sampled.csv"))
