#!/usr/bin/env Rscript
# ==============================================================================
# A18 - CHRONIC INFECTION STATUS AND THE A01 COLUMN-NAME GAPS (R1.6)
# ------------------------------------------------------------------------------
# 1. The registry records chronic infection as "unclear as less than 4 samples
#    tested" when a patient has no history of chronic infection and fewer than
#    four samples were taken in the year (modified Leeds criteria; SACFR Annual
#    Report 2018, Appendix 1). A01 codes this, and blanks, as not chronic.
#    Counts, for the modelled patient-years, of yes / no / unclear / blank, by
#    year, samples and sector.
# 2. A19 found that A01 misses some yearly registry columns whose names have no
#    underscore before the year (2019 P. aeruginosa and MSSA; 2018-2019 MRSA,
#    B. cepacia complex and diabetes/CFRD). This script reads every year with
#    the full set of names and counts the modelled patient-years and patients
#    whose coding would change if A01 were fixed. Nothing is re-fitted.
# 3. CFRD and pancreatic insufficiency: A01 uses each year's recorded value
#    (no carry-forward); counts of later "no" after a "yes".
# Aggregate outputs only. Reads the same inputs as A14.
# Writes: outputs/revision/A18_summary.txt, outputs/revision/row19_chronic_status_raw.csv,
#         outputs/revision/row19_a01_gap_impact.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(readxl); library(janitor) })
out_dir <- "outputs/revision"; summary_file <- file.path(out_dir, "A18_summary.txt"); file.create(summary_file)
say <- function(...) { msg <- paste0(...); message(msg); con <- file(summary_file, open = "a"); writeLines(msg, con); close(con) }
show_tbl <- function(x) { for (l in capture.output(print(as.data.frame(x), row.names = FALSE))) say(l) }
pct <- function(n, d) round(100 * n / d, 1)

# Modelled patient-years, as A02/A14 ------------------------------------------
dat_model <- arrow::read_parquet("outputs/data_derived/analysis_long.parquet") %>% as_tibble() %>%
  mutate(nutrition_cat = as.character(nutrition_cat)) %>%
  filter(nutrition_cat != "Unknown" | is.na(nutrition_cat)) %>%
  group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup() %>%
  filter(!is.na(fev1_pp)) %>% mutate(visit_year = as.integer(visit_year), sector = as.character(sector))

# Raw registry fields, one row per patient-year, all name variants ---------------
s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_df <- s1[grep("^patcode", names(s1), value = TRUE)]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
pick <- function(stems, y) { cands <- c(paste0(stems, "_", y), paste0(stems, y)); x <- cands[cands %in% names(s1)]; if (length(x)) x[1] else NA_character_ }
col_or_na <- function(nm) if (is.na(nm)) rep(NA_character_, nrow(s1)) else tolower(trimws(as.character(s1[[nm]])))
num_cols <- grep("^resp_sample_number_?20[0-9]{2}$", names(s1), value = TRUE)
any_cols <- grep("^resp_samples_?20[0-9]{2}$", names(s1), value = TRUE)
raw <- bind_rows(lapply(2018:2024, function(y) {
  nc <- num_cols[grepl(as.character(y), num_cols)]; ac <- any_cols[grepl(as.character(y), any_cols)]
  nraw <- if (length(nc)) tolower(trimws(as.character(s1[[nc[1]]]))) else rep(NA_character_, nrow(s1))
  araw <- if (length(ac)) tolower(trimws(as.character(s1[[ac[1]]]))) else rep(NA_character_, nrow(s1))
  numv <- suppressWarnings(as.numeric(nraw)); numv[nraw %in% "5 or more"] <- 5
  anyv <- dplyr::case_when(araw %in% c("yes", "1") ~ TRUE, araw %in% c("no", "0") ~ FALSE, TRUE ~ NA)
  tibble(patient_id = pid, visit_year = y,
         pa_raw   = col_or_na(pick(c("chronicpseudomaer", "chronic_pseudomonas"), y)),
         sa_raw   = col_or_na(pick("chronicstaphylaur", y)),
         mrsa_raw = col_or_na(pick("chronic_mrsa", y)),
         bcep_raw = col_or_na(pick("chronic_burkholderia", y)),
         cfrd_raw = col_or_na(pick(c("cfrd", "diabetes"), y)),
         pan_raw  = col_or_na(pick(c("pancreatic_status", "pancreatic"), y)),
         samples  = dplyr::case_when(anyv %in% FALSE ~ 0, !is.na(numv) ~ pmin(numv, 5), TRUE ~ NA_real_))
})) %>% filter(!is.na(patient_id)) %>% distinct(patient_id, visit_year, .keep_all = TRUE)

# Same yes/no rule as A01 (status_from_any_chronic and the CFRD recode)
classify <- function(v) {
  num <- suppressWarnings(as.numeric(v))
  dplyr::case_when(is.na(v) | v == "" ~ "blank",
                   v %in% c("1", "yes", "y", "true", "chronic") | (!is.na(num) & num >= 1) ~ "yes",
                   v %in% c("0", "no", "n", "false") | (!is.na(num) & num == 0) ~ "no",
                   TRUE ~ "unclear")
}
py <- dat_model %>%
  select(patient_id, visit_year, sector, chronic_paeruginosa, chronic_mssa, chronic_mrsa, chronic_bcepacia,
         cfrd_indicator, pancreatic_insufficiency) %>%
  left_join(raw, by = c("patient_id", "visit_year")) %>%
  mutate(pa = classify(pa_raw), sa = classify(sa_raw), mrsa = classify(mrsa_raw), bcep = classify(bcep_raw),
         cfrd_c = classify(cfrd_raw),
         pan_c = dplyr::case_when(pan_raw %in% c("insufficient", "yes", "1", "true", "y") ~ "yes",
                                  is.na(pan_raw) | pan_raw == "" ~ "blank", TRUE ~ "no"),
         samples_cat = dplyr::case_when(is.na(samples) ~ "not recorded", samples >= 4 ~ ">=4", TRUE ~ "<4"))
n_py <- nrow(py)
say("A18 chronic infection status and A01 column-name gaps, run ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Modelled patient-years: ", n_py, " (", n_distinct(py$patient_id), " patients)")

# 1. Unclear / blank chronic status ---------------------------------------------
say(""); say("=== 1. Chronic status as recorded (all name variants) ===")
out <- list()
for (sp in list(list(col = "pa", label = "P. aeruginosa"), list(col = "sa", label = "MSSA"))) {
  d <- py %>% mutate(status = .data[[sp$col]], unk = status %in% c("unclear", "blank"))
  by_year <- d %>% count(visit_year, status) %>% pivot_wider(names_from = status, values_from = n, values_fill = 0)
  for (k in c("yes", "no", "unclear", "blank")) if (!k %in% names(by_year)) by_year[[k]] <- 0L
  by_year <- by_year %>% mutate(total = yes + no + unclear + blank, pct_unclear = pct(unclear, total), pct_blank = pct(blank, total))
  say(""); say(sp$label, ": status by year"); show_tbl(by_year %>% select(visit_year, yes, no, unclear, blank, total, pct_unclear, pct_blank))
  say(sp$label, " overall: unclear ", sum(d$status == "unclear"), " (", pct(sum(d$status == "unclear"), n_py), "%), blank ",
      sum(d$status == "blank"), " (", pct(sum(d$status == "blank"), n_py), "%), yes ", sum(d$status == "yes"), ", no ", sum(d$status == "no"))
  say(sp$label, ": status by respiratory samples in the year")
  show_tbl(d %>% count(samples_cat, status) %>% group_by(samples_cat) %>% mutate(pct = pct(n, sum(n))) %>% ungroup() %>%
             pivot_wider(names_from = status, values_from = c(n, pct), values_fill = 0))
  say(sp$label, ": unclear or blank by sector")
  show_tbl(d %>% group_by(sector) %>% summarise(patient_years = n(), unclear = sum(status == "unclear"), blank = sum(status == "blank"),
                                               pct_unclear_or_blank = pct(sum(unk), n()), .groups = "drop"))
  pl <- d %>% arrange(patient_id, visit_year) %>% group_by(patient_id) %>%
    summarise(ever_yes = any(status == "yes"), first_yes = if (any(status == "yes")) min(visit_year[status == "yes"]) else NA_integer_,
              n_unk_after = if (any(status == "yes")) sum(unk & visit_year > min(visit_year[status == "yes"])) else 0L,
              n_no_after = if (any(status == "yes")) sum(status == "no" & visit_year > min(visit_year[status == "yes"])) else 0L,
              all_unk = all(unk), .groups = "drop")
  say(sp$label, ": patients ever recorded chronic ", sum(pl$ever_yes), "; unclear or blank in every modelled year ", sum(pl$all_unk),
      "; after the first chronic year: ", sum(pl$n_unk_after), " unclear/blank and ", sum(pl$n_no_after), " 'no' patient-years (",
      sum(pl$n_unk_after > 0 | pl$n_no_after > 0), " patients), coded not chronic in the same-year term")
  out[[sp$label]] <- by_year %>% mutate(pathogen = sp$label, .before = 1)
}
write_csv(bind_rows(out), file.path(out_dir, "row19_chronic_status_raw.csv"))

# 2. Patient-years whose coding would change if A01 read every year ----------------
say(""); say("=== 2. Recorded 'yes' but coded 0 in the analysis data (A01 column-name gaps) ===")
imp <- bind_rows(
  py %>% transmute(patient_id, variable = "Chronic P. aeruginosa", visit_year, rec = pa == "yes", coded = chronic_paeruginosa == 1L),
  py %>% transmute(patient_id, variable = "Chronic MSSA", visit_year, rec = sa == "yes", coded = chronic_mssa == 1L),
  py %>% transmute(patient_id, variable = "Chronic MRSA", visit_year, rec = mrsa == "yes", coded = chronic_mrsa == 1L),
  py %>% transmute(patient_id, variable = "Chronic B. cepacia complex", visit_year, rec = bcep == "yes", coded = chronic_bcepacia == 1L),
  py %>% transmute(patient_id, variable = "CFRD", visit_year, rec = cfrd_c == "yes", coded = cfrd_indicator == 1L),
  py %>% transmute(patient_id, variable = "Pancreatic insufficiency", visit_year, rec = pan_c == "yes", coded = pancreatic_insufficiency == 1L))
gap <- imp %>% group_by(variable, visit_year) %>%
  summarise(recorded_yes = sum(rec), coded_1 = sum(coded), yes_coded_0 = sum(rec & !coded), coded_1_not_yes = sum(coded & !rec), .groups = "drop")
show_tbl(gap %>% filter(yes_coded_0 > 0 | coded_1_not_yes > 0))
say("Totals (modelled patient-years recorded yes but coded 0):")
show_tbl(gap %>% group_by(variable) %>% summarise(recorded_yes = sum(recorded_yes), yes_coded_0 = sum(yes_coded_0), .groups = "drop"))
write_csv(gap, file.path(out_dir, "row19_a01_gap_impact.csv"))

# Event classes (first modelled year), as coded now and as recorded
cls <- function(flag, yr) { if (!any(flag)) "never" else if (flag[which.min(yr)]) "prevalent" else "incident" }
say(""); say("Event classes over modelled patient-years (prevalent = chronic in the first modelled year):")
for (v in c("Chronic P. aeruginosa", "Chronic MSSA", "Chronic MRSA", "Chronic B. cepacia complex", "CFRD")) {
  dv <- imp %>% filter(variable == v)
  pc <- dv %>% group_by(patient_id) %>% summarise(now = cls(coded, visit_year), fixed = cls(rec, visit_year),
                                                  onset_now = if (any(coded)) min(visit_year[coded]) else NA_integer_,
                                                  onset_fixed = if (any(rec)) min(visit_year[rec]) else NA_integer_, .groups = "drop")
  say(v, ": ", paste(capture.output(print(table(now = pc$now, fixed = pc$fixed))), collapse = "\n"))
  say(v, ": onset year changes for ", sum(!is.na(pc$onset_fixed) & (is.na(pc$onset_now) | pc$onset_now != pc$onset_fixed)), " patients")
}

# 3. CFRD and pancreatic insufficiency: same-year coding ------------------------------
say(""); say("=== 3. CFRD and pancreatic insufficiency after a recorded 'yes' (all name variants) ===")
for (v in list(list(col = "cfrd_c", label = "CFRD"), list(col = "pan_c", label = "Pancreatic insufficiency"))) {
  pl <- py %>% mutate(s = .data[[v$col]]) %>% arrange(patient_id, visit_year) %>% group_by(patient_id) %>%
    summarise(ever = any(s == "yes"),
              later_no = if (any(s == "yes")) sum(s == "no" & visit_year > min(visit_year[s == "yes"])) else 0L,
              later_blank = if (any(s == "yes")) sum(s == "blank" & visit_year > min(visit_year[s == "yes"])) else 0L, .groups = "drop")
  say(v$label, ": patients ever 'yes' ", sum(pl$ever), "; later 'no' patient-years ", sum(pl$later_no), " (", sum(pl$later_no > 0),
      " patients); later blank patient-years ", sum(pl$later_blank), " (", sum(pl$later_blank > 0), " patients)")
}
say(""); say("Saved row19_chronic_status_raw.csv and row19_a01_gap_impact.csv")
