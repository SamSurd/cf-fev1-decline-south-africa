#!/usr/bin/env Rscript
# ==============================================================================
# A19 - AUDIT OF YEAR-SPECIFIC COLUMN LOOKUPS IN A01 (corrected data)
# ------------------------------------------------------------------------------
# A01's make_cf_longitudinal_from_sheet1() finds each yearly field by exact
# name (first_present). Registry column names are not consistent across years
# (e.g. chronicpseudomaer_2018 but chronicpseudomaer2019), so a lookup can
# return nothing although the field exists. For every year-specific lookup in
# A01 this script reports the column A01 uses, or, when it finds none, columns
# that match once underscores are ignored or with an x<year>_ prefix, and the
# value categories of the fields involved. Headers and category counts only.
# Writes: outputs/revision/A19_summary.txt, outputs/revision/row20_column_audit.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(readxl); library(janitor) })
out_dir <- "outputs/revision"; summary_file <- file.path(out_dir, "A19_summary.txt"); file.create(summary_file)
say <- function(...) { msg <- paste0(...); message(msg); con <- file(summary_file, open = "a"); writeLines(msg, con); close(con) }
show_tbl <- function(x) { for (l in capture.output(print(as.data.frame(x), row.names = FALSE))) say(l) }

s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
hdr <- names(s1); norm <- function(x) gsub("_", "", x)

# Candidate lists copied from A01 (make_cf_longitudinal_from_sheet1)
lookups <- list(
  visit_date = function(y) c(paste0("date_fev1_", y), paste0("date_fev1", y), paste0("date_fev_", y)),
  age        = function(y) c(paste0("age_in_years_", y), paste0("age_years_at_lft_", y), paste0("age_", y), paste0("age_years_", y)),
  sector     = function(y) c(paste0("health_sector_", y), paste0("healthsector_", y), paste0("sector_", y), paste0("public_private_", y)),
  fev1_gli   = function(y) c(paste0("gli_fev1_pp_", y), paste0("gli_percpredfev", y)),
  fev1_local = function(y) c(paste0("fev1_percent_pred_local_ref_", y), paste0("fev1_pp_", y), paste0("fev1_local_", y), paste0("x", y, "_fev1_percent_pred_local_ref")),
  bmi        = function(y) c(paste0("bmi_", y), paste0("bmi", y)),
  bmiz       = function(y) c(paste0("bmiz_", y), paste0("bmiz", y)),
  pancreatic = function(y) c(paste0("pancreatic_status_", y), paste0("pancreatic_", y)),
  cfrd       = function(y) c(paste0("cfrd_", y), paste0("diabetes_", y)),
  pa_any     = function(y) c(paste0("pseudo_any_", y), paste0("ever_pseudomonas_", y)),
  pa_chronic = function(y) c(paste0("chronicpseudomaer_", y), paste0("chronic_pseudomonas_", y)),
  mssa_any   = function(y) c(paste0("any_mssa_", y)),
  mssa_chronic = function(y) c(paste0("chronicstaphylaur_", y)),
  bcep_any   = function(y) c(paste0("any_burkholderia_", y)),
  bcep_chronic = function(y) c(paste0("chronic_burkholderia_", y)),
  mrsa_any   = function(y) c(paste0("any_mrsa_", y)),
  mrsa_chronic = function(y) c(paste0("chronic_mrsa_", y)),
  exacerbations = function(y) c(paste0("exacerbations_per_year_", y), paste0("ivi_number", y), paste0("ivi_number_", y)),
  cftr_drug  = function(y) c(paste0("cftrm_drug_", y), paste0("cftr_modulator_", y))
)
rows <- list()
for (v in names(lookups)) for (y in 2018:2024) {
  cand <- lookups[[v]](y)
  used <- cand[cand %in% hdr][1]
  alt <- if (is.na(used)) {
    a <- hdr[norm(hdr) %in% norm(cand) | hdr %in% paste0("x", y, "_", sub(paste0("_?", y, "$"), "", cand))]
    if (length(a)) paste(a, collapse = ", ") else NA_character_
  } else NA_character_
  rows[[length(rows) + 1]] <- tibble(variable = v, year = y, used = used, missed_but_present = alt,
                                     n_nonblank = if (!is.na(used)) sum(!is.na(s1[[used]]) & trimws(as.character(s1[[used]])) != "") else NA_integer_)
}
audit <- bind_rows(rows)
write_csv(audit, file.path(out_dir, "row20_column_audit.csv"))
say("A19 audit of A01 year-specific lookups, run ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Sheet1: ", nrow(s1), " rows, ", length(hdr), " columns")
say(""); say("=== Lookups that find no column ===")
show_tbl(audit %>% filter(is.na(used)) %>% select(variable, year, missed_but_present))
say(""); say("=== Column used per variable and year (non-blank entries) ===")
show_tbl(audit %>% mutate(cell = ifelse(is.na(used), "-", paste0(used, " (", n_nonblank, ")"))) %>%
           select(variable, year, cell) %>% pivot_wider(names_from = year, values_from = cell))

# Value categories (counts only) for the fields that matter for the gaps
cat_tbl <- function(cols) {
  cols <- cols[cols %in% hdr]
  bind_rows(lapply(cols, function(cn) {
    x <- trimws(tolower(as.character(s1[[cn]]))); x[is.na(x) | x == ""] <- "(blank)"
    tibble(column = cn, value = x) %>% count(column, value)
  }))
}
say(""); say("=== Categories: diabetes / CFRD fields ===")
show_tbl(cat_tbl(grep("^(cfrd|diabetes)_?20[0-9]{2}$|^x20[0-9]{2}_(cfrd|diabetes)", hdr, value = TRUE)))
say(""); say("=== Categories: chronic infection fields, 2018-2019 (names as in the registry) ===")
show_tbl(cat_tbl(grep("^(chronicpseudomaer|chronicstaphylaur|chronic_mrsa|chronic_burkholderia)_?(2018|2019)$", hdr, value = TRUE)))
say(""); say("=== Categories: any-infection fields, 2018-2019 ===")
show_tbl(cat_tbl(grep("^(pseudo_any|ever_ps[e]?udomonas|any_mssa|any_mrsa|any_burkholderia)_?(2018|2019)$", hdr, value = TRUE)))
say(""); say("=== Exacerbation-related columns by year ===")
say(paste(grep("ivi|exacerb", hdr, value = TRUE), collapse = ", "))
