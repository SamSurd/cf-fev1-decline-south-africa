#!/usr/bin/env Rscript

# ==============================================================================
# DATA INTEGRITY & PRE-PROCESSING
# ==============================================================================
# Purpose: 
# 1. Read Master Excel file.
# 2. If 'CF_longitudinal' sheet exists, use it. If not, convert 'Sheet1' (wide) to long.
# 3. Perform data cleaning, column mapping, and feature engineering.
# 4. Generate audit logs and save final parquet file.
# ==============================================================================

set.seed(20251111)

# 1. SETUP & LIBRARIES ---------------------------------------------------------
required_pkgs <- c(
  "tidyverse", "readxl", "lubridate", "janitor",
  "arrow", "digest", "readr", "stringr", "purrr"
)

to_install <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Define Directories
data_path   <- "data_raw/CF_registry_master.xlsx"
out_data    <- "outputs/data_derived"
out_logs    <- "outputs/logs"

dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_logs, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. HELPER FUNCTIONS (GENERAL)
# ==============================================================================

standard_ethnicity_levels <- c(
  "Black African","Mixed Race","Indian","Asian","Caucasian","Other","Unknown"
)

normalize_ethnicity <- function(values) {
  vals <- stringr::str_trim(tolower(as.character(values)))
  vals[is.na(vals)] <- ""
  
  is_white   <- stringr::str_detect(vals, "white|caucas|euro")
  is_black   <- stringr::str_detect(vals, "black|afr")
  is_mixed   <- stringr::str_detect(vals, "mix|colou|color|biracial")
  is_indian  <- stringr::str_detect(vals, "(^|\\b)indian(\\b|$)")
  is_asian   <- stringr::str_detect(vals, "(^|\\b)(asian|pakistan|bangla|chinese|orient|far east)(\\b|$)")
  
  category <- dplyr::case_when(
    vals == "" | vals %in% c("na","n/a","unknown","not known","not recorded","missing") ~ "Unknown",
    is_black ~ "Black African",
    is_mixed ~ "Mixed Race",
    is_white ~ "Caucasian",
    is_indian ~ "Indian",
    is_asian ~ "Asian",
    TRUE ~ "Other"
  )
  factor(category, levels = standard_ethnicity_levels)
}

derive_f508_group <- function(genotype_text, delta_text) {
  clean_vec <- function(x) {
    x <- stringr::str_trim(tolower(as.character(x)))
    x[x %in% c("na","n/a","unknown","not known","")] <- NA_character_
    x
  }
  geno <- clean_vec(genotype_text)
  delta <- clean_vec(delta_text)
  source_vec <- ifelse(!is.na(geno) & geno != "", geno, delta)
  
  classify <- function(x) {
    if (is.na(x) || x == "") return("Other")
    if (stringr::str_detect(x, "homo") ||
        stringr::str_detect(x, "f508[A-z]*[/ ]f508") ||
        stringr::str_detect(x, "2\\s*c(opy|opies)") ||
        stringr::str_detect(x, "two\\s*c(opy|opies)")) {
      return("Homozygous")
    }
    if (stringr::str_detect(x, "hetero") ||
        stringr::str_detect(x, "1\\s*c(opy|opies)") ||
        stringr::str_detect(x, "one\\s*c(opy|opies)") ||
        stringr::str_detect(x, "f508")) {
      return("Heterozygous")
    }
    "Other"
  }
  res <- vapply(source_vec, classify, character(1), USE.NAMES = FALSE)
  factor(res, levels = c("Other","Heterozygous","Homozygous"))
}

center_scale <- function(x){
  mu <- mean(x, na.rm = TRUE); sdv <- stats::sd(x, na.rm = TRUE)
  if (is.finite(sdv) && sdv > 0) (x - mu)/sdv else x - mu
}

# ==============================================================================
# 3. HELPER FUNCTIONS (PRE-PROCESSING WIDE TO LONG)
# ==============================================================================

coerce_date <- function(x) suppressWarnings(as.Date(x))
coerce_num  <- function(x) suppressWarnings(as.numeric(x))
coerce_chr  <- function(x) as.character(x)

first_present <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% names(df)) return(nm)
  }
  NA_character_
}

status_from_any_chronic <- function(any_val, chronic_val) {
  n <- max(length(any_val %||% NA), length(chronic_val %||% NA), na.rm = TRUE)
  any_vec <- rep(NA, n)
  chr_vec <- rep(NA, n)
  
  if (!is.null(chronic_val)) {
    if (is.numeric(chronic_val)) {
      chr_vec <- ifelse(chronic_val >= 1, TRUE, ifelse(chronic_val == 0, FALSE, NA))
    } else {
      cv <- tolower(as.character(chronic_val))
      chr_vec <- ifelse(cv %in% c("1","yes","y","true","chronic"), TRUE,
                        ifelse(cv %in% c("0","no","n","false"), FALSE, NA))
    }
  }
  if (!is.null(any_val)) {
    if (is.numeric(any_val)) {
      any_vec <- ifelse(any_val >= 1, TRUE, ifelse(any_val == 0, FALSE, NA))
    } else {
      av <- tolower(as.character(any_val))
      any_vec <- ifelse(av %in% c("1","yes","y","true","positive","any"), TRUE,
                        ifelse(av %in% c("0","no","n","false","none","negative"), FALSE, NA))
    }
  }
  out <- rep(NA_character_, n)
  out[!is.na(chr_vec) & chr_vec == TRUE] <- "Chronic"
  out[is.na(out) & !is.na(any_vec) & any_vec == TRUE] <- "Intermittent"
  out[is.na(out) & !is.na(any_vec) & any_vec == FALSE] <- "Negative"
  out
}

derive_f508_status_wide <- function(geno, mut1, mut2, other1, other2) {
  vals <- tolower(c(geno, mut1, mut2, other1, other2))
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return(NA_character_)
  has_f508 <- function(s) grepl("f508|\u0394508|deltaf508|\u0394f508|p\\.phe508del|phe508del", s)
  flags <- vapply(vals, has_f508, logical(1))
  if (sum(flags) >= 2) return("Homozygous")
  if (sum(flags) == 1) return("Heterozygous")
  "None"
}

make_cf_longitudinal_from_sheet1 <- function(path) {
  message("Preprocessing: Reading Sheet1 and converting to longitudinal format...")
  df <- readxl::read_excel(path, sheet = "Sheet1") %>% janitor::clean_names()
  years <- 2018:2024
  
  # Patient ID
  pat_cols <- grep("^patcode(_[0-9]{4})?$", names(df), value = TRUE)
  if (!length(pat_cols) && !("patcode" %in% names(df))) stop("No patient ID column found.")
  
  id_vec <- rep(NA_character_, nrow(df))
  for (nm in pat_cols) id_vec <- dplyr::coalesce(id_vec, coerce_chr(df[[nm]]))
  if ("patcode" %in% names(df)) id_vec <- dplyr::coalesce(id_vec, coerce_chr(df$patcode))
  
  # Demographics (Static)
  age_dx_nm <- first_present(df, c("agediagnosis","age_at_diagnosis","age_diag","age_at_dx"))
  age_dx <- if (!is.na(age_dx_nm)) coerce_num(df[[age_dx_nm]]) else NA_real_
  
  sex_nm <- first_present(df, c("sex","gender"))
  sex_vec <- if (!is.na(sex_nm)) coerce_chr(df[[sex_nm]]) else NA_character_
  sex_vec <- case_when(
    tolower(sex_vec) %in% c("m","male") ~ "Male",
    tolower(sex_vec) %in% c("f","female") ~ "Female",
    TRUE ~ NA_character_
  )
  
  eth_nm <- first_present(df, c("ethnicity","race","ethnic_group","ethnicity_raw"))
  ethnicity_vec <- if (!is.na(eth_nm)) coerce_chr(df[[eth_nm]]) else NA_character_
  
  # Genotype Logic
  safe_get <- function(nm) { if (!is.na(nm)) df[[nm]] else rep(NA, nrow(df)) }
  geno_raw <- coerce_chr(safe_get(first_present(df, c("genotype"))))
  
  f508_status <- pmap_chr(list(
    coerce_chr(safe_get(first_present(df, c("f508_status","genotype")))),
    coerce_chr(safe_get(first_present(df, c("mutation1","mut1","allele1","f508_1")))),
    coerce_chr(safe_get(first_present(df, c("mutation2","mut2","allele2","f508_2")))),
    coerce_chr(safe_get(first_present(df, c("othermutation1","other_mutation1")))),
    coerce_chr(safe_get(first_present(df, c("othermutation2","other_mutation2"))))
  ), ~ derive_f508_status_wide(..1, ..2, ..3, ..4, ..5))
  
  # CFTR Start Date (prefer 2024 commencement column when available)
  if ("date_cftrm_commenced_2024" %in% names(df)) {
    cftr_start_date <- coerce_date(df[["date_cftrm_commenced_2024"]])
  } else {
    cftr_date_cols <- c(
      grep("^date_cftrm_commenced_[0-9]{4}$", names(df), value = TRUE),
      grep("^date_cftrm_commenced$", names(df), value = TRUE),
      grep("^treatment_start_date_cftr(_new)?$", names(df), value = TRUE),
      grep("^cftr_start_date$", names(df), value = TRUE)
    )
    cftr_start_per_col <- map(cftr_date_cols, ~ coerce_date(df[[.x]]))
    cftr_start_date <- reduce(cftr_start_per_col, ~ { dplyr::coalesce(.x, .y) }, .init = as.Date(NA))
    
    if (length(cftr_start_per_col)) {
      all_dates_mat <- do.call(cbind, lapply(cftr_start_per_col, as.numeric))
      earliest_num <- suppressWarnings(apply(all_dates_mat, 1, min, na.rm = TRUE))
      earliest_num[is.infinite(earliest_num)] <- NA_real_
      cftr_start_date <- as.Date(earliest_num, origin = "1970-01-01")
    }
  }
  
  # Per-Year Slice Generator
  make_year_slice <- function(y) {
    dob_vec <- suppressWarnings(as.Date(df$dob))
    
    # Dates & Age
    vd_nm <- first_present(df, c(paste0("date_fev1_", y), paste0("date_fev1", y), paste0("date_fev_", y)))
    visit_date <- if (!is.na(vd_nm)) coerce_date(df[[vd_nm]]) else as.Date(NA)
    
    age_nm <- first_present(df, c(paste0("age_in_years_", y), paste0("age_years_at_lft_", y), paste0("age_", y), paste0("age_years_", y)))
    age_at_visit <- if (!is.na(age_nm)) coerce_num(df[[age_nm]]) else NA_real_
    
    # Age fallback
    need_age <- is.na(age_at_visit) & !is.na(dob_vec) & !is.na(visit_date)
    if (any(need_age)) {
      age_at_visit[need_age] <- suppressWarnings(as.numeric(difftime(visit_date[need_age], dob_vec[need_age], units = "days")) / 365.25)
    }
    
    # Sector
    sec_nm <- first_present(df, c(paste0("health_sector_", y), paste0("healthsector_", y), paste0("sector_", y), paste0("public_private_", y)))
    sector <- if (!is.na(sec_nm)) coerce_chr(df[[sec_nm]]) else NA_character_
    
    # FEV1
    get_or_na <- function(nm){ if (!is.na(nm)) coerce_num(df[[nm]]) else rep(NA_real_, nrow(df)) }
    gli_names <- c(paste0("gli_fev1_pp_", y), paste0("gli_percpredfev", y))
    loc_names <- c(paste0("fev1_percent_pred_local_ref_", y), paste0("fev1_pp_", y), paste0("fev1_local_", y), paste0("x", y, "_fev1_percent_pred_local_ref"))  # revision: the 2023 registry column is named x2023_fev1_percent_pred_local_ref
    
    gli_vec <- rep(NA_real_, nrow(df))
    for (nm in gli_names) if (nm %in% names(df)) gli_vec <- dplyr::coalesce(gli_vec, coerce_num(df[[nm]]))
    loc_vec <- rep(NA_real_, nrow(df))
    for (nm in loc_names) if (nm %in% names(df)) loc_vec <- dplyr::coalesce(loc_vec, coerce_num(df[[nm]]))
    fev1_pp <- dplyr::coalesce(gli_vec, loc_vec)
    
    # BMI
    bmi_nm  <- first_present(df, c(paste0("bmi_", y), paste0("bmi", y)))
    bmi     <- if (!is.na(bmi_nm)) coerce_num(df[[bmi_nm]]) else NA_real_
    bmiz_nm <- first_present(df, c(paste0("bmiz_", y), paste0("bmiz", y)))
    bmiz    <- if (!is.na(bmiz_nm)) coerce_num(df[[bmiz_nm]]) else NA_real_
    
    # Pancreatic & CFRD
    pan_nm <- first_present(df, c(paste0("pancreatic_status_", y), paste0("pancreatic_", y)))
    pancreatic <- if (!is.na(pan_nm)) coerce_chr(df[[pan_nm]]) else NA_character_
    pancreatic <- case_when(
      tolower(pancreatic) %in% c("insufficient","yes","1","true","y") ~ "Insufficient",
      tolower(pancreatic) %in% c("sufficient","no","0","false","n") ~ "Sufficient",
      TRUE ~ pancreatic
    )
    
    cfrd_nm <- first_present(df, c(paste0("cfrd_", y), paste0("diabetes_", y), paste0("diabetes", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    cfrd <- if (!is.na(cfrd_nm)) coerce_chr(df[[cfrd_nm]]) else NA_character_
    cfrd <- case_when(
      tolower(cfrd) %in% c("1","yes","true","y") ~ "Yes",
      tolower(cfrd) %in% c("0","no","false","n") ~ "No",
      TRUE ~ NA_character_
    )
    
    # Infections
    pa_any_nm <- first_present(df, c(paste0("pseudo_any_", y), paste0("pseudo_any", y), paste0("ever_pseudomonas_", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    pa_chr_nm <- first_present(df, c(paste0("chronicpseudomaer_", y), paste0("chronicpseudomaer", y), paste0("chronic_pseudomonas_", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    pa_status <- status_from_any_chronic(if (!is.na(pa_any_nm)) df[[pa_any_nm]] else NULL, if (!is.na(pa_chr_nm)) df[[pa_chr_nm]] else NULL)
    
    sa_any_nm <- first_present(df, c(paste0("any_mssa_", y), paste0("any_mssa", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    sa_chr_nm <- first_present(df, c(paste0("chronicstaphylaur_", y), paste0("chronicstaphylaur", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    sa_status <- status_from_any_chronic(if (!is.na(sa_any_nm)) df[[sa_any_nm]] else NULL, if (!is.na(sa_chr_nm)) df[[sa_chr_nm]] else NULL)
    
    bc_any_nm <- first_present(df, c(paste0("any_burkholderia_", y), paste0("any_burkholderia", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    bc_chr_nm <- first_present(df, c(paste0("chronic_burkholderia_", y), paste0("chronic_burkholderia", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    bc_status <- status_from_any_chronic(if (!is.na(bc_any_nm)) df[[bc_any_nm]] else NULL, if (!is.na(bc_chr_nm)) df[[bc_chr_nm]] else NULL)
    
    mr_any_nm <- first_present(df, c(paste0("any_mrsa_", y), paste0("any_mrsa", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    mr_chr_nm <- first_present(df, c(paste0("chronic_mrsa_", y), paste0("chronic_mrsa", y)))  # revision (5 Nov 2026): 2018-2019 registry columns have no underscore before the year
    mrsa_status <- status_from_any_chronic(if (!is.na(mr_any_nm)) df[[mr_any_nm]] else NULL, if (!is.na(mr_chr_nm)) df[[mr_chr_nm]] else NULL)
    
    # Exacerbations
    exac_nm <- first_present(df, c(paste0("exacerbations_per_year_", y), paste0("ivi_number", y), paste0("ivi_number_", y)))
    exacerbations_per_year <- if (!is.na(exac_nm)) coerce_num(df[[exac_nm]]) else NA_real_
    
    # CFTR Modulator (Visit Specific)
    cftr_any_nm <- first_present(df, c(paste0("cftrm_drug_", y), paste0("cftr_modulator_", y), "cftr_modulator_cftr"))
    cftr_raw <- if (!is.na(cftr_any_nm)) df[[cftr_any_nm]] else NA
    cftr_modulator <- case_when(
      tolower(as.character(cftr_raw)) %in% c("1","yes","true","on","y") ~ "Yes",
      !is.na(cftr_raw) & nzchar(as.character(cftr_raw)) ~ "Yes",
      TRUE ~ "No"
    )
    
    tibble(
      patcode = id_vec,
      year = y,
      visit_date = visit_date,
      age_at_visit = age_at_visit,
      age_at_diagnosis = age_dx,
      sex = sex_vec,
      ethnicity = ethnicity_vec,
      sector = sector,
      genotype = geno_raw,
      fev1pp = fev1_pp,
      bmi = bmi,
      bmiz = bmiz,
      pancreatic = pancreatic,
      cfrd = cfrd,
      f508_status = f508_status,
      exacerbations_per_year = exacerbations_per_year,
      paeruginosa_status = pa_status,
      staph_status = sa_status,
      burkholderia_status = bc_status,
      mrsa_status = mrsa_status,
      cftr_modulator = cftr_modulator,
      cftr_start_date = cftr_start_date
    )
  }
  
  slices <- map(years, make_year_slice)
  long <- bind_rows(slices) %>% 
    filter(!is.na(visit_date) | !is.na(fev1pp) | !is.na(bmi) | !is.na(bmiz))
  
  return(long)
}

# ==============================================================================
# 4. MAIN EXECUTION
# ==============================================================================

# 4a. Manifest for raw files ---------------------------------------------------
raw_files <- list.files("data_raw", full.names = TRUE)
manifest <- purrr::map_dfr(raw_files, function(f) {
  info <- file.info(f)
  tibble(
    file = f,
    size_bytes = as.numeric(info$size),
    mtime = as.character(info$mtime),
    md5 = digest::digest(f, file = TRUE, algo = "md5")
  )
})
readr::write_csv(manifest, file.path(out_logs, "data_manifest.csv"))

# 4b. Read Excel ---------------------------------------------------------------
if (!file.exists(data_path)) stop("Expected raw file not found: ", data_path)

# Attempt to read Sheet1 for specific CFTR audits / Patient Level data
raw_sheet1 <- tryCatch(
  readxl::read_excel(data_path, sheet = "Sheet1") %>% janitor::clean_names(),
  error = function(e) NULL
)

# Try reading Long format, otherwise Preprocess
available_sheets <- try(readxl::excel_sheets(data_path), silent = TRUE)

if (!inherits(available_sheets, "try-error") && ("CF_longitudinal" %in% available_sheets)) {
  message("Reading existing 'CF_longitudinal' sheet.")
  raw <- readxl::read_excel(data_path, sheet = "CF_longitudinal") %>%
    janitor::clean_names()
} else {
  message("Sheet 'CF_longitudinal' not found. Using preprocessing logic.")
  raw <- make_cf_longitudinal_from_sheet1(data_path) %>% janitor::clean_names()
}

# 4c. Derive CFTR exposure flags from Sheet1 (Audit Layer) ---------------------
cftr_lookup <- NULL
if (!is.null(raw_sheet1)) {
  id_cols_sheet1 <- grep("^patcode", names(raw_sheet1), value = TRUE)
  if (length(id_cols_sheet1)) {
    id_df <- raw_sheet1[id_cols_sheet1]
    id_df[] <- lapply(id_df, function(x) as.character(x))
    pat_raw <- dplyr::coalesce(!!!id_df)
    
    cftr_cols <- c("cftr_modulator_cftr", "which_treatment_cftr", "cftrm2018", "cftrm_2019",
                   "cftrm_2020", "cftrm_2021...454", "cftrm_2021...455", "cftrm_2022",
                   "cftrm_2023", "cftrm_drug_2024", "treatment_start_date_cftr", "treatment_start_date_cftr_new")
    cftr_cols <- cftr_cols[cftr_cols %in% names(raw_sheet1)]
    
    if(length(cftr_cols) > 0) {
      cftr_mat <- raw_sheet1[, cftr_cols, drop = FALSE]
      cftr_chr <- as.data.frame(lapply(cftr_mat, as.character), stringsAsFactors = FALSE)
      is_nonempty <- !is.na(cftr_chr) & cftr_chr != ""
      ever_cftr_flag <- apply(is_nonempty, 1, any)

      # CFTR start date: prefer the 2024 commencement column when available
      if ("date_cftrm_commenced_2024" %in% names(raw_sheet1)) {
        cftr_start_date_raw <- suppressWarnings(as.Date(raw_sheet1[["date_cftrm_commenced_2024"]]))
      } else {
        date_cols <- grep("^(date_cftrm|date_cftrm_commenced_|treatment_start_date_cftr)", names(raw_sheet1), value = TRUE)
        date_mat <- raw_sheet1[, date_cols, drop = FALSE]
        date_list <- lapply(date_mat, function(x) suppressWarnings(as.Date(x)))
        
        if (length(date_list)) {
          date_num <- do.call(cbind, lapply(date_list, as.numeric))
          earliest_num <- suppressWarnings(apply(date_num, 1, min, na.rm = TRUE))
          earliest_num[!is.finite(earliest_num)] <- NA_real_
          cftr_start_date_raw <- as.Date(earliest_num, origin = "1970-01-01")
        } else {
          cftr_start_date_raw <- rep(as.Date(NA), nrow(raw_sheet1))
        }
      }
      
      cftr_lookup <- tibble(
        patient_id = pat_raw,
        ever_cftr_raw = as.integer(ever_cftr_flag),
        cftr_start_date_raw = cftr_start_date_raw
      ) %>%
        filter(!is.na(patient_id)) %>%
        group_by(patient_id) %>%
        summarise(
          ever_cftr_raw = max(ever_cftr_raw, na.rm = TRUE),
          cftr_start_date_raw = suppressWarnings(min(cftr_start_date_raw, na.rm = TRUE)),
          .groups = "drop"
        )
      cftr_lookup$cftr_start_date_raw[!is.finite(as.numeric(cftr_lookup$cftr_start_date_raw))] <- NA
    }
  }
}

# 4d. Check required columns and Rename ----------------------------------------
# Note: Preprocessing outputs 'fev1pp', but cleaning names might strip punct
# We flexibly check required columns based on what we expect `janitor::clean_names` to produce.

dat <- raw %>%
  rename(
    patient_id                = patcode,
    visit_year                = year,
    visit_date                = visit_date,
    age_years                 = age_at_visit,
    age_diag_years            = age_at_diagnosis,
    sex                       = sex,
    sector                    = sector,
    fev1_pp                   = fev1pp,
    bmi                       = bmi,
    bmiz                      = bmiz,
    pancreatic_status         = pancreatic,
    cfrd                      = cfrd,
    deltaf508_genotype        = f508_status,
    exacerbation_count        = exacerbations_per_year,
    pa_status                 = paeruginosa_status,
    sa_status                 = staph_status,
    bcepacia_status           = burkholderia_status,
    cftr_modulator_use        = cftr_modulator,
    cftr_modulator_start_date = cftr_start_date
  ) %>%
  mutate(
    patient_id = as.character(patient_id),
    visit_date = as.Date(visit_date),
    cftr_modulator_start_date = as.Date(cftr_modulator_start_date),
    deltaf508_genotype = as.character(deltaf508_genotype)
  )

# 4e. Feature Engineering ------------------------------------------------------

# Ethnicity
ethnicity_source <- rep(NA_character_, nrow(dat))
if ("ethnicity" %in% names(raw)) ethnicity_source <- as.character(raw$ethnicity)
dat$ethnicity <- ethnicity_source
dat$ethnicity_cat <- normalize_ethnicity(ethnicity_source)

# Genotype Grouping
if ("genotype" %in% names(raw)) dat$genotype <- as.character(raw$genotype) else dat$genotype <- NA_character_
dat$f508_genotype_group <- derive_f508_group(dat$genotype, dat$deltaf508_genotype)

# Lung Transplant Logic (Exclude Post-Transplant)
if (!is.null(raw_sheet1)) {
  transplant_sheet <- raw_sheet1 # Reuse pre-read sheet
  id_cols <- grep("^patcode", names(transplant_sheet), value = TRUE)
  lung_cols <- grep("^date_lung_tx", names(transplant_sheet), value = TRUE)
  
  if (length(id_cols) > 0 && length(lung_cols) > 0) {
    id_df <- transplant_sheet[id_cols]
    id_df[] <- lapply(id_df, function(x) as.character(x))
    patient_ids <- reduce(id_df, dplyr::coalesce)
    lung_list <- lapply(lung_cols, function(nm) as.Date(transplant_sheet[[nm]]))
    lung_vec <- reduce(lung_list, dplyr::coalesce)
    
    transplant_lookup <- tibble(
      patient_id = patient_ids,
      lung_transplant_date = lung_vec
    ) %>%
      filter(!is.na(patient_id) & !is.na(lung_transplant_date)) %>%
      group_by(patient_id) %>%
      summarise(lung_transplant_date = min(lung_transplant_date, na.rm = TRUE), .groups = "drop")
    
    dat <- dat %>%
      left_join(transplant_lookup, by = "patient_id") %>%
      filter(is.na(lung_transplant_date) | visit_date < lung_transplant_date)
  }
}

# Preserve Copy
dat_all_visits <- dat

# revision (5 Nov 2026): CFRD and pancreatic insufficiency are carried forward from the
# first registry year recorded "Yes"/"Insufficient" (any registry year), as stated in the Methods
# ("coded as permanent binary variables once documented"). Previously each year's own value was used.
cfrd_first <- dat_all_visits %>% filter(cfrd %in% "Yes", !is.na(visit_year)) %>%
  group_by(patient_id) %>% summarise(cfrd_first_year = min(visit_year), .groups = "drop")
pi_first <- dat_all_visits %>% filter(pancreatic_status %in% "Insufficient", !is.na(visit_year)) %>%
  group_by(patient_id) %>% summarise(pi_first_year = min(visit_year), .groups = "drop")
dat <- dat %>%
  left_join(cfrd_first, by = "patient_id") %>% left_join(pi_first, by = "patient_id") %>%
  mutate(
    cfrd = dplyr::if_else(!is.na(cfrd_first_year) & !is.na(visit_year) & visit_year >= cfrd_first_year, "Yes", cfrd),
    pancreatic_status = dplyr::if_else(!is.na(pi_first_year) & !is.na(visit_year) & visit_year >= pi_first_year, "Insufficient", pancreatic_status)
  ) %>%
  select(-cfrd_first_year, -pi_first_year)

# 4f. Inclusion Criteria -------------------------------------------------------
dat <- dat %>% filter(!is.na(age_years) & age_years >= 6)

keep_ids <- dat %>%
  filter(!is.na(fev1_pp) & !is.na(visit_year)) %>%
  group_by(patient_id) %>%
  summarise(n_years = dplyr::n_distinct(visit_year), .groups = "drop") %>%
  filter(n_years >= 2) %>% pull(patient_id)

dat <- dat %>% filter(patient_id %in% keep_ids)

# 4g. Audit Logs ---------------------------------------------------------------
audit_before <- dat %>%
  group_by(patient_id, visit_year) %>%
  summarise(n_visits = n(), any_fev1 = any(!is.na(fev1_pp)), .groups = "drop")
readr::write_csv(audit_before, file.path(out_logs, "extraction_patient_year_before_collapse.csv"))

drop_summary <- tibble(
  step = c("raw_patients", "patients_with_any_fev1", "patients_with_age_ge_6", "patients_with_ge2_years_fev1"),
  n = c(
    dplyr::n_distinct(dat_all_visits$patient_id),
    dplyr::n_distinct(dat_all_visits$patient_id[!is.na(dat_all_visits$fev1_pp)]),
    dplyr::n_distinct(dat_all_visits$patient_id[!is.na(dat_all_visits$age_years) & dat_all_visits$age_years >= 6]),
    dplyr::n_distinct(keep_ids)
  )
)
readr::write_csv(drop_summary, file.path(out_logs, "dropout_summary.csv"))

# 4h. Variable Derivation (Time, Nutrition, Exacerbation, Baseline) ------------
dat <- dat %>%
  group_by(patient_id) %>%
  arrange(visit_date, .by_group = TRUE) %>%
  mutate(
    baseline_date = min(visit_date, na.rm = TRUE),
    time_since_baseline = as.numeric(difftime(visit_date, baseline_date, units = "days")) / 365.25,
    calendar_year = dplyr::if_else(!is.na(visit_year), visit_year, lubridate::year(visit_date))
  ) %>%
  ungroup() %>%
  mutate(calendar_year = factor(calendar_year))

if ("2018" %in% levels(dat$calendar_year)) dat$calendar_year <- stats::relevel(dat$calendar_year, ref = "2018")

dat <- dat %>%
  mutate(
    nutrition_cat = case_when(
      age_years < 18 & !is.na(bmiz) & bmiz < -2          ~ "Underweight",
      age_years < 18 & !is.na(bmiz) & bmiz <= 1          ~ "Normal",
      age_years < 18 & !is.na(bmiz) & bmiz > 1           ~ "Overweight",
      age_years >= 18 & !is.na(bmi) & bmi < 18.5         ~ "Underweight",
      age_years >= 18 & !is.na(bmi) & bmi < 25           ~ "Normal",
      age_years >= 18 & !is.na(bmi) & bmi >= 25          ~ "Overweight",
      TRUE                                               ~ "Unknown"
    ),
    nutrition_cat = factor(nutrition_cat, levels = c("Normal","Underweight","Overweight","Unknown"))
  )

dat <- dat %>%
  mutate(
    exac_count2 = dplyr::coalesce(exacerbation_count, 0),
    exac_cat = case_when(
      exac_count2 == 0 ~ "0",
      exac_count2 == 1 ~ "1",
      exac_count2 >= 2 ~ "2+",
      TRUE             ~ "0"
    ),
    exac_cat = factor(exac_cat, levels = c("0","1","2+","Unknown"))
  ) %>% select(-exac_count2)

# Baseline Extraction
first_visit_df <- dat %>%
  arrange(patient_id, visit_date) %>%
  group_by(patient_id) %>%
  slice_head(n = 1) %>%
  ungroup()

bl_fev <- dat %>%
  arrange(patient_id, visit_date) %>%
  group_by(patient_id) %>%
  filter(!is.na(fev1_pp)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(patient_id, baseline_fev1_pp = fev1_pp)

baseline_df <- first_visit_df %>%
  transmute(
    patient_id,
    baseline_age_years = age_years,
    baseline_age_group = case_when(
      !is.na(age_years) & age_years >=  6 & age_years < 11 ~ "6-10",
      !is.na(age_years) & age_years >= 11 & age_years < 18 ~ "11-17",
      !is.na(age_years) & age_years >= 18 & age_years < 26 ~ "18-25",
      !is.na(age_years) & age_years >= 26 & age_years < 40 ~ "26-39",
      !is.na(age_years) & age_years >= 40                  ~ "40+",
      TRUE                                                  ~ "Unknown"
    ),
    baseline_age_group = factor(baseline_age_group, levels = c("6-10","11-17","18-25","26-39","40+","Unknown")),
    baseline_nutrition_cat        = nutrition_cat,
    baseline_exac_cat             = exac_cat,
    baseline_sector               = sector,
    baseline_pancreatic_insuff    = dplyr::if_else(pancreatic_status == "Insufficient", 1L, 0L, missing = 0L),
    baseline_cfrd_indicator       = dplyr::if_else(cfrd == "Yes", 1L, 0L, missing = 0L),
    baseline_cftr_mod_indicator   = dplyr::if_else(cftr_modulator_use == "Yes", 1L, 0L, missing = 0L),
    baseline_pa_cat               = factor(dplyr::coalesce(pa_status, "Negative"), levels = c("Negative","Intermittent","Chronic")),
    baseline_pa_any               = dplyr::if_else(is.na(pa_status), 0L, if_else(pa_status %in% c("Intermittent","Chronic"), 1L, 0L))
  ) %>%
  left_join(bl_fev, by = "patient_id")

dat <- dat %>% left_join(baseline_df, by = "patient_id")

# Center Baseline FEV1 and Baseline Age
mu_blf <- mean(dat$baseline_fev1_pp, na.rm = TRUE)
mu_age <- mean(dat$baseline_age_years, na.rm = TRUE)
dat <- dat %>% mutate(
  baseline_fev1_pp_c = baseline_fev1_pp - mu_blf,
  baseline_age_years_c = baseline_age_years - mu_age,
  baseline_fev1_cat = case_when(
    is.na(baseline_fev1_pp)          ~ "Unknown",
    baseline_fev1_pp < 40            ~ "Severe",
    baseline_fev1_pp >= 40 & baseline_fev1_pp <= 70 ~ "Moderate",
    baseline_fev1_pp > 70            ~ "Normal",
    TRUE                             ~ "Unknown"
  ),
  baseline_fev1_cat = factor(baseline_fev1_cat, levels = c("Normal","Moderate","Severe","Unknown"))
)

# 4i. Chronic Infections Indicators --------------------------------------------
pa_start <- dat_all_visits %>%
  filter(pa_status == "Chronic") %>%
  group_by(patient_id) %>%
  summarise(pa_chronic_start_date = min(visit_date, na.rm = TRUE), ever_chronic_pa = 1L, .groups = "drop")

dat <- dat %>%
  left_join(pa_start, by = "patient_id") %>%
  mutate(
    ever_chronic_pa = dplyr::if_else(is.na(ever_chronic_pa), 0L, ever_chronic_pa)
  ) %>%
  group_by(patient_id) %>%
  arrange(visit_date, .by_group = TRUE) %>%
  mutate(
    chronic_paeruginosa = dplyr::if_else(pa_status == "Chronic", 1L, 0L, missing = 0L),
    duration_since_chronic_pa = case_when(
      is.na(pa_chronic_start_date) ~ 0,
      visit_date < pa_chronic_start_date ~ 0,
      TRUE ~ as.numeric(difftime(visit_date, pa_chronic_start_date, units = "days")) / 365.25
    )
  ) %>% ungroup()

# B. Cepacia, MSSA, and MRSA Chronic Infection Indicators
bcep_start <- dat_all_visits %>%
  filter(bcepacia_status == "Chronic") %>%
  group_by(patient_id) %>%
  summarise(bcep_chronic_start_date = min(visit_date, na.rm = TRUE), ever_chronic_bcepacia = 1L, .groups = "drop")

mssa_start <- dat_all_visits %>%
  filter(sa_status == "Chronic") %>%
  group_by(patient_id) %>%
  summarise(mssa_chronic_start_date = min(visit_date, na.rm = TRUE), ever_chronic_mssa = 1L, .groups = "drop")

dat <- dat %>%
  left_join(bcep_start, by = "patient_id") %>%
  left_join(mssa_start, by = "patient_id") %>%
  mutate(
    ever_chronic_bcepacia = dplyr::if_else(is.na(ever_chronic_bcepacia), 0L, ever_chronic_bcepacia),
    chronic_bcepacia = dplyr::if_else(bcepacia_status == "Chronic", 1L, 0L, missing = 0L),
    duration_since_chronic_bcep = case_when(
      is.na(bcep_chronic_start_date) ~ 0,
      visit_date < bcep_chronic_start_date ~ 0,
      TRUE ~ as.numeric(difftime(visit_date, bcep_chronic_start_date, units = "days")) / 365.25
    ),
    ever_chronic_mssa = dplyr::if_else(is.na(ever_chronic_mssa), 0L, ever_chronic_mssa),
    chronic_mssa = dplyr::if_else(sa_status == "Chronic", 1L, 0L, missing = 0L),
    duration_since_chronic_mssa = case_when(
      is.na(mssa_chronic_start_date) ~ 0,
      visit_date < mssa_chronic_start_date ~ 0,
      TRUE ~ as.numeric(difftime(visit_date, mssa_chronic_start_date, units = "days")) / 365.25
    )
  )

# MRSA check (variable might not exist in all preprocess variations)
if ("mrsa_status" %in% names(dat_all_visits)) {
  mrsa_start <- dat_all_visits %>%
    filter(mrsa_status == "Chronic") %>%
    group_by(patient_id) %>%
    summarise(mrsa_chronic_start_date = min(visit_date, na.rm = TRUE), ever_chronic_mrsa = 1L, .groups = "drop")
  
  dat <- dat %>%
    left_join(mrsa_start, by = "patient_id") %>%
    mutate(
      ever_chronic_mrsa = dplyr::if_else(is.na(ever_chronic_mrsa), 0L, ever_chronic_mrsa),
      chronic_mrsa = dplyr::if_else(mrsa_status == "Chronic", 1L, 0L, missing = 0L),
      duration_since_chronic_mrsa = case_when(
        is.na(mrsa_chronic_start_date) ~ 0,
        visit_date < mrsa_chronic_start_date ~ 0,
        TRUE ~ as.numeric(difftime(visit_date, mrsa_chronic_start_date, units = "days")) / 365.25
      )
    )
} else {
  dat <- dat %>% mutate(ever_chronic_mrsa = 0L, chronic_mrsa = 0L, duration_since_chronic_mrsa = 0)
}

# 4j. CFTR Modulator Indicators ------------------------------------------------
if (!is.null(cftr_lookup)) {
  dat <- dat %>%
    left_join(cftr_lookup, by = "patient_id") %>%
    mutate(
      ever_cftr_yes = dplyr::coalesce(ever_cftr_raw, 0L),
      cftr_modulator_start_date = cftr_start_date_raw
    ) %>% select(-cftr_start_date_raw)
} else {
  dat <- dat %>% mutate(ever_cftr_raw = 0L, cftr_modulator_start_date = as.Date(NA))
}

dat <- dat %>%
  mutate(
    duration_modulator = dplyr::case_when(
      is.na(cftr_modulator_start_date)       ~ 0,
      visit_date < cftr_modulator_start_date ~ 0,
      TRUE ~ as.numeric(difftime(visit_date, cftr_modulator_start_date, units = "days")) / 365.25
    ),
    cftr_modulator_indicator = dplyr::if_else(
      !is.na(cftr_modulator_start_date) & visit_date >= cftr_modulator_start_date, 1L, 0L
    ),
    post_cftr_indicator = cftr_modulator_indicator
  )

# CFTR Audit Log
cftr_audit <- dat %>%
  group_by(patient_id) %>%
  summarise(
    ever_cftr_yes = max(ever_cftr_yes, na.rm = TRUE),
    cftr_modulator_start_date = dplyr::first(na.omit(cftr_modulator_start_date)),
    n_on_cftr_visits = sum(cftr_modulator_indicator == 1L, na.rm = TRUE),
    n_total_visits = dplyr::n(),
    last_visit_date = max(visit_date, na.rm = TRUE),
    total_exposure_years = suppressWarnings(max(duration_modulator, na.rm = TRUE)),
    .groups = "drop"
  )
readr::write_csv(cftr_audit, file.path(out_logs, "cftr_audit.csv"))

# 4k. Final Clean Up -----------------------------------------------------------
dat <- dat %>%
  mutate(
    pancreatic_insufficiency = dplyr::if_else(pancreatic_status == "Insufficient", 1L, 0L, missing = 0L),
    cfrd_indicator           = dplyr::if_else(cfrd == "Yes", 1L, 0L, missing = 0L),
    # Use time_since_baseline in YEARS (no scaling) so that all time
    # coefficients in downstream models are interpretable as pp/year.
    time_cs = time_since_baseline,
    duration_modulator_cs = center_scale(duration_modulator),
    duration_since_chronic_pa_cs = center_scale(duration_since_chronic_pa),
    duration_since_chronic_bcep_cs = center_scale(duration_since_chronic_bcep),
    duration_since_chronic_mssa_cs = center_scale(duration_since_chronic_mssa),
    duration_since_chronic_mrsa_cs = center_scale(duration_since_chronic_mrsa),
    ethnicity_caucasian = {
      tmp <- ethnicity_cat == "Caucasian"
      tmp[is.na(tmp)] <- FALSE
      as.integer(tmp)
    },
    f508_any = {
      tmp <- f508_genotype_group %in% c("Heterozygous", "Homozygous")
      tmp[is.na(tmp)] <- FALSE
      as.integer(tmp)
    }
  )

# Slope eligibility: require >= 2 post-event FEV1pp measurements for slope terms
slope_counts <- dat %>%
  group_by(patient_id) %>%
  summarise(
    n_post_cftr_fev1      = sum(dplyr::if_else(cftr_modulator_indicator == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    n_post_pa_fev1        = sum(dplyr::if_else(chronic_paeruginosa == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    n_post_bcep_fev1      = sum(dplyr::if_else(chronic_bcepacia == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    n_post_mssa_fev1      = sum(dplyr::if_else(chronic_mssa == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    n_post_mrsa_fev1      = sum(dplyr::if_else(chronic_mrsa == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    n_post_pancreas_fev1  = sum(dplyr::if_else(pancreatic_insufficiency == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    n_post_cfrd_fev1      = sum(dplyr::if_else(cfrd_indicator == 1L & !is.na(fev1_pp), 1L, 0L, missing = 0L)),
    .groups = "drop"
  ) %>%
  mutate(
    slope_cftr_ok     = as.integer(n_post_cftr_fev1     >= 2),
    slope_pa_ok       = as.integer(n_post_pa_fev1       >= 2),
    slope_bcep_ok     = as.integer(n_post_bcep_fev1     >= 2),
    slope_mssa_ok     = as.integer(n_post_mssa_fev1     >= 2),
    slope_mrsa_ok     = as.integer(n_post_mrsa_fev1     >= 2),
    slope_pancreas_ok = as.integer(n_post_pancreas_fev1 >= 2),
    slope_cfrd_ok     = as.integer(n_post_cfrd_fev1     >= 2)
  )

dat <- dat %>%
  left_join(slope_counts, by = "patient_id") %>%
  mutate(
    cftr_modulator_indicator_slope = dplyr::if_else(slope_cftr_ok == 1L, cftr_modulator_indicator, 0L, missing = 0L),
    chronic_paeruginosa_slope      = dplyr::if_else(slope_pa_ok == 1L, chronic_paeruginosa, 0L, missing = 0L),
    chronic_bcepacia_slope         = dplyr::if_else(slope_bcep_ok == 1L, chronic_bcepacia, 0L, missing = 0L),
    chronic_mssa_slope             = dplyr::if_else(slope_mssa_ok == 1L, chronic_mssa, 0L, missing = 0L),
    chronic_mrsa_slope             = dplyr::if_else(slope_mrsa_ok == 1L, chronic_mrsa, 0L, missing = 0L),
    pancreatic_insufficiency_slope = dplyr::if_else(slope_pancreas_ok == 1L, pancreatic_insufficiency, 0L, missing = 0L),
    cfrd_indicator_slope           = dplyr::if_else(slope_cfrd_ok == 1L, cfrd_indicator, 0L, missing = 0L)
  )

# Sector Categorization (Patient Level)
sector_patient <- dat_all_visits %>%
  transmute(patient_id, sector_raw = as.character(sector)) %>%
  mutate(
    sector_clean = dplyr::case_when(
      is.na(sector_raw) ~ NA_character_,
      grepl("\\bpublic\\b|^pub|state|government|\\bgov\\b", tolower(sector_raw))  ~ "Public",
      grepl("\\bprivate\\b|^priv|medical[ -]?aid|insurance|scheme|medaid", tolower(sector_raw)) ~ "Private",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(patient_id) %>%
  summarise(
    has_public  = any(sector_clean == "Public",  na.rm = TRUE),
    has_private = any(sector_clean == "Private", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sector_patient = dplyr::case_when(
      has_public & has_private ~ "Combined",
      has_public              ~ "Public",
      has_private             ~ "Private",
      TRUE                    ~ "Public"
    )
  )

dat <- dat %>%
  left_join(sector_patient, by = "patient_id") %>%
  mutate(sector = factor(sector_patient, levels = c("Public","Private","Combined"))) %>%
  mutate(
    sector_public_only = {
      tmp <- sector_patient == "Public"
      tmp[is.na(tmp)] <- FALSE
      as.integer(tmp)
    }
  ) %>%
  select(-sector_patient) %>%
  filter(!is.na(fev1_pp) & is.finite(fev1_pp))

# 5. OUTPUT --------------------------------------------------------------------
arrow::write_parquet(dat, file.path(out_data, "analysis_long.parquet"))

codebook <- dat %>%
  summarise(across(everything(), ~class(.x)[1])) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "class")

readr::write_csv(codebook, file.path(out_data, "codebook.csv"))

message("Data processing complete: analysis_long.parquet and codebook.csv created.")

