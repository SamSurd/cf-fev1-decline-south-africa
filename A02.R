#!/usr/bin/env Rscript

# MODEL_RUNNER (FINAL)
# Fits longitudinal Linear Mixed Models (M0, M1, M2).
# Logic: Triple-Fallback (Correlated -> Uncorrelated -> Intercept).
# Excludes static 'ever_chronic_pa' in favor of time-varying 'chronic_paeruginosa'.

set.seed(20251111)

# 1. Setup and Dependencies -----------------------------------------------

required_pkgs <- c(
  "tidyverse", "arrow", "lme4", "lmerTest",
  "broom.mixed", "performance", "insight", "tibble", "fs"
)

# Install missing packages if necessary
to_install <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, quiet = TRUE)

invisible(lapply(required_pkgs, library, character.only = TRUE))

center_scale <- function(x){
  mu <- mean(x, na.rm = TRUE); sdv <- stats::sd(x, na.rm = TRUE)
  if (is.finite(sdv) && sdv > 0) (x - mu)/sdv else x - mu
}

# Define directories
data_file   <- "outputs/data_derived/analysis_long.parquet"
log_dir     <- "outputs/logs"
results_dir <- "outputs/models"       # For .rds files (R objects)
tables_dir  <- "outputs/tables"       # For .csv files (Excel readable)

# Create directories safely
fs::dir_create(c(log_dir, results_dir, tables_dir))

# Logging setup
log_file <- file.path(log_dir, "agent02_model_runner.log")
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = " "))
  message(msg)
  cat(msg, file = log_file, append = TRUE, sep = "\n")
}

if (!file.exists(data_file)) {
  log_line("FATAL: analysis_long.parquet not found at: ", data_file)
  stop("Data file missing.")
}

log_line("Starting Agent02. Reading data...")
dat <- arrow::read_parquet(data_file) %>% as_tibble()


# 2. Data Preprocessing ---------------------------------------------------

# Define Required Columns
required_cols <- c(
  "patient_id", "time_cs", "age_years",
  "baseline_fev1_cat", "baseline_pa_any",
  "sex",
  "nutrition_cat", "exac_cat",
  "chronic_paeruginosa", "chronic_bcepacia", "chronic_mssa", "chronic_mrsa",
  "cftr_modulator_indicator", "pancreatic_insufficiency",
  "cfrd_indicator", "sector", "ethnicity_cat", "f508_genotype_group",
  "ethnicity_caucasian", "sector_public_only", "f508_any",
  "cftr_modulator_indicator_slope", "chronic_paeruginosa_slope",
  "chronic_bcepacia_slope", "chronic_mssa_slope", "chronic_mrsa_slope",
  "pancreatic_insufficiency_slope", "cfrd_indicator_slope",
  "slope_cftr_ok"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Refactor columns (Standardizing levels) & Safe Factorization
dat <- dat %>%
  mutate(
    sex          = factor(sex, levels = c("Female", "Male")),
    nutrition_cat = factor(nutrition_cat, levels = c("Normal", "Underweight", "Overweight", "Unknown")),
    exac_cat      = factor(exac_cat, levels = c("0", "1", "2+", "Unknown")),
    sector        = factor(sector, levels = c("Public", "Private", "Combined")),
    ethnicity_cat = factor(ethnicity_cat, levels = c("Black African", "Mixed Race", "Indian", "Asian", "Caucasian", "Other", "Unknown")),
    f508_genotype_group = factor(f508_genotype_group, levels = c("Other", "Heterozygous", "Homozygous")),
    ethnicity_caucasian = dplyr::if_else(ethnicity_cat == "Caucasian", 1L, 0L, missing = 0L),
    sector_public_only = dplyr::if_else(sector == "Public", 1L, 0L, missing = 0L),
    f508_any = dplyr::if_else(f508_genotype_group %in% c("Heterozygous", "Homozygous"), 1L, 0L, missing = 0L)
  )

# Filter Data
n_total <- nrow(dat)
dat_model <- dat %>%
  filter(nutrition_cat != "Unknown" | is.na(nutrition_cat)) %>%
  droplevels() %>%
  group_by(patient_id) %>%
  filter(dplyr::n() >= 2) %>% # Longitudinal requirement
  ungroup()

pancreas_start <- dat_model %>%
  filter(pancreatic_insufficiency == 1L) %>%
  group_by(patient_id) %>%
  summarise(pancreas_start_date = suppressWarnings(min(visit_date, na.rm = TRUE)), .groups = "drop") %>%
  mutate(pancreas_start_date = dplyr::if_else(is.infinite(pancreas_start_date), as.Date(NA), pancreas_start_date))

cfrd_start <- dat_model %>%
  filter(cfrd_indicator == 1L) %>%
  group_by(patient_id) %>%
  summarise(cfrd_start_date = suppressWarnings(min(visit_date, na.rm = TRUE)), .groups = "drop") %>%
  mutate(cfrd_start_date = dplyr::if_else(is.infinite(cfrd_start_date), as.Date(NA), cfrd_start_date))

dat_model <- dat_model %>%
  left_join(pancreas_start, by = "patient_id") %>%
  left_join(cfrd_start, by = "patient_id")

log_line("Filtering complete. Model N: ", nrow(dat_model), 
         ". Patients: ", n_distinct(dat_model$patient_id))

# Baseline age is already derived and centred in A01 (baseline_age_years_c),
# and stored in analysis_long.parquet. We reuse that here directly to avoid
# duplicating logic and accidentally creating suffixes (e.g. *.x, *.y).


# 3. Modeling Logic (Triple Fallback) -------------------------------------

compute_use_flags <- function(df) {
  nlev <- function(x) length(unique(na.omit(x)))
  list(
    use_blf_cat = nlev(df$baseline_fev1_cat) >= 2,
    use_sex   = nlev(df$sex) >= 2,
    use_nut    = nlev(df$nutrition_cat) >= 2,
    use_exac   = nlev(df$exac_cat) >= 2,
    use_sector = nlev(df$sector_public_only) >= 2,
    use_eth    = nlev(df$ethnicity_caucasian) >= 2,
    use_geno   = nlev(df$f508_any) >= 2
  )
}

# The Robust Fitter Function
fit_triple_fallback <- function(rhs_formula, data, use_reml = FALSE) {
  
  # Optimization control: bobyqa is robust for clinical LMMs
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  
  # Formulas
  f_corr   <- as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 + time_cs | patient_id)"))  # Preferred: Correlated
  f_uncorr <- as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 + time_cs || patient_id)")) # Fallback 1: Uncorrelated
  f_int    <- as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 | patient_id)"))            # Fallback 2: Intercept Only

  try_fit <- function(formula, label) {
    log_line("  > Attempting ", label, "...")
    mod_local <- tryCatch({
      lmerTest::lmer(formula, data = data, REML = use_reml, control = ctrl)
    }, error = function(e) {
      message("    -> ", label, " failed: ", conditionMessage(e))
      return(NULL)
    })
    mod_local
  }

  # 1. Correlated random slope
  mod <- try_fit(f_corr, "Correlated Slopes (1 + time | id)")
  if (!is.null(mod)) {
    if (!lme4::isSingular(mod, tol = 1e-4)) {
      return(mod)
    } else {
      message("    -> Correlated model singular; trying simpler random-effects structure.")
    }
  }

  # 2. Uncorrelated random slope
  mod <- try_fit(f_uncorr, "Uncorrelated Slopes (1 + time || id)")
  if (!is.null(mod)) {
    if (!lme4::isSingular(mod, tol = 1e-4)) {
      return(mod)
    } else {
      message("    -> Uncorrelated model singular; falling back to random intercept only.")
    }
  }

  # 3. Random intercept only (return even if singular, as last resort)
  mod <- try_fit(f_int, "Random Intercept Only (1 | id)")
  return(mod)
}

build_rhs <- function(flags, version = c("m1", "m2")) {
  version <- match.arg(version)
  terms <- c("time_cs",
             "pancreatic_insufficiency",
             "cfrd_indicator",
             "baseline_age_years_c",
             "baseline_fev1_cat",
             "baseline_pa_any",
             "sex",
             "ethnicity_caucasian",
             "sector_public_only",
             "f508_any")

  if (flags$use_nut)    terms <- c(terms, "nutrition_cat")
  if (flags$use_exac)   terms <- c(terms, "exac_cat")
  # binary versions already included above for sector/ethnicity/F508

  # INFECTION TERMS (level effects)
  # Use time-varying chronic indicators for B. cepacia, MSSA, MRSA only.
  # Baseline Pseudomonas is handled via baseline_pa_any; we drop the separate
  # chronic_paeruginosa level term to avoid duplicate constructs.
  terms <- c(terms, "chronic_bcepacia", "chronic_mssa", "chronic_mrsa")

  # MODULATOR TERMS
  terms <- c(terms, "cftr_modulator_indicator") # Time Varying for both models

  # INTERACTIONS (Time * X)
  interactions <- c()
  if (version == "m2") {
    interactions <- c("time_cs:pancreatic_insufficiency_slope", "time_cs:cfrd_indicator_slope")
    if (flags$use_nut)     interactions <- c(interactions, "time_cs:nutrition_cat")
    if (flags$use_blf_cat) interactions <- c(interactions, "time_cs:baseline_fev1_cat")
    interactions <- c(interactions, "time_cs:baseline_pa_any")
    if (flags$use_sex)     interactions <- c(interactions, "time_cs:sex")
    if (flags$use_sector)  interactions <- c(interactions, "time_cs:sector_public_only")
    if (flags$use_eth)     interactions <- c(interactions, "time_cs:ethnicity_caucasian")
    if (flags$use_geno)    interactions <- c(interactions, "time_cs:f508_any")
    # Add slope terms for CFTR and infections (sector interactions now included above when available)
    interactions <- c(interactions,
                      "time_cs:cftr_modulator_indicator_slope",
                      "time_cs:chronic_paeruginosa_slope",
                      "time_cs:chronic_mssa_slope")
  }

  paste(c(terms, interactions), collapse = " + ")
}


# 4. Execution ------------------------------------------------------------

flags_main <- compute_use_flags(dat_model)

# CFTR index-anchored time for sensitivity analysis (last pre-CFTR FEV1pp)
cftr_index <- dat_model %>%
  filter(
    ever_cftr_yes == 1L,
    !is.na(cftr_modulator_start_date),
    visit_date < cftr_modulator_start_date,
    !is.na(fev1_pp),
    is.finite(fev1_pp)
  ) %>%
  arrange(patient_id, visit_date) %>%
  group_by(patient_id) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  transmute(patient_id, cftr_index_date = visit_date)

dat_model_cftr <- dat_model %>%
  left_join(cftr_index, by = "patient_id") %>%
  mutate(
    # For never-treated, or if index missing, fall back to baseline_date
    cftr_index_date = dplyr::coalesce(cftr_index_date, baseline_date),
    time_cftr_index = as.numeric(difftime(visit_date, cftr_index_date, units = "days")) / 365.25,
    # For CFTR index-anchored model, treat time in YEARS (centered at 0 = index),
    # so that the CFTR interaction coefficient is interpretable as pp/year.
    time_cftr_index_cs = time_cftr_index,
    # Indicator for post-index period among ever-treated patients
    cftr_post_index = dplyr::if_else(
      slope_cftr_ok == 1L & ever_cftr_yes == 1L & !is.na(time_cftr_index) & time_cftr_index >= 0,
      1L,
      0L,
      missing = 0L
    )
  )

# Helper to build CFTR index-anchored RHS (no CFTR level jump, slope change only)
build_rhs_cftr <- function(flags) {
  terms <- c(
    "time_cftr_index_cs",
    "pancreatic_insufficiency",
    "cfrd_indicator",
    "baseline_age_years_c",
    "baseline_fev1_cat",
    "sex",
    "ethnicity_caucasian",
    "sector_public_only",
    "f508_any"
  )
  if (flags$use_nut)    terms <- c(terms, "nutrition_cat")
  if (flags$use_exac)   terms <- c(terms, "exac_cat")

  # Infection terms
  terms <- c(terms,
             "chronic_paeruginosa", "chronic_bcepacia",
             "chronic_mssa", "chronic_mrsa")

  # Interactions with time (including CFTR post-index slope change)
  interactions <- c(
    "time_cftr_index_cs:pancreatic_insufficiency_slope",
    "time_cftr_index_cs:cfrd_indicator_slope"
  )
  if (flags$use_nut)    interactions <- c(interactions, "time_cftr_index_cs:nutrition_cat")
  if (flags$use_blf_cat) interactions <- c(interactions, "time_cftr_index_cs:baseline_fev1_cat")
  if (flags$use_sex)    interactions <- c(interactions, "time_cftr_index_cs:sex")
  if (flags$use_sector) interactions <- c(interactions, "time_cftr_index_cs:sector_public_only")
  if (flags$use_eth)    interactions <- c(interactions, "time_cftr_index_cs:ethnicity_caucasian")
  if (flags$use_geno)   interactions <- c(interactions, "time_cftr_index_cs:f508_any")
  interactions <- c(
    interactions,
    "time_cftr_index_cs:chronic_paeruginosa_slope",
    "time_cftr_index_cs:chronic_mssa_slope",
    "time_cftr_index_cs:cftr_post_index"
  )

  paste(c(terms, interactions), collapse = " + ")
}

# --- M0: Unconditional Growth ---
log_line("Fitting M0...")
m0 <- fit_triple_fallback("time_cs", dat_model)

# --- M1: Main Effects ---
log_line("Fitting M1...")
rhs_m1 <- build_rhs(flags_main, "m1")
m1 <- fit_triple_fallback(rhs_m1, dat_model)

# --- M2: Interactions ---
log_line("Fitting M2...")
rhs_m2 <- build_rhs(flags_main, "m2")
m2 <- fit_triple_fallback(rhs_m2, dat_model)

# --- M2 REML (For Reporting Variance) ---
log_line("Fitting M2 (REML)...")
m2_reml <- fit_triple_fallback(rhs_m2, dat_model, use_reml = TRUE)

# --- CFTR index-anchored model (slope-only CFTR effect) ---
log_line("Fitting CFTR index-anchored model (M2_cftr)...")
rhs_cftr <- build_rhs_cftr(flags_main)
m2_cftr <- fit_triple_fallback(rhs_cftr, dat_model_cftr, use_reml = TRUE)

# 5. Diagnostics & Exports helpers --------------------------------------------

# Helper: Export RDS and CSV
export_result <- function(name, value, type = c("rds", "csv")) {
  assign(name, value, envir = .GlobalEnv)
  if ("rds" %in% type) saveRDS(value, file.path(results_dir, paste0(name, ".rds")))
  if ("csv" %in% type && (is.data.frame(value) || is_tibble(value))) {
    write_csv(value, file.path(tables_dir, paste0(name, ".csv")))
  }
}

# Helper: Tidy and Add FDR Correction
tidy_and_adjust <- function(mod) {
  if (is.null(mod)) return(NULL)
  broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr")) %>% # Benjamini-Hochberg
    relocate(p.value.adj, .after = p.value)
}

# --- Event-anchored models for other events (post vs pre slope change) -------
event_specs <- list(
  list(name = "pa",    start_col = "pa_chronic_start_date",   slope_ok_col = "slope_pa_ok"),
  list(name = "mssa",  start_col = "mssa_chronic_start_date", slope_ok_col = "slope_mssa_ok"),
  list(name = "pancreas", start_col = "pancreas_start_date", slope_ok_col = "slope_pancreas_ok"),
  list(name = "cfrd",  start_col = "cfrd_start_date",         slope_ok_col = "slope_cfrd_ok")
)

build_event_data <- function(df, spec) {
  time_var <- paste0("time_", spec$name, "_index")
  time_cs_var <- paste0(time_var, "_cs")
  post_var <- paste0(spec$name, "_post_index")
  df %>%
    mutate(
      event_start = .data[[spec$start_col]],
      anchor_date = dplyr::coalesce(event_start, baseline_date),
      !!time_var := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
      # Use event time in YEARS (no scaling) so event-anchored
      # slopes are directly interpretable as pp/year.
      !!time_cs_var := .data[[time_var]],
      !!post_var := dplyr::if_else(
        .data[[spec$slope_ok_col]] == 1L &
          !is.na(event_start) &
          visit_date >= event_start,
        1L, 0L,
        missing = 0L
      )
    ) %>%
    select(-event_start, -anchor_date)
}

build_rhs_event <- function(flags, time_var_cs, post_var) {
  terms <- c(
    time_var_cs, post_var,
    "pancreatic_insufficiency",
    "cfrd_indicator",
    "baseline_fev1_pp_c",
    "baseline_age_years_c"
  )
  if (flags$use_nut)    terms <- c(terms, "nutrition_cat")
  if (flags$use_exac)   terms <- c(terms, "exac_cat")
  if (flags$use_sector) terms <- c(terms, "sector")
  if (flags$use_eth)    terms <- c(terms, "ethnicity_cat")
  if (flags$use_geno)   terms <- c(terms, "f508_genotype_group")

  terms <- c(terms,
             "chronic_paeruginosa", "chronic_bcepacia",
             "chronic_mssa", "chronic_mrsa",
             "cftr_modulator_indicator")

  interactions <- c(paste0(time_var_cs, ":", post_var))

  paste(c(terms, interactions), collapse = " + ")
}

event_model_terms <- list()
for (spec in event_specs) {
  time_var <- paste0("time_", spec$name, "_index")
  time_cs_var <- paste0(time_var, "_cs")
  post_var <- paste0(spec$name, "_post_index")
  rhs_ev <- build_rhs_event(flags_main, time_cs_var, post_var)
  dat_ev <- build_event_data(dat_model, spec)
  log_line("Fitting event-anchored model for ", spec$name, "...")
  mod_ev <- fit_triple_fallback(rhs_ev, dat_ev, use_reml = TRUE)
  export_result(paste0(".model_", spec$name, "_anchored"), mod_ev, "rds")
  export_result(paste0(".model_", spec$name, "_anchored_tidy"), tidy_and_adjust(mod_ev), "csv")
  event_model_terms[[spec$name]] <- list(time_cs_var = time_cs_var, post_var = post_var)
}


# 6. Export Diagnostics & Summaries -------------------------------------------

log_line("Exporting results...")

# 1. Model Objects (RDS)
export_result(".model_m0", m0, "rds")
export_result(".model_m1", m1, "rds")
export_result(".model_m2", m2, "rds")
export_result(".model_m2_reml", m2_reml, "rds")
export_result(".model_m2_cftr", m2_cftr, "rds")

# 2. Tidy Summaries (CSV)
export_result(".model_m0_tidy", tidy_and_adjust(m0), "csv")
export_result(".model_m1_tidy", tidy_and_adjust(m1), "csv")
export_result(".model_m2_tidy", tidy_and_adjust(m2), "csv")
export_result(".model_m2_reml_tidy", tidy_and_adjust(m2_reml), "csv")
export_result(".model_m2_cftr_tidy", tidy_and_adjust(m2_cftr), "csv")

# 3. Variance Components
if(!is.null(m2)) export_result(".model_m2_vc", as_tibble(VarCorr(m2)), "csv")

# 4. R-Squared
calc_r2 <- function(mod) {
  if(is.null(mod)) return(NULL)
  r <- performance::r2_nakagawa(mod)
  icc_val <- NA_real_
  icc_res <- try(performance::icc(mod), silent = TRUE)
  if (!inherits(icc_res, "try-error")) {
    if (is.list(icc_res) && !is.null(icc_res$ICC_adjusted)) {
      icc_val <- icc_res$ICC_adjusted
    } else if (is.numeric(icc_res) && length(icc_res) >= 1) {
      icc_val <- as.numeric(icc_res[1])
    }
  }
  tibble(
    marginal_r2 = r$R2_marginal,
    conditional_r2 = r$R2_conditional,
    icc = icc_val
  )
}
export_result(".model_m2_r2", calc_r2(m2), "csv")

# 5. Collinearity Check (VIF)
log_line("Checking Collinearity on M2...")
col_m2 <- try(performance::check_collinearity(m2), silent = TRUE)
if (!inherits(col_m2, "try-error")) {
  export_result(".model_m2_collinearity", as_tibble(col_m2), "csv")
} else {
  log_line("Collinearity check skipped (Model error).")
}

# 6. ANOVA Comparison
if(!is.null(m0) & !is.null(m1) & !is.null(m2)) {
  export_result(".model_m_anova", broom.mixed::tidy(anova(m0, m1, m2)), "csv")
}

log_line("Agent02 Completed Successfully.")
