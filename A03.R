#!/usr/bin/env Rscript

# TABLE GENERATOR 
# 1. Fixed Effects (Level vs Slope)
# 2. Random Effects (Variance/SD)
# 3. Model Performance (R2, AIC)
# 4. Diagnostics (VIF, t-stats)

set.seed(20251119)

required_pkgs <- c("tidyverse", "broom.mixed", "performance", "gt", "insight", "lme4", "arrow")
to_install <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Define Directories
out_dir <- "outputs/tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load Model (Preferably the REML version for final reporting)
model_path <- "outputs/models/.model_m2_reml.rds"
if (file.exists(model_path)) {
  mod <- readRDS(model_path)
} else {
  # Fallback to searching environment if running interactively
  if (exists(".model_m2_reml", envir = .GlobalEnv)) {
    mod <- get(".model_m2_reml", envir = .GlobalEnv)
  } else {
    stop("Model file not found at ", model_path)
  }
}

# Basic sanity check: ensure we actually loaded a fitted mixed model
if (is.null(mod) || !inherits(mod, "lmerMod")) {
  stop(
    "Loaded model object is NULL or not a valid lmer mixed model.\n",
    "Check that Agent 02 (A02.R) ran successfully and that model M2 (REML) converged.\n",
    "Expected model RDS at: ", model_path
  )
}

# --- HELPER FUNCTIONS --------------------------------------------------------

fmt_num <- function(x, k=2) sprintf(paste0("%.", k, "f"), x)
fmt_p   <- function(p) {
  ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# 1. EXTRACT FIXED EFFECTS ----------------------------------------------------
# We need Estimate, CI, t-stat, p-value
fixed_raw <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
  # Rescale continuous baseline covariates for reporting
  mutate(
    scale_factor = dplyr::case_when(
      term == "baseline_fev1_pp_c" ~ 10,
      term == "baseline_age_years_c" ~ 10,
      TRUE ~ 1
    ),
    estimate = estimate * scale_factor,
    std.error = std.error * scale_factor,
    conf.low = conf.low * scale_factor,
    conf.high = conf.high * scale_factor
  ) %>%
  select(-scale_factor)

# Calculate VIF (if possible) – kept for separate diagnostics via A02, not shown in this table
vif_dat <- try(performance::check_collinearity(mod), silent = TRUE)

# Process Fixed Effects Table
fixed_tab <- fixed_raw %>%
  # Drop baseline-time slope interactions for event-based terms; these will
  # be presented using dedicated event-time (index-anchored) models instead.
  filter(
    !grepl("^time_cs:(cftr_modulator_indicator_slope|chronic_paeruginosa_slope|chronic_mssa_slope)", term),
    # Exclude baseline FEV1 strata (not needed in Table 2 / forest)
    !grepl("^baseline_fev1_cat", term),
    !grepl("^time_cs:baseline_fev1_cat", term)
  ) %>%
  mutate(
    # Classify Terms
    type = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "time_cs" ~ "Time (Slope)",
      grepl(":", term) ~ "Difference in Annual FEV1pp Change (Slope) At Baseline",
      TRUE ~ "Difference in Baseline FEV1pp"
    ),
    type = factor(type, levels = c("Intercept",
                                   "Time (Slope)",
                                   "Difference in Baseline FEV1pp",
                                   "Difference in Annual FEV1pp Change (Slope) At Baseline"),
                  ordered = TRUE),
    # Clean Term Names
    clean_term = term, # Start with raw
    clean_term = str_remove_all(clean_term, "time_cs:"), # Remove interaction prefix for clean reading
    clean_term = case_when(
      term == "(Intercept)" ~ "Intercept (Baseline FEV1pp)",
      term == "time_cs" ~ "Annual Rate of Decline (Global)",
      
      # Manual Mapping for Cleanliness
      term == "baseline_age_years_c" ~ "Baseline age, years (per 10 years above mean)",
      term == "baseline_fev1_pp_c" ~ "Baseline FEV1 % predicted (per 10 pp above mean)",
      term == "baseline_fev1_catModerate" ~ "Baseline FEV1 category: Moderate (vs Normal)",
      term == "baseline_fev1_catSevere" ~ "Baseline FEV1 category: Severe (vs Normal)",
      term == "sexMale" ~ "Male (vs Female)",
      term == "ethnicity_caucasian" ~ "Caucasian vs non-Caucasian",
      term == "sector_public_only" ~ "Public only vs Private/Combined",
      term == "f508_any" ~ "F508 allele (1 or 2) vs other",
      
      term == "cftr_modulator_indicator" ~ "Index FEV1pp on CFTR",
      grepl("^time_cs:cftr_modulator_indicator", term) ~ "CFTR Modulator: Change in annual decline",
      
      term == "nutrition_catUnderweight" ~ "Underweight (vs Normal)",
      term == "nutrition_catOverweight" ~ "Overweight (vs Normal)",
      term == "time_cs:nutrition_catUnderweight" ~ "Slope: Underweight vs Normal",
      term == "time_cs:nutrition_catOverweight" ~ "Slope: Overweight vs Normal",
      
      term == "exac_cat1" ~ "Exacerbations: 1/yr (vs 0)",
      term == "exac_cat2+" ~ "Exacerbations: \u22652/yr (vs 0)",

      term == "baseline_pa_any" ~ "Baseline P. aueruginosa",
      
      term == "cfrd_indicator" ~ "CFRD",
      term == "pancreatic_insufficiency" ~ "Pancreatic Insufficiency",
      grepl("^time_cs:cfrd_indicator", term) ~ "Slope: CFRD",
      grepl("^time_cs:pancreatic_insufficiency", term) ~ "Slope: Pancreatic Insufficiency",
      grepl("^time_cs:sex", term) ~ "Slope: Male vs Female",
      grepl("^time_cs:baseline_fev1_catModerate", term) ~ "Slope: Baseline FEV1 Moderate vs Normal",
      grepl("^time_cs:baseline_fev1_catSevere", term) ~ "Slope: Baseline FEV1 Severe vs Normal",
      grepl("^time_cs:baseline_pa_any", term) ~ "Slope: Baseline P. aeruginosa",
      grepl("^time_cs:ethnicity_caucasian", term) ~ "Slope: Caucasian vs non-Caucasian",
      grepl("^time_cs:sector_public_only", term) ~ "Slope: Public only vs Private/Combined",
      grepl("^time_cs:f508_any", term) ~ "Slope: F508 allele (1 or 2) vs other",
      
      term == "ethnicity_catMixed Race" ~ "Mixed Race vs Black African",
      term == "ethnicity_catIndian" ~ "Indian vs Black African",
      term == "ethnicity_catCaucasian" ~ "Caucasian vs Black African",
      grepl("^time_cs:ethnicity_catMixed Race", term) ~ "Slope: Mixed Race vs Black African",
      grepl("^time_cs:ethnicity_catIndian", term) ~ "Slope: Indian vs Black African",
      grepl("^time_cs:ethnicity_catCaucasian", term) ~ "Slope: Caucasian vs Black African",
      
      term == "sectorPrivate" ~ "Private vs Public",
      term == "sectorCombined" ~ "Combined vs Public",
      term == "time_cs:sectorPrivate" ~ "Slope: Private vs Public",
      term == "time_cs:sectorCombined" ~ "Slope: Combined vs Public",
      term == "chronic_paeruginosa" ~ "Chronic P. aeruginosa",
      term == "chronic_bcepacia" ~ "Chronic B. cepacia",
      term == "chronic_mssa" ~ "Chronic MSSA",
      term == "chronic_mrsa" ~ "Chronic MRSA",
      grepl("^time_cs:chronic_paeruginosa", term) ~ "Slope: Chronic P. aeruginosa",
      grepl("^time_cs:chronic_bcepacia", term) ~ "Slope: Chronic B. cepacia",
      grepl("^time_cs:chronic_mssa", term) ~ "Slope: Chronic MSSA",
      grepl("^time_cs:chronic_mrsa", term) ~ "Slope: Chronic MRSA",
      
      term == "f508_genotype_groupHeterozygous" ~ "F508 Heterozygous vs Other",
      term == "f508_genotype_groupHomozygous" ~ "F508 Homozygous vs Other",
      grepl("^time_cs:f508_genotype_groupHeterozygous", term) ~ "Slope: F508 Heterozygous vs Other",
      grepl("^time_cs:f508_genotype_groupHomozygous", term) ~ "Slope: F508 Homozygous vs Other",
      
      TRUE ~ clean_term # Fallback
    ),
    # Subgrouping within Type (for ordering similar to previous main LME table)
    subgroup = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "time_cs" ~ "Time: Global slope",
      grepl("^time_cs:", term) ~ "Time: Slope interactions",
      
      type == "Baseline Association (Level)" & term == "baseline_age_years_c" ~ "Baseline age (continuous)",
      type == "Baseline Association (Level)" & term == "baseline_fev1_pp_c" ~ "Baseline FEV1 (continuous)",
      type == "Baseline Association (Level)" & grepl("^baseline_fev1_cat", term) ~ "Baseline FEV1 category",
      type == "Baseline Association (Level)" & grepl("^nutrition_cat", term) ~ "Nutrition status",
      type == "Baseline Association (Level)" & term %in% c("cfrd_indicator", "pancreatic_insufficiency") ~ "Metabolic",
      type == "Baseline Association (Level)" & grepl("^sex", term) ~ "Sex",
      type == "Baseline Association (Level)" & grepl("^ethnicity_(cat|caucasian)", term) ~ "Ethnicity",
      type == "Baseline Association (Level)" & grepl("^sector", term) ~ "Sector of care",
      type == "Baseline Association (Level)" & grepl("^exac_cat", term) ~ "Exacerbations",
      type == "Baseline Association (Level)" & (grepl("^chronic_(paeruginosa|bcepacia|mssa|mrsa)$", term) |
                                                grepl("^baseline_pa_", term)) ~ "Chronic infections",
      type == "Baseline Association (Level)" & term == "cftr_modulator_indicator" ~ "CFTR modulator therapy",
      type == "Baseline Association (Level)" & (grepl("^f508_genotype_group", term) | grepl("^f508_any", term)) ~ "F508 genotype",
      
      type == "Interaction (Rate of Change)" & grepl("^time_cs:cftr_modulator_indicator", term) ~ "Slope: CFTR modulator",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:nutrition_cat", term) ~ "Slope: Nutrition",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:baseline_fev1_cat", term) ~ "Slope: Baseline FEV1",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:sex", term) ~ "Slope: Sex",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:(pancreatic_insufficiency|pancreatic_insufficiency_slope|cfrd_indicator|cfrd_indicator_slope)$", term) ~ "Slope: Metabolic",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:chronic_(paeruginosa|bcepacia|mssa|mrsa)(_slope)?$", term) ~ "Slope: Chronic infections",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:sector", term) ~ "Slope: Sector of care",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:(ethnicity_cat|ethnicity_caucasian)", term) ~ "Slope: Ethnicity",
      type == "Interaction (Rate of Change)" & grepl("^time_cs:(f508_genotype_group|f508_any)", term) ~ "Slope: F508 genotype",
      
      TRUE ~ "Other"
    ),
    subgroup_order = case_when(
      subgroup == "Intercept" ~ 1L,
      subgroup == "Time: Global slope" ~ 2L,
      subgroup == "Time: Slope interactions" ~ 3L,
      
      subgroup == "Baseline FEV1 (continuous)" ~ 10L,
      subgroup == "Baseline age (continuous)" ~ 11L,
      subgroup == "Baseline FEV1 category" ~ 12L,
      subgroup == "Sex" ~ 13L,
      subgroup == "Nutrition status" ~ 14L,
      subgroup == "Metabolic" ~ 15L,
      subgroup == "Ethnicity" ~ 16L,
      subgroup == "Sector of care" ~ 17L,
      subgroup == "Exacerbations" ~ 18L,
      subgroup == "Chronic infections" ~ 19L,
      subgroup == "CFTR modulator therapy" ~ 20L,
      subgroup == "F508 genotype" ~ 21L,
      
      subgroup == "Slope: CFTR modulator" ~ 30L,
      subgroup == "Slope: Nutrition" ~ 31L,
      subgroup == "Slope: Metabolic" ~ 32L,
      subgroup == "Slope: Sex" ~ 33L,
      subgroup == "Slope: Chronic infections" ~ 34L,
      subgroup == "Slope: Sector of care" ~ 35L,
      subgroup == "Slope: Ethnicity" ~ 36L,
      subgroup == "Slope: F508 genotype" ~ 37L,
      
      TRUE ~ 99L
    )
  ) %>%
  arrange(type, subgroup_order, term)

# 2. EXTRACT RANDOM EFFECTS ---------------------------------------------------
# Standard: Variance and SD
rand_raw <- as.data.frame(VarCorr(mod))
rand_tab <- rand_raw %>%
  mutate(
    term_label = case_when(
      var1 == "(Intercept)" & is.na(var2) ~ "Random Intercept (Patient)",
      var1 == "time_cs" & is.na(var2) ~ "Random Slope (Time)",
      !is.na(var2) ~ "Correlation (Intercept, Slope)",
      grp == "Residual" ~ "Residual Error"
    ),
    stats = paste0("Var: ", fmt_num(vcov, 2), " (SD: ", fmt_num(sdcor, 2), ")")
  ) %>%
  select(clean_term = term_label, Estimate = stats) %>%
  mutate(type = "Random Effects Structure")

# 3. EXTRACT MODEL FIT --------------------------------------------------------
r2 <- performance::r2_nakagawa(mod)
fit_tab <- tibble(
  clean_term = c("Observations", "Patients", "Marginal R2 (Fixed)", "Conditional R2 (Fixed+Random)", "AIC"),
  Estimate = c(
    as.character(nobs(mod)),
    as.character(ngrps(mod)[1]),
    fmt_num(r2$R2_marginal, 3),
    fmt_num(r2$R2_conditional, 3),
    fmt_num(AIC(mod), 1)
  ),
  type = "Model Performance"
)

# 4. COMBINE AND FORMAT -------------------------------------------------------

# Prepare Main Rows
base_slope_est <- fixed_tab %>%
  filter(term == "time_cs") %>%
  summarise(val = mean(estimate, na.rm = TRUE)) %>%
  pull(val)

main_rows <- fixed_tab %>%
  mutate(
    slope_suffix = ifelse(type %in% c("Time (Slope)", "Interaction (Rate of Change)"), " %/yr", " %"),
    p_label_raw = fmt_p(p.value),
    p_label = ifelse(!is.na(p.value) & p.value <= 0.05,
                     paste0("**", p_label_raw, "**"),
                     p_label_raw),
    t_label = fmt_num(statistic, 2),
    net_effect_raw = dplyr::case_when(
      term == "time_cs" ~ estimate,
      grepl("^time_cs:", term) & is.finite(base_slope_est) ~ base_slope_est + estimate,
      TRUE ~ estimate
    ),
    # est_ci: show estimate + CI for all types (units once)
    est_ci = paste0(
      fmt_num(estimate), slope_suffix,
      " (", fmt_num(conf.low), " to ", fmt_num(conf.high), ")"
    ),
    net_effect = dplyr::case_when(
      is.na(net_effect_raw) ~ "-",
      type %in% c("Time (Slope)",
                  "Difference in Annual FEV1pp Change (Slope) At Baseline",
                  "Difference in Annual FEV1pp Change (Slope) After Event") ~ paste0(
        fmt_num(net_effect_raw), " %/yr",
        " (", fmt_num(conf.low), " to ", fmt_num(conf.high), ")"
      ),
      TRUE ~ ""
    )
  ) %>%
  select(type, clean_term, est_ci, net_effect, t_label, p_label)

cftr_index_row <- NULL
# Event-anchored slopes (post vs pre) for other events
event_specs <- tibble::tribble(
  ~name,      ~label,                                                       ~term,
  "cftr",     "Slope: After CFTR start (event-anchored)",                   "time_cftr_index_cs:cftr_post_index",
  "pa",       "Slope: Chronic P. aeruginosa (event-anchored)",              "time_pa_index_cs:pa_post_index",
  "mssa",     "Slope: Chronic MSSA (event-anchored)",                       "time_mssa_index_cs:mssa_post_index",
  "exac1",    "Slope: Exacerbations 1/year (event-anchored)",               "time_exac_index_cs:exac1_post_index",
  "exac2plus","Slope: Exacerbations ≥2/year (event-anchored)",              "time_exac_index_cs:exac2plus_post_index"
)

event_rows <- purrr::map_dfr(seq_len(nrow(event_specs)), function(i) {
  spec <- event_specs[i, ]
  file <- if (spec$name == "cftr") {
    file.path(out_dir, ".model_m2_cftr_tidy.csv")
  } else {
    file.path(out_dir, paste0(".model_", spec$name, "_anchored_tidy.csv"))
  }
  if (!file.exists(file)) return(NULL)
  df <- readr::read_csv(file, show_col_types = FALSE)
  row <- df %>% dplyr::filter(term == spec$term)
  if (nrow(row) != 1) return(NULL)
  net_val <- if (is.finite(base_slope_est)) base_slope_est + row$estimate else row$estimate
  tibble(
    type = "Difference in Annual FEV1pp Change (Slope) After Event",
    clean_term = spec$label,
    est_ci = paste0(fmt_num(row$estimate), " %/yr (", fmt_num(row$conf.low), " to ", fmt_num(row$conf.high), ")"),
    net_effect = paste0(fmt_num(net_val), " %/yr (", fmt_num(row$conf.low), " to ", fmt_num(row$conf.high), ")"),
    t_label = fmt_num(row$statistic, 2),
    p_label = ifelse(!is.na(row$p.value) && row$p.value <= 0.05,
                     paste0("**", fmt_p(row$p.value), "**"),
                     fmt_p(row$p.value))
  )
})

event_sup <- purrr::map_dfr(seq_len(nrow(event_specs)), function(i) {
  spec <- event_specs[i, ]
  file <- if (spec$name == "cftr") {
    file.path(out_dir, ".model_m2_cftr_tidy.csv")
  } else {
    file.path(out_dir, paste0(".model_", spec$name, "_anchored_tidy.csv"))
  }
  if (!file.exists(file)) return(NULL)
  df <- readr::read_csv(file, show_col_types = FALSE)
  row <- df %>% dplyr::filter(term == spec$term)
  if (nrow(row) != 1) return(NULL)
  tibble(
    event = spec$label,
    estimate = row$estimate,
    ci_low = row$conf.low,
    ci_high = row$conf.high,
    p_value = row$p.value
  )
})

# Supplementary table for event-anchored slopes
if (nrow(event_sup)) {
  event_sup <- event_sup %>%
    mutate(
      estimate_fmt = ifelse(is.finite(estimate), sprintf("%.2f", estimate), "-"),
      ci_fmt = dplyr::case_when(
        is.finite(ci_low) & is.finite(ci_high) ~ paste0(sprintf("%.2f", ci_low), " to ", sprintf("%.2f", ci_high)),
        TRUE ~ "-"
      ),
      p_fmt = dplyr::case_when(
        is.na(p_value) ~ "-",
        p_value < 0.001 ~ "<0.001",
        is.finite(p_value) ~ sprintf("%.3f", p_value),
        TRUE ~ "-"
      )
    )

  event_sup_path_csv <- file.path(out_dir, "Table_Event_Anchored_Slopes.csv")
  event_sup_path_html <- file.path(out_dir, "Table_Event_Anchored_Slopes.html")

  readr::write_csv(event_sup, event_sup_path_csv)

  gt_event <- event_sup %>%
    select(event, estimate_fmt, ci_fmt, p_fmt) %>%
    gt::gt() %>%
    gt::cols_label(
      event = "Event (index-anchored)",
      estimate_fmt = "Estimate (pp/year)",
      ci_fmt = "95% CI",
      p_fmt = "p-value"
    ) %>%
    gt::tab_header(
      title = gt::md("**Event-anchored slope changes (post vs pre)**")
    )

  gt::gtsave(gt_event, event_sup_path_html)
}

# Prepare Footer Rows (Random Effs + Fit)
# We need to align columns. We will put value in 'est_ci' and leave others blank
footer_rows <- bind_rows(rand_tab, fit_tab) %>%
  mutate(
    t_label = "", p_label = "",
    est_ci = Estimate,
    net_effect = ""
  ) %>%
  select(type, clean_term, est_ci, net_effect, t_label, p_label)

# Bind All: main LMM rows, then event-time slopes, then random effects / fit
if (nrow(event_rows)) {
  final_df <- bind_rows(main_rows, event_rows, footer_rows)
} else {
  final_df <- bind_rows(main_rows, footer_rows)
}

# 5. GENERATE GT TABLE --------------------------------------------------------

gt_tbl <- final_df %>%
  gt(groupname_col = "type", rowname_col = "clean_term") %>%
  cols_label(
    clean_term = "Covariate",
    est_ci = "FEV1pp (β)",
    net_effect = "Adjusted Annual Slope",
    t_label = "t-statistic",
    p_label = "p-value"
  ) %>%
  tab_header(
    title = md("**Longitudinal Linear Mixed Model Results (Table 2)**"),
    subtitle = "Outcome: FEV1 % Predicted (pp)"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  cols_align(
    align = c("right"),
    col= p_label
  ) %>%
  fmt_markdown(columns = p_label) %>%
  cols_width(
    clean_term ~ px(350),
    est_ci ~ px(260),
    net_effect ~ px(250),
    t_label ~ px(100),
    p_label ~ px(120)
  ) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    row_group.border.top.color = "black",
    row_group.border.bottom.color = "black",
    table_body.border.bottom.color = "black"
  ) %>%
  tab_source_note(
    source_note = "Abbreviations: CI, Confidence Interval; VIF, Variance Inflation Factor; AIC, Akaike Information Criterion; R2, Nakagawa's R-squared."
  ) %>%
  tab_source_note(
    source_note = "FEV1pp (β) refers to the measured effect of a covariate on FEV1pp. Adjusted Annual Slope is the new change in FEV1pp per year, calculated as the global FEV1pp decline plus the covariate effect."
  )

# 6. EXPORT -------------------------------------------------------------------

# Save HTML
html_path <- file.path(out_dir, "Table_Scientific_Format.html")
gtsave(gt_tbl, html_path)

# Save CSV (Raw)
csv_path <- file.path(out_dir, "Table_Scientific_Format.csv")
write_csv(final_df, csv_path)

message("Success! Tables generated:")
message("  -> HTML (Viewer): ", html_path)
message("  -> CSV (Excel):   ", csv_path)
