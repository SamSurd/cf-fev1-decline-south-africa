#!/usr/bin/env Rscript

# AGENT 04 – TABLE 1 GENERATOR (BASELINE CHARACTERISTICS)
# Creates a comprehensive baseline “Table 1” with:
#  - Demographics
#  - Disease severity and nutrition
#  - Metabolic status
#  - Ethnicity and sector of care
#  - Genotype and chronic infections
#  - CFTR modulator exposure
# Columns are Overall, Never on CFTR, Ever on CFTR.

set.seed(20251119)

required_pkgs <- c("tidyverse", "arrow", "gt")
to_install <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, quiet = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# --------------------------------------------------------------------
# 1. Paths and Data --------------------------------------------------
# --------------------------------------------------------------------

data_file <- "outputs/data_derived/analysis_long.parquet"
out_dir   <- "outputs/tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (!file.exists(data_file)) {
  stop("Expected data file not found at: ", data_file)
}

dat <- arrow::read_parquet(data_file) %>% as_tibble()

# Derive ever-CFRD status across all visits
ever_cfrd_df <- dat %>%
  group_by(patient_id) %>%
  summarise(
    ever_cfrd = dplyr::case_when(
      any(cfrd == "Yes", na.rm = TRUE) ~ "Yes",
      any(cfrd == "No",  na.rm = TRUE) ~ "No",
      TRUE                              ~ NA_character_
    ),
    .groups = "drop"
  )

# One row per patient (baseline record)
baseline <- dat %>%
  arrange(patient_id, visit_date) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  left_join(ever_cfrd_df, by = "patient_id")

# --------------------------------------------------------------------
# 2. Helper Functions ------------------------------------------------
# --------------------------------------------------------------------

fmt_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("-")
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  paste0(sprintf("%.1f", m), " (", sprintf("%.1f", s), ")")
}

fmt_n_pct <- function(n, d) {
  if (!is.finite(d) || d <= 0) return("-")
  paste0(n, " (", sprintf("%.1f", 100 * n / d), "%)")
}

summarise_cont <- function(df, var, label, domain) {
  v_all <- df[[var]]
  overall <- fmt_mean_sd(v_all)

  df_never <- df %>% filter(ever_cftr_yes == 0)
  df_ever  <- df %>% filter(ever_cftr_yes == 1)

  never_val <- if (nrow(df_never)) fmt_mean_sd(df_never[[var]]) else "-"
  ever_val  <- if (nrow(df_ever))  fmt_mean_sd(df_ever[[var]])  else "-"

  tibble(
    domain = domain,
    variable = label,
    level = "Mean (SD)",
    Overall = overall,
    `Never on CFTR` = never_val,
    `Ever on CFTR`  = ever_val
  )
}

summarise_cat <- function(df, var, label, domain) {
  v <- df[[var]]
  fac_all <- if (is.factor(v)) v else factor(v)
  lvls <- levels(fac_all)

  df_never <- df %>% filter(ever_cftr_yes == 0)
  df_ever  <- df %>% filter(ever_cftr_yes == 1)

  fac_never <- if (nrow(df_never)) {
    v0 <- df_never[[var]]
    if (is.factor(v0)) v0 else factor(v0, levels = lvls)
  } else factor(character(0), levels = lvls)

  fac_ever <- if (nrow(df_ever)) {
    v1 <- df_ever[[var]]
    if (is.factor(v1)) v1 else factor(v1, levels = lvls)
  } else factor(character(0), levels = lvls)

  denom_all   <- sum(!is.na(fac_all))
  denom_never <- sum(!is.na(fac_never))
  denom_ever  <- sum(!is.na(fac_ever))

  res <- map_dfr(lvls, function(lvl) {
    n_all   <- sum(fac_all == lvl, na.rm = TRUE)
    n_never <- sum(fac_never == lvl, na.rm = TRUE)
    n_ever  <- sum(fac_ever == lvl, na.rm = TRUE)

    # Drop rows with zero patients overall
    if (n_all == 0) return(NULL)

    tibble(
      domain = domain,
      variable = label,
      level = factor(lvl, levels = lvls),
      Overall = fmt_n_pct(n_all, denom_all),
      `Never on CFTR` = fmt_n_pct(n_never, denom_never),
      `Ever on CFTR`  = fmt_n_pct(n_ever, denom_ever)
    )
  })

  # For CFTR ever (0/1), relabel levels to No / Yes for readability
  if (label == "Ever on CFTR") {
    res <- res %>%
      mutate(
        level = dplyr::recode(
          level,
          "0" = "No",
          "1" = "Yes",
          .default = level
        )
      )
  }

  # For selected binary variables, keep only the "positive" category
  if (label %in% c("Sex", "Ever CFRD", "Pancreatic insufficiency",
                   "Ever on CFTR", "On CFTR at baseline",
                   "Ever chronic P. aeruginosa", "Ever chronic B. cepacia",
                   "Ever chronic MSSA", "Ever chronic MRSA")) {
    res <- res %>%
      filter(
        (label == "Sex" & level == "Male") |
        (label == "Ever CFRD" & level == "Yes") |
        (label == "Pancreatic insufficiency" & level == "Yes") |
        (label == "Ever on CFTR" & level == "Yes") |
        (label == "On CFTR at baseline" & level == "Yes") |
        (label == "Ever chronic P. aeruginosa" & level == "Yes") |
        (label == "Ever chronic B. cepacia" & level == "Yes") |
        (label == "Ever chronic MSSA" & level == "Yes") |
        (label == "Ever chronic MRSA" & level == "Yes")
      )
  }

  res
}

summarise_anchor <- function(df, var, label, domain) {
  v <- df[[var]]
  n_overall <- sum(v == 1, na.rm = TRUE)
  n_never   <- sum(df$ever_cftr_yes == 0 & v == 1, na.rm = TRUE)
  n_ever    <- sum(df$ever_cftr_yes == 1 & v == 1, na.rm = TRUE)
  fmt <- function(x) if (is.finite(x)) as.character(x) else "-"
  tibble(
    domain = domain,
    variable = label,
    level = "Patients with \u2265 2 post-event FEV1pp",
    Overall = fmt(n_overall),
    `Never on CFTR` = fmt(n_never),
    `Ever on CFTR`  = fmt(n_ever)
  )
}

# --------------------------------------------------------------------
# 3. Variable Preparation --------------------------------------------
# --------------------------------------------------------------------

baseline <- baseline %>%
  mutate(
    sex = factor(sex),
    ethnicity_cat = forcats::fct_explicit_na(ethnicity_cat, na_level = "Unknown"),
    sector = forcats::fct_explicit_na(sector, na_level = "Unknown"),
    baseline_age_group = forcats::fct_explicit_na(baseline_age_group, na_level = "Unknown"),
    baseline_age_group = recode(baseline_age_group, `6-10` = "06-10"),
    baseline_age_group = factor(baseline_age_group,
                                levels = c("06-10", "11-17", "18-25", "26-39", "40+", "Unknown")),
    baseline_fev1_cat = forcats::fct_explicit_na(baseline_fev1_cat, na_level = "Unknown"),
    baseline_nutrition_cat = forcats::fct_explicit_na(baseline_nutrition_cat, na_level = "Unknown"),
    baseline_exac_cat = forcats::fct_explicit_na(baseline_exac_cat, na_level = "Unknown"),
    f508_genotype_group = forcats::fct_explicit_na(f508_genotype_group, na_level = "Other"),
    baseline_pa_cat = forcats::fct_explicit_na(baseline_pa_cat, na_level = "Negative"),

    # Binary indicators as labelled factors
    cfrd_indicator = factor(dplyr::if_else(is.na(cfrd_indicator), NA_integer_, cfrd_indicator),
                            levels = c(0L, 1L), labels = c("No", "Yes")),
    pancreatic_insufficiency = factor(dplyr::if_else(is.na(pancreatic_insufficiency), NA_integer_, pancreatic_insufficiency),
                                      levels = c(0L, 1L), labels = c("No", "Yes")),
    ever_chronic_pa = factor(dplyr::if_else(is.na(ever_chronic_pa), NA_integer_, ever_chronic_pa),
                             levels = c(0L, 1L), labels = c("No", "Yes")),
    ever_chronic_bcepacia = factor(dplyr::if_else(is.na(ever_chronic_bcepacia), NA_integer_, ever_chronic_bcepacia),
                                   levels = c(0L, 1L), labels = c("No", "Yes")),
    ever_chronic_mssa = factor(dplyr::if_else(is.na(ever_chronic_mssa), NA_integer_, ever_chronic_mssa),
                               levels = c(0L, 1L), labels = c("No", "Yes")),
    ever_chronic_mrsa = factor(dplyr::if_else(is.na(ever_chronic_mrsa), NA_integer_, ever_chronic_mrsa),
                               levels = c(0L, 1L), labels = c("No", "Yes")),
    ever_cftr_yes = dplyr::if_else(is.na(ever_cftr_yes), 0L, ever_cftr_yes),
    baseline_cftr_mod_indicator = factor(dplyr::if_else(is.na(baseline_cftr_mod_indicator), NA_integer_, baseline_cftr_mod_indicator),
                                         levels = c(0L, 1L), labels = c("No", "Yes")),
    ever_cfrd = factor(ever_cfrd, levels = c("No", "Yes"))
  )

# --------------------------------------------------------------------
# 4. Build Table 1 ---------------------------------------------------
# --------------------------------------------------------------------

table1 <- bind_rows(
  # Demographics
  summarise_cont(baseline, "age_years", "Baseline Age", "Demographics"),
  summarise_cat(baseline, "sex", "Sex", "Demographics"),
  summarise_cat(baseline, "ethnicity_cat", "Ethnicity", "Demographics"),

  # Disease severity
  summarise_cont(baseline, "baseline_fev1_pp", "Baseline FEV1 % predicted", "Disease severity"),
  summarise_cat(baseline, "baseline_fev1_cat", "Baseline FEV1 Stratum", "Disease severity"),
  summarise_cat(baseline, "baseline_age_group", "Baseline Age Stratum", "Disease severity"),

  # Nutrition and exacerbations
  summarise_cat(baseline, "baseline_nutrition_cat", "Baseline Nutritional Status", "Nutrition"),
  summarise_cat(baseline, "baseline_exac_cat", "Exacerbations In Preceding Year", "Exacerbations"),

  # Metabolic
  summarise_cat(baseline, "ever_cfrd", "Ever CFRD", "Metabolic"),
  summarise_cat(baseline, "pancreatic_insufficiency", "Pancreatic insufficiency", "Metabolic"),

  # Sector of care
  summarise_cat(baseline, "sector", "Healthcare Sector", "Service"),

  # Genotype
  summarise_cat(baseline, "f508_genotype_group", "F508 genotype", "Genotype"),

  # Infections
  summarise_cat(baseline, "baseline_pa_cat", "Baseline Pseudomonas status", "Chronic infections"),
  summarise_cat(baseline, "ever_chronic_pa", "Ever chronic P. aeruginosa", "Chronic infections"),
  summarise_cat(baseline, "ever_chronic_bcepacia", "Ever chronic B. cepacia", "Chronic infections"),
  summarise_cat(baseline, "ever_chronic_mssa", "Ever chronic MSSA", "Chronic infections"),
  summarise_cat(baseline, "ever_chronic_mrsa", "Ever chronic MRSA", "Chronic infections"),
)

# Ensure clean ordering within domains
table1 <- table1 %>%
  mutate(
    domain = factor(domain,
                    levels = c("Demographics", "Disease severity", "Nutrition",
                               "Exacerbations", "Metabolic", "Service",
                               "Genotype", "Chronic infections"))
  ) %>%
  arrange(domain, variable, level) %>%
  group_by(domain, variable) %>%
  mutate(variable = if_else(dplyr::row_number() > 1, "", variable)) %>%
  ungroup()

# --------------------------------------------------------------------
# 5. Render with gt and Export ---------------------------------------
# --------------------------------------------------------------------

gt_tbl <- table1 %>%
  select(domain, variable, level, Overall) %>%
  gt(groupname_col = "domain") %>%
  cols_label(
    variable = "Characteristic",
    level = "Category",
    Overall = "Overall"
  ) %>%
  tab_header(
    title = md("**Baseline Characteristics (Table 1)**"),
    subtitle = "Overall cohort"
  ) %>%
  cols_align(
    align = "left",
    columns = c(variable, level)
  ) %>%
  cols_align(
    align = "center",
    columns = c(Overall)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = variable, rows = variable != "")
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
  )

html_path <- file.path(out_dir, "Table1_Baseline_Characteristics.html")
csv_path  <- file.path(out_dir, "Table1_Baseline_Characteristics.csv")

gtsave(gt_tbl, html_path)
readr::write_csv(table1, csv_path)

message("Baseline Table 1 generated:")
message("  -> HTML: ", html_path)
message("  -> CSV:  ", csv_path)
