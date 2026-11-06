#!/usr/bin/env Rscript

# ==============================================================================
# AGENT 03.1 - SPLIT TABLE 2 (2A + 2B) + SUPPLEMENTARY S1
# ------------------------------------------------------------------------------
# Creates:
#   1) Table 2A - Baseline FEV1pp determinants
#   2) Table 2B - FEV1 decline rate determinants
#   3) Supplementary Table S1 - Non-significant baseline slope modifiers
#
# Data source priority:
#   - In-memory model objects (if running interactively after A02)
#   - Saved RDS models in outputs/models
#   - Tidy CSV outputs in outputs/tables
#
# Exports:
#   - outputs/tables/Table2A_Baseline_FEV1pp_Determinants.docx
#   - outputs/tables/Table2B_FEV1_Decline_Rate_Determinants.docx
#   - outputs/tables/TableS1_Nonsignificant_Baseline_Slope_Modifiers.docx
#   - Matching CSV files for each table
# ==============================================================================

required_pkgs <- c("dplyr", "purrr", "readr", "stringr", "tibble", "gt", "broom.mixed")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs)) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    ". Install with install.packages() and rerun."
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

out_dir <- "outputs/tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "-",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

normalize_tidy <- function(df) {
  if ("p_value" %in% names(df) && !("p.value" %in% names(df))) {
    df <- df %>% mutate(p.value = .data$p_value)
  }
  if ("conf_low" %in% names(df) && !("conf.low" %in% names(df))) {
    df <- df %>% mutate(conf.low = .data$conf_low)
  }
  if ("conf_high" %in% names(df) && !("conf.high" %in% names(df))) {
    df <- df %>% mutate(conf.high = .data$conf_high)
  }
  if ("effect" %in% names(df)) {
    df <- df %>% filter(.data$effect == "fixed")
  }
  df
}

tidy_from_object <- function(obj_name) {
  if (!exists(obj_name, envir = .GlobalEnv, inherits = FALSE)) return(NULL)
  obj <- get(obj_name, envir = .GlobalEnv, inherits = FALSE)
  if (!inherits(obj, c("lmerMod", "lmerModLmerTest", "lme"))) return(NULL)
  normalize_tidy(broom.mixed::tidy(obj, effects = "fixed", conf.int = TRUE))
}

tidy_from_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!inherits(obj, c("lmerMod", "lmerModLmerTest", "lme"))) return(NULL)
  normalize_tidy(broom.mixed::tidy(obj, effects = "fixed", conf.int = TRUE))
}

tidy_from_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  normalize_tidy(readr::read_csv(path, show_col_types = FALSE))
}

load_main_tidy <- function() {
  candidates_obj <- c("m2_reml", ".model_m2_reml", "mod", "m2")
  for (nm in candidates_obj) {
    out <- tidy_from_object(nm)
    if (!is.null(out)) {
      message("Using main model from workspace object: ", nm)
      return(out)
    }
  }

  out <- tidy_from_rds("outputs/models/.model_m2_reml.rds")
  if (!is.null(out)) {
    message("Using main model from RDS: outputs/models/.model_m2_reml.rds")
    return(out)
  }

  out <- tidy_from_csv("outputs/tables/.model_m2_reml_tidy.csv")
  if (!is.null(out)) {
    message("Using main model from CSV: outputs/tables/.model_m2_reml_tidy.csv")
    return(out)
  }

  stop("Could not load main mixed model results from workspace, RDS, or tidy CSV.")
}

load_event_tidy <- function(obj_name = NULL, rds_path = NULL, csv_path = NULL) {
  if (!is.null(obj_name)) {
    out <- tidy_from_object(obj_name)
    if (!is.null(out)) return(out)
  }
  if (!is.null(rds_path)) {
    out <- tidy_from_rds(rds_path)
    if (!is.null(out)) return(out)
  }
  if (!is.null(csv_path)) {
    out <- tidy_from_csv(csv_path)
    if (!is.null(out)) return(out)
  }
  NULL
}

extract_term <- function(df, term_key, section, variable, scale = 1) {
  row <- df %>% filter(.data$term == .env$term_key) %>% slice(1)
  if (nrow(row) == 0) {
    warning("Term not found: ", term_key)
    return(tibble(
      section = section,
      variable = variable,
      term = term_key,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_
    ))
  }

  tibble(
    section = section,
    variable = variable,
    term = term_key,
    estimate = row$estimate * scale,
    conf.low = row$conf.low * scale,
    conf.high = row$conf.high * scale,
    p.value = row$p.value
  )
}

add_display_cols <- function(df) {
  df %>%
    mutate(
      sig = is.finite(.data$p.value) & .data$p.value < 0.05,
      est_fmt = ifelse(is.finite(.data$estimate), sprintf("%+.2f", .data$estimate), "-"),
      ci_fmt = ifelse(
        is.finite(.data$conf.low) & is.finite(.data$conf.high),
        paste0(sprintf("%+.2f", .data$conf.low), " to ", sprintf("%+.2f", .data$conf.high)),
        "-"
      ),
      p_fmt = fmt_p(.data$p.value),
      est_fmt = ifelse(.data$sig, paste0("**", .data$est_fmt, "**"), .data$est_fmt),
      ci_fmt = ifelse(.data$sig, paste0("**", .data$ci_fmt, "**"), .data$ci_fmt),
      p_fmt = ifelse(.data$sig, paste0("**", .data$p_fmt, "**"), .data$p_fmt)
    )
}

style_table <- function(gt_tbl) {
  gt_tbl %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_row_groups()
    ) %>%
    gt::tab_options(
      table.border.top.color = "black",
      table.border.bottom.color = "black",
      heading.border.bottom.color = "black",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "black",
      table_body.border.bottom.color = "black"
    )
}

main_tidy <- load_main_tidy()

# Table 2A -----------------------------------------------------------
table2a_rows <- bind_rows(
  extract_term(main_tidy, "baseline_age_years_c", "DEMOGRAPHICS", "Age, years (per 10-year increase)", scale = 10),
  extract_term(main_tidy, "sexMale", "DEMOGRAPHICS", "Male sex (vs female)"),
  extract_term(main_tidy, "ethnicity_caucasian", "DEMOGRAPHICS", "Caucasian ethnicity (vs non-Caucasian)"),

  extract_term(main_tidy, "sector_public_only", "HEALTHCARE ACCESS", "Public only vs Private/Combined"),

  extract_term(main_tidy, "nutrition_catUnderweight", "NUTRITIONAL STATUS", "Underweight (vs normal)"),
  extract_term(main_tidy, "nutrition_catOverweight", "NUTRITIONAL STATUS", "Overweight (vs normal)"),

  extract_term(main_tidy, "cfrd_indicator", "METABOLIC COMPLICATIONS", "CFRD"),
  extract_term(main_tidy, "pancreatic_insufficiency", "METABOLIC COMPLICATIONS", "Pancreatic insufficiency"),

  extract_term(main_tidy, "baseline_pa_any", "CHRONIC INFECTIONS - BASELINE", "Baseline P. aeruginosa"),
  extract_term(main_tidy, "chronic_mrsa", "CHRONIC INFECTIONS - BASELINE", "Chronic MRSA"),
  extract_term(main_tidy, "chronic_mssa", "CHRONIC INFECTIONS - BASELINE", "Chronic MSSA"),
  extract_term(main_tidy, "chronic_bcepacia", "CHRONIC INFECTIONS - BASELINE", "Chronic B. cepacia"),

  extract_term(main_tidy, "exac_cat1", "EXACERBATIONS", "Exacerbations: 1/year (vs 0)"),
  extract_term(main_tidy, "exac_cat2+", "EXACERBATIONS", "Exacerbations: >=2/year (vs 0)"),

  extract_term(main_tidy, "f508_any", "GENOTYPE", "F508del (any vs none)"),

  extract_term(main_tidy, "cftr_modulator_indicator", "CFTR MODULATOR THERAPY", "Index FEV1pp on CFTR")
) %>%
  mutate(
    section = factor(
      .data$section,
      levels = c(
        "DEMOGRAPHICS",
        "HEALTHCARE ACCESS",
        "NUTRITIONAL STATUS",
        "METABOLIC COMPLICATIONS",
        "CHRONIC INFECTIONS - BASELINE",
        "EXACERBATIONS",
        "GENOTYPE",
        "CFTR MODULATOR THERAPY"
      )
    )
  ) %>%
  add_display_cols()

key_rows_2a <- c(
  "Public only vs Private/Combined",
  "Underweight (vs normal)",
  "CFRD",
  "Index FEV1pp on CFTR"
)

table2a_gt <- table2a_rows %>%
  select(section, variable, est_fmt, ci_fmt, p_fmt) %>%
  gt::gt(groupname_col = "section") %>%
  gt::cols_label(
    variable = "Variable",
    est_fmt = "Difference (%)",
    ci_fmt = "95% CI",
    p_fmt = "p-value"
  ) %>%
  gt::tab_header(
    title = gt::md("**Table 2A. Baseline FEV1pp Determinants**"),
    subtitle = "Outcome: FEV1 % Predicted (pp)"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(variable)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(est_fmt, ci_fmt)
  ) %>%
  gt::cols_align(align = "right", columns = c(p_fmt)) %>%
  gt::fmt_markdown(columns = c(est_fmt, ci_fmt, p_fmt)) %>%
  gt::cols_width(
    variable ~ gt::px(360),
    est_fmt ~ gt::px(170),
    ci_fmt ~ gt::px(240),
    p_fmt ~ gt::px(120)
  ) %>%
  style_table() %>%
  gt::tab_source_note(
    source_note = "Bold indicates p < 0.05."
  ) %>%
  gt::tab_source_note(
    source_note = paste0(
      "Key findings: public sector (",
      sprintf("%.2f", table2a_rows$estimate[table2a_rows$variable == key_rows_2a[1]]), "%), ",
      "underweight (",
      sprintf("%.2f", table2a_rows$estimate[table2a_rows$variable == key_rows_2a[2]]), "%), ",
      "CFRD (",
      sprintf("%.2f", table2a_rows$estimate[table2a_rows$variable == key_rows_2a[3]]), "%), ",
      "CFTR modulator (",
      sprintf("%.2f", table2a_rows$estimate[table2a_rows$variable == key_rows_2a[4]]), "%)."
    )
  )

table2a_csv <- file.path(out_dir, "Table2A_Baseline_FEV1pp_Determinants.csv")
table2a_docx <- file.path(out_dir, "Table2A_Baseline_FEV1pp_Determinants.docx")
table2a_html <- file.path(out_dir, "Table2A_Baseline_FEV1pp_Determinants.html")
readr::write_csv(table2a_rows, table2a_csv)
gt::gtsave(table2a_gt, table2a_html)
gt::gtsave(table2a_gt, table2a_docx)

# Table 2B -----------------------------------------------------------
cftr_tidy <- load_event_tidy(
  obj_name = "m2_cftr",
  rds_path = "outputs/models/.model_m2_cftr.rds",
  csv_path = "outputs/tables/.model_m2_cftr_tidy.csv"
)
pa_tidy <- load_event_tidy(csv_path = "outputs/tables/.model_pa_anchored_tidy.csv")
mssa_tidy <- load_event_tidy(csv_path = "outputs/tables/.model_mssa_anchored_tidy.csv")
exac1_tidy <- load_event_tidy(csv_path = "outputs/tables/.model_exac1_anchored_tidy.csv")
exac2_tidy <- load_event_tidy(csv_path = "outputs/tables/.model_exac2plus_anchored_tidy.csv")

table2b_rows <- bind_rows(
  extract_term(main_tidy, "time_cs", "GLOBAL ANNUAL DECLINE", "Global annual slope"),

  extract_term(cftr_tidy %||% tibble(term = character()), "time_cftr_index_cs:cftr_post_index", "EVENT-ANCHORED TRAJECTORY CHANGES", "CFTR modulator initiation"),
  extract_term(pa_tidy %||% tibble(term = character()), "time_pa_index_cs:pa_post_index", "EVENT-ANCHORED TRAJECTORY CHANGES", "Chronic P. aeruginosa acquisition"),
  extract_term(mssa_tidy %||% tibble(term = character()), "time_mssa_index_cs:mssa_post_index", "EVENT-ANCHORED TRAJECTORY CHANGES", "Chronic MSSA acquisition"),
  extract_term(exac1_tidy %||% tibble(term = character()), "time_exac_index_cs:exac1_post_index", "EVENT-ANCHORED TRAJECTORY CHANGES", "Exacerbations: 1/year"),
  extract_term(exac2_tidy %||% tibble(term = character()), "time_exac_index_cs:exac2plus_post_index", "EVENT-ANCHORED TRAJECTORY CHANGES", "Exacerbations: >=2/year"),

  extract_term(main_tidy, "time_cs:sector_public_only", "BASELINE SLOPE MODIFIERS", "Public only vs Private/Combined")
) %>%
  mutate(
    section = factor(
      .data$section,
      levels = c(
        "GLOBAL ANNUAL DECLINE",
        "EVENT-ANCHORED TRAJECTORY CHANGES",
        "BASELINE SLOPE MODIFIERS"
      )
    )
  ) %>%
  add_display_cols()

table2b_gt <- table2b_rows %>%
  select(section, variable, est_fmt, ci_fmt, p_fmt) %>%
  gt::gt(groupname_col = "section") %>%
  gt::cols_label(
    variable = "Variable",
    est_fmt = "Annual Rate of Change (%/yr)",
    ci_fmt = "95% CI",
    p_fmt = "p-value"
  ) %>%
  gt::tab_header(
    title = gt::md("**Table 2B. FEV1 Decline Rate Determinants**"),
    subtitle = "Outcome: FEV1 % Predicted (pp)"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(variable)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(est_fmt, ci_fmt)
  ) %>%
  gt::cols_align(align = "right", columns = c(p_fmt)) %>%
  gt::fmt_markdown(columns = c(est_fmt, ci_fmt, p_fmt)) %>%
  gt::cols_width(
    variable ~ gt::px(380),
    est_fmt ~ gt::px(220),
    ci_fmt ~ gt::px(220),
    p_fmt ~ gt::px(120)
  ) %>%
  style_table() %>%
  gt::tab_source_note(
    source_note = "Event-anchored effects represent the change in slope. Net slope = global slope + event effect."
  ) %>%
  gt::tab_source_note(
    source_note = "Bold indicates p < 0.05."
  )

table2b_csv <- file.path(out_dir, "Table2B_FEV1_Decline_Rate_Determinants.csv")
table2b_docx <- file.path(out_dir, "Table2B_FEV1_Decline_Rate_Determinants.docx")
table2b_html <- file.path(out_dir, "Table2B_FEV1_Decline_Rate_Determinants.html")
readr::write_csv(table2b_rows, table2b_csv)
gt::gtsave(table2b_gt, table2b_html)
gt::gtsave(table2b_gt, table2b_docx)

# Supplementary Table S1 ------------------------------------------------------
supp_s1_rows <- bind_rows(
  extract_term(main_tidy, "time_cs:baseline_pa_any", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "Baseline P. aeruginosa"),
  extract_term(main_tidy, "time_cs:cfrd_indicator_slope", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "CFRD"),
  extract_term(main_tidy, "time_cs:ethnicity_caucasian", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "Caucasian ethnicity"),
  extract_term(main_tidy, "time_cs:f508_any", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "F508del"),
  extract_term(main_tidy, "time_cs:nutrition_catOverweight", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "Overweight"),
  extract_term(main_tidy, "time_cs:nutrition_catUnderweight", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "Underweight"),
  extract_term(main_tidy, "time_cs:pancreatic_insufficiency_slope", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "Pancreatic insufficiency"),
  extract_term(main_tidy, "time_cs:sexMale", "BASELINE SLOPE MODIFIERS (NON-SIGNIFICANT)", "Male sex")
) %>%
  add_display_cols() %>%
  mutate(non_sig_expected = is.na(.data$p.value) | .data$p.value > 0.05)

if (any(!supp_s1_rows$non_sig_expected, na.rm = TRUE)) {
  warning("One or more S1 rows are statistically significant at p<0.05 in current model output.")
}

supp_s1_gt <- supp_s1_rows %>%
  select(variable, est_fmt, ci_fmt, p_fmt) %>%
  gt::gt() %>%
  gt::cols_label(
    variable = "Variable",
    est_fmt = "Annual Rate Modifier (%/yr)",
    ci_fmt = "95% CI",
    p_fmt = "p-value"
  ) %>%
  gt::tab_header(
    title = gt::md("**Supplementary Table S1. Non-significant Baseline Slope Modifiers**"),
    subtitle = "Outcome: FEV1 % Predicted (pp)"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(variable)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(est_fmt, ci_fmt)
  ) %>%
  gt::cols_align(align = "right", columns = c(p_fmt)) %>%
  gt::fmt_markdown(columns = c(est_fmt, ci_fmt, p_fmt)) %>%
  gt::cols_width(
    variable ~ gt::px(360),
    est_fmt ~ gt::px(220),
    ci_fmt ~ gt::px(220),
    p_fmt ~ gt::px(120)
  ) %>%
  style_table() %>%
  gt::tab_source_note(
    source_note = "Rows listed were prespecified baseline slope modifiers with p > 0.05 in the primary model."
  )

supp_s1_csv <- file.path(out_dir, "TableS1_Nonsignificant_Baseline_Slope_Modifiers.csv")
supp_s1_docx <- file.path(out_dir, "TableS1_Nonsignificant_Baseline_Slope_Modifiers.docx")
supp_s1_html <- file.path(out_dir, "TableS1_Nonsignificant_Baseline_Slope_Modifiers.html")
readr::write_csv(supp_s1_rows, supp_s1_csv)
gt::gtsave(supp_s1_gt, supp_s1_html)
gt::gtsave(supp_s1_gt, supp_s1_docx)

message("Split Table 2 outputs generated:")
message("  - ", table2a_docx)
message("  - ", table2b_docx)
message("  - ", supp_s1_docx)
