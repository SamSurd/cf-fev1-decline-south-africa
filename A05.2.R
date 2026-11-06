#!/usr/bin/env Rscript

# ==============================================================================
# AGENT 05.2 - TABLE 1 (PUBLIC ONLY VS PRIVATE/COMBINED)
# ------------------------------------------------------------------------------
# Builds a baseline characteristics table stratified by healthcare sector:
#   - Public Only
#   - Private/Combined
#   - Overall
# Includes p-values for Public Only vs Private/Combined:
#   - t-test for continuous variables
#   - chi-square test (or Fisher's exact when expected counts are small)
# Exports:
#   - outputs/tables/Table1_Baseline_Public_vs_PrivateCombined.docx
#   - outputs/tables/Table1_Baseline_Public_vs_PrivateCombined.csv
# ==============================================================================

required_pkgs <- c(
  "dplyr", "tibble", "purrr", "stringr", "readr", "gt"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs)) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them with install.packages() and rerun."
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

pick_var <- function(df, candidates) {
  nms <- names(df)

  # Exact match first
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[[1]])

  # Case-insensitive fallback
  idx <- match(tolower(candidates), tolower(nms), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx) > 0) return(nms[idx[[1]]])

  NA_character_
}

as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

as_factor_clean <- function(x) {
  out <- trimws(as.character(x))
  out[out == ""] <- NA_character_
  factor(out)
}

extract_col <- function(df, var_name, what = "column") {
  out <- tryCatch(df[[var_name]], error = function(e) NULL)
  if (!is.null(out) && length(out) == nrow(df)) return(out)

  nms <- names(df)
  idx <- which(tolower(trimws(nms)) == tolower(trimws(var_name)))
  if (length(idx) >= 1) {
    out2 <- tryCatch(df[[idx[[1]]]], error = function(e) NULL)
    if (!is.null(out2) && length(out2) == nrow(df)) return(out2)
  }

  stop(
    "Could not extract ", what, " '", var_name, "'. ",
    "This usually means the object is not a plain data frame or the column is malformed."
  )
}

as_yes_no <- function(x, yes_values = c("yes", "1", "true")) {
  if (is.factor(x)) x <- as.character(x)

  if (is.logical(x)) {
    out <- ifelse(is.na(x), NA_character_, ifelse(x, "Yes", "No"))
    return(factor(out, levels = c("No", "Yes")))
  }

  if (is.numeric(x)) {
    out <- ifelse(is.na(x), NA_character_, ifelse(x == 1, "Yes", "No"))
    return(factor(out, levels = c("No", "Yes")))
  }

  x_chr <- tolower(trimws(as.character(x)))
  x_chr[x_chr == ""] <- NA_character_

  yes <- x_chr %in% tolower(yes_values)
  no <- x_chr %in% c("no", "0", "false", "negative", "normal", "sufficient")

  out <- ifelse(yes, "Yes", ifelse(no, "No", NA_character_))
  factor(out, levels = c("No", "Yes"))
}

as_sex <- function(x) {
  if (is.numeric(x)) {
    out <- ifelse(is.na(x), NA_character_, ifelse(x == 1, "Male", "Female"))
    return(factor(out))
  }

  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  low <- tolower(x_chr)

  out <- dplyr::case_when(
    low %in% c("m", "male", "1") ~ "Male",
    low %in% c("f", "female", "0") ~ "Female",
    TRUE ~ x_chr
  )
  factor(out)
}

as_f508 <- function(x, source_var) {
  if (is.numeric(x) || source_var == "f508_any") {
    out <- ifelse(is.na(x), NA_character_, ifelse(as_numeric_safe(x) == 1, "Any F508del", "No F508del"))
    return(factor(out, levels = c("No F508del", "Any F508del")))
  }
  as_factor_clean(x)
}

fmt_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("-")
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
}

fmt_n_pct <- function(n, d) {
  if (!is.finite(d) || d <= 0) return("-")
  sprintf("%d (%.1f%%)", n, 100 * n / d)
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("-")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

p_cont <- function(x, grp) {
  keep <- !is.na(x) & !is.na(grp)
  x <- x[keep]
  grp <- droplevels(factor(grp[keep]))
  if (length(unique(grp)) != 2) return(NA_real_)
  split_x <- split(x, grp)
  if (any(vapply(split_x, length, integer(1)) < 2)) return(NA_real_)
  tryCatch(stats::t.test(x ~ grp)$p.value, error = function(e) NA_real_)
}

p_cat <- function(x, grp) {
  keep <- !is.na(x) & !is.na(grp)
  x <- x[keep]
  grp <- droplevels(factor(grp[keep]))
  if (length(unique(grp)) != 2 || length(unique(x)) < 2) return(NA_real_)

  tab <- table(x, grp)
  if (any(dim(tab) < 2)) return(NA_real_)

  tryCatch({
    chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    if (any(chi$expected < 5)) stats::fisher.test(tab)$p.value else chi$p.value
  }, error = function(e) {
    tryCatch(stats::fisher.test(tab)$p.value, error = function(e2) NA_real_)
  })
}

cont_row <- function(df, var, label, domain) {
  x <- df[[var]]
  tibble(
    domain = domain,
    variable = label,
    level = "Mean (SD)",
    Public_Only = fmt_mean_sd(x[df$sector_group == "Public Only"]),
    Private_Combined = fmt_mean_sd(x[df$sector_group == "Private/Combined"]),
    Overall = fmt_mean_sd(x),
    p_value = fmt_p(p_cont(x, df$sector_group))
  )
}

binary_row <- function(df, var, label, domain, positive_level = "Yes") {
  x <- factor(df[[var]])
  lvls <- levels(x)
  lvl_match <- lvls[tolower(lvls) == tolower(positive_level)]
  pos <- if (length(lvl_match)) lvl_match[[1]] else positive_level

  tibble(
    domain = domain,
    variable = label,
    level = as.character(pos),
    Public_Only = fmt_n_pct(
      sum(x[df$sector_group == "Public Only"] == pos, na.rm = TRUE),
      sum(!is.na(x[df$sector_group == "Public Only"]))
    ),
    Private_Combined = fmt_n_pct(
      sum(x[df$sector_group == "Private/Combined"] == pos, na.rm = TRUE),
      sum(!is.na(x[df$sector_group == "Private/Combined"]))
    ),
    Overall = fmt_n_pct(
      sum(x == pos, na.rm = TRUE),
      sum(!is.na(x))
    ),
    p_value = fmt_p(p_cat(x, df$sector_group))
  )
}

multicat_rows <- function(df, var, label, domain) {
  x <- factor(df[[var]])
  lvls <- levels(droplevels(x))
  lvls <- lvls[vapply(lvls, function(lvl) any(x == lvl, na.rm = TRUE), logical(1))]

  p_fmt <- fmt_p(p_cat(x, df$sector_group))
  body <- purrr::map_dfr(seq_along(lvls), function(i) {
    lvl <- lvls[[i]]
    tibble(
      domain = domain,
      variable = ifelse(i == 1, label, ""),
      level = as.character(lvl),
      Public_Only = fmt_n_pct(
        sum(x[df$sector_group == "Public Only"] == lvl, na.rm = TRUE),
        sum(!is.na(x[df$sector_group == "Public Only"]))
      ),
      Private_Combined = fmt_n_pct(
        sum(x[df$sector_group == "Private/Combined"] == lvl, na.rm = TRUE),
        sum(!is.na(x[df$sector_group == "Private/Combined"]))
      ),
      Overall = fmt_n_pct(
        sum(x == lvl, na.rm = TRUE),
        sum(!is.na(x))
      ),
      p_value = ifelse(i == 1, p_fmt, "")
    )
  })
  body
}

get_input_data <- function() {
  candidates <- c("cf_data", "sacfr_data", "dat", "analysis_long", "dat_model", "df_model")
  hits <- candidates[candidates %in% ls(envir = .GlobalEnv)]

  if (length(hits) > 0) {
    obj_name <- hits[[1]]
    obj <- get(obj_name, envir = .GlobalEnv)
    if (!inherits(obj, c("data.frame", "tbl_df", "tbl"))) {
      stop("Detected object '", obj_name, "' but it is not a data frame.")
    }
    return(list(data = tibble::as_tibble(obj), source = paste0("workspace object: ", obj_name)))
  }

  parquet_path <- "outputs/data_derived/analysis_long.parquet"
  if (file.exists(parquet_path)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("No workspace data object found and package 'arrow' is unavailable for parquet fallback.")
    }
    obj <- arrow::read_parquet(parquet_path) %>% tibble::as_tibble()
    return(list(data = obj, source = parquet_path))
  }

  stop(
    "Could not find an analysis data object. Checked workspace names: ",
    paste(candidates, collapse = ", "),
    "; and fallback file: ", parquet_path
  )
}

input <- get_input_data()
dat <- input$data %>% tibble::as_tibble() %>% dplyr::ungroup()
message("Using dataset from ", input$source)

id_var <- pick_var(dat, c("patient_id", "id", "patientid"))
visit_date_var <- pick_var(dat, c("visit_date", "encounter_date", "date"))

if (is.na(id_var)) {
  stop("Could not find a patient identifier column (e.g., patient_id).")
}

if (!(id_var %in% names(dat))) {
  stop(
    "Detected id column '", id_var, "' is not in data. Available columns include: ",
    paste(head(names(dat), 30), collapse = ", "),
    if (ncol(dat) > 30) ", ..."
  )
}

if (!is.na(visit_date_var) && !(visit_date_var %in% names(dat))) {
  warning("Detected visit date column not found in data. Proceeding without visit-date ordering.")
  visit_date_var <- NA_character_
}

baseline <- as.data.frame(dat, stringsAsFactors = FALSE)
id_vec <- extract_col(baseline, id_var, what = "ID column")

if (!is.na(visit_date_var)) {
  visit_vec <- extract_col(baseline, visit_date_var, what = "visit date column")
  ord <- order(
    id_vec,
    visit_vec,
    na.last = TRUE
  )
  baseline <- baseline[ord, , drop = FALSE]
  id_vec <- id_vec[ord]
}

baseline <- baseline[!duplicated(id_vec), , drop = FALSE]
baseline <- tibble::as_tibble(baseline)

if (!("ever_cfrd" %in% names(baseline)) && ("cfrd" %in% names(dat))) {
  ever_cfrd <- dat %>%
    group_by(.data[[id_var]]) %>%
    summarise(
      ever_cfrd = dplyr::case_when(
        any(tolower(as.character(cfrd)) == "yes", na.rm = TRUE) ~ "Yes",
        any(tolower(as.character(cfrd)) == "no", na.rm = TRUE) ~ "No",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    )
  names(ever_cfrd)[1] <- id_var
  baseline <- baseline %>% left_join(ever_cfrd, by = id_var)
}

var_map <- list(
  sector = c("sector", "baseline_sector", "healthcare_sector", "health_sector"),
  sector_public_only = c("sector_public_only"),
  age = c("age_years", "baseline_age_years", "age"),
  sex = c("sex", "gender"),
  ethnicity = c("ethnicity_cat", "ethnicity", "race_ethnicity"),
  fev1pp = c("baseline_fev1_pp", "fev1_pp", "fev1pp"),
  fev1_stratum = c("baseline_fev1_cat", "fev1_stratum", "fev1_cat"),
  nutrition = c("baseline_nutrition_cat", "nutrition_cat", "nutritional_status"),
  cfrd = c("ever_cfrd", "cfrd_indicator", "baseline_cfrd_indicator", "cfrd"),
  pancreatic = c("pancreatic_insufficiency", "baseline_pancreatic_insuff", "pancreatic_status"),
  pa = c("ever_chronic_pa", "chronic_paeruginosa", "baseline_pa_any", "pa_status"),
  mssa = c("ever_chronic_mssa", "chronic_mssa", "sa_status"),
  mrsa = c("ever_chronic_mrsa", "chronic_mrsa", "mrsa_status"),
  bcep = c("ever_chronic_bcepacia", "chronic_bcepacia", "bcepacia_status"),
  f508 = c("f508_genotype_group", "f508_any", "f508_status"),
  mod = c("ever_cftr_yes", "ever_modulator", "cftr_ever", "cftr_modulator_indicator")
)

vars <- purrr::map_chr(var_map, pick_var, df = baseline)

required_vars <- c(
  "age", "sex", "ethnicity", "fev1pp", "fev1_stratum", "nutrition",
  "cfrd", "pancreatic", "pa", "mssa", "mrsa", "bcep", "f508", "mod"
)

missing_required <- required_vars[is.na(vars[required_vars])]
if (length(missing_required) > 0) {
  stop("Missing required variables for Table 1: ", paste(missing_required, collapse = ", "))
}

if (is.na(vars["sector"]) && is.na(vars["sector_public_only"])) {
  stop("Could not find sector variable. Checked: sector / healthcare_sector / sector_public_only.")
}

sector_group <- NULL
if (!is.na(vars["sector"])) {
  sec <- tolower(trimws(as.character(baseline[[vars["sector"]]])))
  sector_group <- dplyr::case_when(
    is.na(sec) ~ NA_character_,
    stringr::str_detect(sec, "public") ~ "Public Only",
    stringr::str_detect(sec, "private|combined") ~ "Private/Combined",
    TRUE ~ NA_character_
  )
} else {
  sec_bin <- as_numeric_safe(baseline[[vars["sector_public_only"]]])
  sector_group <- ifelse(is.na(sec_bin), NA_character_, ifelse(sec_bin == 1, "Public Only", "Private/Combined"))
}

tbl_dat <- tibble(
  sector_group = factor(sector_group, levels = c("Public Only", "Private/Combined")),
  age = as_numeric_safe(baseline[[vars["age"]]]),
  sex = as_sex(baseline[[vars["sex"]]]),
  ethnicity = as_factor_clean(baseline[[vars["ethnicity"]]]),
  fev1pp = as_numeric_safe(baseline[[vars["fev1pp"]]]),
  fev1_stratum = as_factor_clean(baseline[[vars["fev1_stratum"]]]),
  nutrition = as_factor_clean(baseline[[vars["nutrition"]]]),
  cfrd = as_yes_no(baseline[[vars["cfrd"]]], yes_values = c("yes", "1", "true")),
  pancreatic = as_yes_no(baseline[[vars["pancreatic"]]], yes_values = c("yes", "1", "true", "insufficient")),
  pa = as_yes_no(baseline[[vars["pa"]]], yes_values = c("yes", "1", "true", "chronic", "intermittent", "positive")),
  mssa = as_yes_no(baseline[[vars["mssa"]]], yes_values = c("yes", "1", "true", "chronic", "positive")),
  mrsa = as_yes_no(baseline[[vars["mrsa"]]], yes_values = c("yes", "1", "true", "chronic", "positive")),
  bcep = as_yes_no(baseline[[vars["bcep"]]], yes_values = c("yes", "1", "true", "chronic", "positive")),
  f508 = as_f508(baseline[[vars["f508"]]], vars["f508"]),
  modulator = as_yes_no(baseline[[vars["mod"]]], yes_values = c("yes", "1", "true"))
) %>%
  filter(!is.na(sector_group))

n_public <- sum(tbl_dat$sector_group == "Public Only", na.rm = TRUE)
n_private_combined <- sum(tbl_dat$sector_group == "Private/Combined", na.rm = TRUE)
n_overall <- nrow(tbl_dat)

message("Sector counts after collapse: Public Only=", n_public,
        ", Private/Combined=", n_private_combined,
        ", Overall=", n_overall)

f508_levels <- levels(droplevels(tbl_dat$f508))
f508_is_multilevel <- length(f508_levels) > 2

table1 <- bind_rows(
  cont_row(tbl_dat, "age", "Age at baseline, years", "Demographics"),
  binary_row(tbl_dat, "sex", "Sex", "Demographics", positive_level = "Male"),
  multicat_rows(tbl_dat, "ethnicity", "Ethnicity", "Demographics"),

  cont_row(tbl_dat, "fev1pp", "Baseline FEV1 % predicted", "Disease severity"),
  multicat_rows(tbl_dat, "fev1_stratum", "FEV1 stratum", "Disease severity"),

  multicat_rows(tbl_dat, "nutrition", "Nutritional status", "Nutrition"),

  binary_row(tbl_dat, "cfrd", "CFRD", "Metabolic", positive_level = "Yes"),
  binary_row(tbl_dat, "pancreatic", "Pancreatic insufficiency", "Metabolic", positive_level = "Yes"),

  binary_row(tbl_dat, "pa", "P. aeruginosa", "Chronic infections", positive_level = "Yes"),
  binary_row(tbl_dat, "mssa", "MSSA", "Chronic infections", positive_level = "Yes"),
  binary_row(tbl_dat, "mrsa", "MRSA", "Chronic infections", positive_level = "Yes"),
  binary_row(tbl_dat, "bcep", "B. cepacia", "Chronic infections", positive_level = "Yes"),

  if (f508_is_multilevel) multicat_rows(tbl_dat, "f508", "F508del status", "Genotype")
  else binary_row(tbl_dat, "f508", "F508del status", "Genotype", positive_level = "Any F508del"),

  binary_row(tbl_dat, "modulator", "Ever on modulator", "CFTR modulator therapy", positive_level = "Yes")
)

table1 <- table1 %>%
  mutate(
    domain = factor(
      domain,
      levels = c(
        "Demographics",
        "Disease severity",
        "Nutrition",
        "Metabolic",
        "Chronic infections",
        "Genotype",
        "CFTR modulator therapy"
      )
    )
  )

out_dir <- "outputs/tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

csv_path <- file.path(out_dir, "Table1_Baseline_Public_vs_PrivateCombined.csv")
html_path <- file.path(out_dir, "Table1_Baseline_Public_vs_PrivateCombined.html")
docx_path <- file.path(out_dir, "Table1_Baseline_Public_vs_PrivateCombined.docx")

readr::write_csv(table1, csv_path)

gt_tbl <- table1 %>%
  gt::gt(groupname_col = "domain") %>%
  gt::cols_label(
    variable = "Characteristic",
    level = "Category",
    Public_Only = paste0("Public Only (n=", n_public, ")"),
    Private_Combined = paste0("Private/Combined (n=", n_private_combined, ")"),
    Overall = paste0("Overall (n=", n_overall, ")"),
    p_value = "p-value"
  ) %>%
  gt::tab_header(
    title = gt::md("**Baseline Characteristics (Table 1)**"),
    subtitle = "Stratified by healthcare sector"
  ) %>%
  gt::cols_align(
    align = "left",
    columns = c(variable, level)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(Public_Only, Private_Combined, Overall, p_value)
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_body(columns = variable, rows = variable != "")
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
  ) %>%
  gt::tab_source_note(
    source_note = "Values are mean (SD) or n (%). p-values compare Public Only vs Private/Combined using t-test (continuous) and chi-square/Fisher's exact (categorical)."
  )

gt::gtsave(gt_tbl, filename = html_path)
gt::gtsave(gt_tbl, filename = docx_path)

message("Table 1 exported:")
message("  - ", html_path)
message("  - ", csv_path)
message("  - ", docx_path)
