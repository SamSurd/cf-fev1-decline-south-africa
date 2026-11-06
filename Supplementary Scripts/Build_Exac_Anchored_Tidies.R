#!/usr/bin/env Rscript
library(tidyverse)
library(arrow)
library(lmerTest)
library(broom.mixed)

out_dir <- "outputs/tables"
fs::dir_create(out_dir)

data_file <- "outputs/data_derived/analysis_long.parquet"
if (!file.exists(data_file)) stop("Missing data file: ", data_file)

dat <- read_parquet(data_file) %>%
  as_tibble() %>%
  select(patient_id, visit_date, fev1_pp, baseline_date, exac_cat) %>%
  mutate(
    visit_date = as.Date(visit_date),
    baseline_date = as.Date(baseline_date)
  ) %>%
  drop_na(patient_id, visit_date, fev1_pp, baseline_date, exac_cat)

# Helper to fit event-anchored model for a given exac category label ("1" or "2+")
fit_exac_model <- function(cat_value, label_suffix) {
  event_df <- dat %>%
    mutate(has_event = if_else(exac_cat == cat_value, 1L, 0L, missing = 0L)) %>%
    group_by(patient_id) %>%
    mutate(event_date = suppressWarnings(min(visit_date[has_event == 1L], na.rm = TRUE))) %>%
    mutate(event_date = ifelse(is.infinite(event_date), NA, event_date)) %>%
    ungroup() %>%
    mutate(
      anchor_date = coalesce(as.Date(event_date, origin = "1970-01-01"), baseline_date),
      time_exac_index_cs = as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
      post_index = if_else(!is.na(event_date) & visit_date >= event_date, 1L, 0L, missing = 0L)
    )

  # Require at least some variability in time to fit a slope
  if (n_distinct(event_df$time_exac_index_cs) < 3) return(NULL)

  mod <- try(
    lmer(fev1_pp ~ time_exac_index_cs * post_index + (1 | patient_id),
         data = event_df,
         REML = TRUE),
    silent = TRUE
  )
  if (inherits(mod, "try-error")) return(NULL)

  tidy_mod <- tidy(mod, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      # Ensure p.value exists (lmerTest should supply; else fill NA)
      p.value = if (!"p.value" %in% names(.)) NA_real_ else p.value
    ) %>%
    mutate(
      term = gsub(
        "time_exac_index_cs:post_index",
        paste0("time_exac_index_cs:exac", label_suffix, "_post_index"),
        term,
        fixed = TRUE
      )
    )
  out_path <- file.path(out_dir, paste0(".model_exac", label_suffix, "_anchored_tidy.csv"))
  readr::write_csv(tidy_mod, out_path)
  message("Saved: ", out_path)
  tidy_mod
}

# Exacerbations: 1/year vs 0
fit_exac_model("1", "1")

# Exacerbations: >=2/year vs 0
fit_exac_model("2+", "2plus")
