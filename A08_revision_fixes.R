#!/usr/bin/env Rscript
# ==============================================================================
# A08 - REVISION FIXES FOR JCF-D-26-00377 (Reviewer 1: model specification)
# ------------------------------------------------------------------------------
# Fix 1  CFTR modulator: segmented model with a LEVEL change and a SLOPE change
#        at initiation (A02's M2_cftr estimated a slope change only), anchored
#        at the start date, with every on-treatment visit coded as exposed.
# Fix 2  Population-averaged adjusted annual slope from M2. A02 reported the
#        time_cs coefficient, which is the slope at the reference level of
#        every covariate that interacts with time.
# Fix 3  Net post-event slopes computed WITHIN each model (pre-event slope +
#        slope change) with Satterthwaite 95% CIs, replacing
#        "global slope + event effect".
# Fix 4  Event-anchored models for chronic P. aeruginosa and MSSA restricted
#        to incident events (>=1 non-chronic registry year before onset).
#
# Run from the project folder:  Rscript A08_revision_fixes.R
# Reads : outputs/data_derived/analysis_long.parquet (built by A01.R)
#         data_raw/CF_registry_master.xlsx (Sheet1: CFTR start dates + flags)
# Writes: outputs/revision/ only. Existing outputs are not touched.
# ==============================================================================

set.seed(20251111)

suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(lme4); library(lmerTest)
  library(readxl); library(janitor)
})

data_file <- "outputs/data_derived/analysis_long.parquet"
raw_file  <- "data_raw/CF_registry_master.xlsx"
out_dir   <- "outputs/revision"
dir.create(file.path(out_dir, "models"), showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path(out_dir, "A08_summary.txt")
file.create(summary_file)
say <- function(...) {
  msg <- paste0(...)
  message(msg)
  con <- file(summary_file, open = "a"); writeLines(msg, con); close(con)
}
show_tbl <- function(x, digits = 3) {
  txt <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  for (l in txt) say(l)
}
med_iqr <- function(x, d = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  q <- quantile(x, c(0.5, 0.25, 0.75), names = FALSE)
  sprintf(paste0("%.", d, "f (IQR %.", d, "f-%.", d, "f)"), q[1], q[2], q[3])
}

# ------------------------------------------------------------------------------
# 1. Data preparation: identical to A02.R
# ------------------------------------------------------------------------------
dat <- arrow::read_parquet(data_file) %>% as_tibble()

dat <- dat %>%
  mutate(
    sex           = factor(sex, levels = c("Female", "Male")),
    nutrition_cat = factor(nutrition_cat, levels = c("Normal", "Underweight", "Overweight", "Unknown")),
    exac_cat      = factor(exac_cat, levels = c("0", "1", "2+", "Unknown")),
    sector        = factor(sector, levels = c("Public", "Private", "Combined")),
    ethnicity_cat = factor(ethnicity_cat, levels = c("Black African", "Mixed Race", "Indian", "Asian", "Caucasian", "Other", "Unknown")),
    f508_genotype_group = factor(f508_genotype_group, levels = c("Other", "Heterozygous", "Homozygous")),
    ethnicity_caucasian = dplyr::if_else(ethnicity_cat == "Caucasian", 1L, 0L, missing = 0L),
    sector_public_only  = dplyr::if_else(sector == "Public", 1L, 0L, missing = 0L),
    f508_any = dplyr::if_else(f508_genotype_group %in% c("Heterozygous", "Homozygous"), 1L, 0L, missing = 0L)
  )

dat_model <- dat %>%
  filter(nutrition_cat != "Unknown" | is.na(nutrition_cat)) %>%
  droplevels() %>%
  group_by(patient_id) %>%
  filter(dplyr::n() >= 2) %>%
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

say("A08 revision fixes, run ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Model data: ", nrow(dat_model), " rows, ", n_distinct(dat_model$patient_id), " patients, ",
    sum(!is.na(dat_model$fev1_pp)), " rows with FEV1pp")

# ------------------------------------------------------------------------------
# 2. Model builders: copied from A02.R so the original models reproduce exactly
# ------------------------------------------------------------------------------
compute_use_flags <- function(df) {
  nlev <- function(x) length(unique(na.omit(x)))
  list(
    use_blf_cat = nlev(df$baseline_fev1_cat) >= 2,
    use_sex     = nlev(df$sex) >= 2,
    use_nut     = nlev(df$nutrition_cat) >= 2,
    use_exac    = nlev(df$exac_cat) >= 2,
    use_sector  = nlev(df$sector_public_only) >= 2,
    use_eth     = nlev(df$ethnicity_caucasian) >= 2,
    use_geno    = nlev(df$f508_any) >= 2
  )
}

# Same fallback logic as A02 (correlated -> uncorrelated -> intercept-only),
# but it also records which structure was used and any fitting warnings.
fit_log <- list()
fit_triple_fallback <- function(rhs_formula, data, use_reml = FALSE, label) {
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  forms <- list(
    correlated   = as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 + time_cs | patient_id)")),
    uncorrelated = as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 + time_cs || patient_id)")),
    intercept    = as.formula(paste("fev1_pp ~", rhs_formula, "+ (1 | patient_id)"))
  )
  for (nm in names(forms)) {
    w <- character(0)
    mod <- tryCatch(
      withCallingHandlers(
        lmerTest::lmer(forms[[nm]], data = data, REML = use_reml, control = ctrl),
        warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
      ),
      error = function(e) { say("    ", label, ": ", nm, " fit failed: ", conditionMessage(e)); NULL }
    )
    if (is.null(mod)) next
    sing <- lme4::isSingular(mod, tol = 1e-4)
    if (!sing || nm == "intercept") {
      fit_log[[label]] <<- tibble(
        model = label,
        re_structure = paste0(nm, if (sing) " (singular)" else ""),
        n_patients = unname(lme4::ngrps(mod)[["patient_id"]]),
        n_obs = nobs(mod),
        n_warnings = length(w),
        first_warning = if (length(w)) w[1] else ""
      )
      return(mod)
    }
  }
  NULL
}

build_rhs <- function(flags, version = c("m1", "m2")) {
  version <- match.arg(version)
  terms <- c("time_cs", "pancreatic_insufficiency", "cfrd_indicator", "baseline_age_years_c",
             "baseline_fev1_cat", "baseline_pa_any", "sex", "ethnicity_caucasian",
             "sector_public_only", "f508_any")
  if (flags$use_nut)  terms <- c(terms, "nutrition_cat")
  if (flags$use_exac) terms <- c(terms, "exac_cat")
  terms <- c(terms, "chronic_bcepacia", "chronic_mssa", "chronic_mrsa")
  terms <- c(terms, "cftr_modulator_indicator")
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
    interactions <- c(interactions, "time_cs:cftr_modulator_indicator_slope",
                      "time_cs:chronic_paeruginosa_slope", "time_cs:chronic_mssa_slope")
  }
  paste(c(terms, interactions), collapse = " + ")
}

build_rhs_cftr <- function(flags) {
  terms <- c("time_cftr_index_cs", "pancreatic_insufficiency", "cfrd_indicator", "baseline_age_years_c",
             "baseline_fev1_cat", "sex", "ethnicity_caucasian", "sector_public_only", "f508_any")
  if (flags$use_nut)  terms <- c(terms, "nutrition_cat")
  if (flags$use_exac) terms <- c(terms, "exac_cat")
  terms <- c(terms, "chronic_paeruginosa", "chronic_bcepacia", "chronic_mssa", "chronic_mrsa")
  interactions <- c("time_cftr_index_cs:pancreatic_insufficiency_slope", "time_cftr_index_cs:cfrd_indicator_slope")
  if (flags$use_nut)     interactions <- c(interactions, "time_cftr_index_cs:nutrition_cat")
  if (flags$use_blf_cat) interactions <- c(interactions, "time_cftr_index_cs:baseline_fev1_cat")
  if (flags$use_sex)     interactions <- c(interactions, "time_cftr_index_cs:sex")
  if (flags$use_sector)  interactions <- c(interactions, "time_cftr_index_cs:sector_public_only")
  if (flags$use_eth)     interactions <- c(interactions, "time_cftr_index_cs:ethnicity_caucasian")
  if (flags$use_geno)    interactions <- c(interactions, "time_cftr_index_cs:f508_any")
  interactions <- c(interactions, "time_cftr_index_cs:chronic_paeruginosa_slope",
                    "time_cftr_index_cs:chronic_mssa_slope", "time_cftr_index_cs:cftr_post_index")
  paste(c(terms, interactions), collapse = " + ")
}

build_event_data <- function(df, spec) {
  time_var <- paste0("time_", spec$name, "_index")
  time_cs_var <- paste0(time_var, "_cs")
  post_var <- paste0(spec$name, "_post_index")
  df %>%
    mutate(
      event_start = .data[[spec$start_col]],
      anchor_date = dplyr::coalesce(event_start, baseline_date),
      !!time_var := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
      !!time_cs_var := .data[[time_var]],
      !!post_var := dplyr::if_else(
        .data[[spec$slope_ok_col]] == 1L & !is.na(event_start) & visit_date >= event_start,
        1L, 0L, missing = 0L)
    ) %>%
    select(-event_start, -anchor_date)
}

# Level covariates used in A02's event-anchored models (build_rhs_event),
# without the CFTR indicator so they can be reused for the CFTR model.
event_level_covariates <- function(flags) {
  terms <- c("pancreatic_insufficiency", "cfrd_indicator", "baseline_fev1_pp_c", "baseline_age_years_c")
  if (flags$use_nut)    terms <- c(terms, "nutrition_cat")
  if (flags$use_exac)   terms <- c(terms, "exac_cat")
  if (flags$use_sector) terms <- c(terms, "sector")
  if (flags$use_eth)    terms <- c(terms, "ethnicity_cat")
  if (flags$use_geno)   terms <- c(terms, "f508_genotype_group")
  c(terms, "chronic_paeruginosa", "chronic_bcepacia", "chronic_mssa", "chronic_mrsa")
}

build_rhs_event <- function(flags, time_var_cs, post_var) {
  paste(c(time_var_cs, post_var, event_level_covariates(flags), "cftr_modulator_indicator",
          paste0(time_var_cs, ":", post_var)), collapse = " + ")
}

rhs_drop <- function(rhs, drop) {
  paste(setdiff(strsplit(rhs, " + ", fixed = TRUE)[[1]], drop), collapse = " + ")
}

flags_main <- compute_use_flags(dat_model)

# ------------------------------------------------------------------------------
# 3. Contrast helpers (Satterthwaite df via lmerTest::contest1D)
# ------------------------------------------------------------------------------
lc <- function(mod, w, label) {
  fe <- lme4::fixef(mod)
  if (!all(names(w) %in% names(fe))) stop("Unknown term(s) in contrast: ", paste(setdiff(names(w), names(fe)), collapse = ", "))
  L <- setNames(numeric(length(fe)), names(fe)); L[names(w)] <- w
  r <- tryCatch(as.data.frame(lmerTest::contest1D(mod, unname(L), confint = TRUE)), error = function(e) NULL)
  need <- c("Estimate", "Std. Error", "df", "lower", "upper", "Pr(>|t|)")
  if (!is.null(r) && all(need %in% names(r))) {
    tibble(contrast = label, estimate = r[["Estimate"]], se = r[["Std. Error"]], df = r[["df"]],
           conf.low = r[["lower"]], conf.high = r[["upper"]], p.value = r[["Pr(>|t|)"]])
  } else {
    V <- as.matrix(vcov(mod)); est <- sum(L * fe); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    tibble(contrast = label, estimate = est, se = se, df = NA_real_,
           conf.low = est - 1.96 * se, conf.high = est + 1.96 * se, p.value = 2 * pnorm(-abs(est / se)))
  }
}

seg_rows <- function(mod, label, t_term, slope_term, post_term = NULL) {
  one <- function(a) setNames(1, a)
  rows <- list(
    lc(mod, one(t_term), "Pre-event slope (pp/yr)"),
    if (!is.null(post_term)) lc(mod, one(post_term), "Level change at event (pp)"),
    lc(mod, one(slope_term), "Slope change after event (pp/yr)"),
    lc(mod, setNames(c(1, 1), c(t_term, slope_term)), "Post-event slope, net (pp/yr)"),
    lc(mod, if (is.null(post_term)) one(slope_term) else setNames(c(1, 1), c(post_term, slope_term)),
       "Difference from pre-event trend at +1 yr (pp)"),
    lc(mod, if (is.null(post_term)) setNames(2, slope_term) else setNames(c(1, 2), c(post_term, slope_term)),
       "Difference from pre-event trend at +2 yr (pp)")
  )
  bind_rows(rows) %>% mutate(model = label, .before = 1)
}

ci_txt <- function(tb, model, contrast) {
  r <- tb %>% filter(.data$model == !!model, .data$contrast == !!contrast)
  if (!nrow(r)) return("not estimated")
  sprintf("%+.2f (95%% CI %+.2f to %+.2f; p=%.3f)", r$estimate, r$conf.low, r$conf.high, r$p.value)
}

# ------------------------------------------------------------------------------
# 4. FIX 2: global annual slope
# ------------------------------------------------------------------------------
say("")
say("=== FIX 2: global annual slope ===")
m0      <- fit_triple_fallback("time_cs", dat_model, use_reml = FALSE, label = "M0 unconditional (ML)")
m1      <- fit_triple_fallback(build_rhs(flags_main, "m1"), dat_model, use_reml = FALSE, label = "M1 main effects (ML)")
m2_reml <- fit_triple_fallback(build_rhs(flags_main, "m2"), dat_model, use_reml = TRUE, label = "M2 interactions (REML)")

# Derivative of each row's linear predictor with respect to time (exact, since
# the model is linear in time_cs), averaged over rows or over patients.
slope_design <- function(mod, time_var = "time_cs") {
  mf <- model.frame(mod)
  ff <- lme4::nobars(formula(mod))
  d1 <- mf; d1[[time_var]] <- d1[[time_var]] + 1
  D <- model.matrix(ff, data = d1) - model.matrix(ff, data = mf)
  list(D = D[, names(lme4::fixef(mod)), drop = FALSE], id = mf$patient_id)
}
sd2 <- slope_design(m2_reml)
w_obs <- colMeans(sd2$D)
w_pat <- colMeans(rowsum(sd2$D, sd2$id) / as.vector(table(sd2$id)[unique(sort(sd2$id))]))

fix2 <- bind_rows(
  lc(m0, c(time_cs = 1), "Unadjusted slope (M0)") %>% mutate(model = "M0 unconditional (ML)", .before = 1),
  lc(m1, c(time_cs = 1), "Adjusted slope, main-effects model (M1)") %>% mutate(model = "M1 main effects (ML)", .before = 1),
  lc(m2_reml, c(time_cs = 1), "M2 time_cs coefficient = reference-profile slope (as reported)") %>% mutate(model = "M2 interactions (REML)", .before = 1),
  lc(m2_reml, w_obs, "M2 population-averaged slope (visit-weighted)") %>% mutate(model = "M2 interactions (REML)", .before = 1),
  lc(m2_reml, w_pat, "M2 population-averaged slope (patient-weighted)") %>% mutate(model = "M2 interactions (REML)", .before = 1)
)
show_tbl(fix2)
say("Slope-modifier means used for the visit-weighted average (share of visits with each modifier):")
show_tbl(tibble(term = names(w_obs), mean = round(w_obs, 3)) %>% filter(term != "time_cs", mean != 0))
write_csv(fix2, file.path(out_dir, "fix2_global_slope.csv"))

# ------------------------------------------------------------------------------
# 5. FIX 1: CFTR modulator segmented model
# ------------------------------------------------------------------------------
say("")
say("=== FIX 1: CFTR modulator start dates ===")
s1 <- suppressMessages(readxl::read_excel(raw_file, sheet = "Sheet1")) %>% janitor::clean_names()
id_cols <- grep("^patcode", names(s1), value = TRUE)
id_df <- s1[id_cols]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)

to_date <- function(x) {
  if (inherits(x, c("POSIXct", "POSIXt", "Date"))) return(as.Date(x))
  suppressWarnings(lubridate::dmy(as.character(x), quiet = TRUE))
}
date_cols <- intersect(c("date_cftrm2018", "date_cftrm2019", "date_cftrm2020", "date_cftrm2021", "date_cftrm2022",
                         "date_cftrm2023", "date_cftrm_commenced_2024", "treatment_start_date_cftr",
                         "treatment_start_date_cftr_new"), names(s1))
lo <- as.Date("2010-01-01"); hi <- as.Date("2024-12-31")
dm <- sapply(date_cols, function(v) {
  d <- to_date(s1[[v]]); d[!is.na(d) & (d < lo | d > hi)] <- NA; as.numeric(d)
})
earliest_num <- suppressWarnings(apply(dm, 1, min, na.rm = TRUE))
earliest_num[!is.finite(earliest_num)] <- NA
earliest <- as.Date(earliest_num, origin = "1970-01-01")

flag_cols <- grep("^cftr_modifier_?20[0-9]{2}$", names(s1), value = TRUE)
flag_year <- as.integer(stringr::str_extract(flag_cols, "20[0-9]{2}"))
Y <- sapply(flag_cols, function(v) as.character(s1[[v]]) %in% "Yes")
first_yes <- apply(Y, 1, function(r) if (any(r)) min(flag_year[r]) else NA_integer_)
impute <- !is.na(first_yes) & (is.na(earliest) | as.integer(format(earliest, "%Y")) > first_yes)
start_rev <- earliest
start_rev[impute] <- as.Date(paste0(first_yes[impute], "-07-01"))

starts <- tibble(patient_id = pid, cftr_start_rev = start_rev, start_imputed = impute) %>%
  filter(!is.na(patient_id)) %>%
  group_by(patient_id) %>%
  summarise(cftr_start_rev = suppressWarnings(min(cftr_start_rev, na.rm = TRUE)),
            start_imputed = any(start_imputed), .groups = "drop") %>%
  mutate(cftr_start_rev = dplyr::if_else(is.finite(as.numeric(cftr_start_rev)), cftr_start_rev, as.Date(NA)))

pt <- dat_model %>%
  group_by(patient_id) %>%
  summarise(cftr_start_a01 = dplyr::first(na.omit(cftr_modulator_start_date)),
            ever_cftr_a01 = max(ever_cftr_yes, na.rm = TRUE), .groups = "drop") %>%
  left_join(starts, by = "patient_id")

date_counts <- tibble(
  item = c("Patients in model data", "Matched to Sheet1", "Ever on modulator (A01 flag)",
           "With A01 start date", "With revised start date (2010-2024)",
           "Revised start imputed from annual Yes flag (1 July of first flagged year)",
           "Revised start >30 days earlier than A01 start", "Revised start >30 days later than A01 start",
           "A01 start missing, revised start present", "Ever (A01) but no start date or flag"),
  n = c(nrow(pt), sum(pt$patient_id %in% starts$patient_id), sum(pt$ever_cftr_a01 == 1),
        sum(!is.na(pt$cftr_start_a01)), sum(!is.na(pt$cftr_start_rev)),
        sum(pt$start_imputed %in% TRUE & !is.na(pt$cftr_start_rev)),
        sum(as.numeric(pt$cftr_start_a01 - pt$cftr_start_rev) > 30, na.rm = TRUE),
        sum(as.numeric(pt$cftr_start_rev - pt$cftr_start_a01) > 30, na.rm = TRUE),
        sum(is.na(pt$cftr_start_a01) & !is.na(pt$cftr_start_rev)),
        sum(pt$ever_cftr_a01 == 1 & is.na(pt$cftr_start_rev)))
)
show_tbl(date_counts)
write_csv(date_counts, file.path(out_dir, "fix1_start_date_reconciliation.csv"))

dat_c <- dat_model %>%
  left_join(pt %>% select(patient_id, cftr_start_rev, cftr_start_a01, ever_cftr_a01), by = "patient_id") %>%
  mutate(
    cftr_post     = as.integer(!is.na(cftr_start_rev) & visit_date >= cftr_start_rev),
    cftr_tsi      = dplyr::if_else(cftr_post == 1L, as.numeric(visit_date - cftr_start_rev) / 365.25, 0),
    cftr_post_a01 = as.integer(!is.na(cftr_start_a01) & visit_date >= cftr_start_a01),
    cftr_tsi_a01  = dplyr::if_else(cftr_post_a01 == 1L, as.numeric(visit_date - cftr_start_a01) / 365.25, 0)
  )

cftr_pt <- dat_c %>%
  group_by(patient_id) %>%
  summarise(
    has_start   = any(!is.na(cftr_start_rev)),
    ever_a01    = max(ever_cftr_a01),
    n_pre_fev1  = sum(!is.na(fev1_pp) & cftr_post == 0L & !is.na(cftr_start_rev)),
    n_post_fev1 = sum(!is.na(fev1_pp) & cftr_post == 1L),
    post_years  = if (any(cftr_post == 1L & !is.na(fev1_pp))) max(cftr_tsi[cftr_post == 1L & !is.na(fev1_pp)]) else NA_real_,
    .groups = "drop") %>%
  mutate(group = case_when(
    !has_start & ever_a01 == 1 ~ "Ever (A01) but start unknown",
    !has_start                 ~ "Never treated",
    n_post_fev1 == 0           ~ "Started after last FEV1",
    n_pre_fev1 == 0            ~ "Treated from first FEV1 (prevalent user)",
    TRUE                       ~ "New user (FEV1 before and after start)"))

grp_counts <- cftr_pt %>% count(group, name = "n_patients")
say("CFTR exposure groups (revised start dates):")
show_tbl(grp_counts)
nu <- cftr_pt %>% filter(group == "New user (FEV1 before and after start)")
say("New users: ", nrow(nu), "; pre-start FEV1 per patient ", med_iqr(nu$n_pre_fev1, 0),
    "; post-start FEV1 ", med_iqr(nu$n_post_fev1, 0), "; >=2 post-start FEV1: ", sum(nu$n_post_fev1 >= 2),
    "; post-start follow-up (years) ", med_iqr(nu$post_years, 1))
write_csv(grp_counts, file.path(out_dir, "fix1_cftr_exposure_groups.csv"))

ids_primary <- cftr_pt %>% filter(group %in% c("Never treated", "Started after last FEV1", "New user (FEV1 before and after start)")) %>% pull(patient_id)
ids_2post   <- cftr_pt %>% filter(group %in% c("Never treated", "Started after last FEV1") |
                                    (group == "New user (FEV1 before and after start)" & n_post_fev1 >= 2)) %>% pull(patient_id)

# C0: A02's M2_cftr, reproduced (slope change only, anchored at last pre-CFTR FEV1)
cftr_index <- dat_model %>%
  filter(ever_cftr_yes == 1L, !is.na(cftr_modulator_start_date), visit_date < cftr_modulator_start_date,
         !is.na(fev1_pp), is.finite(fev1_pp)) %>%
  arrange(patient_id, visit_date) %>%
  group_by(patient_id) %>% slice_tail(n = 1) %>% ungroup() %>%
  transmute(patient_id, cftr_index_date = visit_date)
dat_model_cftr <- dat_model %>%
  left_join(cftr_index, by = "patient_id") %>%
  mutate(
    cftr_index_date = dplyr::coalesce(cftr_index_date, baseline_date),
    time_cftr_index = as.numeric(difftime(visit_date, cftr_index_date, units = "days")) / 365.25,
    time_cftr_index_cs = time_cftr_index,
    cftr_post_index = dplyr::if_else(
      slope_cftr_ok == 1L & ever_cftr_yes == 1L & !is.na(time_cftr_index) & time_cftr_index >= 0, 1L, 0L, missing = 0L)
  )
rhs_c0 <- build_rhs_cftr(flags_main)
m_c0 <- fit_triple_fallback(rhs_c0, dat_model_cftr, use_reml = TRUE, label = "C0 original (slope only)")
m_c1 <- fit_triple_fallback(paste(rhs_c0, "+ cftr_post_index"), dat_model_cftr, use_reml = TRUE, label = "C1 original + level term")

# C2-C6: segmented at the start date, level + slope change, level covariates only
cov_tv <- event_level_covariates(flags_main)
cov_fixed <- c("baseline_fev1_pp_c", "baseline_age_years_c", "sex", "sector", "ethnicity_cat",
               "f508_genotype_group", "baseline_pancreatic_insuff")
rhs_seg  <- paste(c("time_cs", "cftr_post", "cftr_tsi", cov_tv), collapse = " + ")
rhs_segf <- paste(c("time_cs", "cftr_post", "cftr_tsi", cov_fixed), collapse = " + ")
rhs_a01  <- paste(c("time_cs", "cftr_post_a01", "cftr_tsi_a01", cov_tv), collapse = " + ")

m_c2 <- fit_triple_fallback(rhs_seg,  dat_c %>% filter(patient_id %in% ids_primary), TRUE, "C2 corrected primary")
m_c3 <- fit_triple_fallback(rhs_segf, dat_c %>% filter(patient_id %in% ids_primary), TRUE, "C3 baseline covariates only")
m_c4 <- fit_triple_fallback(rhs_a01,  dat_c %>% filter(patient_id %in% ids_primary), TRUE, "C4 A01 start dates")
m_c5 <- fit_triple_fallback(rhs_seg,  dat_c %>% filter(patient_id %in% ids_2post),   TRUE, "C5 new users with >=2 post FEV1")
m_c6 <- fit_triple_fallback(rhs_seg,  dat_c,                                         TRUE, "C6 all patients")

fix1 <- bind_rows(
  seg_rows(m_c0, "C0 original (slope only)", "time_cftr_index_cs", "time_cftr_index_cs:cftr_post_index"),
  seg_rows(m_c1, "C1 original + level term", "time_cftr_index_cs", "time_cftr_index_cs:cftr_post_index", "cftr_post_index"),
  seg_rows(m_c2, "C2 corrected primary", "time_cs", "cftr_tsi", "cftr_post"),
  seg_rows(m_c3, "C3 baseline covariates only", "time_cs", "cftr_tsi", "cftr_post"),
  seg_rows(m_c4, "C4 A01 start dates", "time_cs", "cftr_tsi_a01", "cftr_post_a01"),
  seg_rows(m_c5, "C5 new users with >=2 post FEV1", "time_cs", "cftr_tsi", "cftr_post"),
  seg_rows(m_c6, "C6 all patients", "time_cs", "cftr_tsi", "cftr_post")
)
say("CFTR models (C0/C1 pre-event slope is the reference-profile slope, as in A02):")
show_tbl(fix1 %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
write_csv(fix1, file.path(out_dir, "fix1_cftr_models.csv"))
saveRDS(m_c2, file.path(out_dir, "models", "C2_cftr_segmented.rds"))

# ------------------------------------------------------------------------------
# 6. FIX 3 + FIX 4: infection event models, net slopes within model, incident-only
# ------------------------------------------------------------------------------
say("")
say("=== FIX 3 + 4: infection event-anchored models ===")
event_specs <- list(
  list(name = "pa",   start_col = "pa_chronic_start_date",   slope_ok_col = "slope_pa_ok",
       status_col = "pa_status", own = "chronic_paeruginosa", label = "P. aeruginosa"),
  list(name = "mssa", start_col = "mssa_chronic_start_date", slope_ok_col = "slope_mssa_ok",
       status_col = "sa_status", own = "chronic_mssa", label = "MSSA")
)

classify_events <- function(df, start_col, status_col) {
  df %>%
    group_by(patient_id) %>%
    summarise(
      onset = dplyr::first(.data[[start_col]]),
      n_pre_nonchronic = sum(!is.na(onset) & visit_date < onset & !is.na(.data[[status_col]]) & .data[[status_col]] != "Chronic"),
      n_pre_fev1  = sum(!is.na(onset) & visit_date <  onset & !is.na(fev1_pp)),
      n_post_fev1 = sum(!is.na(onset) & visit_date >= onset & !is.na(fev1_pp)),
      .groups = "drop") %>%
    mutate(event_class = case_when(is.na(onset) ~ "never", n_pre_nonchronic >= 1 ~ "incident", TRUE ~ "prevalent"))
}

event_rows <- list(); event_counts <- list(); repro <- list()
for (spec in event_specs) {
  tv <- paste0("time_", spec$name, "_index_cs"); pv <- paste0(spec$name, "_post_index"); iv <- paste0(tv, ":", pv)

  ev <- classify_events(dat_model, spec$start_col, spec$status_col)
  inc <- ev %>% filter(event_class == "incident")
  event_counts[[spec$label]] <- tibble(
    organism = spec$label,
    never = sum(ev$event_class == "never"), prevalent = sum(ev$event_class == "prevalent"),
    incident = nrow(inc), incident_with_pre_and_2post = sum(inc$n_pre_fev1 >= 1 & inc$n_post_fev1 >= 2),
    incident_pre_fev1 = med_iqr(inc$n_pre_fev1, 0), incident_post_fev1 = med_iqr(inc$n_post_fev1, 0))

  # E0: A02 original (post requires >=2 post-onset FEV1; own chronic indicator kept)
  d0 <- build_event_data(dat_model, spec)
  m_e0 <- fit_triple_fallback(build_rhs_event(flags_main, tv, pv), d0, TRUE, paste(spec$label, "E0 original"))
  repro[[spec$label]] <- unname(lme4::fixef(m_e0)[[iv]])

  # E1: post coded from onset for every patient; own time-varying chronic indicator dropped
  d1 <- dat_model %>%
    mutate(event_start = .data[[spec$start_col]],
           anchor_date = dplyr::coalesce(event_start, baseline_date),
           !!tv := as.numeric(difftime(visit_date, anchor_date, units = "days")) / 365.25,
           !!pv := as.integer(!is.na(event_start) & visit_date >= event_start))
  rhs1 <- rhs_drop(build_rhs_event(flags_main, tv, pv), spec$own)
  m_e1 <- fit_triple_fallback(rhs1, d1, TRUE, paste(spec$label, "E1 all events"))

  # E2: incident events only (prevalent patients excluded; never-infected kept)
  keep2 <- ev %>% filter(event_class %in% c("never", "incident")) %>% pull(patient_id)
  m_e2 <- fit_triple_fallback(rhs1, d1 %>% filter(patient_id %in% keep2), TRUE, paste(spec$label, "E2 incident only"))

  # E3: incident events with >=1 pre-onset and >=2 post-onset FEV1
  keep3 <- ev %>% filter(event_class == "never" | (event_class == "incident" & n_pre_fev1 >= 1 & n_post_fev1 >= 2)) %>% pull(patient_id)
  m_e3 <- fit_triple_fallback(rhs1, d1 %>% filter(patient_id %in% keep3), TRUE, paste(spec$label, "E3 incident, >=1 pre and >=2 post"))

  event_rows[[spec$label]] <- bind_rows(
    seg_rows(m_e0, paste(spec$label, "E0 original"), tv, iv, pv),
    seg_rows(m_e1, paste(spec$label, "E1 all events"), tv, iv, pv),
    seg_rows(m_e2, paste(spec$label, "E2 incident only"), tv, iv, pv),
    seg_rows(m_e3, paste(spec$label, "E3 incident, >=1 pre and >=2 post"), tv, iv, pv)
  )
  saveRDS(m_e2, file.path(out_dir, "models", paste0("E2_", spec$name, "_incident.rds")))
}
ev_counts <- bind_rows(event_counts)
ev_tab <- bind_rows(event_rows)
say("Event classes (patients):")
show_tbl(ev_counts)
show_tbl(ev_tab %>% select(model, contrast, estimate, conf.low, conf.high, p.value))
write_csv(ev_counts, file.path(out_dir, "fix4_event_classes.csv"))
write_csv(ev_tab, file.path(out_dir, "fix34_event_models.csv"))

# ------------------------------------------------------------------------------
# 7. Reproduction check, fit log and headline comparison
# ------------------------------------------------------------------------------
say("")
say("=== Reproduction check (original estimates should match the manuscript) ===")
rc <- tibble(
  estimate = c("M0 time_cs", "M2 REML time_cs (reported global slope)", "C0 CFTR slope change",
               "P. aeruginosa E0 slope change", "MSSA E0 slope change"),
  reproduced = c(unname(lme4::fixef(m0)[["time_cs"]]), unname(lme4::fixef(m2_reml)[["time_cs"]]),
                 unname(lme4::fixef(m_c0)[["time_cftr_index_cs:cftr_post_index"]]),
                 repro[["P. aeruginosa"]], repro[["MSSA"]]),
  manuscript_or_saved = c(-0.602, -1.553, 2.949, -0.966, -1.669)
) %>% mutate(difference = round(reproduced - manuscript_or_saved, 3))
show_tbl(rc)
write_csv(rc, file.path(out_dir, "reproduction_check.csv"))

fl <- bind_rows(fit_log)
say("")
say("=== Random-effects structure used and fitting warnings ===")
show_tbl(fl %>% select(model, re_structure, n_patients, n_obs, n_warnings))
write_csv(fl, file.path(out_dir, "fit_log.csv"))

say("")
say("=== Headline comparison ===")
say("Global annual slope: reported -1.55 (M2 reference profile)")
say("   unadjusted M0:          ", ci_txt(fix2, "M0 unconditional (ML)", "Unadjusted slope (M0)"))
say("   adjusted M1:            ", ci_txt(fix2, "M1 main effects (ML)", "Adjusted slope, main-effects model (M1)"))
say("   M2 population-averaged: ", ci_txt(fix2, "M2 interactions (REML)", "M2 population-averaged slope (visit-weighted)"))
for (o in c("P. aeruginosa", "MSSA")) {
  say(o, " net post-onset slope (reported ", if (o == "MSSA") "-3.22" else "-2.52", "):")
  for (m in c("E0 original", "E1 all events", "E2 incident only", "E3 incident, >=1 pre and >=2 post")) {
    say("   ", m, ": ", ci_txt(ev_tab, paste(o, m), "Post-event slope, net (pp/yr)"),
        " | level change ", ci_txt(ev_tab, paste(o, m), "Level change at event (pp)"),
        " | slope change ", ci_txt(ev_tab, paste(o, m), "Slope change after event (pp/yr)"))
  }
}
say("CFTR (reported slope change +2.95, net +1.40):")
for (m in unique(fix1$model)) {
  say("   ", m, ": level ", ci_txt(fix1, m, "Level change at event (pp)"),
      " | slope change ", ci_txt(fix1, m, "Slope change after event (pp/yr)"),
      " | pre slope ", ci_txt(fix1, m, "Pre-event slope (pp/yr)"),
      " | post slope ", ci_txt(fix1, m, "Post-event slope, net (pp/yr)"))
}
say("")
say("Done. Outputs in ", normalizePath(out_dir))
