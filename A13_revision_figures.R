#!/usr/bin/env Rscript
# ==============================================================================
# A13 - REVISION FIGURES (JCF-D-26-00377, corrected data)
# ------------------------------------------------------------------------------
# Figure 2 (revised): CFTR modulator panel only. New users aligned at the start
#   date; observed mean change from the last pre-treatment FEV1pp in half-year
#   windows, with the C2 segmented-model estimate and the pre-treatment trend
#   continued (counterfactual).
# Figure S1: cohort flowchart (STROBE) with longitudinal data availability.
# Run after A08 (needs outputs/revision/models/C2_cftr_segmented.rds) and A09.
# Writes: outputs/revision/figures/ only.
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(readxl); library(janitor); library(lme4) })
fig_dir <- "outputs/revision/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
use_ragg <- requireNamespace("ragg", quietly = TRUE)
save_fig <- function(p, name, w, h) {
  if (use_ragg) ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, device = ragg::agg_png)
  else ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, device = pdf_dev)
  message("Saved ", name)
}

# ------------------------------------------------------------------------------
# Figure 2: CFTR modulator initiation
# ------------------------------------------------------------------------------
dat <- arrow::read_parquet("outputs/data_derived/analysis_long.parquet") %>% as_tibble()
s1 <- suppressMessages(readxl::read_excel("data_raw/CF_registry_master.xlsx", sheet = "Sheet1")) %>% janitor::clean_names()
id_df <- s1[grep("^patcode", names(s1), value = TRUE)]; id_df[] <- lapply(id_df, as.character)
pid <- dplyr::coalesce(!!!id_df)
to_date2 <- function(x) { if (inherits(x, c("POSIXct", "POSIXt", "Date"))) return(as.Date(x)); suppressWarnings(lubridate::dmy(as.character(x), quiet = TRUE)) }
dcols <- intersect(c("date_cftrm2018", "date_cftrm2019", "date_cftrm2020", "date_cftrm2021", "date_cftrm2022", "date_cftrm2023",
                     "date_cftrm_commenced_2024", "treatment_start_date_cftr", "treatment_start_date_cftr_new"), names(s1))
dmat <- sapply(dcols, function(v) { d <- to_date2(s1[[v]]); d[!is.na(d) & (d < as.Date("2010-01-01") | d > as.Date("2024-12-31"))] <- NA; as.numeric(d) })
en <- suppressWarnings(apply(dmat, 1, min, na.rm = TRUE)); en[!is.finite(en)] <- NA
earliest <- as.Date(en, origin = "1970-01-01")
fcols <- grep("^cftr_modifier_?20[0-9]{2}$", names(s1), value = TRUE); fy <- as.integer(stringr::str_extract(fcols, "20[0-9]{2}"))
Yf <- sapply(fcols, function(v) as.character(s1[[v]]) %in% "Yes")
fyes <- apply(Yf, 1, function(r) if (any(r)) min(fy[r]) else NA_integer_)
imp <- !is.na(fyes) & (is.na(earliest) | as.integer(format(earliest, "%Y")) > fyes)
srev <- earliest; srev[imp] <- as.Date(paste0(fyes[imp], "-07-01"))
starts <- tibble(patient_id = pid, start = srev) %>% filter(!is.na(patient_id), !is.na(start)) %>%
  group_by(patient_id) %>% summarise(start = min(start), .groups = "drop")

d <- dat %>% filter(!is.na(fev1_pp)) %>% inner_join(starts, by = "patient_id") %>%
  mutate(t = as.numeric(visit_date - start) / 365.25)
anchor <- d %>% filter(t < 0) %>% group_by(patient_id) %>% slice_max(visit_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(patient_id, anchor_date = visit_date, anchor_fev1 = fev1_pp, t_anchor = t)
new_users <- d %>% filter(t >= 0) %>% distinct(patient_id) %>% semi_join(anchor, by = "patient_id") %>% pull(patient_id)
obs <- d %>% filter(patient_id %in% new_users) %>% inner_join(anchor, by = "patient_id") %>%
  filter(visit_date != anchor_date, t >= -3, t < 2.5) %>%
  mutate(delta = fev1_pp - anchor_fev1, bin = floor(t / 0.5) * 0.5 + 0.25)
binned <- obs %>% group_by(bin) %>%
  summarise(n = n(), mean = mean(delta), se = sd(delta) / sqrt(n()), .groups = "drop") %>%
  filter(n >= 20) %>%
  mutate(lo = mean - qt(0.975, n - 1) * se, hi = mean + qt(0.975, n - 1) * se)
t_a <- median(anchor$t_anchor[anchor$patient_id %in% new_users])
message("New users plotted: ", length(new_users), "; median time of last pre-treatment FEV1: ", round(t_a, 2), " years")

c2 <- readRDS("outputs/revision/models/C2_cftr_segmented.rds")
fe <- lme4::fixef(c2); V <- as.matrix(vcov(c2))
pred_line <- function(tt, post, counterfactual = FALSE) {
  L <- matrix(0, nrow = length(tt), ncol = length(fe), dimnames = list(NULL, names(fe)))
  L[, "time_cs"] <- tt - t_a
  if (post && !counterfactual) { L[, "cftr_post"] <- 1; L[, "cftr_tsi"] <- tt }
  est <- as.vector(L %*% fe); se <- sqrt(rowSums((L %*% V) * L))
  tibble(t = tt, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}
pre  <- pred_line(seq(-3, 0, by = 0.02), post = FALSE) %>% mutate(segment = "pre")
post <- pred_line(seq(0, 2.5, by = 0.02), post = TRUE) %>% mutate(segment = "post")
cf   <- pred_line(seq(0, 2.5, by = 0.02), post = TRUE, counterfactual = TRUE)

col_model <- "#1F4E79"; col_obs <- "grey15"; col_cf <- "grey45"
y_lo <- min(c(binned$lo, pre$lo, cf$lo), na.rm = TRUE) - 3
gap_lo <- cf$est[1]; gap_hi <- post$est[1]; gap_x <- 0.1
p2 <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(data = bind_rows(pre, post), aes(x = t, ymin = lo, ymax = hi, group = segment), fill = col_model, alpha = 0.15) +
  geom_line(data = cf, aes(x = t, y = est, linetype = "Pre-treatment trend continued"), colour = col_cf, linewidth = 0.7) +
  geom_line(data = bind_rows(pre, post), aes(x = t, y = est, group = segment, linetype = "Model estimate (95% CI)"), colour = col_model, linewidth = 0.9) +
  geom_errorbar(data = binned, aes(x = bin, ymin = lo, ymax = hi), width = 0.08, colour = col_obs, linewidth = 0.5) +
  geom_point(data = binned, aes(x = bin, y = mean, shape = "Observed mean (95% CI)"), colour = col_obs, size = 2.2) +
  geom_text(data = binned, aes(x = bin, y = y_lo, label = n), size = 2.8, colour = "grey35") +
  annotate("text", x = -3.05, y = y_lo, label = "n", hjust = 1, size = 2.8, colour = "grey35", fontface = "italic") +
  annotate("text", x = 0.05, y = max(c(binned$hi, post$hi), na.rm = TRUE) + 1, label = "Modulator start",
           hjust = 0, size = 3, colour = "grey30") +
  annotate("segment", x = gap_x, xend = gap_x, y = gap_lo + 0.3, yend = gap_hi - 0.3, colour = col_model,
           linewidth = 0.5, arrow = arrow(ends = "both", length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = gap_x + 0.07, y = (gap_lo + gap_hi) / 2, label = sprintf("+%.1f points", fe[["cftr_post"]]),
           hjust = 0, size = 3.2, colour = col_model, fontface = "bold") +
  scale_linetype_manual(NULL, values = c("Model estimate (95% CI)" = "solid", "Pre-treatment trend continued" = "dashed")) +
  scale_shape_manual(NULL, values = c("Observed mean (95% CI)" = 16)) +
  scale_x_continuous(breaks = seq(-3, 2.5, by = 0.5), limits = c(-3.2, 2.6)) +
  labs(x = "Years from CFTR modulator initiation",
       y = expression(atop("Change in" ~ FEV[1] ~ "% predicted", "from last pre-treatment test (points)"))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "left",
        legend.key.width = unit(1.6, "lines"), legend.box = "horizontal", legend.margin = margin(0, 0, 0, 0),
        axis.title.y = element_text(size = 10), plot.margin = margin(4, 10, 6, 6)) +
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2))
save_fig(p2, "Figure2_CFTR_initiation", 7, 4.6)
write_csv(binned, file.path(fig_dir, "Figure2_observed_bins.csv"))

# ------------------------------------------------------------------------------
# Figure S1: cohort flowchart
# ------------------------------------------------------------------------------
flow <- read_csv("outputs/revision/row7_cohort_flow.csv", show_col_types = FALSE)
nfev <- dat %>% filter(!is.na(fev1_pp)) %>% count(patient_id, name = "k") %>% count(k, name = "patients")
fu <- dat %>% filter(!is.na(fev1_pp)) %>% group_by(patient_id) %>%
  summarise(k = n(), span = as.numeric(max(visit_date) - min(visit_date)) / 365.25, .groups = "drop")
q <- function(x) sprintf("%.1f (IQR %.1f–%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
qi <- function(x) sprintf("%d (IQR %d–%d)", as.integer(median(x)), as.integer(quantile(x, 0.25)), as.integer(quantile(x, 0.75)))
fmt <- function(x) format(x, big.mark = ",")
main_lab <- c(
  paste0("Registered in the SACFR extract\nn = ", fmt(flow$patients[1])),
  paste0("At least one annual registry record\nn = ", fmt(flow$patients[2])),
  paste0("After removing records on or after lung transplantation\nn = ", fmt(flow$patients[3])),
  paste0("A registry record at age ≥6 years\nn = ", fmt(flow$patients[4])),
  paste0("An FEV1pp value at age ≥6 years\nn = ", fmt(flow$patients[5])),
  paste0("Analytic cohort: FEV1pp in ≥2 calendar years\nn = ", fmt(flow$patients[6]), " patients; ", fmt(sum(fu$k)), " FEV1pp values"),
  paste0("Included in mixed-effects models\nn = ", fmt(flow$patients[7]), " patients; ", fmt(flow$fev1_values[7]), " FEV1pp values"))
excl_lab <- c(NA,
  paste0("No annual registry record (n = ", flow$excluded_at_step[2], ")"),
  paste0("Records only after lung transplantation (n = ", flow$excluded_at_step[3], ")"),
  paste0("No record at age ≥6 years (n = ", flow$excluded_at_step[4], ")"),
  paste0("No FEV1pp at age ≥6 years (n = ", flow$excluded_at_step[5], ")"),
  paste0("FEV1pp in fewer than 2 calendar years (n = ", flow$excluded_at_step[6], ")"),
  paste0("Fewer than 2 records with known\nnutritional status (n = ", flow$excluded_at_step[7], ")"))
avail_lab <- paste0("Longitudinal data in the analytic cohort\n",
                    "FEV1pp values per patient: median ", qi(fu$k), "\n",
                    paste(paste0(nfev$k, " values: ", nfev$patients), collapse = "  ·  "), "\n",
                    "Follow-up from first to last FEV1pp: median ", q(fu$span), " years")
ys <- seq(0, by = -1.15, length.out = 7)
boxes <- tibble(x = 0, y = ys, label = main_lab)
# exclusion i (i = 2..7) happens between main box i-1 and main box i
ex <- tibble(x = 1.25, y = ys + 0.575, label = excl_lab) %>% filter(!is.na(label))
pf <- ggplot() +
  geom_segment(data = tibble(y0 = ys[-7] - 0.3, y1 = ys[-1] + 0.3), aes(x = 0, xend = 0, y = y0, yend = y1),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"), colour = "grey30") +
  geom_segment(data = ex, aes(x = 0, xend = 1.22, y = y, yend = y), arrow = arrow(length = unit(0.12, "cm"), type = "closed"), colour = "grey30") +
  geom_label(data = boxes, aes(x = x, y = y, label = label), size = 3, label.r = unit(0.08, "lines"),
             label.padding = unit(0.35, "lines"), fill = "white", lineheight = 1.05) +
  geom_label(data = ex, aes(x = x, y = y, label = label), hjust = 0, size = 2.8, label.r = unit(0.08, "lines"),
             label.padding = unit(0.3, "lines"), fill = "grey96", colour = "grey20", lineheight = 1.05) +
  geom_label(aes(x = 0, y = ys[7] - 1.15, label = avail_lab), size = 2.8, label.r = unit(0.08, "lines"),
             label.padding = unit(0.35, "lines"), fill = "#EEF3F8", lineheight = 1.1) +
  coord_cartesian(xlim = c(-2.2, 3.6), ylim = c(ys[7] - 1.65, 0.35), clip = "off") +
  theme_void() + theme(plot.margin = margin(10, 10, 10, 10))
save_fig(pf, "FigureS1_cohort_flowchart", 8, 9)
message("Done.")
