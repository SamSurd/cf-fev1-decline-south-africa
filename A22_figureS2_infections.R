#!/usr/bin/env Rscript
# ==============================================================================
# A22 - SUPPLEMENTARY FIGURE S2 (JCF-D-26-00377, corrected data after fix 2)
# ------------------------------------------------------------------------------
# Redraws the former Figure 2B-C on the same basis as revised Figure 2:
# patients who acquired chronic P. aeruginosa (A) or chronic MSSA (B) during
# follow-up (incident events, A08 rule: >=1 earlier registry year with recorded
# non-chronic status), aligned at onset (first registry year recorded as
# chronic). Points: observed mean change from each patient's last FEV1pp before
# onset, in half-year windows with >=10 patients. Line: incident-only
# event-anchored model (A08 E2), set to zero at the median time of the last
# pre-onset test; dashed line: pre-onset trend continued.
# Replaces the panels from Supplementary Scripts/Normalised * ever_vs_never_nosmooth.R.
# Run after A08 (needs outputs/revision/models/E2_*_incident.rds).
# Writes: outputs/revision/figures/FigureS2_infection_onset.{png,pdf}, FigureS2_observed_bins.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(arrow); library(lme4); library(grid) })
fig_dir <- "outputs/revision/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

dat <- arrow::read_parquet("outputs/data_derived/analysis_long.parquet") %>% as_tibble() %>%
  mutate(nutrition_cat = as.character(nutrition_cat)) %>%
  filter(nutrition_cat != "Unknown" | is.na(nutrition_cat)) %>%
  group_by(patient_id) %>% filter(dplyr::n() >= 2) %>% ungroup() %>%
  mutate(visit_date = as.Date(visit_date))

specs <- list(
  list(name = "pa",   start_col = "pa_chronic_start_date",   status_col = "pa_status", label = "Chronic P. aeruginosa", tag = "A", sp_it = "P. aeruginosa"),
  list(name = "mssa", start_col = "mssa_chronic_start_date", status_col = "sa_status", label = "Chronic MSSA",          tag = "B", sp_it = NULL))
col_model <- "#1F4E79"; col_obs <- "grey15"; col_cf <- "grey45"
panels <- list(); bins_out <- list()
for (sp in specs) {
  cls <- dat %>% mutate(onset = as.Date(.data[[sp$start_col]])) %>% group_by(patient_id) %>%
    summarise(onset = dplyr::first(onset),
              n_pre_nonchronic = sum(!is.na(onset) & visit_date < onset & !is.na(.data[[sp$status_col]]) & .data[[sp$status_col]] != "Chronic"),
              .groups = "drop") %>%
    filter(!is.na(onset), n_pre_nonchronic >= 1)
  d <- dat %>% filter(!is.na(fev1_pp)) %>% inner_join(cls %>% select(patient_id, onset), by = "patient_id") %>%
    mutate(t = as.numeric(visit_date - onset) / 365.25)
  anchor <- d %>% filter(t < 0) %>% group_by(patient_id) %>% slice_max(visit_date, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(patient_id, anchor_date = visit_date, anchor_fev1 = fev1_pp, t_anchor = t)
  obs <- d %>% inner_join(anchor, by = "patient_id") %>% filter(visit_date != anchor_date, t >= -3, t < 3) %>%
    mutate(delta = fev1_pp - anchor_fev1, bin = floor(t / 0.5) * 0.5 + 0.25)
  binned <- obs %>% group_by(bin) %>% summarise(n = n_distinct(patient_id), mean = mean(delta), se = sd(delta) / sqrt(n()), .groups = "drop") %>%
    filter(n >= 10) %>% mutate(lo = mean - qt(0.975, n - 1) * se, hi = mean + qt(0.975, n - 1) * se)
  t_a <- median(anchor$t_anchor)
  message(sp$label, ": incident patients ", nrow(cls), ", with a pre-onset FEV1 ", nrow(anchor),
          "; median time of last pre-onset FEV1 ", round(t_a, 2), " years")
  m <- readRDS(file.path("outputs/revision/models", paste0("E2_", sp$name, "_incident.rds")))
  fe <- lme4::fixef(m); V <- as.matrix(vcov(m))
  tv <- paste0("time_", sp$name, "_index_cs"); pv <- paste0(sp$name, "_post_index"); iv <- paste0(tv, ":", pv)
  if (!iv %in% names(fe)) iv <- paste0(pv, ":", tv)
  pred_line <- function(tt, post, counterfactual = FALSE) {
    L <- matrix(0, nrow = length(tt), ncol = length(fe), dimnames = list(NULL, names(fe)))
    L[, tv] <- tt - t_a
    if (post && !counterfactual) { L[, pv] <- 1; L[, iv] <- tt }
    est <- as.vector(L %*% fe); se <- sqrt(rowSums((L %*% V) * L))
    tibble(t = tt, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }
  pre  <- pred_line(seq(-3, 0, by = 0.02), post = FALSE) %>% mutate(segment = "pre")
  post <- pred_line(seq(0, 3, by = 0.02), post = TRUE) %>% mutate(segment = "post")
  cf   <- pred_line(seq(0, 3, by = 0.02), post = TRUE, counterfactual = TRUE)
  y_lo <- -21
  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5, linetype = "dashed") +
    geom_ribbon(data = bind_rows(pre, post), aes(x = t, ymin = lo, ymax = hi, group = segment), fill = col_model, alpha = 0.15) +
    geom_line(data = cf, aes(x = t, y = est, linetype = "Pre-onset trend continued"), colour = col_cf, linewidth = 0.7) +
    geom_line(data = bind_rows(pre, post), aes(x = t, y = est, group = segment, linetype = "Model estimate (95% CI)"), colour = col_model, linewidth = 0.9) +
    geom_errorbar(data = binned, aes(x = bin, ymin = lo, ymax = hi), width = 0.08, colour = col_obs, linewidth = 0.5) +
    geom_point(data = binned, aes(x = bin, y = mean, shape = "Observed mean (95% CI)"), colour = col_obs, size = 2.2) +
    geom_text(data = binned, aes(x = bin, y = y_lo, label = n), size = 2.8, colour = "grey35") +
    annotate("text", x = -3.05, y = y_lo, label = "n", hjust = 1, size = 2.8, colour = "grey35", fontface = "italic") +
    scale_linetype_manual(NULL, values = c("Model estimate (95% CI)" = "solid", "Pre-onset trend continued" = "dashed")) +
    scale_shape_manual(NULL, values = c("Observed mean (95% CI)" = 16)) +
    scale_x_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.2, 3.1)) +
    coord_cartesian(ylim = c(y_lo - 1, 14)) +
    labs(title = if (!is.null(sp$sp_it)) bquote(bold(.(paste0(sp$tag, "    Chronic"))) ~ bolditalic(.(sp$sp_it)) ~ bold(.(paste0("(", nrow(anchor), " patients)"))))
                 else bquote(bold(.(paste0(sp$tag, "    ", sp$label, " (", nrow(anchor), " patients)")))),
         x = "Years from onset of chronic infection",
         y = expression(atop("Change in" ~ FEV[1] ~ "% predicted", "from last test before onset (points)"))) +
    theme_classic(base_size = 10.5) +
    theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "left",
          legend.key.width = unit(1.6, "lines"), legend.margin = margin(0, 0, 0, 0),
          plot.title = element_text(size = 10.5, face = "bold"), plot.title.position = "plot",
          axis.title.y = element_text(size = 9.5), plot.margin = margin(4, 10, 6, 6)) +
    guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2))
  if (sp$tag == "B") p <- p + theme(legend.position = "none")
  panels[[sp$name]] <- p
  bins_out[[sp$name]] <- binned %>% mutate(organism = sp$label, patients = nrow(anchor), t_anchor_median = t_a, .before = 1)
}
g <- rbind(ggplotGrob(panels$pa), ggplotGrob(panels$mssa), size = "max")
w <- 7; h <- 8.4
if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(file.path(fig_dir, "FigureS2_infection_onset.png"), width = w, height = h, units = "in", res = 300) else
  png(file.path(fig_dir, "FigureS2_infection_onset.png"), width = w, height = h, units = "in", res = 300)
grid.newpage(); grid.draw(g); invisible(dev.off())
pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_dev(file.path(fig_dir, "FigureS2_infection_onset.pdf"), width = w, height = h); grid.newpage(); grid.draw(g); invisible(dev.off())
write_csv(bind_rows(bins_out), file.path(fig_dir, "FigureS2_observed_bins.csv"))
message("Saved FigureS2_infection_onset")
