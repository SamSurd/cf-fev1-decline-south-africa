#!/usr/bin/env Rscript
# ==============================================================================
# A16 - REVISED FIGURE 1 (JCF-D-26-00377, corrected data)
# ------------------------------------------------------------------------------
# Forest plot in two panels, built only from saved model outputs:
#   A  level of FEV1pp: M2 level effects (fixed at entry; time-varying, same
#      registry year) and the step at each event from the event-anchored models
#   B  annual change: M2 population-averaged slope, M2 slope modifiers, and the
#      change in slope after each event from the event-anchored models
# Event-anchored sources: C2 (CFTR modulators), E1 (chronic P. aeruginosa and
# MSSA, all patients), adjusted exacerbation model (A15, first year with >=2).
# Inputs: outputs/tables/Table2A*, Table2B*, TableS1*; outputs/revision/fix1*,
#         fix2*, fix34*, row18_exacerbations.csv. Run after A08 and A15.
# Writes: outputs/revision/figures/Figure1_determinants.{png,pdf} and
#         outputs/revision/figure1_values.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(grid) })
fig_dir <- "outputs/revision/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read_csv(f, show_col_types = FALSE)
t2a <- rd("outputs/tables/Table2A_Baseline_FEV1pp_Determinants.csv")
t2b <- rd("outputs/tables/Table2B_FEV1_Decline_Rate_Determinants.csv")
ts1 <- rd("outputs/tables/TableS1_Nonsignificant_Baseline_Slope_Modifiers.csv")
fx1 <- rd("outputs/revision/fix1_cftr_models.csv")
fx2 <- rd("outputs/revision/fix2_global_slope.csv")
f34 <- rd("outputs/revision/fix34_event_models.csv")
r18 <- rd("outputs/revision/row18_exacerbations.csv")
slopes <- bind_rows(t2b, ts1)

m2_term <- function(tab, term) {
  r <- tab %>% filter(.data$term == !!term)
  if (nrow(r) != 1) stop("Term not found once: ", term)
  tibble(estimate = r$estimate, lo = r$conf.low, hi = r$conf.high, p = r$p.value)
}
rev_row <- function(tab, model, contrast) {
  r <- tab %>% filter(.data$model == !!model, .data$contrast == !!contrast)
  if (nrow(r) != 1) stop("Row not found once: ", model, " / ", contrast)
  tibble(estimate = r$estimate, lo = r$conf.low, hi = r$conf.high, p = r$p.value)
}
# Labels are plotmath strings so that headers are bold and species names italic.
hdr <- function(label) tibble(label = label, pm = paste0('bold("', label, '")'), header = TRUE, source = NA_character_,
                              estimate = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
itm <- function(label, src, est, pm = NULL) {
  tibble(label = label, pm = if (is.null(pm)) paste0('"   ', label, '"') else pm, header = FALSE, source = src) %>% bind_cols(est)
}
pm_it <- function(before, species, after) paste0('paste("   ', before, '", italic("', species, '"), "', after, '")')
ge2 <- "≥"
lv <- "Level change at event (pp)"; sl <- "Slope change after event (pp/yr)"
exac2 <- "First year with >=2 exacerbations (adjusted)"

panel_a <- bind_rows(
  hdr("Fixed at entry"),
  itm("Age at entry (per 10 years)", "M2", m2_term(t2a, "baseline_age_years_c")),
  itm("Male (vs female)", "M2", m2_term(t2a, "sexMale")),
  itm("Caucasian (vs other)", "M2", m2_term(t2a, "ethnicity_caucasian")),
  itm("Public sector only (vs private/combined)", "M2", m2_term(t2a, "sector_public_only")),
  itm("F508del, 1 or 2 alleles (vs none)", "M2", m2_term(t2a, "f508_any")),
  itm("P. aeruginosa at entry", "M2", m2_term(t2a, "baseline_pa_any"), pm_it("", "P. aeruginosa", " at entry")),
  hdr("Time-varying, same registry year"),
  itm("Underweight (vs normal)", "M2", m2_term(t2a, "nutrition_catUnderweight")),
  itm("Overweight (vs normal)", "M2", m2_term(t2a, "nutrition_catOverweight")),
  itm("CFRD", "M2", m2_term(t2a, "cfrd_indicator")),
  itm("Pancreatic insufficiency", "M2", m2_term(t2a, "pancreatic_insufficiency")),
  itm("Chronic MSSA", "M2", m2_term(t2a, "chronic_mssa")),
  itm("Chronic MRSA", "M2", m2_term(t2a, "chronic_mrsa")),
  itm("Chronic B. cepacia complex", "M2", m2_term(t2a, "chronic_bcepacia"), pm_it("Chronic ", "B. cepacia", " complex")),
  itm("1 exacerbation in the year (vs none)", "M2", m2_term(t2a, "exac_cat1")),
  itm(paste0(ge2, "2 exacerbations in the year (vs none)"), "M2", m2_term(t2a, "exac_cat2+")),
  hdr("At the event (event-anchored models)"),
  itm("CFTR modulator initiation", "event", rev_row(fx1, "C2 corrected primary", lv)),
  itm("Chronic P. aeruginosa onset", "event", rev_row(f34, "P. aeruginosa E1 all events", lv), pm_it("Chronic ", "P. aeruginosa", " onset")),
  itm("Chronic MSSA onset", "event", rev_row(f34, "MSSA E1 all events", lv)),
  itm(paste0("First year with ", ge2, "2 exacerbations"), "event", rev_row(r18, exac2, lv))
)
panel_b <- bind_rows(
  hdr("Adjusted annual change"),
  itm("All patients (population-averaged)", "ref", rev_row(fx2, "M2 interactions (REML)", "M2 population-averaged slope (visit-weighted)")),
  hdr("Fixed at entry"),
  itm("Male (vs female)", "M2", m2_term(slopes, "time_cs:sexMale")),
  itm("Caucasian (vs other)", "M2", m2_term(slopes, "time_cs:ethnicity_caucasian")),
  itm("Public sector only (vs private/combined)", "M2", m2_term(slopes, "time_cs:sector_public_only")),
  itm("F508del, 1 or 2 alleles (vs none)", "M2", m2_term(slopes, "time_cs:f508_any")),
  itm("P. aeruginosa at entry", "M2", m2_term(slopes, "time_cs:baseline_pa_any"), pm_it("", "P. aeruginosa", " at entry")),
  hdr("Time-varying"),
  itm("Underweight (vs normal)", "M2", m2_term(slopes, "time_cs:nutrition_catUnderweight")),
  itm("Overweight (vs normal)", "M2", m2_term(slopes, "time_cs:nutrition_catOverweight")),
  itm("CFRD", "M2", m2_term(slopes, "time_cs:cfrd_indicator_slope")),
  itm("Pancreatic insufficiency", "M2", m2_term(slopes, "time_cs:pancreatic_insufficiency_slope")),
  hdr("After the event (event-anchored models)"),
  itm("CFTR modulator initiation", "event", rev_row(fx1, "C2 corrected primary", sl)),
  itm("Chronic P. aeruginosa onset", "event", rev_row(f34, "P. aeruginosa E1 all events", sl), pm_it("Chronic ", "P. aeruginosa", " onset")),
  itm("Chronic MSSA onset", "event", rev_row(f34, "MSSA E1 all events", sl)),
  itm(paste0("First year with ", ge2, "2 exacerbations"), "event", rev_row(r18, exac2, sl))
)

fmt <- function(x) { x <- ifelse(abs(x) < 0.005, 0, x); sub("^-", "−", sprintf("%.2f", x)) }
prep <- function(d) {
  d %>% mutate(y = rev(seq_len(n())),
               ci = ifelse(header, "", paste0(fmt(estimate), " (", fmt(lo), " to ", fmt(hi), ")")),
               excl0 = !header & (lo > 0 | hi < 0),
               col = dplyr::case_when(header ~ "grey10", source == "event" ~ "#1F4E79", TRUE ~ "grey20"))
}
pa <- prep(panel_a); pb <- prep(panel_b)
pa$ci[1] <- "Estimate (95% CI)"; pb$ci[1] <- "Estimate (95% CI)"
write_csv(bind_rows(pa %>% mutate(panel = "A"), pb %>% mutate(panel = "B")) %>% filter(!header) %>%
            transmute(panel, label, source, estimate, conf.low = lo, conf.high = hi, p.value = p),
          "outputs/revision/figure1_values.csv")

col_event <- "#1F4E79"
forest <- function(d, title, xlab, xlim, xbreaks) {
  ax <- d %>% arrange(y)
  items <- d %>% filter(!header)
  ggplot(items) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = source), linewidth = 0.55) +
    geom_point(aes(x = estimate, y = y, colour = source, fill = interaction(source, excl0), shape = source), size = 2.1, stroke = 0.7) +
    scale_colour_manual(values = c(M2 = "grey20", event = col_event, ref = "black"), guide = "none") +
    scale_shape_manual(values = c(M2 = 21, event = 21, ref = 23), guide = "none") +
    scale_fill_manual(values = c(M2.TRUE = "grey20", M2.FALSE = "white", event.TRUE = col_event, event.FALSE = "white",
                                 ref.TRUE = "black", ref.FALSE = "white"), guide = "none") +
    scale_x_continuous(limits = xlim, breaks = xbreaks, labels = function(b) sub("^-", "−", b), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.4, nrow(d) + 0.6), expand = c(0, 0), breaks = ax$y, labels = parse(text = ax$pm),
                       sec.axis = dup_axis(labels = ax$ci, name = NULL)) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 8.5) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y.left  = element_text(hjust = 0, colour = ax$col, size = 8),
          axis.text.y.right = element_text(hjust = 0, colour = ax$col, face = ifelse(ax$header, "bold", "plain"), size = 8),
          axis.text.x = element_text(colour = "grey20"), axis.title.x = element_text(size = 8.5, margin = margin(t = 4)),
          plot.title = element_text(size = 9.5, margin = margin(b = 6)), plot.title.position = "plot",
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          plot.margin = margin(6, 8, 6, 4))
}
p_a <- forest(pa, expression(bold("A    Level of") ~ bold(FEV[1]) ~ bold("% predicted")),
              expression("Difference in" ~ FEV[1] ~ "% predicted (percentage points)"), c(-9, 15), seq(-5, 15, 5))
p_b <- forest(pb, expression(bold("B    Annual change in") ~ bold(FEV[1]) ~ bold("% predicted")),
              "Annual change, or difference in annual change (percentage points per year)", c(-2.5, 5), seq(-2, 4, 2))
ga <- suppressWarnings(ggplotGrob(p_a)); gb <- suppressWarnings(ggplotGrob(p_b))
g <- rbind(ga, gb, size = "max")
pan <- which(g$layout$name == "panel")
g$heights[g$layout$t[pan[1]]] <- unit(nrow(pa), "null")
g$heights[g$layout$t[pan[2]]] <- unit(nrow(pb), "null")
w <- 7.2; h <- 9.6
ragg::agg_png(file.path(fig_dir, "Figure1_determinants.png"), width = w, height = h, units = "in", res = 300)
grid.newpage(); grid.draw(g); invisible(dev.off())
grDevices::cairo_pdf(file.path(fig_dir, "Figure1_determinants.pdf"), width = w, height = h)
grid.newpage(); grid.draw(g); invisible(dev.off())
message("Saved Figure1_determinants (", nrow(pa) + nrow(pb), " lines)")
