#!/usr/bin/env Rscript
# Main figure: proportion of transmitted rifampicin resistance across sampling
# fractions, with the observed estimate carrying the bootstrap confidence
# interval and the extrapolated estimate carrying the GAM confidence interval.
#
# Usage:  SIMEX_ROOT=/path/to/project Rscript 09_figure.R
# Writes: output/results/figure_simex.pdf and .png
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(mgcv); library(scales)
})
script_dir <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(script_dir(), "common.R"))

# On macOS, cairo_pdf() requires XQuartz and the base pdf() device cannot
# encode the typographic minus sign, so text is rendered through showtext when
# it is available.
if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 72)
}

SCENARIO <- "primary (denominator = resistant cases)"
COLOURS <- c("N" = "steelblue", "N-1" = "darkgreen")

simulations <- read_csv(file.path(OUTPUT, "results", "simulations.csv"),
                        show_col_types = FALSE) %>%
  filter(scenario == SCENARIO) %>%
  mutate(transmitted_N  = 1 - acquired_N  / trials,
         transmitted_N1 = 1 - acquired_N1 / trials)

observed <- filter(simulations, run_id == "observed")

fit_curve <- function(column, label) {
  d <- mutate(simulations, successes = trials - .data[[column]])
  model <- gam(cbind(successes, trials - successes) ~ s(lambda, bs = "tp", k = 10),
               family = binomial, data = d, method = "REML")
  grid <- tibble(lambda = seq(min(simulations$lambda), 1, length.out = 200))
  prediction <- predict(model, grid, type = "link", se.fit = TRUE)
  # Kept on the link scale under distinct names: assigning a column called
  # "fit" inside mutate() would shadow this vector for the lines that follow.
  link <- as.numeric(prediction$fit)
  link_se <- as.numeric(prediction$se.fit)
  last <- length(link)
  list(curve = mutate(grid,
                      fit   = plogis(link),
                      lower = plogis(link - 1.96 * link_se),
                      upper = plogis(link + 1.96 * link_se),
                      method = label),
       estimate = plogis(link[last]),
       lower = plogis(link[last] - 1.96 * link_se[last]),
       upper = plogis(link[last] + 1.96 * link_se[last]))
}
curve_N  <- fit_curve("acquired_N",  "N")
curve_N1 <- fit_curve("acquired_N1", "N-1")

curves <- bind_rows(curve_N$curve, curve_N1$curve) %>%
  mutate(method = factor(method, levels = c("N", "N-1")))

points_by_lambda <- simulations %>%
  pivot_longer(c(transmitted_N, transmitted_N1),
               names_to = "method", values_to = "transmitted") %>%
  mutate(method = factor(ifelse(method == "transmitted_N", "N", "N-1"),
                         levels = c("N", "N-1"))) %>%
  group_by(lambda, method) %>%
  summarise(mean_transmitted = mean(transmitted, na.rm = TRUE), .groups = "drop")

bootstrap <- read_csv(file.path(OUTPUT, "results", "bootstrap.csv"),
                      show_col_types = FALSE)
observed_points <- tibble(
  method = factor(c("N", "N-1"), levels = c("N", "N-1")),
  value = c(observed$transmitted_N[1], observed$transmitted_N1[1]),
  lower = c(quantile(bootstrap$transmitted_N,  0.025),
            quantile(bootstrap$transmitted_N1, 0.025)) / 100,
  upper = c(quantile(bootstrap$transmitted_N,  0.975),
            quantile(bootstrap$transmitted_N1, 0.975)) / 100,
  lambda = OBS_LAMBDA + c(-0.012, 0.012))

extrapolated_points <- tibble(
  method = factor(c("N", "N-1"), levels = c("N", "N-1")),
  value = c(curve_N$estimate, curve_N1$estimate),
  lower = c(curve_N$lower, curve_N1$lower),
  upper = c(curve_N$upper, curve_N1$upper),
  lambda = 1.0 + c(-0.012, 0.012))

figure <- ggplot() +
  geom_ribbon(data = curves,
              aes(x = lambda, ymin = lower, ymax = upper, fill = method),
              alpha = 0.25) +
  geom_line(data = curves, aes(x = lambda, y = fit, colour = method),
            linewidth = 1) +
  geom_point(data = points_by_lambda,
             aes(x = lambda, y = mean_transmitted, colour = method),
             size = 2, alpha = 0.75) +
  geom_vline(xintercept = OBS_LAMBDA, linetype = "dashed", colour = "gray40") +
  annotate("text", x = OBS_LAMBDA, y = 1.0,
           label = paste0("Observed (", round(OBS_LAMBDA * 100), "%)"),
           hjust = 1.1, vjust = 1, size = 3, colour = "gray40") +
  geom_vline(xintercept = 1.0, linetype = "dashed", colour = "firebrick") +
  annotate("text", x = 1.0, y = 1.0, label = "Extrapolated",
           hjust = 1.1, vjust = 1, size = 3, colour = "firebrick") +
  geom_errorbar(data = observed_points,
                aes(x = lambda, ymin = lower, ymax = upper),
                width = 0.016, linewidth = 0.9, colour = "gray25") +
  geom_point(data = observed_points,
             aes(x = lambda, y = value, shape = method),
             size = 4, colour = "gray25") +
  geom_errorbar(data = extrapolated_points,
                aes(x = lambda, ymin = lower, ymax = upper),
                width = 0.016, linewidth = 0.9, colour = "black") +
  geom_point(data = extrapolated_points,
             aes(x = lambda, y = value, shape = method),
             size = 4, colour = "black") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.6, 1)) +
  scale_colour_manual(name = "Method", values = COLOURS) +
  scale_fill_manual(name = "Method", values = COLOURS) +
  scale_shape_manual(name = "Method", values = c("N" = 16, "N-1" = 17)) +
  labs(x = "Sampling fraction (lambda)", y = "% Transmitted",
       title = "SIMEX Extrapolation: Transmitted vs Acquired Rifampicin Resistance",
       subtitle = paste("Observed estimate with bootstrap 95% CI;",
                        "extrapolated estimate with GAM 95% CI")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT, "results", "figure_simex.pdf"), figure,
       width = 10, height = 7)
ggsave(file.path(OUTPUT, "results", "figure_simex.png"), figure,
       width = 10, height = 7, dpi = 150)

cat(sprintf("observed     N = %.1f%% (bootstrap CI %.1f-%.1f)   N-1 = %.1f%% (%.1f-%.1f)\n",
            100 * observed_points$value[1], 100 * observed_points$lower[1],
            100 * observed_points$upper[1], 100 * observed_points$value[2],
            100 * observed_points$lower[2], 100 * observed_points$upper[2]))
cat(sprintf("extrapolated N = %.1f%% (%.1f-%.1f)   N-1 = %.1f%% (%.1f-%.1f)\n",
            100 * curve_N$estimate, 100 * curve_N$lower, 100 * curve_N$upper,
            100 * curve_N1$estimate, 100 * curve_N1$lower, 100 * curve_N1$upper))
