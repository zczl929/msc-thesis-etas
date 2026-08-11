#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

output_dir <- "results/submission_v1/figures/thesis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

colours <- c("Naive" = "#666666", "Plug-in" = "#0072B2")
chain_colours <- c(
  "Chain 1" = "#0072B2", "Chain 2" = "#D55E00",
  "Chain 3" = "#009E73", "Chain 4" = "#CC79A7"
)
parameter_levels <- c("mu", "K", "alpha", "c", "p")
parameter_labels <- c(
  mu = "mu", K = "K", alpha = "alpha", c = "c", p = "p"
)

theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(colour = "grey30")
    )
}

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 320, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = grDevices::pdf,
    useDingbats = FALSE, bg = "white"
  )
}

as_chain_long <- function(fit, model, keep_every = 1L) {
  bind_rows(lapply(seq_along(fit$chains), function(chain_id) {
    # Explicit matrix conversion preserves the MCMC parameter names even when
    # tibble's as.data.frame methods are attached.
    x <- as.data.frame(unclass(fit$chains[[chain_id]]), check.names = FALSE)
    index <- seq(1L, nrow(x), by = keep_every)
    missing_parameters <- setdiff(parameter_levels, names(x))
    if (length(missing_parameters) > 0L) {
      stop(
        "Missing MCMC columns for ", model, " chain ", chain_id, ": ",
        paste(missing_parameters, collapse = ", "),
        "; available columns: ", paste(names(x), collapse = ", "),
        call. = FALSE
      )
    }
    x <- x[index, match(parameter_levels, names(x)), drop = FALSE]
    names(x) <- parameter_levels
    x$draw <- index
    x$chain <- paste("Chain", chain_id)
    pivot_longer(
      x, cols = all_of(parameter_levels),
      names_to = "parameter", values_to = "value"
    )
  })) |>
    mutate(
      model = model,
      chain = factor(chain, levels = names(chain_colours)),
      parameter = factor(parameter, levels = parameter_levels)
    )
}

# -------------------------------------------------------------------------
# 0. One realised synthetic catalogue illustrating the observation process.
# -------------------------------------------------------------------------
synthetic_observed <- readRDS("data/submission_v1/observed_001.rds")
synthetic_audit <- attr(synthetic_observed, "detection_audit") |>
  mutate(
    days_after_mainshock = ts - 500,
    status = ifelse(detected, "Observed", "Missed"),
    status = factor(status, levels = c("Missed", "Observed"))
  ) |>
  filter(days_after_mainshock >= 0, days_after_mainshock <= 0.35)

synthetic_curve_time <- c(
  10^seq(-5, log10(0.35), length.out = 900), 0.35
)
synthetic_curve <- data.frame(
  days_after_mainshock = synthetic_curve_time,
  mc = pmax(
    2.5,
    6.5 - 4.5 - 0.75 * log10(synthetic_curve_time)
  )
)
synthetic_recovery <- 10^((6.5 - 4.5 - 2.5) / 0.75)

synthetic_observation_plot <- ggplot(
  synthetic_audit,
  aes(days_after_mainshock, magnitudes, colour = status, shape = status)
) +
  geom_point(alpha = 0.72, size = 1.8) +
  geom_line(
    data = synthetic_curve,
    aes(days_after_mainshock, mc), inherit.aes = FALSE,
    colour = "#B2182B", linewidth = 0.9
  ) +
  geom_vline(
    xintercept = synthetic_recovery, linetype = "dashed",
    colour = "#B2182B"
  ) +
  annotate(
    "text", x = synthetic_recovery, y = 6.15,
    label = sprintf("Return to baseline: %.3f d", synthetic_recovery),
    hjust = -0.04, colour = "#B2182B", size = 3.2
  ) +
  scale_colour_manual(values = c("Missed" = "#D55E00", "Observed" = "#0072B2")) +
  scale_shape_manual(values = c("Missed" = 4, "Observed" = 16)) +
  coord_cartesian(xlim = c(0, 0.35), ylim = c(2.4, 6.7)) +
  labs(
    x = "Days after the seeded M6.5 mainshock", y = "Magnitude",
    colour = NULL, shape = NULL,
    title = "Illustration of simulated short-term incompleteness",
    subtitle = paste(
      "Catalogue 1; red curve is the data-generating completeness threshold",
      "and crossed points are latent missed events"
    )
  ) +
  theme_thesis()

save_figure(
  synthetic_observation_plot,
  "00_synthetic_observation_mechanism", 8.2, 5.2
)
write.csv(
  synthetic_audit,
  file.path(output_dir, "00_synthetic_observation_mechanism_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 1. Synthetic MCMC diagnostics across all catalogues, models and parameters.
# -------------------------------------------------------------------------
diagnostics <- read.csv(
  "results/submission_v1/mcmc_primary/summary/mcmc_diagnostics.csv"
) |>
  mutate(
    model = recode(model, naive = "Naive", plugin = "Plug-in"),
    model = factor(model, levels = names(colours)),
    parameter = factor(parameter, levels = parameter_levels)
  )

rhat_plot <- ggplot(
  diagnostics,
  aes(parameter, rhat, colour = model)
) +
  geom_hline(yintercept = 1.05, linetype = "dashed", colour = "#B2182B") +
  geom_boxplot(
    aes(group = interaction(parameter, model)),
    position = position_dodge(width = 0.55), width = 0.44,
    outlier.shape = NA, linewidth = 0.45
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12, dodge.width = 0.55, seed = 20260810
    ),
    alpha = 0.22, size = 0.8
  ) +
  scale_colour_manual(values = colours) +
  labs(
    x = NULL, y = expression(hat(R)), colour = NULL,
    title = "(a) Between-chain convergence",
    subtitle = "Dashed line: pre-specified threshold"
  ) +
  theme_thesis()

ess_plot <- ggplot(
  diagnostics,
  aes(parameter, effective_size, colour = model)
) +
  geom_hline(yintercept = 400, linetype = "dashed", colour = "#B2182B") +
  geom_boxplot(
    aes(group = interaction(parameter, model)),
    position = position_dodge(width = 0.55), width = 0.44,
    outlier.shape = NA, linewidth = 0.45
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12, dodge.width = 0.55, seed = 20260810
    ),
    alpha = 0.22, size = 0.8
  ) +
  scale_colour_manual(values = colours) +
  scale_y_log10() +
  labs(
    x = "ETAS parameter", y = "Effective sample size (log scale)",
    colour = NULL, title = "(b) Effective posterior information",
    subtitle = "Dashed line: pre-specified threshold"
  ) +
  theme_thesis()

synthetic_diagnostics_plot <- rhat_plot / ess_plot +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Synthetic-study MCMC diagnostics",
    subtitle = paste(
      "All 1,000 parameter-specific diagnostics from 100 paired catalogues;",
      "max R-hat = 1.038 and min ESS = 432"
    )
  ) &
  theme(legend.position = "top")

save_figure(
  synthetic_diagnostics_plot,
  "01_synthetic_mcmc_diagnostics", 8.2, 7.1
)
write.csv(
  diagnostics,
  file.path(output_dir, "01_synthetic_mcmc_diagnostics_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 2. Trace plot for the fit containing the largest R-hat.
# -------------------------------------------------------------------------
worst_row <- diagnostics[which.max(diagnostics$rhat), ]
worst_path <- sprintf(
  "results/submission_v1/mcmc_primary/%s_%03d.rds",
  tolower(as.character(worst_row$model)), worst_row$sim_id
)
worst_fit <- readRDS(worst_path)
worst_trace <- as_chain_long(
  worst_fit, as.character(worst_row$model), keep_every = 2L
)

worst_trace_plot <- ggplot(
  worst_trace,
  aes(draw, value, colour = chain)
) +
  geom_line(linewidth = 0.25, alpha = 0.72) +
  facet_wrap(~parameter, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = chain_colours) +
  labs(
    x = "Retained draw index", y = "Parameter value", colour = NULL,
    title = "Trace plots for the synthetic fit with the largest R-hat",
    subtitle = sprintf(
      "Catalogue %d, %s ETAS; largest R-hat = %.3f for %s",
      worst_row$sim_id, worst_row$model, worst_row$rhat,
      as.character(worst_row$parameter)
    )
  ) +
  theme_thesis()

save_figure(
  worst_trace_plot,
  "02_synthetic_worst_fit_trace", 8.2, 9.5
)

# -------------------------------------------------------------------------
# 3. Paired posterior-median errors across simulated catalogues.
# -------------------------------------------------------------------------
posterior <- read.csv(
  "results/submission_v1/mcmc_primary/summary/posterior_long.csv"
) |>
  mutate(
    model = recode(model, naive = "Naive", plugin = "Plug-in"),
    model = factor(model, levels = names(colours)),
    parameter = factor(parameter, levels = parameter_levels),
    error = median - true
  )

recovery_plot <- ggplot(
  posterior,
  aes(model, error, colour = model)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35") +
  geom_line(
    aes(group = sim_id), colour = "grey72", alpha = 0.28,
    linewidth = 0.3
  ) +
  geom_point(alpha = 0.45, size = 1.1, position = position_jitter(width = 0.035)) +
  geom_boxplot(
    width = 0.45, outlier.shape = NA, fill = "white",
    linewidth = 0.55
  ) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = colours) +
  labs(
    x = NULL, y = "Posterior median minus true value", colour = NULL,
    title = "Paired parameter-recovery errors",
    subtitle = paste(
      "Each line joins Naive and Plug-in fits to the same incomplete",
      "catalogue; dashed line indicates zero error"
    )
  ) +
  theme_thesis()

save_figure(recovery_plot, "03_synthetic_parameter_recovery", 8.2, 6.2)
write.csv(
  posterior,
  file.path(output_dir, "03_synthetic_parameter_recovery_data.csv"),
  row.names = FALSE
)

# Copy the already validated coverage figure into the thesis figure set.
coverage_source_png <- "results/submission_v1/figures/portfolio_coverage_visualisation.png"
coverage_source_pdf <- "results/submission_v1/figures/portfolio_coverage_visualisation.pdf"
file.copy(
  coverage_source_png,
  file.path(output_dir, "04_synthetic_empirical_coverage.png"),
  overwrite = TRUE
)
file.copy(
  coverage_source_pdf,
  file.path(output_dir, "04_synthetic_empirical_coverage.pdf"),
  overwrite = TRUE
)

# -------------------------------------------------------------------------
# 4. Ridgecrest catalogue, train/test split, and fitted completeness curve.
# -------------------------------------------------------------------------
source("config/submission_ridgecrest.R")
cfg <- submission_ridgecrest_config()
ridgecrest <- readRDS(cfg$full_path)
ridgecrest$period <- ifelse(
  ridgecrest$ts <= cfg$train_end, "Training", "Held-out test"
)
ridgecrest$period <- factor(
  ridgecrest$period, levels = c("Training", "Held-out test")
)
train <- readRDS(cfg$train_path)
mainshock <- train[train$ID == cfg$mainshock_id, , drop = FALSE]
gh <- readRDS(
  file.path(cfg$results_dir, "plugin_GH_single.rds")
)$estimate

positive_grid <- unique(c(
  10^seq(-5, log10(cfg$forecast_end), length.out = 1600),
  seq(0.001, cfg$forecast_end, length.out = 1600)
))
positive_grid <- sort(positive_grid[positive_grid <= cfg$forecast_end])
mc_curve <- data.frame(
  ts = positive_grid,
  mc = pmax(
    cfg$M0,
    mainshock$magnitudes[[1]] - gh[["G"]] - gh[["H"]] * log10(positive_grid)
  )
)
recovery_end <- 10^(
  (mainshock$magnitudes[[1]] - gh[["G"]] - cfg$M0) / gh[["H"]]
)

catalogue_full <- ggplot(
  ridgecrest,
  aes(ts, magnitudes, colour = period)
) +
  geom_point(alpha = 0.62, size = 1.15) +
  geom_line(
    data = mc_curve, aes(ts, mc), inherit.aes = FALSE,
    colour = "#B2182B", linewidth = 0.9
  ) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "grey25") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey25") +
  scale_colour_manual(values = c("Training" = "#666666", "Held-out test" = "#0072B2")) +
  coord_cartesian(xlim = range(ridgecrest$ts), ylim = c(2.85, 7.25)) +
  labs(
    x = "Days relative to the M7.1 mainshock", y = "Magnitude",
    colour = NULL, title = "(a) Catalogue and held-out split"
  ) +
  theme_thesis()

catalogue_zoom <- ggplot(
  ridgecrest |> filter(ts >= 0, ts <= 0.06),
  aes(ts, magnitudes)
) +
  geom_point(alpha = 0.62, size = 1.25, colour = "#666666") +
  geom_line(
    data = mc_curve |> filter(ts <= 0.06),
    aes(ts, mc), colour = "#B2182B", linewidth = 0.9
  ) +
  geom_vline(
    xintercept = recovery_end, linetype = "dashed", colour = "#B2182B"
  ) +
  annotate(
    "text", x = recovery_end, y = 6.8,
    label = sprintf(
      "Return to baseline: %.3f days (%.0f min)",
      recovery_end, recovery_end * 24 * 60
    ),
    hjust = -0.05, size = 3.2, colour = "#B2182B"
  ) +
  coord_cartesian(xlim = c(0, 0.06), ylim = c(2.85, 7.25)) +
  labs(
    x = "Days after mainshock", y = "Magnitude",
    title = "(b) Early-period completeness recovery"
  ) +
  theme_thesis()

catalogue_plot <- catalogue_full / catalogue_zoom +
  plot_layout(guides = "collect", heights = c(0.85, 1.15)) +
  plot_annotation(
    subtitle = sprintf(
      "Fitted completeness curve: G = %.3f and H = %.3f",
      gh[["G"]], gh[["H"]]
    )
  ) &
  theme(legend.position = "top")

save_figure(catalogue_plot, "05_ridgecrest_catalogue_completeness", 8.2, 5.8)
write.csv(
  mc_curve,
  file.path(output_dir, "05_ridgecrest_completeness_curve.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 5. Ridgecrest MCMC traces and marginal posterior distributions.
# -------------------------------------------------------------------------
ridgecrest_fit_dir <- file.path(cfg$results_dir, "mcmc_conditioned_single")
ridgecrest_fits <- list(
  "Naive" = readRDS(file.path(ridgecrest_fit_dir, "naive.rds")),
  "Plug-in" = readRDS(file.path(ridgecrest_fit_dir, "plugin.rds"))
)
ridgecrest_chains <- bind_rows(
  as_chain_long(ridgecrest_fits[["Naive"]], "Naive", keep_every = 4L),
  as_chain_long(ridgecrest_fits[["Plug-in"]], "Plug-in", keep_every = 4L)
) |>
  mutate(model = factor(model, levels = names(colours)))

ridgecrest_trace_plot <- ggplot(
  ridgecrest_chains,
  aes(draw, value, colour = chain)
) +
  geom_line(linewidth = 0.22, alpha = 0.66) +
  facet_grid(parameter ~ model, scales = "free_y") +
  scale_colour_manual(values = chain_colours) +
  labs(
    x = "Retained draw index", y = "Parameter value", colour = NULL,
    title = "Ridgecrest MCMC trace plots",
    subtitle = "Four final chains per model; every fourth retained draw shown"
  ) +
  theme_thesis()

save_figure(ridgecrest_trace_plot, "06_ridgecrest_mcmc_trace", 8.2, 10.2)

ridgecrest_posterior <- ridgecrest_chains |>
  filter(draw %% 16L == 1L | draw == 1L)

ridgecrest_posterior_plot <- ggplot(
  ridgecrest_posterior,
  aes(value, colour = model, fill = model)
) +
  geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.05) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(
    x = "Parameter value", y = "Posterior density",
    colour = NULL, fill = NULL,
    title = "Ridgecrest marginal posterior distributions",
    subtitle = "Posterior location and uncertainty under Naive and Plug-in ETAS"
  ) +
  theme_thesis()

save_figure(
  ridgecrest_posterior_plot,
  "07_ridgecrest_posterior_marginals", 8.2, 6.2
)

# -------------------------------------------------------------------------
# 6. Ridgecrest branching-ratio stability.
# -------------------------------------------------------------------------
branching <- bind_rows(lapply(names(ridgecrest_fits), function(model) {
  draws <- as.data.frame(
    unclass(do.call(rbind, ridgecrest_fits[[model]]$chains)),
    check.names = FALSE
  )
  n <- ifelse(draws$alpha < cfg$b, draws$K * cfg$b / (cfg$b - draws$alpha), Inf)
  data.frame(
    model = model,
    n = n,
    n_plot = pmin(n, 5),
    supercritical = n >= 1,
    alpha_ge_b = draws$alpha >= cfg$b
  )
})) |>
  mutate(model = factor(model, levels = names(colours)))

branching_intervals <- branching |>
  group_by(model) |>
  summarise(
    median = median(n),
    q25 = quantile(n, 0.25),
    q75 = quantile(n, 0.75),
    .groups = "drop"
  )

branching_interval_plot <- ggplot(
  branching_intervals,
  aes(median, model, colour = model)
) +
  annotate(
    "rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = "#B2182B", alpha = 0.06
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#B2182B") +
  geom_errorbar(
    aes(xmin = q25, xmax = q75), orientation = "y",
    width = 0.12, linewidth = 0.8
  ) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = sprintf("median %.3f", median)),
    nudge_y = 0.22, size = 3.1, show.legend = FALSE
  ) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.14))) +
  labs(
    x = "Branching ratio (line shows interquartile interval)", y = NULL,
    colour = NULL, title = "(a) Posterior branching ratio"
  ) +
  guides(colour = "none") +
  theme_thesis()

branching_summary <- branching |>
  group_by(model) |>
  summarise(
    `P(n >= 1)` = mean(supercritical),
    `P(alpha >= b)` = mean(alpha_ge_b),
    .groups = "drop"
  ) |>
  pivot_longer(-model, names_to = "measure", values_to = "probability")

branching_bars <- ggplot(
  branching_summary,
  aes(probability, measure, colour = model)
) +
  geom_point(position = position_dodge(width = 0.45), size = 3.1) +
  geom_text(
    aes(label = sprintf("%.3f", probability)),
    position = position_dodge(width = 0.45), hjust = -0.35, size = 3.1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = colours) +
  scale_y_discrete(
    labels = c(
      "P(n >= 1)" = expression(P(n >= 1)),
      "P(alpha >= b)" = expression(P(alpha >= b))
    )
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Posterior probability", y = NULL, colour = NULL,
    title = "(b) Posterior instability probabilities"
  ) +
  theme_thesis()

branching_plot <- branching_interval_plot / branching_bars +
  plot_layout(guides = "collect") +
  plot_annotation() &
  theme(legend.position = "top")

save_figure(branching_plot, "08_ridgecrest_branching_stability", 8.2, 5.2)
write.csv(
  branching_summary,
  file.path(output_dir, "08_ridgecrest_branching_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 7. Ridgecrest posterior predictive event counts and upper tail.
# -------------------------------------------------------------------------
predictive <- readRDS(
  file.path(ridgecrest_fit_dir, "forecast", "forecast_draws.rds")
)
count_draws <- bind_rows(
  data.frame(
    model = "Naive", count = predictive$naive$counts,
    overflow = predictive$naive$overflow
  ),
  data.frame(
    model = "Plug-in", count = predictive$plugin$counts,
    overflow = predictive$plugin$overflow
  )
) |>
  mutate(model = factor(model, levels = names(colours)))

observed_count <- sum(
  ridgecrest$ts > cfg$train_end & ridgecrest$ts <= cfg$forecast_end
)

count_summary <- count_draws |>
  group_by(model) |>
  summarise(
    q025 = quantile(count, 0.025),
    q25 = quantile(count, 0.25),
    median = median(count),
    q75 = quantile(count, 0.75),
    q975 = quantile(count, 0.975),
    overflow_probability = mean(overflow),
    .groups = "drop"
  )

count_distribution <- ggplot(
  count_summary,
  aes(median, model, colour = model)
) +
  geom_vline(
    xintercept = observed_count, linetype = "dashed", colour = "#B2182B"
  ) +
  geom_errorbar(
    aes(xmin = q025, xmax = q975), orientation = "y",
    width = 0, linewidth = 0.75, alpha = 0.65
  ) +
  geom_errorbar(
    aes(xmin = q25, xmax = q75), orientation = "y",
    width = 0.14, linewidth = 2.2
  ) +
  geom_point(size = 3.4) +
  annotate(
    "text", x = observed_count, y = 2.38,
    label = paste("Observed =", observed_count), hjust = 0,
    vjust = 0, size = 3.1, colour = "#B2182B"
  ) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(limits = c(0, 1400), breaks = seq(0, 1400, 200)) +
  labs(
    x = "Test-period event count", y = NULL, colour = NULL,
    title = "(a) Posterior predictive count intervals",
    subtitle = "Thick and thin lines show 50% and 95% intervals"
  ) +
  guides(colour = "none") +
  theme_thesis()

count_tail <- ggplot(
  count_summary,
  aes(overflow_probability, model, colour = model)
) +
  geom_segment(aes(x = 0, xend = overflow_probability, yend = model), linewidth = 0.8) +
  geom_point(size = 3.4) +
  geom_text(
    aes(label = scales::percent(overflow_probability, accuracy = 0.1)),
    hjust = -0.45, size = 3.2
  ) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(
    limits = c(0, 0.022), breaks = c(0, 0.005, 0.01, 0.015, 0.02),
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  labs(
    x = "Probability of reaching the 2,000-event simulation limit",
    y = NULL, colour = NULL, title = "(b) Predictive overflow probability"
  ) +
  guides(colour = "none") +
  theme_thesis()

count_plot <- count_distribution / count_tail +
  plot_layout(guides = "collect", heights = c(1.2, 0.8)) +
  plot_annotation() &
  theme(legend.position = "top")

save_figure(count_plot, "09_ridgecrest_predictive_counts", 8.2, 5.5)
write.csv(
  count_draws,
  file.path(output_dir, "09_ridgecrest_predictive_count_draws.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8. Ridgecrest temporal predictive performance.
# -------------------------------------------------------------------------
forecast_summary <- read.csv(
  file.path(ridgecrest_fit_dir, "forecast", "forecast_summary.csv")
) |>
  transmute(
    model = factor(
      recode(model, naive = "Naive", plugin = "Plug-in"),
      levels = names(colours)
    ),
    lpd_per_event = log_score_per_event
  )
mc_error <- read.csv(
  file.path(ridgecrest_fit_dir, "forecast", "forecast_mc_error.csv")
)

ig_data <- data.frame(
  estimate = mc_error$information_gain_per_event,
  lower = mc_error$mc_q025,
  upper = mc_error$mc_q975
)
ig_panel <- ggplot(ig_data, aes(estimate, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +
  geom_errorbar(
    aes(xmin = lower, xmax = upper), orientation = "y", width = 0.10,
    colour = colours[["Plug-in"]], linewidth = 0.8
  ) +
  geom_point(colour = colours[["Plug-in"]], size = 3.6) +
  geom_text(
    aes(label = sprintf("%.4f", estimate)),
    nudge_y = 0.12, size = 3.4
  ) +
  scale_y_continuous(breaks = 1, labels = "Plug-in minus Naive") +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    x = "Information gain per test event", y = NULL,
    title = "Held-out temporal predictive performance",
    subtitle = "Interval shows posterior Monte Carlo integration error"
  ) +
  theme_thesis()

temporal_prediction_plot <- ig_panel

save_figure(
  temporal_prediction_plot,
  "10_ridgecrest_temporal_prediction", 6.8, 2.8
)
write.csv(
  forecast_summary,
  file.path(output_dir, "10_ridgecrest_temporal_prediction_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 8b. Compact Ridgecrest prediction summary for the dissertation body.
# -------------------------------------------------------------------------
prediction_lpd_panel <- ggplot(
  forecast_summary,
  aes(lpd_per_event, model, colour = model)
) +
  geom_segment(
    aes(
      x = min(forecast_summary$lpd_per_event) - 0.004,
      xend = lpd_per_event, yend = model
    ),
    linewidth = 0.8, alpha = 0.8
  ) +
  geom_point(size = 3.5) +
  geom_text(
    aes(label = sprintf("%.4f", lpd_per_event)),
    nudge_x = 0.0015, hjust = 0, size = 3.2, show.legend = FALSE
  ) +
  annotate(
    "text",
    x = mean(forecast_summary$lpd_per_event), y = 2.42,
    label = sprintf("Plug-in gain: +%.4f per event", ig_data$estimate),
    colour = colours[["Plug-in"]], size = 3.1
  ) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(
    limits = c(2.397, 2.435),
    breaks = c(2.40, 2.41, 2.42, 2.43)
  ) +
  labs(
    x = "Temporal log predictive density per test event (higher is better)",
    y = NULL, colour = NULL, title = "(a) Temporal prediction"
  ) +
  guides(colour = "none") +
  theme_thesis()

prediction_count_panel <- ggplot(
  count_summary,
  aes(median, model, colour = model)
) +
  geom_vline(
    xintercept = observed_count, linetype = "dashed", colour = "#B2182B"
  ) +
  geom_errorbar(
    aes(xmin = q025, xmax = q975), orientation = "y",
    width = 0.12, linewidth = 0.9
  ) +
  geom_point(size = 3.5) +
  geom_text(
    aes(
      label = sprintf(
        "median %.1f; overflow %.1f%%",
        median, 100 * overflow_probability
      )
    ),
    nudge_y = 0.20, size = 3.0, show.legend = FALSE
  ) +
  annotate(
    "text", x = observed_count, y = 2.42,
    label = sprintf("Observed = %d", observed_count),
    hjust = -0.06, colour = "#B2182B", size = 3.1
  ) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(limits = c(0, 1350), breaks = seq(0, 1200, 200)) +
  labs(
    x = "Held-out event count (line shows 95% predictive interval)",
    y = NULL, colour = NULL, title = "(b) Event-count prediction"
  ) +
  guides(colour = "none") +
  theme_thesis()

ridgecrest_prediction_summary <-
  prediction_lpd_panel / prediction_count_panel +
  plot_layout(heights = c(1, 1.05))

save_figure(
  ridgecrest_prediction_summary,
  "12_ridgecrest_prediction_summary", 8.2, 5.5
)

# -------------------------------------------------------------------------
# 9. Catalogue-specific plug-in completeness-parameter recovery.
# -------------------------------------------------------------------------
gh_estimates <- bind_rows(lapply(1:100, function(sim_id) {
  fit <- readRDS(sprintf(
    "results/submission_v1/initialisation/plugin_GH_%03d.rds", sim_id
  ))
  data.frame(
    sim_id = sim_id,
    parameter = c("G", "H"),
    estimate = unname(fit$estimate[c("G", "H")]),
    truth = c(4.5, 0.75)
  )
})) |>
  mutate(
    parameter = factor(parameter, levels = c("G", "H")),
    error = estimate - truth
  )

gh_recovery_plot <- ggplot(
  gh_estimates,
  aes(parameter, error, colour = parameter)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35") +
  geom_violin(alpha = 0.12, width = 0.75, linewidth = 0.6) +
  geom_boxplot(width = 0.20, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.08, alpha = 0.42, size = 1.1) +
  scale_colour_manual(values = c("G" = "#009E73", "H" = "#CC79A7")) +
  labs(
    x = NULL, y = "Estimate minus true value", colour = NULL,
    title = "Recovery of plug-in completeness parameters",
    subtitle = paste(
      "Across 100 catalogues: bias (RMSE) = 0.001 (0.135) for G",
      "and 0.026 (0.127) for H"
    )
  ) +
  guides(colour = "none") +
  theme_thesis()

save_figure(
  gh_recovery_plot,
  "11_synthetic_completeness_parameter_recovery", 7.2, 4.8
)
write.csv(
  gh_estimates,
  file.path(output_dir, "11_synthetic_completeness_parameter_recovery_data.csv"),
  row.names = FALSE
)

message("Wrote thesis figures to ", output_dir)
