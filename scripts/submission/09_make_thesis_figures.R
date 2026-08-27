#!/usr/bin/env Rscript

# Recreate the five dissertation-body figures and two appendix figures from
# the compact submission.
# The script uses frozen source data and authoritative plotting CSV files. It
# does not require retained MCMC objects; Figure 5 uses exported predictive
# count draws so that the event-count distribution is shown directly.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

source("config/submission_experiment.R")
source("config/submission_ridgecrest.R")
source("R/etas_likelihood.R")
source("R/simulate_etas.R")
source("R/loglinear_incompleteness.R")

output_dir <- "results/submission_v1/figures/thesis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_colours <- c("Naive" = "#666666", "Plug-in" = "#0072B2")
benchmark_model_colours <- c(
  "Naive" = "#666666",
  "Estimated Plug-in" = "#0072B2",
  "Oracle" = "#D55E00",
  "Complete-Data" = "#009E73"
)
mechanism_colours <- c(
  "Observed" = "#111111",
  "Missed" = "#D73027",
  "Data-generating completeness threshold" = "#B2182B",
  "Baseline threshold" = "#4D4D4D"
)
red <- "#B2182B"
blue <- "#0072B2"
orange <- "#D55E00"

theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = rel(1.03)),
      plot.title = element_text(face = "bold", size = rel(1.02)),
      plot.subtitle = element_text(colour = "grey30"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, device = grDevices::pdf,
    useDingbats = FALSE, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 320, bg = "white"
  )
}

# -------------------------------------------------------------------------
# Figure 1. Representative synthetic incomplete catalogue.
# -------------------------------------------------------------------------
synthetic_cfg <- submission_experiment_config()
synthetic_complete <- simulate_etas(
  theta = as.list(synthetic_cfg$theta_true),
  T1 = synthetic_cfg$T1,
  T2 = synthetic_cfg$T2,
  M_cut = synthetic_cfg$M0,
  beta = synthetic_cfg$beta,
  Ht = synthetic_cfg$mainshock,
  seed = synthetic_cfg$seed_start + 1L,
  max_events = synthetic_cfg$max_events
)
synthetic_observed <- apply_loglinear_incompleteness(
  synthetic_complete,
  M0 = synthetic_cfg$M0,
  G = synthetic_cfg$incompleteness$G,
  H = synthetic_cfg$incompleteness$H,
  b = synthetic_cfg$b,
  trigger_mode = synthetic_cfg$incompleteness$trigger_mode,
  M_trigger = synthetic_cfg$incompleteness$M_trigger,
  cap_at_trigger = synthetic_cfg$incompleteness$cap_at_trigger,
  retain_seeded = TRUE
)

synthetic_audit_full <- attr(synthetic_observed, "detection_audit") |>
  mutate(
    days_after_mainshock = ts - synthetic_cfg$mainshock$ts[[1]],
    status = ifelse(detected, "Observed", "Missed"),
    status = factor(status, levels = c("Observed", "Missed"))
  )

synthetic_audit <- synthetic_audit_full |>
  filter(days_after_mainshock > 0, days_after_mainshock <= 0.35)

synthetic_curve_time <- 10^seq(-5, log10(0.35), length.out = 900)
synthetic_curve <- data.frame(
  days_after_mainshock = synthetic_curve_time,
  mc = pmax(
    synthetic_cfg$M0,
    synthetic_cfg$mainshock$magnitudes[[1]] -
      synthetic_cfg$incompleteness$G -
      synthetic_cfg$incompleteness$H * log10(synthetic_curve_time)
  )
) |>
  mutate(
    psi = 10^(-synthetic_cfg$b * pmax(mc - synthetic_cfg$M0, 0))
  )

synthetic_lines <- bind_rows(
  synthetic_curve |>
    transmute(
      days_after_mainshock, value = mc,
      series = "Data-generating completeness threshold"
    ),
  data.frame(
    days_after_mainshock = c(0, 0.35),
    value = synthetic_cfg$M0,
    series = "Baseline threshold"
  )
)

synthetic_a <- ggplot(
  synthetic_audit_full, aes(ts, magnitudes)
) +
  geom_hline(
    yintercept = synthetic_cfg$M0, linetype = "dashed",
    colour = "grey40"
  ) +
  geom_point(colour = "#111111", shape = 16, alpha = 0.84, size = 1.15) +
  geom_vline(
    xintercept = synthetic_cfg$mainshock$ts[[1]],
    linetype = "dotted", colour = "grey30", linewidth = 0.75
  ) +
  annotate(
    "point", x = synthetic_cfg$mainshock$ts[[1]],
    y = synthetic_cfg$mainshock$magnitudes[[1]],
    shape = 16, size = 3.8, colour = "black"
  ) +
  annotate(
    "text", x = synthetic_cfg$mainshock$ts[[1]] + 24,
    y = synthetic_cfg$mainshock$magnitudes[[1]],
    label = "Seeded M6.5 mainshock", hjust = 0, size = 3.0
  ) +
  coord_cartesian(
    xlim = c(synthetic_cfg$T1, synthetic_cfg$T2),
    ylim = c(2.4, 6.8)
  ) +
  labs(
    x = "Simulation time (days)", y = "Magnitude",
    title = "(a) Representative complete simulated ETAS catalogue"
  ) +
  theme_thesis() +
  theme(
    legend.position = "none",
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    )
  )

synthetic_b <- ggplot() +
  geom_line(
    data = synthetic_lines,
    aes(days_after_mainshock, value, colour = series, linetype = series),
    linewidth = 1.0
  ) +
  geom_point(
    data = synthetic_audit,
    aes(
      days_after_mainshock, magnitudes,
      colour = status, shape = status
    ),
    alpha = 0.82, size = 2.1, stroke = 0.9
  ) +
  annotate(
    "point", x = 0, y = synthetic_cfg$mainshock$magnitudes[[1]],
    shape = 16, size = 3.8, colour = "black"
  ) +
  annotate(
    "text", x = 0.008, y = 6.43, label = "Seeded mainshock",
    hjust = 0, size = 3.1
  ) +
  scale_colour_manual(
    values = mechanism_colours,
    breaks = c(
      "Observed", "Missed", "Data-generating completeness threshold",
      "Baseline threshold"
    ),
    labels = parse(text = c(
      "'Observed'", "'Missed'",
      "'Data-generating'~m[c](t)",
      "'Baseline threshold'~m[0]"
    ))
  ) +
  scale_shape_manual(values = c("Observed" = 16, "Missed" = 16)) +
  scale_linetype_manual(values = c(
    "Data-generating completeness threshold" = "solid",
    "Baseline threshold" = "dashed"
  )) +
  guides(
    shape = "none", linetype = "none",
    colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(
        shape = c(16, 16, NA, NA),
        linetype = c(0, 0, 1, 2),
        linewidth = c(0, 0, 1.0, 0.8)
      )
    )
  ) +
  coord_cartesian(xlim = c(0, 0.35), ylim = c(2.4, 6.7)) +
  labs(
    x = "Days after the seeded M6.5 mainshock",
    y = "Magnitude", colour = NULL,
    title = "(b) Imposed short-term catalogue incompleteness"
  ) +
  theme_thesis() +
  theme(
    legend.position = c(0.73, 0.84),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.88), colour = "grey75",
      linewidth = 0.35
    ),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.margin = margin(3, 5, 3, 5),
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    )
  )

synthetic_figure <- synthetic_a / synthetic_b +
  plot_layout(heights = c(0.90, 1.15))
save_figure(synthetic_figure, "02_synthetic_observation_mechanism", 8.2, 6.8)
write.csv(
  synthetic_audit,
  file.path(output_dir, "02_synthetic_observation_mechanism_data.csv"),
  row.names = FALSE
)
write.csv(
  synthetic_curve,
  file.path(output_dir, "02_synthetic_observation_curve.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 2. Paired parameter-recovery errors.
# -------------------------------------------------------------------------
parameter_order <- c("mu", "K", "alpha", "c", "p")
parameter_labels <- c(
  K = "K", alpha = "alpha", c = "c", mu = "mu", p = "p"
)

posterior <- read.csv(
  "results/submission_v1/mcmc_primary/summary/posterior_long.csv"
) |>
  mutate(
    model = recode(model, naive = "Naive", plugin = "Plug-in"),
    model = factor(model, levels = names(model_colours)),
    parameter = factor(parameter, levels = parameter_order),
    error = median - true
  )

model_summary <- read.csv(
  "results/submission_v1/mcmc_primary/summary/model_summary.csv"
)
paired_recovery <- posterior |>
  transmute(
    sim_id, parameter, model = as.character(model),
    absolute_error = abs(error)
  ) |>
  pivot_wider(names_from = model, values_from = absolute_error) |>
  mutate(
    improved = `Plug-in` < Naive,
    comparison = factor(
      ifelse(
        improved, "Lower error under Plug-in",
        "Equal or higher error under Plug-in"
      ),
      levels = c(
        "Lower error under Plug-in", "Equal or higher error under Plug-in"
      )
    )
  )

paired_recovery_annotations <- paired_recovery |>
  group_by(parameter) |>
  summarise(
    panel_max = max(Naive, `Plug-in`),
    improved_percent = 100 * mean(improved),
    label = sprintf("Plug-in improved: %.0f%%", improved_percent),
    .groups = "drop"
  )

paired_recovery_plot <- ggplot(
  paired_recovery,
  aes(Naive, `Plug-in`, colour = comparison)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    colour = "grey35", linewidth = 0.65
  ) +
  geom_point(size = 1.8, alpha = 0.70) +
  geom_blank(
    data = paired_recovery_annotations,
    aes(x = 0, y = 0), inherit.aes = FALSE
  ) +
  geom_blank(
    data = paired_recovery_annotations,
    aes(x = panel_max, y = panel_max), inherit.aes = FALSE
  ) +
  geom_label(
    data = paired_recovery_annotations,
    aes(x = 0.04 * panel_max, y = 0.94 * panel_max, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1,
    size = 2.8, linewidth = 0, fill = "white", colour = "grey20"
  ) +
  facet_wrap(
    ~parameter, scales = "free", ncol = 3,
    labeller = as_labeller(parameter_labels, label_parsed)
  ) +
  scale_colour_manual(values = c(
    "Lower error under Plug-in" = blue,
    "Equal or higher error under Plug-in" = orange
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.07))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.07))) +
  labs(
    x = "Naive posterior-median absolute error",
    y = "Plug-in posterior-median absolute error",
    colour = NULL
  ) +
  theme_thesis() +
  theme(
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    ),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.35),
    legend.position = "top",
    aspect.ratio = 1,
    panel.spacing.x = grid::unit(0.55, "lines"),
    panel.spacing.y = grid::unit(0.65, "lines"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(3, 3, 3, 3)
  )

save_figure(
  paired_recovery_plot,
  "03_synthetic_parameter_recovery",
  8.8, 5.6
)
write.csv(
  paired_recovery,
  file.path(output_dir, "03_synthetic_parameter_recovery_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 3. Empirical coverage with paired model changes.
# -------------------------------------------------------------------------
coverage_order <- c("mu", "K", "alpha", "c", "p")
coverage <- model_summary |>
  filter(parameter %in% coverage_order) |>
  mutate(
    model = recode(model, naive = "Naive", plugin = "Plug-in"),
    model = factor(model, levels = names(model_colours)),
    parameter = factor(parameter, levels = coverage_order),
    y_base = length(coverage_order) + 1 - match(
      as.character(parameter), coverage_order
    ),
    y = y_base + ifelse(model == "Naive", 0.14, -0.14),
    label_y = y + ifelse(model == "Naive", 0.15, -0.15)
  )

coverage_plot <- ggplot() +
  geom_vline(xintercept = 0.95, linetype = "dashed", colour = red) +
  geom_segment(
    data = coverage,
    aes(
      x = coverage_ci_low, xend = coverage_ci_high, y = y, yend = y,
      colour = model
    ),
    linewidth = 1.05
  ) +
  geom_segment(
    data = coverage,
    aes(
      x = coverage_ci_low, xend = coverage_ci_low,
      y = y - 0.06, yend = y + 0.06, colour = model
    ),
    linewidth = 0.95
  ) +
  geom_segment(
    data = coverage,
    aes(
      x = coverage_ci_high, xend = coverage_ci_high,
      y = y - 0.06, yend = y + 0.06, colour = model
    ),
    linewidth = 0.95
  ) +
  geom_point(
    data = coverage, aes(coverage, y, colour = model, shape = model),
    size = 4.2
  ) +
  geom_text(
    data = coverage,
    aes(coverage, label_y, label = sprintf("%.2f", coverage), colour = model),
    vjust = 0.5, size = 3.15, fontface = "bold", show.legend = FALSE
  ) +
  annotate(
    "text", x = 0.77, y = 5.42,
    label = "Nominal coverage = 0.95", hjust = 0,
    colour = red, size = 3.0
  ) +
  scale_colour_manual(values = model_colours) +
  scale_shape_manual(values = c("Naive" = 16, "Plug-in" = 17)) +
  scale_x_continuous(
    limits = c(0, 1.01), breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(coverage_order),
    labels = parse(text = rev(parameter_labels[coverage_order])),
    limits = c(0.55, 5.52), expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Empirical coverage", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_thesis() +
  theme(
    legend.position = "top",
    legend.key.width = grid::unit(1.15, "lines"),
    legend.spacing.x = grid::unit(0.2, "lines"),
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    ),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

save_figure(coverage_plot, "04_synthetic_empirical_coverage", 8.2, 4.6)
write.csv(
  coverage,
  file.path(output_dir, "04_synthetic_empirical_coverage_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 4. RMSE comparison across primary and simulation-only benchmarks.
# -------------------------------------------------------------------------
benchmark_rmse <- read.csv(
  "results/submission_v1/mcmc_primary/summary/benchmark_decomposition.csv"
) |>
  select(
    parameter, rmse_naive, rmse_plugin, rmse_oracle, rmse_complete
  ) |>
  pivot_longer(
    starts_with("rmse_"), names_to = "model", values_to = "rmse"
  ) |>
  mutate(
    model = recode(
      model,
      rmse_naive = "Naive",
      rmse_plugin = "Estimated Plug-in",
      rmse_oracle = "Oracle",
      rmse_complete = "Complete-Data"
    ),
    model = factor(model, levels = names(benchmark_model_colours)),
    parameter = factor(parameter, levels = rev(parameter_order))
  ) |>
  group_by(parameter) |>
  mutate(relative_rmse = rmse / rmse[model == "Naive"]) |>
  ungroup()

benchmark_ranges <- benchmark_rmse |>
  group_by(parameter) |>
  summarise(
    xmin = min(relative_rmse), xmax = max(relative_rmse), .groups = "drop"
  )

benchmark_rmse_plot <- ggplot(
  benchmark_rmse,
  aes(relative_rmse, parameter, colour = model, shape = model)
) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_segment(
    data = benchmark_ranges,
    aes(x = xmin, xend = xmax, y = parameter, yend = parameter),
    inherit.aes = FALSE, colour = "grey76", linewidth = 0.95
  ) +
  geom_point(
    position = position_dodge(width = 0.48), size = 3.5, alpha = 1
  ) +
  scale_colour_manual(values = benchmark_model_colours) +
  scale_shape_manual(values = c(
    "Naive" = 16,
    "Estimated Plug-in" = 17,
    "Oracle" = 15,
    "Complete-Data" = 18
  )) +
  scale_x_continuous(
    limits = c(0.20, 1.12), breaks = seq(0.2, 1.0, by = 0.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = "RMSE relative to Naive ETAS", y = "Parameter",
    colour = NULL, shape = NULL
  ) +
  theme_thesis() +
  theme(
    legend.position = "top",
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    ),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35)
  )

save_figure(benchmark_rmse_plot, "05_synthetic_benchmark_rmse", 8.2, 4.6)
write.csv(
  benchmark_rmse,
  file.path(output_dir, "05_synthetic_benchmark_rmse_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Supplementary Figure S1. Worst-case diagnostics for all 400 synthetic fits.
# -------------------------------------------------------------------------
primary_diagnostics <- read.csv(
  "results/submission_v1/mcmc_primary/summary/mcmc_diagnostics.csv"
) |>
  mutate(model = recode(
    model, naive = "Naive", plugin = "Estimated Plug-in"
  ))
benchmark_diagnostics <- read.csv(
  "results/submission_v1/mcmc_primary/summary/benchmark_mcmc_diagnostics.csv"
) |>
  mutate(model = recode(
    model, oracle = "Oracle", complete = "Complete-Data"
  ))
synthetic_diagnostics <- bind_rows(
  primary_diagnostics, benchmark_diagnostics
)
fit_diagnostics <- synthetic_diagnostics |>
  group_by(sim_id, model) |>
  summarise(
    max_rhat = max(rhat), min_ess = min(effective_size), .groups = "drop"
  ) |>
  mutate(model = factor(model, levels = names(benchmark_model_colours)))

diagnostic_plot <- ggplot(
  fit_diagnostics, aes(max_rhat, min_ess, colour = model, shape = model)
) +
  geom_vline(xintercept = 1.05, linetype = "dashed", colour = red) +
  geom_hline(yintercept = 400, linetype = "dashed", colour = red) +
  geom_point(size = 2.0, alpha = 0.70) +
  annotate(
    "text", x = 1.049, y = 1370, label = expression(hat(R) == 1.05),
    hjust = 1, colour = red, size = 3.0
  ) +
  annotate(
    "text", x = 1.001, y = 420, label = "ESS = 400",
    hjust = 0, colour = red, size = 3.0
  ) +
  scale_colour_manual(values = benchmark_model_colours) +
  scale_shape_manual(values = c(
    "Naive" = 16,
    "Estimated Plug-in" = 17,
    "Oracle" = 15,
    "Complete-Data" = 18
  )) +
  scale_x_continuous(
    limits = c(1.000, 1.052), breaks = seq(1.00, 1.05, by = 0.01)
  ) +
  scale_y_continuous(limits = c(350, 1425), breaks = seq(400, 1400, 200)) +
  labs(
    x = expression("Maximum " * hat(R) * " across parameters"),
    y = "Minimum ESS across parameters", colour = NULL, shape = NULL
  ) +
  theme_thesis() +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA)
  )
save_figure(diagnostic_plot, "S1_synthetic_mcmc_diagnostics", 8.2, 4.8)
write.csv(
  fit_diagnostics,
  file.path(output_dir, "S1_synthetic_mcmc_diagnostics_data.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 4. Ridgecrest catalogue, held-out split, and completeness curve.
# -------------------------------------------------------------------------
ridgecrest_cfg <- submission_ridgecrest_config()
ridgecrest_raw <- read.csv(ridgecrest_cfg$raw_path, stringsAsFactors = FALSE)
ridgecrest_raw$time_posix <- as.POSIXct(
  ridgecrest_raw$time, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"
)
ridgecrest_main <- ridgecrest_raw[
  ridgecrest_raw$id == ridgecrest_cfg$mainshock_id, , drop = FALSE
]
rad <- pi / 180
dlat <- (ridgecrest_raw$latitude - ridgecrest_cfg$centre_latitude) * rad
dlon <- (ridgecrest_raw$longitude - ridgecrest_cfg$centre_longitude) * rad
a <- sin(dlat / 2)^2 +
  cos(ridgecrest_cfg$centre_latitude * rad) *
  cos(ridgecrest_raw$latitude * rad) * sin(dlon / 2)^2
ridgecrest_raw$distance_km <- 6371.0088 * 2 * atan2(sqrt(a), sqrt(1 - a))
ridgecrest_raw$ts <- as.numeric(difftime(
  ridgecrest_raw$time_posix, ridgecrest_main$time_posix, units = "days"
))
ridgecrest <- ridgecrest_raw |>
  filter(
    distance_km <= ridgecrest_cfg$radius_km,
    mag >= ridgecrest_cfg$M0,
    ts <= ridgecrest_cfg$forecast_end
  ) |>
  transmute(
    ID = id, ts, magnitudes = mag,
    period = factor(
      ifelse(ts <= ridgecrest_cfg$train_end, "Training", "Held-out test"),
      levels = c("Training", "Held-out test")
    )
  ) |>
  arrange(ts, ID)

completeness <- read.csv(
  "results/submission_v1/ridgecrest/completeness_summary.csv"
)
ridgecrest_recovery <- completeness$return_to_baseline_days[[1]]
ridgecrest_time <- 10^seq(-6, log10(ridgecrest_cfg$forecast_end), length.out = 1800)
ridgecrest_curve <- data.frame(
  ts = ridgecrest_time,
  mc = pmax(
    ridgecrest_cfg$M0,
    ridgecrest_main$mag[[1]] - completeness$G[[1]] -
      completeness$H[[1]] * log10(ridgecrest_time)
  )
)

ridgecrest_a <- ggplot(ridgecrest, aes(ts, magnitudes, colour = period)) +
  geom_hline(
    yintercept = ridgecrest_cfg$M0, linetype = "dashed",
    colour = "grey35"
  ) +
  geom_point(alpha = 0.64, size = 1.1) +
  geom_line(
    data = ridgecrest_curve |> filter(ts <= ridgecrest_recovery),
    aes(ts, mc), inherit.aes = FALSE,
    colour = red, linewidth = 0.95
  ) +
  annotate(
    "point", x = 0, y = ridgecrest_main$mag[[1]],
    shape = 16, size = 2.2, colour = "black"
  ) +
  geom_vline(
    xintercept = ridgecrest_cfg$train_end,
    linetype = "dotted", colour = "grey20", linewidth = 0.8
  ) +
  annotate(
    "text", x = ridgecrest_cfg$train_end + 0.10, y = 6.75,
    label = "Train/test split", hjust = 0, size = 3.0
  ) +
  scale_colour_manual(values = c(
    "Training" = "#666666", "Held-out test" = blue
  )) +
  coord_cartesian(xlim = range(ridgecrest$ts), ylim = c(2.85, 7.25)) +
  labs(
    x = "Days relative to the M7.1 mainshock", y = "Magnitude",
    colour = NULL, title = "(a) Training and held-out periods"
  ) +
  theme_thesis() +
  theme(
    legend.position = c(0.985, 0.965),
    legend.justification = c(1, 1),
    legend.direction = "horizontal",
    legend.background = element_rect(
      fill = scales::alpha("white", 0.88), colour = "grey70", linewidth = 0.35
    ),
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    )
  )

ridgecrest_b <- ggplot(
  ridgecrest |> filter(ts >= 0, ts <= 0.06),
  aes(ts, magnitudes)
) +
  geom_hline(
    yintercept = ridgecrest_cfg$M0, linetype = "dashed",
    colour = "grey35"
  ) +
  geom_point(alpha = 0.66, size = 1.2, colour = "#666666") +
  geom_line(
    data = ridgecrest_curve |> filter(ts <= 0.06),
    aes(ts, mc), colour = red, linewidth = 1.0
  ) +
  annotate(
    "text", x = 0.003, y = 5.55,
    label = expression(hat(m)[c](t)), colour = red,
    hjust = 0, size = 3.0
  ) +
  annotate(
    "point", x = 0, y = ridgecrest_main$mag[[1]],
    shape = 16, size = 2.2, colour = "black"
  ) +
  geom_vline(
    xintercept = ridgecrest_recovery,
    linetype = "dotted", colour = red, linewidth = 0.8
  ) +
  annotate(
    "text", x = ridgecrest_recovery + 0.0015, y = 6.78,
    label = sprintf(
      "Return to baseline\n%.3f days (%.0f min)",
      ridgecrest_recovery, 24 * 60 * ridgecrest_recovery
    ),
    hjust = 0, colour = red, size = 3.0
  ) +
  annotate(
    "text", x = 0.058, y = 2.76,
    label = expression(m[0] == 2.95), colour = "grey30", size = 3.0,
    hjust = 1
  ) +
  coord_cartesian(xlim = c(0, 0.06), ylim = c(2.65, 7.25)) +
  labs(
    x = "Days after mainshock", y = "Magnitude",
    title = "(b) Early completeness recovery"
  ) +
  theme_thesis() +
  theme(
    legend.position = "none",
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    )
  )

ridgecrest_figure <- ridgecrest_a / ridgecrest_b +
  plot_layout(heights = c(0.9, 1.1))
save_figure(
  ridgecrest_figure,
  "05_ridgecrest_catalogue_completeness", 8.2, 6.0
)
write.csv(
  ridgecrest_curve,
  file.path(output_dir, "05_ridgecrest_completeness_curve.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Figure 5. Ridgecrest held-out predictive performance.
# -------------------------------------------------------------------------
forecast_dir <- file.path(
  "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast"
)
forecast <- read.csv(file.path(forecast_dir, "forecast_summary.csv")) |>
  mutate(
    model = factor(
      recode(model, naive = "Naive", plugin = "Plug-in"),
      levels = names(model_colours)
    )
  )
forecast_draws <- read.csv(
  file.path(forecast_dir, "forecast_count_draws.csv")
) |>
  mutate(
    model = recode(model, naive = "Naive", plugin = "Plug-in"),
    model = factor(model, levels = names(model_colours))
  )
observed_count <- unique(forecast$n_test)

tail_markers <- forecast |>
  transmute(
    model,
    q975 = count_q975,
    label = paste0(
      "97.5th percentile = ",
      format(round(count_q975), big.mark = ",", scientific = FALSE)
    )
  )

prediction_figure <- ggplot(
  forecast_draws |> filter(!overflow),
  aes(x = count, fill = model, colour = model)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 35, alpha = 0.25, linewidth = 0.25
  ) +
  geom_density(linewidth = 0.85, adjust = 1.1, fill = NA) +
  geom_vline(
    xintercept = observed_count, linetype = "dashed",
    linewidth = 0.85, colour = red
  ) +
  geom_vline(
    data = tail_markers,
    aes(xintercept = q975, colour = model),
    inherit.aes = FALSE, linetype = "dotdash", linewidth = 0.75,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = observed_count + 18, y = Inf,
    label = paste("Observed =", observed_count),
    hjust = 0, vjust = 1.2, colour = red, size = 3.0
  ) +
  geom_text(
    data = tail_markers,
    aes(
      x = q975, y = Inf, label = label,
      colour = model
    ),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.25, size = 3.0,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = model_colours) +
  scale_fill_manual(values = model_colours) +
  facet_grid(model ~ .) +
  labs(
    x = "Held-out event count", y = "Predictive density"
  ) +
  guides(colour = "none", fill = "none") +
  theme_thesis() +
  theme(
    legend.position = "none",
    panel.border = element_rect(
      colour = "grey45", fill = NA, linewidth = 0.55
    )
  )
save_figure(prediction_figure, "06_ridgecrest_prediction_summary", 8.2, 4.5)
write.csv(
  forecast_draws,
  file.path(output_dir, "06_ridgecrest_prediction_summary_data.csv"),
  row.names = FALSE
)

message(
  "Created six dissertation-body figures and two appendix figures in ",
  output_dir
)
