#!/usr/bin/env Rscript

# Supplementary robustness figures for fixed-b sensitivity.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

output_dir <- "results/submission_v1/figures/thesis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

blue <- "#0072B2"
grey <- "#666666"
red <- "#B2182B"
parameter_colours <- c(
  mu = "#666666", K = "#0072B2", alpha = "#D55E00",
  c = "#009E73", p = "#CC79A7"
)
model_colours <- c("Naive" = grey, "Plug-in" = blue)
parameter_labels <- c(mu = expression(mu), K = expression(K),
                      alpha = expression(alpha), c = expression(c),
                      p = expression(p))

theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.5),
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

# Figure S2: synthetic fixed-b sensitivity.
b_summary <- read.csv(
  "results/submission_v1/b_sensitivity/synthetic/summary/model_summary.csv"
) |>
  mutate(parameter = factor(parameter, levels = c("mu", "K", "alpha", "c", "p")))
naive_reference <- read.csv(
  "results/submission_v1/mcmc_primary/summary/model_summary.csv"
) |>
  filter(model == "naive") |>
  transmute(
    parameter = factor(parameter, levels = c("mu", "K", "alpha", "c", "p")),
    rmse_naive = rmse
  )
b_plot_data <- b_summary |>
  left_join(naive_reference, by = "parameter") |>
  mutate(relative_rmse = rmse / rmse_naive)

s3a <- ggplot(
  b_plot_data,
  aes(assumed_b, relative_rmse, colour = parameter, group = parameter)
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = parameter_colours, labels = parameter_labels,
    name = "Parameter"
  ) +
  scale_x_continuous(breaks = c(0.8, 1, 1.2)) +
  labs(
    x = expression("Assumed fixed " * b),
    y = "Plug-in RMSE / Naive RMSE",
    title = "(a) RMSE relative to Naive ETAS"
  ) +
  theme_thesis()

s3b <- ggplot(
  b_plot_data,
  aes(assumed_b, coverage, colour = parameter, group = parameter)
) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey45") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = parameter_colours, labels = parameter_labels,
    name = "Parameter"
  ) +
  scale_x_continuous(breaks = c(0.8, 1, 1.2)) +
  scale_y_continuous(limits = c(0.75, 1), breaks = seq(0.75, 1, 0.05)) +
  labs(
    x = expression("Assumed fixed " * b),
    y = "Empirical coverage",
    title = "(b) Credible-interval coverage"
  ) +
  theme_thesis()

s3 <- (s3a | s3b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
save_figure(s3, "S2_fixed_b_synthetic_sensitivity", 8.2, 4.2)
write.csv(
  b_plot_data,
  file.path(output_dir, "S2_fixed_b_synthetic_sensitivity_data.csv"),
  row.names = FALSE
)

# Figure S3: Ridgecrest count stability under fixed b.
read_forecast <- function(path, assumed_b) {
  read.csv(path) |>
    mutate(
      assumed_b = assumed_b,
      model = recode(model, naive = "Naive", plugin = "Plug-in")
    )
}
ridge_b <- bind_rows(
  read_forecast(
    "results/submission_v1/ridgecrest/b_sensitivity/b_0p8/forecast/forecast_summary.csv",
    0.8
  ),
  read_forecast(
    "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/forecast_summary.csv",
    1
  ),
  read_forecast(
    "results/submission_v1/ridgecrest/b_sensitivity/b_1p2/forecast/forecast_summary.csv",
    1.2
  )
) |>
  mutate(
    model = factor(model, levels = c("Naive", "Plug-in")),
    cap_dominated = count_median >= 2000,
    overflow_label = sprintf("%.1f%%", 100 * overflow_rate),
    overflow_label_y = case_when(
      assumed_b == 0.8 ~ 62,
      assumed_b == 1 & model == "Naive" ~ 12,
      TRUE ~ 100 * overflow_rate + 4
    )
  )

read_mc_error <- function(path, assumed_b) {
  read.csv(path) |>
    transmute(
      assumed_b = assumed_b,
      ig_mc_q025 = mc_q025,
      ig_mc_q975 = mc_q975,
      probability_positive = probability_positive
    )
}
ridge_score_sensitivity <- ridge_b |>
  select(assumed_b, model, log_score_per_event) |>
  mutate(model = recode(as.character(model), "Naive" = "naive", "Plug-in" = "plugin")) |>
  pivot_wider(
    names_from = model, values_from = log_score_per_event,
    names_glue = "{model}_lpd_per_event"
  ) |>
  left_join(
    ridge_b |>
      filter(model == "Plug-in") |>
      select(assumed_b, information_gain_per_event = information_gain_vs_naive),
    by = "assumed_b"
  ) |>
  left_join(
    bind_rows(
      read_mc_error(
        "results/submission_v1/ridgecrest/b_sensitivity/b_0p8/forecast/forecast_mc_error.csv",
        0.8
      ),
      read_mc_error(
        "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/forecast_mc_error.csv",
        1
      ),
      read_mc_error(
        "results/submission_v1/ridgecrest/b_sensitivity/b_1p2/forecast/forecast_mc_error.csv",
        1.2
      )
    ),
    by = "assumed_b"
  ) |>
  arrange(assumed_b)

s4a <- ggplot(
  ridge_b,
  aes(assumed_b, count_median, colour = model, shape = model)
) +
  geom_hline(yintercept = 224, linetype = "dashed", colour = red) +
  geom_errorbar(
    aes(ymin = count_q025, ymax = pmin(count_q975, 2000)),
    width = 0.025, linewidth = 0.7,
    position = position_dodge(width = 0.055)
  ) +
  geom_point(
    aes(fill = ifelse(cap_dominated, "Cap-dominated", "Identifiable")),
    size = 3, stroke = 0.9, position = position_dodge(width = 0.055)
  ) +
  annotate(
    "segment", x = 0.78625, xend = 0.78625, y = 2025, yend = 2100,
    colour = grey, linewidth = 0.6,
    arrow = arrow(length = grid::unit(0.07, "inches"), type = "closed")
  ) +
  annotate(
    "segment", x = 0.81375, xend = 0.81375, y = 2025, yend = 2100,
    colour = blue, linewidth = 0.6,
    arrow = arrow(length = grid::unit(0.07, "inches"), type = "closed")
  ) +
  annotate(
    "text", x = 0.82, y = 2185, hjust = 0,
    label = "59.4% Naive; 59.5% Plug-in reached cap",
    colour = "grey20", size = 2.9
  ) +
  annotate("text", x = 1.10, y = 340, label = "Observed = 224",
           colour = red, size = 3.1) +
  scale_colour_manual(values = model_colours, name = "Model") +
  scale_shape_manual(values = c("Naive" = 21, "Plug-in" = 22), name = "Model") +
  scale_fill_manual(
    values = c("Cap-dominated" = "white", "Identifiable" = "grey70"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = c(0.8, 1, 1.2)) +
  scale_y_continuous(limits = c(0, 2250), breaks = c(0, 500, 1000, 1500, 2000)) +
  labs(
    x = expression("Assumed fixed " * b),
    y = "Posterior predictive event count"
  ) +
  theme_thesis()
save_figure(s4a, "S3_ridgecrest_fixed_b_forecast", 8.2, 4.4)
write.csv(
  ridge_b,
  file.path(output_dir, "S3_ridgecrest_fixed_b_forecast_data.csv"),
  row.names = FALSE
)
write.csv(
  ridge_score_sensitivity,
  file.path(output_dir, "ridgecrest_fixed_b_predictive_scores_data.csv"),
  row.names = FALSE
)

# Figure S4: Ridgecrest posterior sensitivity to fixed b.
read_posterior_summary <- function(path, assumed_b) {
  read.csv(path) |>
    filter(model == "plugin", parameter %in% c("K", "alpha")) |>
    mutate(assumed_b = assumed_b)
}
ridge_posterior_b <- bind_rows(
  read_posterior_summary(
    "results/submission_v1/ridgecrest/b_sensitivity/b_0p8/forecast/posterior_summary.csv",
    0.8
  ),
  read_posterior_summary(
    "results/submission_v1/ridgecrest/mcmc_conditioned_single/forecast/posterior_summary.csv",
    1
  ),
  read_posterior_summary(
    "results/submission_v1/ridgecrest/b_sensitivity/b_1p2/forecast/posterior_summary.csv",
    1.2
  )
) |>
  mutate(
    parameter = factor(parameter, levels = c("K", "alpha")),
    assumed_b_factor = factor(
      assumed_b, levels = c(1.2, 1, 0.8), labels = c("1.2", "1.0", "0.8")
    )
  )
write.csv(
  ridge_posterior_b,
  file.path(output_dir, "S4_ridgecrest_fixed_b_posterior_data.csv"),
  row.names = FALSE
)

cat("Created Supplementary Figures S2--S3 in", output_dir, "\n")
