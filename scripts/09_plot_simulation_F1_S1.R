# Produce ONLY F1 (pooled MAE) and S1 (3 x 3 setting-specific MAE).
# Match Figure 2: monochrome graphics, 6.5-inch width, ~8-point text.
stopifnot(requireNamespace("here", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
source(here::here("R", "simulation.R"))
cfg <- simulation_config()
folder <- here::here("results", "simulation", cfg$name)
fig_dir <- here::here("figures", cfg$name)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
scenario <- utils::read.csv(file.path(folder, "scenario_summary.csv"))
pooled <- utils::read.csv(file.path(folder, "pooled_MAE_by_K.csv"))
methods <- c("NW", "CW", "EQ", "ER")
K_values <- cfg$K

# Black-and-white distinctions remain visible when printed in grayscale.
lines <- c(NW = "dotted", CW = "dotdash", EQ = "solid", ER = "dashed")
shapes <- c(NW = 1L, CW = 0L, EQ = 16L, ER = 17L)
fig_width <- 6.5
fig_fontsize <- 8

format_data <- function(d) {
  d$method <- factor(d$method, levels = methods)
  d$K_label <- factor(d$K, levels = K_values)
  if (anyNA(d$K_label) || anyNA(d$method)) {
    stop("Unexpected figure dimensions")
  }
  d
}

plot_base <- function(d, outcome) {
  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = K_label, y = .data[[outcome]],
      linetype = method, shape = method, group = method
    )
  ) +
    ggplot2::geom_line(color = "black", linewidth = 0.4) +
    ggplot2::geom_point(color = "black", size = 1.7, stroke = 0.4) +
    ggplot2::scale_linetype_manual(values = lines, breaks = methods) +
    ggplot2::scale_shape_manual(values = shapes, breaks = methods) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      x = "Rank-grid resolution (K)", y = "Mean absolute error",
      linetype = NULL, shape = NULL
    ) +
    ggplot2::theme_bw(base_size = fig_fontsize) +
    #ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 7.2, color = "black"),
      axis.title = ggplot2::element_text(size = 7.2, color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.35),
      strip.text = ggplot2::element_text(size = 7.2, face = "plain"),
      strip.background = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.65, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 7.2),
      legend.key.width = grid::unit(1, "cm"),
      legend.spacing.x = grid::unit(5, "pt"),
      plot.margin = ggplot2::margin(4, 4, 4, 4, unit = "pt"),
      axis.text.x = ggplot2::element_text(
        size = 7.2, color = "black",
        margin = ggplot2::margin(t = 2)
      ),
      axis.text.y = ggplot2::element_text(
        size = 7.2, color = "black",
        margin = ggplot2::margin(r = 2)
      ),
      axis.title.x = ggplot2::element_text(
        size = 7.2,
        margin = ggplot2::margin(t = 4)
      ),
      axis.title.y = ggplot2::element_text(
        size = 7.2,
        margin = ggplot2::margin(r = 4)
      ),
      panel.grid.major = ggplot2::element_line(
        color = "grey70",
        linewidth = 0.1
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.spacing = grid::unit(1, "pt")
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    )
}


save_figure <- function(plot, stem, width, height) {
  for (ext in c("pdf", "png")) {
    ggplot2::ggsave(
      filename = file.path(fig_dir, paste0(stem, ".", ext)),
      plot = plot, width = width, height = height,
      units = "in", dpi = 300, bg = "white"
    )
  }
}

p <- plot_base(format_data(pooled), "pooled_mae")
save_figure(p, "F1_pooled_MAE", width = 3.2, height = 2.4)

scenario$eta_label <- factor(
  scenario$eta,
  levels = c(0, 2, 4),
  labels = c("eta == 0", "eta == 2", "eta == 4")
)
scenario$rho_label <- factor(
  scenario$rho_uv,
  levels = c(0, 0.25, 0.50),
  labels = c("rho[uv] == 0.00", "rho[uv] == 0.25", "rho[uv] == 0.50")
)
q <- plot_base(format_data(scenario), "mae") +
  ggplot2::facet_grid(eta_label ~ rho_label, labeller = ggplot2::label_parsed)

save_figure(q, "S1_MAE_3x3_grid", width = 6.5, height = 6.2)

message("Saved F1 and S1 to: ", fig_dir)
