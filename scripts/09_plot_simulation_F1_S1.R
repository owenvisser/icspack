
# Produce ONLY F1 (pooled MAE) and S1 (3 x 3 setting-specific MAE).

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

palette <- c(
  NW = "#777777",
  CW = "#0072B2",
  EQ = "#009E73",
  ER = "#7B4FA3"
)

lines <- c(
  NW = "dotted",
  CW = "solid",
  EQ = "solid",
  ER = "dashed"
)

shapes <- c(
  NW = 1L,
  CW = 16L,
  EQ = 15L,
  ER = 18L
)

format_data <- function(d, y) {
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
      x = K_label,
      y = .data[[outcome]],
      color = method,
      linetype = method,
      shape = method,
      group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_linetype_manual(values = lines) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(
      x = "Rank-grid resolution (K)",
      y = "Mean absolute error",
      color = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    )
}

p <- plot_base(format_data(pooled, "pooled_mae"), "pooled_mae")

for (ext in c("pdf", "png")) {
  ggplot2::ggsave(
    file.path(fig_dir, paste0("F1_pooled_MAE.", ext)),
    p,
    width = 5.3,
    height = 3.8,
    dpi = 300
  )
}

scenario$eta_label <- factor(
  scenario$eta,
  levels = c(0, 2, 4),
  labels = c("eta == 0", "eta == 2", "eta == 4")
)

scenario$rho_label <- factor(
  scenario$rho_uv,
  levels = c(0, 0.25, 0.50),
  labels = c(
    "rho[uv] == 0.00",
    "rho[uv] == 0.25",
    "rho[uv] == 0.50"
  )
)

q <- plot_base(format_data(scenario, "mae"), "mae") +
  ggplot2::facet_grid(
    eta_label ~ rho_label,
    labeller = ggplot2::label_parsed
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white")
  )

for (ext in c("pdf", "png")) {
  ggplot2::ggsave(
    file.path(fig_dir, paste0("S1_MAE_3x3_grid.", ext)),
    q,
    width = 7.0,
    height = 6.2,
    dpi = 300
  )
}

message("Saved only F1 and S1 to: ", fig_dir)
