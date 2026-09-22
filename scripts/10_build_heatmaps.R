####################################################################################################
#
# SCRIPT 10: PEARSON M = 20 HEATMAPS
#
# Uses the compiled Pearson M = 20 summary table from Script 9
# and creates smaller heatmaps for selected simulation scenarios.
#
# Input:
#   icspack/results/sim/correlations/Pearson/summaries/compiled_pearson_M20_summary.rds
#
# Output:
#   icspack/figures/sim/heatmaps/Pearson_M20/
#
####################################################################################################


# ==================================================================================================
# 1. LOAD LIBRARIES
# ==================================================================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)


# ==================================================================================================
# 2. DIRECTORIES
# ==================================================================================================

summary_file <- file.path(
  "results",
  "sim",
  "correlations",
  "Pearson",
  "summaries",
  "compiled_pearson_M20_summary.rds"
)

figure_dir <- file.path(
  "figures",
  "sim",
  "heatmaps",
  "Pearson_M20"
)

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==================================================================================================
# 3. LOAD SUMMARY DATA
# ==================================================================================================

compiled_m20 <- readRDS(summary_file)


# ==================================================================================================
# 4. ENSURE METRICS EXIST
# ==================================================================================================

if (!"bias_true" %in% names(compiled_m20)) {
  compiled_m20 <- compiled_m20 %>%
    mutate(
      bias_true = mean_estimate - rho_true,
      abs_bias_true = abs(bias_true),
      bias_obs = mean_estimate - rho_obs_hat,
      abs_bias_obs = abs(bias_obs),
      coverage_true_error = abs(coverage_true - 0.95),
      coverage_obs_error = abs(coverage_obs - 0.95)
    )
}


# ==================================================================================================
# 5. KEEP ONLY CURRENT WEIGHTS
# ==================================================================================================

weight_order <- c(
  "no_weight",
  "CW",
  "PPW",
  "OPW"
)

compiled_m20 <- compiled_m20 %>%
  filter(weight_type %in% weight_order) %>%
  mutate(
    weight_type = factor(weight_type, levels = weight_order)
  )


# ==================================================================================================
# 6. PARSE SUBGROUPING LABELS
# ==================================================================================================

parse_subgrouping <- function(x) {
  
  parts <- str_match(
    x,
    "^L_Y_(.*?)(?:__K_X_|_K_X_)(.*)$"
  )
  
  if (any(is.na(parts))) {
    stop("Could not parse subgrouping label: ", x)
  }
  
  clean_label <- function(z) {
    case_when(
      z == "sev" ~ "Severity",
      z == "dich_1_5" ~ "1 | 2-5",
      z == "dich_2_5" ~ "1-2 | 3-5",
      z == "dich_3_5" ~ "1-3 | 4-5",
      z == "dich_4_5" ~ "1-4 | 5",
      TRUE ~ z
    )
  }
  
  tibble(
    Y_subgrouping = clean_label(parts[, 2]),
    X_subgrouping = clean_label(parts[, 3])
  )
}

subgroup_lookup <- tibble(
  subgrouping = unique(compiled_m20$subgrouping)
) %>%
  rowwise() %>%
  mutate(
    parsed = list(parse_subgrouping(subgrouping))
  ) %>%
  tidyr::unnest(parsed) %>%
  ungroup()

subgroup_order <- c(
  "1 | 2-5",
  "1-2 | 3-5",
  "1-3 | 4-5",
  "1-4 | 5",
  "Severity"
)

compiled_m20 <- compiled_m20 %>%
  left_join(subgroup_lookup, by = "subgrouping") %>%
  mutate(
    X_subgrouping = factor(X_subgrouping, levels = subgroup_order),
    Y_subgrouping = factor(Y_subgrouping, levels = subgroup_order)
  )


# ==================================================================================================
# 7. SELECT SCENARIOS
#
# These are the same representative scenarios as before, but now only for Pearson / M = 20.
# ==================================================================================================

representative_scenarios <- tribble(
  ~rho_xy, ~rho_uv, ~eta_x, ~eta_y,
  0,       0,       0,      0,
  0,       0,       4,      4,
  0,       0.5,     0,      4,
  0,       0.5,     4,      0,
  0,       0.5,     4,      4,
  0.5,     0,       4,      0,
  0.5,     0,       0,      4,
  0.5,     0,       4,      4,
  0.5,     0.5,     4,      0,
  0.5,     0.5,     0,      4,
  0.5,     0.5,     4,      4
)


# ==================================================================================================
# 8. METRICS TO PLOT
# ==================================================================================================

metric_specs <- tribble(
  ~value_col,              ~value_label,
  "abs_bias_true",         "Absolute Bias to True Correlation",
  "abs_bias_obs",          "Absolute Bias to Observed Correlation",
  "coverage_true_error",   "Coverage Error to True Correlation",
  "coverage_obs_error",    "Coverage Error to Observed Correlation"
)


# ==================================================================================================
# 9. HEATMAP FUNCTION
# ==================================================================================================

make_heatmap_pearson_m20 <- function(dat, value_col, value_label, out_file) {
  
  subgroup_order <- c(
    "1 | 2-5",
    "1-2 | 3-5",
    "1-3 | 4-5",
    "1-4 | 5",
    "Severity"
  )
  
  weight_order_plot <- c("NW", "CW", "PPW", "OPW")
  
  dat_plot <- dat %>%
    mutate(
      weight_type_plot = case_when(
        as.character(weight_type) == "no_weight" ~ "NW",
        TRUE ~ as.character(weight_type)
      ),
      weight_type_plot = factor(weight_type_plot, levels = weight_order_plot),
      X_subgrouping = factor(X_subgrouping, levels = subgroup_order),
      Y_subgrouping = factor(Y_subgrouping, levels = subgroup_order),
      X_pos = as.numeric(X_subgrouping),
      Y_pos = as.numeric(Y_subgrouping)
    )
  
  fill_limits <- range(dat_plot[[value_col]], na.rm = TRUE)
  
  # NW and CW should be constant across all subgroupings for a fixed scenario.
  bar_dat <- dat_plot %>%
    filter(weight_type_plot %in% c("NW", "CW")) %>%
    group_by(weight_type_plot) %>%
    summarise(
      value = unique(round(.data[[value_col]], 10))[1],
      n_unique = n_distinct(round(.data[[value_col]], 10)),
      .groups = "drop"
    ) %>%
    mutate(
      xmin = 0.5,
      xmax = 5.5,
      ymin = 0.00,
      ymax = 0.70,
      x_text = 3,
      y_text = 0.35,
      label = sprintf("%.3f", value)
    )
  
  if (any(bar_dat$n_unique > 1)) {
    warning("NW or CW is not constant within this scenario.")
  }
  
  heat_dat <- dat_plot %>%
    filter(weight_type_plot %in% c("PPW", "OPW"))
  
  p <- ggplot() +
    
    geom_rect(
      data = bar_dat,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = value
      ),
      color = "white"
    ) +
    
    geom_text(
      data = bar_dat,
      aes(
        x = x_text,
        y = y_text,
        label = label
      ),
      size = 4
    ) +
    
    geom_tile(
      data = heat_dat,
      aes(
        x = X_pos,
        y = Y_pos,
        fill = .data[[value_col]]
      ),
      color = "white"
    ) +
    
    geom_text(
      data = heat_dat,
      aes(
        x = X_pos,
        y = Y_pos,
        label = sprintf("%.3f", .data[[value_col]])
      ),
      size = 3.5
    ) +
    
    scale_x_continuous(
      breaks = 1:5,
      labels = subgroup_order,
      limits = c(0.5, 5.5),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      breaks = 1:5,
      labels = subgroup_order,
      expand = c(0, 0)
    ) +
    
    scale_fill_gradient(
      low = "white",
      high = "red",
      limits = fill_limits,
      na.value = "grey90"
    ) +
    
    facet_grid(
      weight_type_plot ~ .,
      scales = "free_y",
      space = "free_y",
      drop = FALSE
    ) +
    
    labs(
      x = "X subgrouping",
      y = "Y subgrouping",
      fill = value_label,
      title = unique(dat$plot_title)
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 12),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 7,
    height = 10,
    dpi = 300
  )
  
  invisible(p)
}


# ==================================================================================================
# 10. MAKE HEATMAPS
# ==================================================================================================

for (s in seq_len(nrow(representative_scenarios))) {
  
  scen <- representative_scenarios[s, ]
  
  dat_s <- compiled_m20 %>%
    filter(
      rho_xy == scen$rho_xy,
      rho_uv == scen$rho_uv,
      eta_x == scen$eta_x,
      eta_y == scen$eta_y
    ) %>%
    mutate(
      plot_title = paste0(
        "Pearson, M = 20; ",
        "rho_xy = ", rho_xy,
        ", rho_uv = ", rho_uv,
        ", eta_x = ", eta_x,
        ", eta_y = ", eta_y
      )
    )
  
  if (nrow(dat_s) == 0) {
    warning(
      "No data found for scenario: ",
      paste(
        c(
          scen$rho_xy,
          scen$rho_uv,
          scen$eta_x,
          scen$eta_y
        ),
        collapse = ", "
      )
    )
    next
  }
  
  for (i in seq_len(nrow(metric_specs))) {
    
    value_col <- metric_specs$value_col[i]
    value_label <- metric_specs$value_label[i]
    
    out_file <- file.path(
      figure_dir,
      paste0(
        "pearson_m20_",
        value_col,
        "_rhoxy_", scen$rho_xy,
        "_rhouv_", scen$rho_uv,
        "_etax_", scen$eta_x,
        "_etay_", scen$eta_y,
        ".png"
      )
    )
    
    make_heatmap_pearson_m20(
      dat = dat_s,
      value_col = value_col,
      value_label = value_label,
      out_file = out_file
    )
  }
}


# ==================================================================================================
# 11. MESSAGE
# ==================================================================================================

message("Finished creating Pearson M = 20 heatmaps.")
message("Saved to: ", figure_dir)