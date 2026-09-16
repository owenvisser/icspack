library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(readr)

base_dir <- r"(C:\Users\owvis\OneDrive - University of Florida\nHANES\Sim_Res)"
compiled_summary <- readRDS("C:/Users/owvis/OneDrive - University of Florida/nHANES/Sim_Res/compiled_correlation_summary_true_and_observed_long.rds")

display_dir <- file.path(base_dir, "Display_Checks_M100")
dir.create(display_dir, showWarnings = FALSE, recursive = TRUE)

table_dir <- file.path(display_dir, "Tables")
figure_dir <- file.path(display_dir, "Figures")

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

weight_order <- c("no_weight", "CW", "PPW", "OPW", "MOPW")

compiled_m100 <- compiled_summary %>%
  filter(M == 100) %>%
  mutate(
    weight_type = factor(weight_type, levels = weight_order),
    corr_setting = paste0(
      "rho_xy=", rho_xy,
      ", rho_uv=", rho_uv,
      ", rho_true=", rho_true
    ),
    eta_setting = paste0(
      "eta_x=", eta_x,
      ", eta_y=", eta_y
    ),
    scenario_label = paste0(
      "rho_xy=", rho_xy,
      ", rho_uv=", rho_uv,
      "; eta_x=", eta_x,
      ", eta_y=", eta_y
    ),
    bias_true = mean_estimate - rho_true,
    abs_bias_true = abs(bias_true),
    bias_obs = mean_estimate - rho_obs_hat,
    abs_bias_obs = abs(bias_obs),
    coverage_true_dist = abs(coverage_true - 0.95),
    coverage_obs_dist = abs(coverage_obs - 0.95)
  )


parse_subgrouping <- function(x) {
  
  y_part <- str_extract(x, "L_Y_(.*?)__K_X") %>%
    str_remove("^L_Y_") %>%
    str_remove("__K_X$")
  
  x_part <- str_extract(x, "__K_X_(.*)$") %>%
    str_remove("^__K_X_")
  
  clean_label <- function(z) {
    case_when(
      z == "sev" ~ "Severity",
      str_detect(z, "dich_1_5") ~ "1 | 2-5",
      str_detect(z, "dich_2_5") ~ "1-2 | 3-5",
      str_detect(z, "dich_3_5") ~ "1-3 | 4-5",
      str_detect(z, "dich_4_5") ~ "1-4 | 5",
      TRUE ~ z
    )
  }
  
  tibble(
    Y_subgrouping = clean_label(y_part),
    X_subgrouping = clean_label(x_part)
  )
}

subgrouping_lookup <- tibble(
  subgrouping = unique(compiled_m100$subgrouping)
) %>%
  mutate(parsed = map(subgrouping, parse_subgrouping)) %>%
  unnest(parsed)

subgroup_order <- c(
  "1 | 2-5",
  "1-2 | 3-5",
  "1-3 | 4-5",
  "1-4 | 5",
  "Severity"
)

compiled_m100 <- compiled_m100 %>%
  left_join(subgrouping_lookup, by = "subgrouping") %>%
  mutate(
    Y_subgrouping = factor(Y_subgrouping, levels = subgroup_order),
    X_subgrouping = factor(X_subgrouping, levels = subgroup_order)
  )




#=======================================================================================
# MAIN SEVERITY TABLE
#=======================================================================================

main_severity_long <- compiled_m100 %>%
  filter(subgrouping == "L_Y_sev__K_X_sev") %>%
  mutate(
    est_cov = sprintf(
      "%.3f (%.3f, %.3f)",
      mean_estimate,
      coverage_true,
      coverage_obs
    )
  ) %>%
  arrange(
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    weight_type
  )

write_csv(
  main_severity_long,
  file.path(table_dir, "main_severity_long_M100.csv")
)

main_severity_wide <- main_severity_long %>%
  select(
    correlation_measure,
    rho_xy,
    rho_uv,
    rho_true,
    eta_x,
    eta_y,
    rho_obs_hat,
    weight_type,
    est_cov
  ) %>%
  pivot_wider(
    names_from = weight_type,
    values_from = est_cov
  ) %>%
  arrange(
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  )

write_csv(
  main_severity_wide,
  file.path(table_dir, "main_severity_wide_M100.csv")
)

print(main_severity_wide, n = 100)



#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================


for (cm in unique(main_severity_wide$correlation_measure)) {
  
  tab_cm <- main_severity_wide %>%
    filter(correlation_measure == cm)
  
  write_csv(
    tab_cm,
    file.path(
      table_dir,
      paste0("main_severity_wide_M100_", cm, ".csv")
    )
  )
}



#=======================================================================================
# Heatmaps over all 25 subgroupings
#=======================================================================================

make_heatmap <- function(dat, value_col, value_label, out_file) {
  
  p <- ggplot(
    dat,
    aes(
      x = X_subgrouping,
      y = Y_subgrouping,
      fill = .data[[value_col]]
    )
  ) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = sprintf("%.3f", .data[[value_col]])),
      size = 3
    ) +
    scale_fill_gradient(
      low = "white",
      high = "red"
    ) +
    facet_grid(weight_type ~ correlation_measure) +
    labs(
      x = "X subgrouping",
      y = "Y subgrouping",
      fill = value_label,
      title = unique(dat$plot_title)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 11)
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 13,
    height = 9,
    dpi = 300
  )
  
  invisible(p)
}

# 4A. Representative heatmaps only

representative_scenarios <- tribble(
  ~rho_xy, ~rho_uv, ~eta_x, ~eta_y,
  0,       0,       0,      0,
  0,       0,       4,      4,
  0.5,     0,       4,      0,
  0,       0.5,     0,      4,
  0.5,     0.5,     4,      4
)

heatmap_metrics <- tribble(
  ~value_col,        ~value_label,
  "mean_estimate",   "Mean estimate",
  "bias_true",       "Bias vs true rho",
  "abs_bias_true",   "Abs. bias vs true rho",
  "coverage_true",   "Coverage true",
  "coverage_obs",    "Coverage observed"
)

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
        "M=20; rho_xy=", rho_xy,
        ", rho_uv=", rho_uv,
        ", eta_x=", eta_x,
        ", eta_y=", eta_y
      )
    )
  
  for (m in seq_len(nrow(heatmap_metrics))) {
    
    value_col <- heatmap_metrics$value_col[m]
    value_label <- heatmap_metrics$value_label[m]
    
    out_file <- file.path(
      figure_dir,
      paste0(
        "heatmap_",
        value_col,
        "_rhoxy_", scen$rho_xy,
        "_rhouv_", scen$rho_uv,
        "_etax_", scen$eta_x,
        "_etay_", scen$eta_y,
        ".png"
      )
    )
    
    make_heatmap(
      dat = dat_s,
      value_col = value_col,
      value_label = value_label,
      out_file = out_file
    )
  }
}



#=======================================================================================
# Heatmaps for selected weights only
#=======================================================================================
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

make_heatmap_with_cw_bar <- function(dat, value_col, value_label, out_file,
                                     red_is_better = TRUE) {
  
  subgroup_order <- c(
    "1 | 2-5",
    "1-2 | 3-5",
    "1-3 | 4-5",
    "1-4 | 5",
    "Severity"
  )
  
  corr_order <- c("Pearson", "Spearman", "Phi")
  weight_order <- c("NW", "CW", "PPW", "OPW", "MOPW")
  
  dat_plot <- dat %>%
    mutate(
      correlation_measure = factor(correlation_measure, levels = corr_order),
      weight_type_plot = case_when(
        as.character(weight_type) == "no_weight" ~ "NW",
        TRUE ~ as.character(weight_type)
      ),
      weight_type_plot = factor(weight_type_plot, levels = weight_order),
      X_subgrouping = factor(X_subgrouping, levels = subgroup_order),
      Y_subgrouping = factor(Y_subgrouping, levels = subgroup_order),
      X_pos = as.numeric(X_subgrouping),
      Y_pos = as.numeric(Y_subgrouping)
    )
  
  fill_limits <- range(dat_plot[[value_col]], na.rm = TRUE)
  
  # Set palette direction depending on what "better" means
  low_col  <- if (red_is_better) "white"   else "red"
  high_col <- if (red_is_better) "red" else "white"
  
  # NW and CW: short numeric bar rows
  bar_dat <- dat_plot %>%
    filter(weight_type_plot %in% c("NW", "CW")) %>%
    group_by(correlation_measure, weight_type_plot) %>%
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
    warning("NW or CW is not exactly constant for at least one correlation measure.")
  }
  
  heat_dat <- dat_plot %>%
    filter(weight_type_plot %in% c("PPW", "OPW", "MOPW"))
  
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
      size = 12 / .pt
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
      size = 12 / .pt
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
      low = low_col,
      high = high_col,
      limits = fill_limits,
      na.value = "grey90"
    ) +
    
    facet_grid(
      weight_type_plot ~ correlation_measure,
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
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 12),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  invisible(p)
}

representative_scenarios <- tribble(
  ~rho_xy, ~rho_uv, ~eta_x, ~eta_y,
  0,       0,       0,      0
  ,
  0,       0,       4,      4,
  0,       0.5,     0,      4,
  0,       0.5,     4,      0,
  0,       0.5,     4,      4,
  0.5,       0,     4,      0,
  0.5,       0,     0,      4,
  0.5,       0,     4,      4,
  0.5,       0.5,     4,      0,
  0.5,       0.5,     0,      4,
  0.5,       0.5,     4,      4
)


metric_specs <- tibble::tribble(
  ~value_col,              ~value_label,                    ~red_is_better,
  "abs_bias_true",         "Absolute\nBias\nTrue",          TRUE,
  "abs_bias_obs",          "Absolute\nBias\nObserved",      TRUE,
  "coverage_true_error",   "Coverage\nTrue\nError",         TRUE,
  "coverage_obs_error",    "Coverage\nObserved\nError",     TRUE
)

metric_labels <- c(
  abs_bias_true       = "Absolute\nBias\nTrue",
  abs_bias_obs        = "Absolute\nBias\nObserved",
  coverage_true_error = "Coverage\nTrue\nError",
  coverage_obs_error  = "Coverage\nObserved\nError"
)
selected_weights <- c("no_weight", "CW", "PPW", "OPW", "MOPW")

for (s in seq_len(nrow(representative_scenarios))) {
  
  scen <- representative_scenarios[s, ]
  
  dat_s <- compiled_m100 %>%
    filter(
      rho_xy == scen$rho_xy,
      rho_uv == scen$rho_uv,
      eta_x == scen$eta_x,
      eta_y == scen$eta_y,
      weight_type %in% selected_weights
    ) %>%
    mutate(
      coverage_true_error = abs(coverage_true - 0.95),
      coverage_obs_error  = abs(coverage_obs  - 0.95),
      plot_title = paste0(
        "M=100; rho_xy=", rho_xy,
        ", rho_uv=", rho_uv,
        ", eta_x=", eta_x,
        ", eta_y=", eta_y
      )
    )
  
  for (i in seq_len(nrow(metric_specs))) {
    
    value_col <- metric_specs$value_col[i]
    value_label <- metric_specs$value_label[i]
    red_is_better <- metric_specs$red_is_better[i]
    
    out_file <- file.path(
      figure_dir,
      paste0(
        "nw_cw_bar_heatmap_",
        value_col,
        "_rhoxy_", scen$rho_xy,
        "_rhouv_", scen$rho_uv,
        "_etax_", scen$eta_x,
        "_etay_", scen$eta_y,
        ".png"
      )
    )
    
    make_heatmap_with_cw_bar(
      dat = dat_s,
      value_col = value_col,
      value_label = value_label,
      out_file = out_file,
      red_is_better = red_is_better
    )
  }
}

#=======================================================================================
# Winner summaries
#=======================================================================================

winner_true_bias <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(abs_bias_true, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  select(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    rho_true,
    rho_obs_hat,
    best_weight_true_bias = weight_type,
    mean_estimate,
    abs_bias_true,
    coverage_true,
    coverage_obs
  )

write_csv(
  winner_true_bias,
  file.path(table_dir, "winner_by_abs_bias_true_M20.csv")
)


winner_cov_true <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(coverage_true_dist, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  select(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    rho_true,
    rho_obs_hat,
    best_weight_coverage_true = weight_type,
    mean_estimate,
    coverage_true,
    coverage_true_dist,
    coverage_obs
  )

write_csv(
  winner_cov_true,
  file.path(table_dir, "winner_by_coverage_true_M20.csv")
)

winner_cov_obs <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(coverage_obs_dist, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  select(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    rho_true,
    rho_obs_hat,
    best_weight_coverage_obs = weight_type,
    mean_estimate,
    coverage_obs,
    coverage_obs_dist,
    coverage_true
  )

write_csv(
  winner_cov_obs,
  file.path(table_dir, "winner_by_coverage_observed_M20.csv")
)


#=======================================================================================
# Compact winner table for full severity subgrouping
#=======================================================================================
winner_summary_severity <- compiled_m20 %>%
  filter(subgrouping == "L_Y_sev__K_X_sev") %>%
  group_by(
    correlation_measure,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    rho_true,
    rho_obs_hat
  ) %>%
  summarise(
    best_true_bias = weight_type[which.min(abs_bias_true)],
    min_abs_bias_true = min(abs_bias_true, na.rm = TRUE),
    
    best_cov_true = weight_type[which.min(coverage_true_dist)],
    best_coverage_true = coverage_true[which.min(coverage_true_dist)],
    
    best_cov_obs = weight_type[which.min(coverage_obs_dist)],
    best_coverage_obs = coverage_obs[which.min(coverage_obs_dist)],
    
    .groups = "drop"
  ) %>%
  arrange(
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  )

write_csv(
  winner_summary_severity,
  file.path(table_dir, "winner_summary_full_severity_M20.csv")
)

print(winner_summary_severity, n = 100)



#=======================================================================================
# Count how often each weight wins
#=======================================================================================

library(dplyr)
library(readr)

# assumes compiled_m20 already exists and has:
# abs_bias_true, abs_bias_obs, coverage_true_dist, coverage_obs_dist

winner_dir <- file.path(table_dir, "Winner_Counts_By_Scenario")
dir.create(winner_dir, showWarnings = FALSE, recursive = TRUE)

winner_true_bias_by_scenario <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(abs_bias_true, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  rename(winning_weight = weight_type) %>%
  mutate(win_criterion = "Smallest absolute bias to true rho")

winner_obs_bias_by_scenario <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(abs_bias_obs, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  rename(winning_weight = weight_type) %>%
  mutate(win_criterion = "Smallest absolute bias to observed rho")

winner_cov_true_by_scenario <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(coverage_true_dist, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  rename(winning_weight = weight_type) %>%
  mutate(win_criterion = "Coverage closest to 0.95 for true rho")

winner_cov_obs_by_scenario <- compiled_m20 %>%
  group_by(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  ) %>%
  slice_min(coverage_obs_dist, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  rename(winning_weight = weight_type) %>%
  mutate(win_criterion = "Coverage closest to 0.95 for observed rho")

count_wins_by_scenario <- function(winner_dat) {
  
  winner_dat %>%
    count(
      correlation_measure,
      rho_xy,
      rho_uv,
      eta_x,
      eta_y,
      winning_weight,
      name = "n_wins"
    ) %>%
    group_by(
      correlation_measure,
      rho_xy,
      rho_uv,
      eta_x,
      eta_y
    ) %>%
    mutate(
      total_wins = sum(n_wins),
      prop_wins = n_wins / total_wins
    ) %>%
    ungroup() %>%
    arrange(
      correlation_measure,
      rho_xy,
      rho_uv,
      eta_x,
      eta_y,
      desc(n_wins),
      winning_weight
    )
}

win_counts_true_bias_by_scenario <- count_wins_by_scenario(
  winner_true_bias_by_scenario
)

win_counts_obs_bias_by_scenario <- count_wins_by_scenario(
  winner_obs_bias_by_scenario
)

win_counts_cov_true_by_scenario <- count_wins_by_scenario(
  winner_cov_true_by_scenario
)

win_counts_cov_obs_by_scenario <- count_wins_by_scenario(
  winner_cov_obs_by_scenario
)

write_csv(
  win_counts_true_bias_by_scenario,
  file.path(winner_dir, "win_counts_true_bias_by_scenario_M20.csv")
)

write_csv(
  win_counts_obs_bias_by_scenario,
  file.path(winner_dir, "win_counts_observed_bias_by_scenario_M20.csv")
)

write_csv(
  win_counts_cov_true_by_scenario,
  file.path(winner_dir, "win_counts_coverage_true_by_scenario_M20.csv")
)

write_csv(
  win_counts_cov_obs_by_scenario,
  file.path(winner_dir, "win_counts_coverage_observed_by_scenario_M20.csv")
)


all_win_counts_by_scenario <- bind_rows(
  win_counts_true_bias_by_scenario %>%
    mutate(win_criterion = "Smallest absolute bias to true rho"),
  
  win_counts_obs_bias_by_scenario %>%
    mutate(win_criterion = "Smallest absolute bias to observed rho"),
  
  win_counts_cov_true_by_scenario %>%
    mutate(win_criterion = "Coverage closest to 0.95 for true rho"),
  
  win_counts_cov_obs_by_scenario %>%
    mutate(win_criterion = "Coverage closest to 0.95 for observed rho")
) %>%
  select(
    win_criterion,
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    winning_weight,
    n_wins,
    total_wins,
    prop_wins
  ) %>%
  arrange(
    win_criterion,
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    desc(n_wins)
  )

write_csv(
  all_win_counts_by_scenario,
  file.path(winner_dir, "all_win_counts_by_scenario_M20.csv")
)

print(all_win_counts_by_scenario, n = 200)


all_win_counts_by_scenario_wide <- all_win_counts_by_scenario %>%
  select(
    win_criterion,
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    winning_weight,
    n_wins
  ) %>%
  pivot_wider(
    names_from = winning_weight,
    values_from = n_wins,
    values_fill = 0
  ) %>%
  arrange(
    win_criterion,
    correlation_measure,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  )

write_csv(
  all_win_counts_by_scenario_wide,
  file.path(winner_dir, "all_win_counts_by_scenario_wide_M20.csv")
)

print(all_win_counts_by_scenario_wide, n = 200)

#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================




#=======================================================================================
# Pearson, Spearman, and Phi separate full-severity tables
#=======================================================================================

