####################################################################################################
#
# 05c_nhanes_lineplots.R
#
# Create publication-style line plots of age-specific NHANES Pearson correlations
# calculated in 05b_nhanes_correlation_calc_by_age_group.R.
#
# Two master figures are created:
#   1. CAL figure
#      Columns: FS, DS
#      Rows:    CAL, CAL >= 3, CAL >= 4, CAL >= 5
#
#   2. PD figure
#      Columns: FS, DS
#      Rows:    PD, PD >= 4, PD >= 5, PD >= 6
#
# Within each panel:
#   x-axis: Age group (30-49, 50-64, 65+)
#   y-axis: Pearson correlation (0 to 0.20)
#   lines:  No Weight, CW, OPW, PPW
#
# IMPORTANT:
#   The 05b results contain threshold-specific L variables only:
#       CAL >= 3, CAL >= 4, CAL >= 5
#       PD  >= 4, PD  >= 5, PD  >= 6
#
#   Therefore, the unthresholded CAL and PD rows can be defined directly only for
#   No Weight and CW, because those two weighting schemes do not require an L threshold.
#   PPW and OPW are threshold-specific and are consequently not drawn in the
#   unthresholded CAL/PD row.
#
####################################################################################################


# ==================================================================================================
# 1. LOAD PACKAGES
# ==================================================================================================

library(dplyr)
library(ggplot2)
library(here)


# ==================================================================================================
# 2. LOAD AGE-SPECIFIC PEARSON RESULTS FROM SCRIPT 05b
# ==================================================================================================

results_file <- here::here(
  "results",
  "correlations",
  "by_age",
  "nhanes_pearson_results_by_age.rds"
)

all_pearson_age_results <- readRDS(results_file)


# ==================================================================================================
# 3. DEFINE DISPLAY ORDER AND APPEARANCE
# ==================================================================================================

age_order <- c(
  "30-49",
  "50-64",
  "65+"
)

caries_order <- c(
  "FS",
  "FSI",
  "DS",
  "DSI"
)

weight_order <- c(
  "no_weight",
  "CW",
  "OPW",
  "PPW"
)

weight_labels <- c(
  no_weight = "NW",
  CW        = "CW",
  OPW       = "OPW",
  PPW       = "PPW"
)

# Color-blind-friendly colors.
# Linetypes provide a second visual encoding for black-and-white printing.
weight_colors <- c(
  no_weight = "#000000",
  CW        = "#0072B2",
  OPW       = "#D55E00",
  PPW       = "#009E73"
)

weight_linetypes <- c(
  no_weight = "solid",
  CW        = "dashed",
  OPW       = "dotdash",
  PPW       = "dotted"
)

weight_shapes <- c(
  no_weight = 16,
  CW        = 17,
  OPW       = 15,
  PPW       = 18
)


# ==================================================================================================
# 4. PREPARE THRESHOLD-SPECIFIC RESULTS
# ==================================================================================================

threshold_data <- all_pearson_age_results %>%
  filter(
    K %in% caries_order,
    weight_type %in% weight_order,
    Age %in% age_order
  ) %>%
  mutate(
    Age = factor(
      Age,
      levels = age_order
    ),

    K = factor(
      K,
      levels = caries_order
    ),

    weight_type = factor(
      weight_type,
      levels = weight_order
    ),

    perio_family = case_when(
      Y == "maxCAL" ~ "CAL",
      Y == "maxPD"  ~ "PD",
      TRUE          ~ NA_character_
    ),

    perio_measure = case_when(
      L == "Cal>=3" ~ "CAL >= 3",
      L == "Cal>=4" ~ "CAL >= 4",
      L == "Cal>=5" ~ "CAL >= 5",
      L == "Pd>=4"  ~ "PD >= 4",
      L == "Pd>=5"  ~ "PD >= 5",
      L == "Pd>=6"  ~ "PD >= 6",
      TRUE          ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(perio_family),
    !is.na(perio_measure)
  )


# ==================================================================================================
# 5. CREATE UNTHRESHOLDED CAL AND PD ROWS
# ==================================================================================================

# No Weight and CW do not depend on the periodontal threshold L.
# Script 05b therefore contains the same estimate repeatedly for these methods
# across the three L thresholds belonging to a given periodontal outcome.
#
# We collapse those repeated estimates into one unthresholded CAL or PD estimate.
#
# PPW and OPW are not included here because their definitions depend on L.

baseline_data <- threshold_data %>%
  filter(
    as.character(weight_type) %in% c(
      "no_weight",
      "CW"
    )
  ) %>%
  group_by(
    Age,
    K,
    perio_family,
    weight_type
  ) %>%
  summarise(
    estimate = first(estimate),
    .groups = "drop"
  ) %>%
  mutate(
    perio_measure = perio_family
  )


# ==================================================================================================
# 6. COMBINE DATA USED FOR PLOTTING
# ==================================================================================================

plot_data <- bind_rows(
  threshold_data %>%
    select(
      Age,
      K,
      weight_type,
      estimate,
      perio_family,
      perio_measure
    ),
  baseline_data %>%
    select(
      Age,
      K,
      weight_type,
      estimate,
      perio_family,
      perio_measure
    )
)


# ==================================================================================================
# 7. HELPER FUNCTION FOR MASTER PLOTS
# ==================================================================================================

make_periodontal_plot <- function(
    dat,
    periodontal_family,
    row_order,
    plot_title
) {

  plot_dat <- dat %>%
    filter(
      perio_family == periodontal_family,
      perio_measure %in% row_order
    ) %>%
    mutate(
      perio_measure = factor(
        perio_measure,
        levels = row_order
      )
    )

  ggplot(
    plot_dat,
    aes(
      x = Age,
      y = estimate,
      group = weight_type,
      color = weight_type,
      linetype = weight_type,
      shape = weight_type
    )
  ) +
    geom_hline(
      yintercept = 0,
      color = "grey70",
      linewidth = 0.3
    ) +

    geom_line(
      linewidth = 0.8,
      na.rm = TRUE
    ) +

    geom_point(
      size = 2.1,
      na.rm = TRUE
    ) +

    facet_grid(
      rows = vars(perio_measure),
      cols = vars(K),
      scales = "fixed"
    ) +

    scale_color_manual(
      name = "Weight: ",
      values = weight_colors,
      breaks = weight_order,
      labels = weight_labels
    ) +

    scale_linetype_manual(
      name = "Weight: ",
      values = weight_linetypes,
      breaks = weight_order,
      labels = weight_labels
    ) +

    scale_shape_manual(
      name = "Weight: ",
      values = weight_shapes,
      breaks = weight_order,
      labels = weight_labels
    ) +

    scale_y_continuous(
      breaks = seq(
        -0.05,
        0.30,
        by = 0.1
      ),
      labels = scales::label_number(
        accuracy = 0.01
      )
    ) +

    coord_cartesian(
      ylim = c(
        -0.05,
        0.30
      )
    ) +

    labs(
      x = "Age Group",
      y = "Pearson Correlation",
      title = plot_title
    ) +

    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE
      ),
      linetype = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    ) +

    theme_bw(
      base_size = 12
    ) +

    theme(
      axis.text.x = element_text(
        size = 12,
        color = "black"
      ),

      axis.text.y = element_text(
        size = 12,
        color = "black"
      ),

      axis.title.x = element_text(
        size = 13
      ),

      axis.title.y = element_text(
        size = 13
      ),

      strip.text = element_text(
        size = 12
      ),

      # Keep row facet labels on the right, but horizontal.
      strip.text.y.right = element_text(
        angle = 270,
        size = 12
      ),

      legend.title = element_text(
        size = 13
      ),

      legend.key.width = grid::unit(0.75, units = "in"),

      legend.text = element_text(
        size = 12
      ),

      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),

      legend.position = "bottom",

      # Same clean panel treatment as the custom heatmaps.
      panel.grid = element_blank()
    )
}


# ==================================================================================================
# 8. CREATE CAL MASTER PLOT
# ==================================================================================================

cal_plot <- make_periodontal_plot(
  dat = plot_data,
  periodontal_family = "CAL",
  row_order = c(
    "CAL >= 3",
    "CAL >= 4",
    "CAL >= 5"
  ),
  plot_title = "Pearson Correlations Between Caries and Clinical Attachment Loss"
)

cal_plot


# ==================================================================================================
# 9. CREATE PD MASTER PLOT
# ==================================================================================================

pd_plot <- make_periodontal_plot(
  dat = plot_data,
  periodontal_family = "PD",
  row_order = c(
    "PD >= 4",
    "PD >= 5",
    "PD >= 6"
  ),
  plot_title = "Pearson Correlations Between Caries and Pocket Depth"
)

pd_plot


# ==================================================================================================
# 10. SAVE FIGURES
# ==================================================================================================

figures_dir <- here::here(
  "figures"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# PNG files for quick viewing / manuscript drafting.

ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_CAL_300dpi.png"
  ),
  plot = cal_plot,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_PD_300dpi.png"
  ),
  plot = pd_plot,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)


# PDF files preserve vector graphics for publication.

ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_CAL.pdf"
  ),
  plot = cal_plot,
  width = 8,
  height = 5,
  units = "in"
)

ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_PD.pdf"
  ),
  plot = pd_plot,
  width = 8,
  height = 5,
  units = "in"
)
