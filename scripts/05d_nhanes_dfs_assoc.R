####################################################################################################
#
# 05d_nhanes_dfs_assoc.R
#
# Calculate age-specific weighted Pearson correlations between periodontal
# measures and DFS / DFSI measures using the original periodontal dichotomizations.
#
# Continuous periodontal outcomes:
#   - maxCAL
#   - maxPD
#
# Continuous caries outcomes:
#   - DFS_score
#   - DFSI_score
#
# Binary caries subgroup variables:
#   - DFS
#   - DFSI
#
# Periodontal threshold variables used in the weighting equations:
#   - CAL >= 3, CAL >= 4, CAL >= 5
#   - PD  >= 4, PD  >= 5, PD  >= 6
#
# Age groups:
#   - 30-49
#   - 50-64
#   - 65+
#
# Weighting methods:
#   - No weighting
#   - Cluster weighting (CW)
#   - Pairwise probability weighting (PPW)
#   - Outcome probability weighting (OPW)
#
####################################################################################################


# ==================================================================================================
# 1. LOAD PACKAGES
# ==================================================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(here)


# ==================================================================================================
# 2. LOAD FUNCTIONS
# ==================================================================================================

source(
  here::here("R", "Weighting_Equations.R")
)

source(
  here::here("R", "Correlation_Equations.R")
)


# ==================================================================================================
# 3. LOAD NHANES DATA
# ==================================================================================================

load(
  here::here("data", "Tooth_Level_Data.RData")
)

load(
  here::here("data", "Person_Level_Data.RData")
)


# Age is stored at the person level, so attach it to each tooth-level observation.
age_data <- df_person_level %>%
  select(
    SEQN,
    Age
  ) %>%
  distinct()


df_tooth_level <- df_tooth_level %>%
  left_join(
    age_data,
    by = "SEQN"
  )


# ==================================================================================================
# 4. CLEAN ANALYSIS DATA
# ==================================================================================================

# Keep only observations complete for variables needed in this analysis.
required_vars <- c(
  "SEQN",
  "Age",
  "Tooth_Count",
  "WTMEC4YR",
  "maxCAL",
  "maxPD",
  "DFS_score",
  "DFSI_score",
  "DFS",
  "DFSI"
)


df_clean <- df_tooth_level %>%
  filter(
    if_all(
      all_of(required_vars),
      ~ !is.na(.x)
    )
  ) %>%
  group_by(SEQN) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(Tooth_Count > 9)


# ==================================================================================================
# 5. CREATE ORIGINAL PERIODONTAL THRESHOLD INDICATORS
# ==================================================================================================

# These are the same periodontal dichotomizations used in 05b.
# The continuous Pearson outcome remains maxCAL or maxPD; these binary variables
# are used only as L in the weighting equations.

df_clean <- df_clean %>%
  mutate(
    `Cal>=3` = as.integer(maxCAL >= 3),
    `Cal>=4` = as.integer(maxCAL >= 4),
    `Cal>=5` = as.integer(maxCAL >= 5),
    `Pd>=4`  = as.integer(maxPD  >= 4),
    `Pd>=5`  = as.integer(maxPD  >= 5),
    `Pd>=6`  = as.integer(maxPD  >= 6)
  )


# ==================================================================================================
# 6. DEFINE ANALYSIS SETTINGS
# ==================================================================================================

# Continuous DFS / DFSI scores used in the Pearson association.
X_vars <- c(
  "DFS_score",
  "DFSI_score"
)


# Binary DFS / DFSI indicators used as K in the weighting equations.
K_vars <- c(
  "DFS",
  "DFSI"
)


# Continuous periodontal outcomes.
pearson_Y_vars <- c(
  "maxCAL",
  "maxCAL",
  "maxCAL",
  "maxPD",
  "maxPD",
  "maxPD"
)


# Original periodontal dichotomizations used as L in the weighting equations.
pearson_L_vars <- c(
  "Cal>=3",
  "Cal>=4",
  "Cal>=5",
  "Pd>=4",
  "Pd>=5",
  "Pd>=6"
)


weight_types <- c(
  "no_weight",
  "CW",
  "PPW",
  "OPW"
)


age_groups <- c(
  "30-49",
  "50-64",
  "65+"
)


# ==================================================================================================
# 7. HELPER FUNCTION: RUN WEIGHTED PEARSON GRID
# ==================================================================================================

run_nhanes_dfs_grid <- function(
    dat,
    Y_vars,
    X_vars,
    clusterID,
    clusterSize,
    K_vars,
    L_vars,
    surveyWeight,
    weight_type = "no_weight"
) {

  cid <- dat[[clusterID]]
  cis <- dat[[clusterSize]]
  sw  <- dat[[surveyWeight]]

  results <- list()
  counter <- 1

  if (length(Y_vars) != length(L_vars)) {
    stop("Y_vars and L_vars must match length")
  }

  if (length(X_vars) != length(K_vars)) {
    stop("X_vars and K_vars must match length")
  }

  for (l in seq_along(Y_vars)) {

    for (k in seq_along(X_vars)) {

      ics_weight <- make_weights(
        clusterID = cid,
        clusterSize = cis,
        K = dat[[K_vars[k]]],
        L = dat[[L_vars[l]]],
        weight_type = weight_type
      )

      # Combine NHANES sampling weights with the informative cluster/subgroup weight.
      omega <- sw * ics_weight

      fit <- pearson_association(
        clusterID = cid,
        Y = dat[[Y_vars[l]]],
        X = dat[[X_vars[k]]],
        omega = omega
      )

      results[[counter]] <- data.frame(
        Y = Y_vars[l],
        L = L_vars[l],
        X = X_vars[k],
        K = K_vars[k],
        weight_type = weight_type,
        association_type = "pearson",
        estimate = fit$rho,
        variance = fit$variance,
        se = fit$se,
        lower_95 = unname(fit$ci[1]),
        upper_95 = unname(fit$ci[2])
      )

      counter <- counter + 1
    }
  }

  bind_rows(results)
}


run_all_weights <- function(dat) {

  results <- lapply(
    weight_types,
    function(wtype) {

      run_nhanes_dfs_grid(
        dat = dat,
        Y_vars = pearson_Y_vars,
        X_vars = X_vars,
        clusterID = "SEQN",
        clusterSize = "Tooth_Count",
        K_vars = K_vars,
        L_vars = pearson_L_vars,
        surveyWeight = "WTMEC4YR",
        weight_type = wtype
      )
    }
  )

  bind_rows(results)
}


# ==================================================================================================
# 8. RUN ANALYSIS WITHIN EACH AGE GROUP
# ==================================================================================================

dfs_results_by_age <- list()


for (age_group in age_groups) {

  cat(
    "\nRunning age group:",
    age_group,
    "\n"
  )

  df_age <- df_clean %>%
    filter(Age == age_group)

  age_results <- run_all_weights(
    dat = df_age
  )

  age_results$Age <- age_group

  dfs_results_by_age[[age_group]] <- age_results
}


# ==================================================================================================
# 9. COMBINE AND LABEL RESULTS
# ==================================================================================================

all_dfs_age_results <- bind_rows(
  dfs_results_by_age
) %>%
  mutate(
    Periodontal_Family = case_when(
      Y == "maxCAL" ~ "CAL",
      Y == "maxPD"  ~ "PD",
      TRUE          ~ NA_character_
    ),

    Periodontal_Measure = case_when(
      L == "Cal>=3" ~ "CAL >= 3",
      L == "Cal>=4" ~ "CAL >= 4",
      L == "Cal>=5" ~ "CAL >= 5",
      L == "Pd>=4"  ~ "PD >= 4",
      L == "Pd>=5"  ~ "PD >= 5",
      L == "Pd>=6"  ~ "PD >= 6",
      TRUE          ~ NA_character_
    ),

    Caries_Measure = case_when(
      K == "DFS"  ~ "DFS",
      K == "DFSI" ~ "DFSI",
      TRUE        ~ K
    )
  )


print(
  all_dfs_age_results,
  n = 200
)


# ==================================================================================================
# 10. CREATE AGE-SPECIFIC WIDE TABLES
# ==================================================================================================

weight_order <- c(
  "no_weight",
  "CW",
  "OPW",
  "PPW"
)

caries_order <- c(
  "DFS",
  "DFSI"
)

perio_measure_order <- c(
  "CAL >= 3",
  "CAL >= 4",
  "CAL >= 5",
  "PD >= 4",
  "PD >= 5",
  "PD >= 6"
)


make_dfs_table <- function(results) {

  results %>%
    select(
      Caries_Measure,
      Periodontal_Family,
      Periodontal_Measure,
      weight_type,
      estimate,
      se
    ) %>%
    mutate(
      weight_type = factor(
        weight_type,
        levels = weight_order
      ),
      Caries_Measure = factor(
        Caries_Measure,
        levels = caries_order
      ),
      Periodontal_Measure = factor(
        Periodontal_Measure,
        levels = perio_measure_order
      ),
      result = ifelse(
        is.na(estimate) | is.na(se),
        "-",
        sprintf(
          "%.3f (%.3f)",
          estimate,
          se
        )
      )
    ) %>%
    select(
      Caries_Measure,
      Periodontal_Family,
      Periodontal_Measure,
      weight_type,
      result
    ) %>%
    pivot_wider(
      names_from = weight_type,
      values_from = result
    ) %>%
    arrange(
      Caries_Measure,
      Periodontal_Measure
    )
}


dfs_tables_by_age <- lapply(
  dfs_results_by_age,
  function(x) {
    x %>%
      mutate(
        Periodontal_Family = case_when(
          Y == "maxCAL" ~ "CAL",
          Y == "maxPD"  ~ "PD",
          TRUE          ~ NA_character_
        ),
        Periodontal_Measure = case_when(
          L == "Cal>=3" ~ "CAL >= 3",
          L == "Cal>=4" ~ "CAL >= 4",
          L == "Cal>=5" ~ "CAL >= 5",
          L == "Pd>=4"  ~ "PD >= 4",
          L == "Pd>=5"  ~ "PD >= 5",
          L == "Pd>=6"  ~ "PD >= 6",
          TRUE          ~ NA_character_
        ),
        Caries_Measure = case_when(
          K == "DFS"  ~ "DFS",
          K == "DFSI" ~ "DFSI",
          TRUE        ~ K
        )
      ) %>%
      make_dfs_table()
  }
)


for (age_group in age_groups) {

  cat(
    "\n\n============================================================\n"
  )

  cat(
    "AGE GROUP:",
    age_group,
    "\n"
  )

  cat(
    "============================================================\n\n"
  )

  print(
    dfs_tables_by_age[[age_group]],
    n = 100
  )
}


# ==================================================================================================
# 11. SAVE NUMERIC RESULTS
# ==================================================================================================

results_dir <- here::here(
  "results",
  "correlations",
  "dfs_by_age"
)


dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  all_dfs_age_results,
  file = file.path(
    results_dir,
    "nhanes_dfs_periodontal_threshold_results_by_age.rds"
  )
)


write.csv(
  all_dfs_age_results,
  file = file.path(
    results_dir,
    "nhanes_dfs_periodontal_threshold_results_by_age.csv"
  ),
  row.names = FALSE
)


for (age_group in age_groups) {

  age_file_label <- str_replace_all(
    age_group,
    "\\+",
    "plus"
  )

  write.csv(
    dfs_tables_by_age[[age_group]],
    file = file.path(
      results_dir,
      paste0(
        "nhanes_dfs_periodontal_threshold_table_age_",
        age_file_label,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# ==================================================================================================
# 12. PREPARE DATA FOR AGE-GROUP LINE PLOTS
# ==================================================================================================

plot_data <- all_dfs_age_results %>%
  filter(
    !is.na(Periodontal_Family),
    !is.na(Periodontal_Measure),
    Caries_Measure %in% caries_order,
    weight_type %in% weight_order,
    Age %in% age_groups
  ) %>%
  mutate(
    Age = factor(
      Age,
      levels = age_groups
    ),
    Caries_Measure = factor(
      Caries_Measure,
      levels = caries_order
    ),
    weight_type = factor(
      weight_type,
      levels = weight_order
    )
  )


weight_labels <- c(
  no_weight = "NW",
  CW        = "CW",
  OPW       = "OPW",
  PPW       = "PPW"
)

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
# 13. HELPER FUNCTION FOR 3 x 2 MASTER PLOTS
# ==================================================================================================

make_dfs_periodontal_plot <- function(
    dat,
    periodontal_family,
    row_order,
    plot_title
) {

  plot_dat <- dat %>%
    filter(
      Periodontal_Family == periodontal_family,
      Periodontal_Measure %in% row_order
    ) %>%
    mutate(
      Periodontal_Measure = factor(
        Periodontal_Measure,
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
      rows = vars(Periodontal_Measure),
      cols = vars(Caries_Measure),
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
      limits = c(-0.05, 0.30),
      breaks = c(-0.05, 0.05, 0.15, 0.25),
      expand = c(0, 0),
      labels = scales::label_number(
        accuracy = 0.01
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
      strip.text.y.right = element_text(
        angle = 270,
        size = 12
      ),
      legend.title = element_text(
        size = 13
      ),
      legend.key.width = grid::unit(
        0.75,
        units = "in"
      ),
      legend.text = element_text(
        size = 12
      ),
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}


# ==================================================================================================
# 14. CREATE CAL AND PD MASTER PLOTS
# ==================================================================================================

cal_plot <- make_dfs_periodontal_plot(
  dat = plot_data,
  periodontal_family = "CAL",
  row_order = c(
    "CAL >= 3",
    "CAL >= 4",
    "CAL >= 5"
  ),
  plot_title = "Pearson Correlations Between DFS and Clinical Attachment Loss"
)


pd_plot <- make_dfs_periodontal_plot(
  dat = plot_data,
  periodontal_family = "PD",
  row_order = c(
    "PD >= 4",
    "PD >= 5",
    "PD >= 6"
  ),
  plot_title = "Pearson Correlations Between DFS and Pocket Depth"
)


cal_plot
pd_plot


# ==================================================================================================
# 15. SAVE MASTER PLOTS
# ==================================================================================================

figures_dir <- here::here(
  "figures"
)


dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_DFS_CAL_300dpi.png"
  ),
  plot = cal_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_DFS_PD_300dpi.png"
  ),
  plot = pd_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)


ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_DFS_CAL.pdf"
  ),
  plot = cal_plot,
  width = 7,
  height = 5,
  units = "in"
)


ggsave(
  filename = file.path(
    figures_dir,
    "nhanes_age_correlations_DFS_PD.pdf"
  ),
  plot = pd_plot,
  width = 7,
  height = 5,
  units = "in"
)
