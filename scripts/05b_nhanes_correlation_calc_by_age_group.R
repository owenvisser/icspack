####################################################################################################
#
# 05b_nhanes_correlation_calc_by_age_group.R
#
# Calculate weighted Pearson correlations between periodontal
# measures and caries measures within NHANES age groups.
#
# Age groups:
#   - 30-49
#   - 50-64
#   - 65+
#
# Continuous periodontal scores:
#   - maxCAL
#   - maxPD
#
# Association measure:
#   - Pearson correlation
#
# Weighting methods:
#   - No weighting
#   - Cluster weighting (CW)
#   - Pairwise probability weighting (PPW)
#   - Outcome probability weighting (OPW)
#   - Modified outcome probability weighting (MOPW)
#
####################################################################################################


# ==================================================================================================
# 1. LOAD PACKAGES
# ==================================================================================================

library(dplyr)
library(tidyr)
library(stringr)


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

# Tooth-level analysis data
load(
  here::here("data", "Tooth_Level_Data.RData")
)

# Person-level demographic data
load(
  here::here("data", "Person_Level_Data.RData")
)

# Age is stored at the person level, so we attach it to each tooth-level
# observation using the participant identifier SEQN.
#
# Keep only one SEQN-Age record per participant before joining.

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

df_clean <- df_tooth_level %>%
  na.omit() %>%
  group_by(SEQN) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(Tooth_Count > 9)


# ==================================================================================================
# 5. CREATE PERIODONTAL THRESHOLD INDICATORS
# ==================================================================================================

# These are used as L variables in the weighting equations.
# Pearson itself is calculated using continuous maxCAL and maxPD.

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

# Caries variables

X_vars <- c(
  "DS_score",
  "FS_score",
  "DSI_score",
  "FSI_score"
)

K_vars <- c(
  "DS",
  "FS",
  "DSI",
  "FSI"
)


# Continuous periodontal outcomes

pearson_Y_vars <- c(
  "maxCAL",
  "maxCAL",
  "maxCAL",
  "maxPD",
  "maxPD",
  "maxPD"
)


# Threshold variables used in the weighting equations

pearson_L_vars <- c(
  "Cal>=3",
  "Cal>=4",
  "Cal>=5",
  "Pd>=4",
  "Pd>=5",
  "Pd>=6"
)


# Weighting methods

weight_types <- c(
  "no_weight",
  "CW",
  "PPW",
  "OPW"
)


# Age groups

age_groups <- c(
  "30-49",
  "50-64",
  "65+"
)


# ==================================================================================================
# 7. HELPER FUNCTION: RUN ALL WEIGHTS
# ==================================================================================================

run_all_weights <- function(dat) {

  results <- lapply(
    weight_types,
    function(wtype) {

      run_association_grid(
        dat = dat,
        Y_vars = pearson_Y_vars,
        X_vars = X_vars,
        clusterID = "SEQN",
        clusterSize = "Tooth_Count",
        K_vars = K_vars,
        L_vars = pearson_L_vars,
        weight_type = wtype,
        association_type = "pearson"
      )

    }
  )

  do.call(rbind, results)
}


# ==================================================================================================
# 8. RUN PEARSON ANALYSIS WITHIN EACH AGE GROUP
# ==================================================================================================

# Store one results data frame per age group.

pearson_results_by_age <- list()


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

  pearson_results_by_age[[age_group]] <- age_results
}


# ==================================================================================================
# 9. COMBINE RESULTS
# ==================================================================================================

all_pearson_age_results <- bind_rows(
  pearson_results_by_age
)


head(
  all_pearson_age_results,
  n = 200
)


# ==================================================================================================
# 10. CREATE NICE WIDE TABLE FOR EACH AGE GROUP
# ==================================================================================================

weight_order <- c(
  "no_weight",
  "CW",
  "PPW",
  "OPW",
  "MOPW"
)

K_order <- c(
  "FS",
  "FSI",
  "DS",
  "DSI"
)


make_pearson_table <- function(results) {

  results %>%

    select(
      Y,
      L,
      K,
      weight_type,
      estimate,
      se
    ) %>%

    mutate(

      weight_type = factor(
        weight_type,
        levels = weight_order
      ),

      K = factor(
        K,
        levels = K_order
      ),

      perio_group = case_when(
        Y == "maxCAL" ~ "CAL",
        Y == "maxPD"  ~ "PD",
        TRUE          ~ Y
      ),

      threshold = as.numeric(
        str_extract(L, "\\d+")
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
      K,
      perio_group,
      threshold,
      weight_type,
      result
    ) %>%

    pivot_wider(
      names_from = weight_type,
      values_from = result
    ) %>%

    arrange(
      K,
      factor(
        perio_group,
        levels = c("CAL", "PD")
      ),
      threshold
    )
}


# Create one table per age group.

pearson_tables_by_age <- lapply(
  pearson_results_by_age,
  make_pearson_table
)


# ==================================================================================================
# 11. PRINT EACH AGE-SPECIFIC TABLE
# ==================================================================================================

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
    pearson_tables_by_age[[age_group]],
    n = 100
  )
}


# ==================================================================================================
# 12. SAVE RESULTS
# ==================================================================================================

age_results_dir <- here::here(
  "results",
  "correlations",
  "by_age"
)

dir.create(
  age_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Save the combined result object.

saveRDS(
  all_pearson_age_results,
  file = file.path(
    age_results_dir,
    "nhanes_pearson_results_by_age.rds"
  )
)


# Save the combined results as CSV.

write.csv(
  all_pearson_age_results,
  file = file.path(
    age_results_dir,
    "nhanes_pearson_results_by_age.csv"
  ),
  row.names = FALSE
)


# Save one wide table per age group.

for (age_group in age_groups) {

  # Replace characters that are awkward in filenames.
  age_file_label <- str_replace_all(
    age_group,
    "\\+",
    "plus"
  )

  write.csv(
    pearson_tables_by_age[[age_group]],
    file = file.path(
      age_results_dir,
      paste0(
        "nhanes_pearson_table_age_",
        age_file_label,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# ==================================================================================================
# 13. SAVE AGE-SPECIFIC PEARSON TABLES AS IMAGES
# ==================================================================================================

library(gt)

# Directory for figure output
figures_dir <- here::here("figures")

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Helper function to turn a Pearson table into a gt table
make_pearson_gt_table <- function(tab, age_group) {

  tab_display <- tab %>%
    mutate(
      `Caries Measure` = as.character(K),
      `Periodontal Measure` = paste0(perio_group, " \u2265 ", threshold)
    ) %>%
    select(
      `Caries Measure`,
      `Periodontal Measure`,
      no_weight,
      CW,
      PPW,
      OPW,
      MOPW
    )

  gt(tab_display) %>%
    tab_header(
      title = md("**NHANES Pearson Correlations**"),
      subtitle = paste("Age group:", age_group)
    ) %>%
    cols_label(
      no_weight = "No Weight",
      CW = "CW",
      PPW = "PPW",
      OPW = "OPW",
      MOPW = "MOPW"
    ) %>%
    tab_spanner(
      label = "Weighting Method",
      columns = c(no_weight, CW, PPW, OPW, MOPW)
    ) %>%
    cols_align(
      align = "center",
      columns = everything()
    ) %>%
    tab_options(
      table.font.size = px(12),
      data_row.padding = px(4),
      heading.align = "center"
    )
}


# Save one image per age group
for (age_group in age_groups) {

  age_file_label <- age_group %>%
    stringr::str_replace_all("\\+", "plus") %>%
    stringr::str_replace_all("-", "_")

  gt_table <- make_pearson_gt_table(
    tab = pearson_tables_by_age[[age_group]],
    age_group = age_group
  )

  gtsave(
    data = gt_table,
    filename = file.path(
      figures_dir,
      paste0("nhanes_pearson_table_age_", age_file_label, ".png")
    )
  )
}
