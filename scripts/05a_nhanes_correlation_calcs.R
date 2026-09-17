####################################################################################################
#
# 05_nhanes_correlation_calcs.R
#
# Calculate weighted tooth-level associations between periodontal
# measures and caries measures in NHANES.
#
# Continuous periodontal scores:
#   - maxCAL
#   - maxPD
#
# Binary periodontal indicators:
#   - CAL >= 3, 4, 5
#   - PD  >= 4, 5, 6
#
# Association measures:
#   - Pearson correlation for continuous scores
#   - Phi coefficient for binary indicators
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

# Weighting functions
source(
  here::here("R", "Weighting_Equations.R")
)

# Association functions
source(
  here::here("R", "Correlation_Equations.R")
)


# ==================================================================================================
# 3. LOAD TOOTH-LEVEL NHANES DATA
# ==================================================================================================

load(
  here::here("data", "Tooth_Level_Data.RData")
)


# ==================================================================================================
# 4. CLEAN ANALYSIS DATA
# ==================================================================================================

# Remove rows containing missing values.
#
# This preserves the same complete-case analysis used in the
# original version of the script.

df_clean <- na.omit(df_tooth_level)


# Keep:
#
#   1. Subjects contributing more than one tooth-level observation.
#   2. Subjects with more than 9 teeth.
#
# SEQN is the subject/cluster identifier.

df_clean <- df_clean %>%
  group_by(SEQN) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(Tooth_Count > 9)


# ==================================================================================================
# 5. CREATE BINARY PERIODONTAL INDICATORS
# ==================================================================================================

# These indicators are used for the phi coefficient analysis.
#
# The original continuous maxCAL and maxPD variables remain unchanged
# and are used for the Pearson correlation analysis below.

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
# 6. DEFINE COMMON ANALYSIS SETTINGS
# ==================================================================================================

# Caries variables.
#
# These variables are used both as the observed X variables and
# as K in the weighting equations.

X_vars <- c(
  "DS",
  "FS",
  "DSI",
  "FSI"
)

K_vars <- X_vars


# Weighting methods.

weight_types <- c(
  "no_weight",
  "CW",
  "PPW",
  "OPW",
  "MOPW"
)


# ==================================================================================================
# 7. HELPER FUNCTION: RUN ALL FIVE WEIGHTING METHODS
# ==================================================================================================

# run_all_weights() is simply a convenience wrapper around
# run_association_grid().
#
# It does NOT change any of the association or weighting mathematics.
# It only prevents us from repeating the same lapply()/rbind code for
# Pearson and phi analyses.

run_all_weights <- function(dat,
                            Y_vars,
                            X_vars,
                            K_vars,
                            L_vars,
                            association_type) {

  results <- lapply(
    weight_types,
    function(wtype) {

      run_association_grid(
        dat = dat,
        Y_vars = Y_vars,
        X_vars = X_vars,
        clusterID = "SEQN",
        clusterSize = "Tooth_Count",
        K_vars = K_vars,
        L_vars = L_vars,
        weight_type = wtype,
        association_type = association_type
      )

    }
  )

  do.call(rbind, results)
}


# ==================================================================================================
# 8. PEARSON CORRELATION: CONTINUOUS PERIODONTAL SCORES
# ==================================================================================================

# We now use the continuous periodontal scores directly.
#
# maxCAL is repeated three times because each version corresponds to a
# different L variable used by the weighting equations:
#
#   CAL >= 3
#   CAL >= 4
#   CAL >= 5
#
# Likewise, maxPD is paired with the three PD threshold indicators.
#
# The association itself is calculated using the CONTINUOUS variables
# maxCAL and maxPD. The binary L variables are only used in constructing
# the weighting scheme.

pearson_Y_vars <- c(
  "maxCAL",
  "maxCAL",
  "maxCAL",
  "maxPD",
  "maxPD",
  "maxPD"
)

pearson_L_vars <- c(
  "Cal>=3",
  "Cal>=4",
  "Cal>=5",
  "Pd>=4",
  "Pd>=5",
  "Pd>=6"
)


# Run Pearson correlation under all five weighting schemes.

all_pearson_results <- run_all_weights(
  dat = df_clean,
  Y_vars = pearson_Y_vars,
  X_vars = X_vars,
  K_vars = K_vars,
  L_vars = pearson_L_vars,
  association_type = "pearson"
)


# View results.

print(all_pearson_results)


# ==================================================================================================
# 9. SAVE ASSOCIATION RESULTS
# ==================================================================================================

# Create a results directory specifically for the correlation analyses.

correlation_results_dir <- here::here(
  "results",
  "correlations"
)

dir.create(
  correlation_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Save complete R objects.
saveRDS(
  all_pearson_results,
  file = file.path(
    correlation_results_dir,
    "nhanes_pearson_results.rds"
  )
)

# Also save human-readable CSV versions.
write.csv(
  all_pearson_results,
  file = file.path(
    correlation_results_dir,
    "nhanes_pearson_results.csv"
  ),
  row.names = FALSE
)

# ==================================================================================================
# 10. CREATE WIDE PEARSON RESULTS TABLE
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


pearson_wide_tab <- all_pearson_results %>%
  
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
    
    # Combine estimate and standard error into one display column
    result = ifelse(
      is.na(estimate) | is.na(se),
      "-",
      sprintf("%.3f (%.3f)", estimate, se)
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
    factor(perio_group, levels = c("CAL", "PD")),
    threshold
  )


print(
  pearson_wide_tab,
  n = 100
)


# ==================================================================================================
# 11. CREATE PEARSON LATEX TABLE
# ==================================================================================================

# Helper for LaTeX line breaks
latex_break <- paste0("\\", "\\")


pearson_latex_rows <- pearson_wide_tab %>%
  
  group_by(K) %>%
  
  mutate(
    
    # Only display the caries variable on the first row
    K_tex = ifelse(
      row_number() == 1,
      as.character(K),
      ""
    ),
    
    # Create periodontal labels
    perio_tex = case_when(
      perio_group == "CAL" ~ paste0("CAL$\\geq$", threshold),
      perio_group == "PD"  ~ paste0("PD$\\geq$", threshold),
      TRUE                 ~ perio_group
    ),
    
    # Construct each LaTeX row
    row_tex = paste(
      K_tex,
      "&",
      perio_tex,
      "&&",
      no_weight,
      "&",
      CW,
      "&",
      PPW,
      "&",
      OPW,
      "&",
      MOPW,
      latex_break
    )
  ) %>%
  
  ungroup()


# Add a blank row before the first PD row within each caries variable
pearson_latex_rows <- pearson_latex_rows %>%
  
  group_by(K) %>%
  
  mutate(
    
    first_PD = perio_group == "PD" & !duplicated(perio_group),
    
    row_tex = ifelse(
      first_PD,
      paste(
        latex_break,
        row_tex,
        sep = "\n"
      ),
      row_tex
    )
  ) %>%
  
  ungroup()


# Collapse rows within each caries-variable block
pearson_latex_rows <- pearson_latex_rows %>%
  
  group_by(K) %>%
  
  summarise(
    block_tex = paste(
      row_tex,
      collapse = "\n"
    ),
    .groups = "drop"
  ) %>%
  
  pull(block_tex)


# Combine all caries-variable blocks
pearson_latex_table <- paste(
  pearson_latex_rows,
  collapse = paste0(
    "\n",
    latex_break,
    "\n"
  )
)


# Print LaTeX rows to console
cat(pearson_latex_table)


# Save to file
writeLines(
  pearson_latex_table,
  con = file.path(
    correlation_results_dir,
    "nhanes_pearson_table.tex"
  )
)
