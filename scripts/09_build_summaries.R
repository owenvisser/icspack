####################################################################################################
#
# SIMULATION CORRELATION SUMMARIES
#
# Summarizes Pearson correlation simulation results for M = 20.
#
# Input:
#   icspack/results/sim/correlations/Pearson/
#
# Output:
#   icspack/results/sim/correlations/Pearson/summaries/
#
####################################################################################################


library(dplyr)
library(purrr)
library(stringr)
library(tidyr)


# ==================================================================================================
# 1. DIRECTORIES
# ==================================================================================================

pearson_dir <- file.path(
  "results",
  "sim",
  "correlations",
  "Pearson"
)

summary_dir <- file.path(
  pearson_dir,
  "summaries"
)

dir.create(
  summary_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==================================================================================================
# 2. SETTINGS
# ==================================================================================================

M_keep <- 20

weight_order <- c(
  "no_weight",
  "CW",
  "PPW",
  "OPW"
)


# ==================================================================================================
# 3. TRUE CORRELATION LOOKUP
# ==================================================================================================

get_true_corr <- function(rho_xy, rho_uv) {
  
  dplyr::case_when(
    rho_xy == 0   & rho_uv == 0   ~ 0,
    rho_xy == 0   & rho_uv == 0.5 ~ 0.4,
    rho_xy == 0.5 & rho_uv == 0   ~ 0.1,
    rho_xy == 0.5 & rho_uv == 0.5 ~ 0.5,
    TRUE ~ NA_real_
  )
}


# ==================================================================================================
# 4. SUMMARIZE ONE SUBGROUPING
# ==================================================================================================

summarize_one_file <- function(file_path) {
  
  message("Reading: ", file_path)
  
  dat <- readRDS(file_path)
  
  
  # -----------------------------------------------------------------------------------------------
  # Check required columns
  # -----------------------------------------------------------------------------------------------
  
  needed_cols <- c(
    "weight_type",
    "estimate",
    "se",
    "M",
    "rho_xy",
    "rho_uv",
    "eta_x",
    "eta_y"
  )
  
  missing_cols <- setdiff(
    needed_cols,
    names(dat)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ",
      file_path
    )
  }
  
  
  # -----------------------------------------------------------------------------------------------
  # Keep M = 20 and the four current weighting methods
  # -----------------------------------------------------------------------------------------------
  
  dat <- dat %>%
    filter(
      M == M_keep,
      weight_type %in% weight_order
    )
  
  
  # -----------------------------------------------------------------------------------------------
  # Identify subgrouping
  #
  # Script 8 already stores dich_label, so use that rather than trying to reconstruct it
  # from the file name.
  # -----------------------------------------------------------------------------------------------
  
  if ("dich_label" %in% names(dat)) {
    
    subgrouping <- unique(dat$dich_label)
    
    if (length(subgrouping) != 1) {
      stop(
        "More than one dich_label found in:\n",
        file_path
      )
    }
    
  } else {
    
    subgrouping <- basename(dirname(file_path))
  }
  
  
  # -----------------------------------------------------------------------------------------------
  # Simulation scenario variables
  # -----------------------------------------------------------------------------------------------
  
  scenario_vars <- c(
    "M",
    "rho_xy",
    "rho_uv",
    "eta_x",
    "eta_y"
  )
  
  
  # -----------------------------------------------------------------------------------------------
  # Add true correlation and confidence intervals
  # -----------------------------------------------------------------------------------------------
  
  dat2 <- dat %>%
    mutate(
      
      rho_true = get_true_corr(
        rho_xy,
        rho_uv
      ),
      
      lower_95 = estimate - qnorm(0.975) * se,
      
      upper_95 = estimate + qnorm(0.975) * se
    )
  
  
  if (any(is.na(dat2$rho_true))) {
    
    bad_vals <- dat2 %>%
      filter(is.na(rho_true)) %>%
      distinct(
        rho_xy,
        rho_uv
      )
    
    print(bad_vals)
    
    stop(
      "Some rho_xy/rho_uv combinations do not have a defined true correlation."
    )
  }
  
  
  # -----------------------------------------------------------------------------------------------
  # Observed-data target
  #
  # For each simulation scenario, use the mean no-weight estimate as the observed-data
  # correlation target.
  # -----------------------------------------------------------------------------------------------
  
  obs_targets <- dat2 %>%
    filter(
      weight_type == "no_weight"
    ) %>%
    group_by(
      across(
        all_of(scenario_vars)
      )
    ) %>%
    summarise(
      
      rho_obs_hat = mean(
        estimate,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  dat2 <- dat2 %>%
    left_join(
      obs_targets,
      by = scenario_vars
    )
  
  
  if (any(is.na(dat2$rho_obs_hat))) {
    
    stop(
      "Some scenarios did not receive rho_obs_hat. ",
      "Check that no_weight exists for every scenario."
    )
  }
  
  
  # -----------------------------------------------------------------------------------------------
  # Determine coverage for each replicate
  # -----------------------------------------------------------------------------------------------
  
  dat2 <- dat2 %>%
    mutate(
      
      covered_true =
        lower_95 <= rho_true &
        rho_true <= upper_95,
      
      covered_obs =
        lower_95 <= rho_obs_hat &
        rho_obs_hat <= upper_95
    )
  
  
  # -----------------------------------------------------------------------------------------------
  # Summarize simulation replicates
  # -----------------------------------------------------------------------------------------------
  
  out <- dat2 %>%
    group_by(
      across(
        all_of(scenario_vars)
      ),
      weight_type
    ) %>%
    summarise(
      
      rho_true = first(rho_true),
      
      rho_obs_hat = first(rho_obs_hat),
      
      mean_estimate = mean(
        estimate,
        na.rm = TRUE
      ),
      
      mean_se = mean(
        se,
        na.rm = TRUE
      ),
      
      coverage_true = mean(
        covered_true,
        na.rm = TRUE
      ),
      
      coverage_obs = mean(
        covered_obs,
        na.rm = TRUE
      ),
      
      n_reps = n(),
      
      .groups = "drop"
    ) %>%
    mutate(
      
      correlation_measure = "Pearson",
      
      subgrouping = subgrouping,
      
      weight_type = factor(
        weight_type,
        levels = weight_order
      )
    ) %>%
    arrange(
      rho_xy,
      rho_uv,
      eta_x,
      eta_y,
      weight_type
    ) %>%
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
      weight_type,
      mean_estimate,
      mean_se,
      coverage_true,
      coverage_obs,
      n_reps
    )
  
  
  # -----------------------------------------------------------------------------------------------
  # Save subgroup summary
  # -----------------------------------------------------------------------------------------------
  
  safe_name <- str_replace_all(
    subgrouping,
    "[^A-Za-z0-9_-]",
    "_"
  )
  
  saveRDS(
    out,
    file.path(
      summary_dir,
      paste0(
        safe_name,
        "_summary_M20.rds"
      )
    )
  )
  
  write.csv(
    out,
    file.path(
      summary_dir,
      paste0(
        safe_name,
        "_summary_M20.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  message(
    "Finished: ",
    subgrouping,
    " | rows = ",
    nrow(out)
  )
  
  invisible(out)
}


# ==================================================================================================
# 5. FIND ALL COMBINED PEARSON RESULT FILES
# ==================================================================================================

result_files <- list.files(
  pearson_dir,
  pattern = "^all_pearson_.*_results\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)


message(
  "Number of Pearson subgrouping files found: ",
  length(result_files)
)


if (length(result_files) == 0) {
  
  stop(
    "No combined Pearson result files were found under:\n",
    pearson_dir
  )
}


# ==================================================================================================
# 6. SUMMARIZE ALL SUBGROUPINGS
# ==================================================================================================

compiled_pearson_m20 <- map_dfr(
  result_files,
  summarize_one_file
)


# ==================================================================================================
# 7. ADD HEATMAP METRICS
# ==================================================================================================

compiled_pearson_m20 <- compiled_pearson_m20 %>%
  mutate(
    
    bias_true =
      mean_estimate - rho_true,
    
    abs_bias_true =
      abs(bias_true),
    
    bias_obs =
      mean_estimate - rho_obs_hat,
    
    abs_bias_obs =
      abs(bias_obs),
    
    coverage_true_error =
      abs(coverage_true - 0.95),
    
    coverage_obs_error =
      abs(coverage_obs - 0.95)
  )


# ==================================================================================================
# 8. SAVE COMPILED SUMMARY
# ==================================================================================================

saveRDS(
  compiled_pearson_m20,
  file.path(
    summary_dir,
    "compiled_pearson_M20_summary.rds"
  )
)

write.csv(
  compiled_pearson_m20,
  file.path(
    summary_dir,
    "compiled_pearson_M20_summary.csv"
  ),
  row.names = FALSE
)


# ==================================================================================================
# 9. DISPLAY
# ==================================================================================================

print(
  compiled_pearson_m20,
  n = 100
)