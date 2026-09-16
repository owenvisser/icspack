library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

# Main directory
base_dir <- r"(C:\Users\owvis\OneDrive - University of Florida\nHANES\Sim_Res)"

# Correlation folders
corr_dirs <- c("Phi", "Pearson", "Spearman")

# Weight order
weight_order <- c("no_weight", "CW", "PPW", "OPW", "MOPW")

# Output directory
summary_dir <- file.path(base_dir, "Scenario_Summaries_True_and_Observed")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

# True full-data correlation lookup
get_true_corr <- function(rho_xy, rho_uv) {
  dplyr::case_when(
    rho_xy == 0   & rho_uv == 0   ~ 0,
    rho_xy == 0   & rho_uv == 0.5 ~ 0.4,
    rho_xy == 0.5 & rho_uv == 0   ~ 0.1,
    rho_xy == 0.5 & rho_uv == 0.5 ~ 0.5,
    TRUE ~ NA_real_
  )
}

summarize_one_file <- function(file_path, corr_measure) {
  
  file_name <- basename(file_path)
  
  message("Reading file: ", file_name)
  
  dat <- readRDS(file_path)
  
  needed_cols <- c(
    "weight_type",
    "estimate",
    "se",
    "rho_xy",
    "rho_uv"
  )
  
  missing_cols <- setdiff(needed_cols, names(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "File is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFile: ", file_path
    )
  }
  
  # Identify the subgrouping/scenario name from the file name.
  # This is the L/K subgrouping label, not the simulation scenario.
  subgrouping <- file_name %>%
    str_remove("^all_") %>%
    str_remove("_results\\.rds$") %>%
    str_remove(paste0("^", tolower(corr_measure), "_"))
  
  # Use only scenario variables that actually exist in the file.
  candidate_scenario_vars <- c(
    "M",
    "rho_xy",
    "rho_uv",
    "eta_x",
    "eta_y"
  )
  
  scenario_vars <- intersect(candidate_scenario_vars, names(dat))
  
  if (!all(c("rho_xy", "rho_uv") %in% scenario_vars)) {
    stop(
      "rho_xy and rho_uv must be available as scenario variables.\n",
      "File: ", file_path
    )
  }
  
  message("Scenario variables used: ", paste(scenario_vars, collapse = ", "))
  
  # Add true rho and confidence intervals row-wise
  dat2 <- dat %>%
    mutate(
      rho_true = get_true_corr(rho_xy, rho_uv),
      lower_95 = estimate - qnorm(0.975) * se,
      upper_95 = estimate + qnorm(0.975) * se
    )
  
  if (any(is.na(dat2$rho_true))) {
    bad_vals <- dat2 %>%
      filter(is.na(rho_true)) %>%
      distinct(rho_xy, rho_uv)
    
    print(bad_vals)
    
    stop(
      "Some rho_xy/rho_uv combinations do not have a defined true correlation.\n",
      "File: ", file_path
    )
  }
  
  # Compute observed-data target within each actual simulation scenario:
  # average no_weight estimate for the same scenario.
  obs_targets <- dat2 %>%
    filter(weight_type == "no_weight") %>%
    group_by(across(all_of(scenario_vars))) %>%
    summarise(
      rho_obs_hat = mean(estimate, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join observed target back to all rows
  dat2 <- dat2 %>%
    left_join(obs_targets, by = scenario_vars)
  
  if (any(is.na(dat2$rho_obs_hat))) {
    stop(
      "Some rows did not receive rho_obs_hat. ",
      "Check whether no_weight exists for every scenario.\n",
      "File: ", file_path
    )
  }
  
  # Summarize by scenario and weight type
  out <- dat2 %>%
    mutate(
      covered_true = lower_95 <= rho_true    & rho_true    <= upper_95,
      covered_obs  = lower_95 <= rho_obs_hat & rho_obs_hat <= upper_95
    ) %>%
    group_by(across(all_of(scenario_vars)), weight_type) %>%
    summarise(
      rho_true      = first(rho_true),
      rho_obs_hat   = first(rho_obs_hat),
      mean_estimate = mean(estimate, na.rm = TRUE),
      mean_se       = mean(se, na.rm = TRUE),
      coverage_true = mean(covered_true, na.rm = TRUE),
      coverage_obs  = mean(covered_obs, na.rm = TRUE),
      n_reps        = n(),
      .groups = "drop"
    ) %>%
    mutate(
      correlation_measure = corr_measure,
      subgrouping = subgrouping,
      file_name = file_name,
      weight_type = factor(weight_type, levels = weight_order)
    ) %>%
    arrange(across(all_of(scenario_vars)), weight_type) %>%
    select(
      correlation_measure,
      subgrouping,
      file_name,
      all_of(scenario_vars),
      rho_true,
      rho_obs_hat,
      weight_type,
      mean_estimate,
      mean_se,
      coverage_true,
      coverage_obs,
      n_reps
    )
  
  save_name_base <- file_name %>%
    str_remove("\\.rds$")
  
  rds_out <- file.path(
    summary_dir,
    paste0(save_name_base, "_summary.rds")
  )
  
  csv_out <- file.path(
    summary_dir,
    paste0(save_name_base, "_summary.csv")
  )
  
  saveRDS(out, rds_out)
  write.csv(out, csv_out, row.names = FALSE)
  
  message("Saved summary: ", basename(rds_out))
  message("Number of summarized rows: ", nrow(out))
  message("Finished subgrouping: ", subgrouping)
  message("----------------------------------------")
  
  invisible(out)
}


for (corr_measure in corr_dirs) {
  
  this_dir <- file.path(base_dir, corr_measure)
  
  files <- list.files(
    this_dir,
    pattern = "_results\\.rds$",
    full.names = TRUE
  )
  
  message("========================================")
  message("Starting correlation measure: ", corr_measure)
  message("Number of files found: ", length(files))
  message("========================================")
  
  for (file_path in files) {
    
    tryCatch(
      {
        summarize_one_file(
          file_path = file_path,
          corr_measure = corr_measure
        )
      },
      error = function(e) {
        message("ERROR in file: ", basename(file_path))
        message("Error message: ", e$message)
        message("Skipping to next file.")
        message("----------------------------------------")
      }
    )
  }
  
  message("Finished correlation measure: ", corr_measure)
  message("========================================")
}

summary_files <- list.files(
  summary_dir,
  pattern = "_summary\\.rds$",
  full.names = TRUE
)

compiled_summary <- map_dfr(summary_files, readRDS) %>%
  mutate(
    weight_type = factor(weight_type, levels = weight_order)
  ) %>%
  arrange(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y,
    weight_type
  )

saveRDS(
  compiled_summary,
  file.path(base_dir, "compiled_correlation_summary_true_and_observed_long.rds")
)

write.csv(
  compiled_summary,
  file.path(base_dir, "compiled_correlation_summary_true_and_observed_long.csv"),
  row.names = FALSE
)

print(compiled_summary, n = 100)


compiled_summary_wide <- compiled_summary %>%
  mutate(
    est_cov_true_obs = sprintf(
      "%.3f (%.3f, %.3f)",
      mean_estimate,
      coverage_true,
      coverage_obs
    )
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
    est_cov_true_obs
  ) %>%
  pivot_wider(
    names_from = weight_type,
    values_from = est_cov_true_obs
  ) %>%
  arrange(
    correlation_measure,
    subgrouping,
    M,
    rho_xy,
    rho_uv,
    eta_x,
    eta_y
  )

saveRDS(
  compiled_summary_wide,
  file.path(base_dir, "compiled_correlation_summary_true_and_observed_wide.rds")
)

write.csv(
  compiled_summary_wide,
  file.path(base_dir, "compiled_correlation_summary_true_and_observed_wide.csv"),
  row.names = FALSE
)

print(compiled_summary_wide, n = 100)




