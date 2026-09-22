####################################################################################################
#
# SCRIPT 08: CALCULATE SIMULATION CORRELATION VALUES
#
# Calculates Pearson, Spearman, and Phi association estimates
# across the simulation grid and severity subgrouping choices.
#
# Simulation data:
#   data/sim/
#
# Correlation results:
#   results/sim/correlations/
#
####################################################################################################


# ==================================================================================================
# 1. LOAD PACKAGES AND HELPER FUNCTIONS
# ==================================================================================================

library(dplyr)

source(
  file.path(
    "R",
    "Weighting_Equations.R"
  )
)

source(
  file.path(
    "R",
    "Correlation_Equations.R"
  )
)


# ==================================================================================================
# 2. DEFINE FILE LOCATIONS
# ==================================================================================================

# Simulation data created in Script 06.
sim_dir <- file.path(
  "data",
  "sim"
)

# Correlation results created by this script.
results_dir <- file.path(
  "results",
  "sim",
  "correlations"
)

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Load the simulation grid.
grid <- readRDS(
  file.path(
    sim_dir,
    "simulation_grid_index.rds"
  )
)

cat(
  "Loaded simulation grid with",
  nrow(grid),
  "saved simulation chunks.\n"
)


# ==================================================================================================
# 3. HELPER: CREATE SAFE FILE/FOLDER NAMES
# ==================================================================================================

safe_name <- function(x) {

  x <- gsub(
    "\\.",
    "_",
    x
  )

  x <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    x
  )

  x
}


# ==================================================================================================
# 4. HELPER: CREATE DICHOTOMIZED SEVERITY VARIABLES
# ==================================================================================================

add_severity_dichotomies <- function(
  dat,
  cutpoints = c(
    1.5,
    2.5,
    3.5,
    4.5
  )
) {

  # Convert severity variables to numeric in case they
  # were stored as factors or characters.
  X_sev_num <- as.numeric(
    as.character(
      dat$X_sev
    )
  )

  Y_sev_num <- as.numeric(
    as.character(
      dat$Y_sev
    )
  )


  # Create one binary severity variable for each cutpoint.
  for (cp in cutpoints) {

    cp_lab <- safe_name(
      as.character(cp)
    )

    x_name <- paste0(
      "X_dich_",
      cp_lab
    )

    y_name <- paste0(
      "Y_dich_",
      cp_lab
    )

    dat[[x_name]] <- as.integer(
      X_sev_num > cp
    )

    dat[[y_name]] <- as.integer(
      Y_sev_num > cp
    )
  }


  return(dat)
}


# ==================================================================================================
# 5. HELPER: BUILD THE SUBGROUPING GRID
# ==================================================================================================

make_subgroup_grid <- function(
  cutpoints = c(
    1.5,
    2.5,
    3.5,
    4.5
  )
) {

  cp_labs <- safe_name(
    as.character(
      cutpoints
    )
  )


  # L subgrouping variables.
  L_vars <- c(
    "Y_sev",
    paste0(
      "Y_dich_",
      cp_labs
    )
  )


  # K subgrouping variables.
  K_vars <- c(
    "X_sev",
    paste0(
      "X_dich_",
      cp_labs
    )
  )


  # Store the corresponding cutpoints.
  #
  # NA corresponds to the original ordinal severity variable.
  L_cut <- c(
    NA,
    cutpoints
  )

  K_cut <- c(
    NA,
    cutpoints
  )


  L_info <- data.frame(
    L_var = L_vars,
    L_cut = L_cut,
    stringsAsFactors = FALSE
  )

  K_info <- data.frame(
    K_var = K_vars,
    K_cut = K_cut,
    stringsAsFactors = FALSE
  )


  # All combinations of L and K subgroup definitions.
  subgroup_grid <- merge(
    L_info,
    K_info,
    all = TRUE
  )


  # Create a label that can safely be used in filenames.
  subgroup_grid$dich_label <- paste0(
    "L_",
    subgroup_grid$L_var,
    "__K_",
    subgroup_grid$K_var
  )

  subgroup_grid$dich_label <- safe_name(
    subgroup_grid$dich_label
  )


  return(subgroup_grid)
}


# ==================================================================================================
# 6. HELPER: CHOOSE VARIABLES FOR EACH ASSOCIATION TYPE
# ==================================================================================================

get_association_variables <- function(
  association_type
) {

  if (association_type == "pearson") {

    # Pearson uses the original continuous outcomes.
    out <- list(
      Y_var = "Y",
      X_var = "X"
    )

  } else if (association_type == "spearman") {

    # Spearman uses the ordinal severity outcomes.
    out <- list(
      Y_var = "Y_sev",
      X_var = "X_sev"
    )

  } else if (association_type == "phi") {

    # Phi uses binary versions of the continuous outcomes.
    out <- list(
      Y_var = "Y_bin",
      X_var = "X_bin"
    )

  } else {

    stop(
      "association_type must be one of: pearson, spearman, phi"
    )
  }


  return(out)
}


# ==================================================================================================
# 7. PROCESS ONE SAVED SIMULATION CHUNK
# ==================================================================================================

process_one_sim_file <- function(
  row_id,
  association_type,
  L_var,
  K_var,
  L_cut,
  K_cut,
  dich_label,
  res_dir,
  grid,
  cutpoints = c(
    1.5,
    2.5,
    3.5,
    4.5
  )
) {

  # Get the simulation file directly from the saved grid.
  sim_file <- grid$out_file[row_id]


  cat("\n")
  cat(
    "Processing simulation chunk",
    row_id,
    "of",
    nrow(grid),
    "\n"
  )

  cat(
    "Reading:",
    sim_file,
    "\n"
  )

  cat(
    "Association type:",
    association_type,
    "\n"
  )

  cat(
    "Subgrouping:",
    dich_label,
    "\n"
  )


  # Load simulation chunk.
  sim <- readRDS(
    sim_file
  )


  # ----------------------------------------------------------------------------------------------
  # Prepare each replicate in the simulation chunk.
  # ----------------------------------------------------------------------------------------------

  mdat_obs_list <- lapply(
    sim$data,
    function(dat) {

      # Restrict to observed observations.
      dat <- dat[
        dat$R == 1,
        ,
        drop = FALSE
      ]


      # Remove clusters with only one observation.
      dat <- dat %>%
        filter(
          size > 1
        )


      # Create unique cluster ID.
      dat$cID <- as.integer(
        factor(
          paste0(
            dat$replicate,
            "_",
            dat$cluster
          )
        )
      )


      # Binary variables used for the Phi coefficient.
      #
      # 1 = positive outcome
      # 0 = non-positive outcome
      dat$X_bin <- as.integer(
        dat$X > 0
      )

      dat$Y_bin <- as.integer(
        dat$Y > 0
      )


      # Create severity dichotomies used for subgroup definitions.
      dat <- add_severity_dichotomies(
        dat = dat,
        cutpoints = cutpoints
      )


      return(dat)
    }
  )


  # ----------------------------------------------------------------------------------------------
  # Determine association variables.
  # ----------------------------------------------------------------------------------------------

  assoc_vars <- get_association_variables(
    association_type
  )

  Y_vars <- assoc_vars$Y_var
  X_vars <- assoc_vars$X_var


  # Subgrouping variables.
  L_vars <- L_var
  K_vars <- K_var


  # Weighting methods.
  weight_types <- c(
    "no_weight",
    "CW",
    "PPW",
    "OPW"#,
    #"MOPW"
  )


  # ----------------------------------------------------------------------------------------------
  # Calculate correlations for every replicate in this chunk.
  # ----------------------------------------------------------------------------------------------

  all_res <- do.call(
    rbind,
    lapply(
      seq_along(mdat_obs_list),
      function(rep_id) {

        mdat_obs <- mdat_obs_list[[rep_id]]


        # Run all weighting methods for this replicate.
        assoc_results <- do.call(
          rbind,
          lapply(
            weight_types,
            function(wtype) {

              run_association_grid(
                dat = mdat_obs,
                Y_vars = Y_vars,
                X_vars = X_vars,
                clusterID = "cID",
                clusterSize = "size",
                K_vars = K_vars,
                L_vars = L_vars,
                weight_type = wtype,
                association_type = association_type
              )
            }
          )
        )


        # ------------------------------------------------------------------------------------------
        # Add simulation information.
        # ------------------------------------------------------------------------------------------

        assoc_results$replicate <- (
          sim$rep_start +
          rep_id -
          1
        )

        assoc_results$chunk_id <- sim$chunk_id

        assoc_results$label <- sim$label

        assoc_results$M <- sim$M

        assoc_results$eta_x <- sim$params$eta_x

        assoc_results$eta_y <- sim$params$eta_y

        assoc_results$rho_xy <- sim$params$rho_xy

        assoc_results$rho_uv <- sim$params$rho_uv


        # ------------------------------------------------------------------------------------------
        # Add subgrouping information.
        # ------------------------------------------------------------------------------------------

        assoc_results$association_type <- association_type

        assoc_results$dich_label <- dich_label

        assoc_results$L_subgroup_var <- L_var

        assoc_results$K_subgroup_var <- K_var

        assoc_results$L_cut <- L_cut

        assoc_results$K_cut <- K_cut


        return(assoc_results)
      }
    )
  )


  # ----------------------------------------------------------------------------------------------
  # Save results for this simulation chunk.
  # ----------------------------------------------------------------------------------------------

  out_file <- file.path(
    res_dir,
    paste0(
      association_type,
      "_results_",
      dich_label,
      "_",
      sim$label,
      "_",
      sim$chunk_id,
      ".rds"
    )
  )


  saveRDS(
    all_res,
    out_file
  )


  cat(
    "Saved:",
    out_file,
    "\n"
  )


  # Remove large objects before moving to the next chunk.
  rm(
    sim,
    mdat_obs_list,
    all_res
  )

  gc()


  return(
    invisible(
      out_file
    )
  )
}


# ==================================================================================================
# 8. RUN ONE ASSOCIATION TYPE FOR ONE SUBGROUPING CHOICE
# ==================================================================================================

run_one_association_subgrouping <- function(
  association_type,
  subgroup_row,
  results_dir,
  grid,
  cutpoints = c(
    1.5,
    2.5,
    3.5,
    4.5
  )
) {

  L_var <- subgroup_row$L_var

  K_var <- subgroup_row$K_var

  L_cut <- subgroup_row$L_cut

  K_cut <- subgroup_row$K_cut

  dich_label <- subgroup_row$dich_label


  # Create association-specific folder name.
  association_folder <- tools::toTitleCase(
    association_type
  )


  # Results for this exact association/subgrouping combination.
  res_dir <- file.path(
    results_dir,
    association_folder,
    dich_label
  )


  dir.create(
    res_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )


  cat("\n")
  cat("================================================================================\n")
  cat(
    "Running association type:",
    association_type,
    "\n"
  )
  cat(
    "L subgroup variable:",
    L_var,
    "\n"
  )
  cat(
    "K subgroup variable:",
    K_var,
    "\n"
  )
  cat(
    "Saving to:",
    res_dir,
    "\n"
  )
  cat("================================================================================\n")


  # ----------------------------------------------------------------------------------------------
  # Run each saved simulation chunk sequentially.
  # ----------------------------------------------------------------------------------------------

  res_files <- lapply(
    seq_len(
      nrow(grid)
    ),
    function(row_id) {

      process_one_sim_file(
        row_id = row_id,
        association_type = association_type,
        L_var = L_var,
        K_var = K_var,
        L_cut = L_cut,
        K_cut = K_cut,
        dich_label = dich_label,
        res_dir = res_dir,
        grid = grid,
        cutpoints = cutpoints
      )
    }
  )


  res_files <- unlist(
    res_files,
    use.names = FALSE
  )


  # ----------------------------------------------------------------------------------------------
  # Save an index of the result files.
  # ----------------------------------------------------------------------------------------------

  grid_out <- grid

  grid_out$association_type <- association_type

  grid_out$dich_label <- dich_label

  grid_out$L_subgroup_var <- L_var

  grid_out$K_subgroup_var <- K_var

  grid_out$L_cut <- L_cut

  grid_out$K_cut <- K_cut

  grid_out$result_file <- res_files


  index_file <- file.path(
    res_dir,
    paste0(
      association_type,
      "_",
      dich_label,
      "_results_index.rds"
    )
  )


  saveRDS(
    grid_out,
    index_file
  )


  cat("\n")
  cat(
    "Saved index file:",
    index_file,
    "\n"
  )


  # ----------------------------------------------------------------------------------------------
  # Combine all simulation chunks.
  # ----------------------------------------------------------------------------------------------

  cat(
    "Combining result files for",
    association_type,
    "-",
    dich_label,
    "\n"
  )


  all_results <- do.call(
    rbind,
    lapply(
      res_files,
      readRDS
    )
  )


  all_results_file <- file.path(
    res_dir,
    paste0(
      "all_",
      association_type,
      "_",
      dich_label,
      "_results.rds"
    )
  )


  saveRDS(
    all_results,
    all_results_file
  )


  cat(
    "Saved combined results:",
    all_results_file,
    "\n"
  )


  rm(
    all_results
  )

  gc()


  return(
    invisible(
      list(
        association_type = association_type,
        dich_label = dich_label,
        L_subgroup_var = L_var,
        K_subgroup_var = K_var,
        L_cut = L_cut,
        K_cut = K_cut,
        res_dir = res_dir,
        result_files = res_files,
        index_file = index_file,
        all_results_file = all_results_file
      )
    )
  )
}


# ==================================================================================================
# 9. RUN ONE ASSOCIATION TYPE OVER ALL SUBGROUPING CHOICES
# ==================================================================================================

run_simulation_association <- function(
  association_type,
  results_dir,
  grid,
  subgroup_grid,
  cutpoints = c(
    1.5,
    2.5,
    3.5,
    4.5
  )
) {

  if (
    !association_type %in%
    c(
      "pearson",
      "spearman",
      "phi"
    )
  ) {

    stop(
      "association_type must be one of: pearson, spearman, phi"
    )
  }


  cat("\n")
  cat("################################################################################\n")
  cat(
    "Starting association type:",
    association_type,
    "\n"
  )
  cat(
    "Number of subgrouping choices:",
    nrow(subgroup_grid),
    "\n"
  )
  cat("################################################################################\n")


  # Run each subgrouping choice sequentially.
  subgroup_outputs <- lapply(
    seq_len(
      nrow(subgroup_grid)
    ),
    function(sg_id) {

      run_one_association_subgrouping(
        association_type = association_type,
        subgroup_row = subgroup_grid[sg_id, ],
        results_dir = results_dir,
        grid = grid,
        cutpoints = cutpoints
      )
    }
  )


  names(subgroup_outputs) <- subgroup_grid$dich_label


  # ----------------------------------------------------------------------------------------------
  # Save association-level index.
  # ----------------------------------------------------------------------------------------------

  association_folder <- tools::toTitleCase(
    association_type
  )


  association_dir <- file.path(
    results_dir,
    association_folder
  )


  dir.create(
    association_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )


  association_index_file <- file.path(
    association_dir,
    paste0(
      association_type,
      "_all_subgroupings_index.rds"
    )
  )


  saveRDS(
    subgroup_outputs,
    association_index_file
  )


  cat("\n")
  cat(
    "Saved association-level subgrouping index:",
    association_index_file,
    "\n"
  )


  return(
    invisible(
      list(
        association_type = association_type,
        association_dir = association_dir,
        subgroup_outputs = subgroup_outputs,
        association_index_file = association_index_file
      )
    )
  )
}


# ==================================================================================================
# 10. DEFINE SEVERITY CUTPOINTS AND SUBGROUPING GRID
# ==================================================================================================

cutpoints <- c(
  1.5,
  2.5,
  3.5,
  4.5
)


subgroup_grid <- make_subgroup_grid(
  cutpoints = cutpoints
)


cat(
  "Number of subgrouping combinations:",
  nrow(subgroup_grid),
  "\n"
)


# ==================================================================================================
# 11. RUN PEARSON, SPEARMAN, AND PHI
# ==================================================================================================

# association_types <- c(
#   "pearson",
#   "spearman",
#   "phi"
# )

association_types <- c(
  "pearson"
)

association_outputs <- lapply(
  association_types,
  function(a_type) {

    run_simulation_association(
      association_type = a_type,
      results_dir = results_dir,
      grid = grid,
      subgroup_grid = subgroup_grid,
      cutpoints = cutpoints
    )
  }
)

names(
  association_outputs
) <- association_types


# ==================================================================================================
# 12. SAVE MASTER CORRELATION RESULTS INDEX
# ==================================================================================================

master_index_file <- file.path(
  results_dir,
  "correlation_results_master_index.rds"
)


saveRDS(
  association_outputs,
  master_index_file
)

