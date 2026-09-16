# Calculating the simulation correlation values with severity cutpoint subgrouping

HiperGator <- TRUE

if (HiperGator) {
  
  save_dir <- "/orange/somnath.datta/NHANES/Simulation_RDS"
  
  grid <- readRDS(file.path(save_dir, "simulation_grid_index.rds"))
  
  library(dplyr)
  
  source("/orange/somnath.datta/NHANES/Weighting_Equations.R")
  source("/orange/somnath.datta/NHANES/Correlation_Equations.R")
  
} else {
  
  save_dir <- "C:/Users/owvis/OneDrive - University of Florida/nHANES/Simulation_RDS"
  
  grid <- readRDS(file.path(save_dir, "simulation_grid_index.rds"))
  
  library(dplyr)
  
  source("C:/Users/owvis/OneDrive - University of Florida/nHANES/Code/Weighting_Equations.R")
  source("C:/Users/owvis/OneDrive - University of Florida/nHANES/Code/Correlation_Equations.R")
}

## ------------------------------------------------------------
## Helper: make safe folder names
## ------------------------------------------------------------

safe_name <- function(x) {
  x <- gsub("\\.", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x
}

## ------------------------------------------------------------
## Helper: create dichotomized severity columns
## ------------------------------------------------------------

add_severity_dichotomies <- function(dat,
                                     cutpoints = c(1.5, 2.5, 3.5, 4.5)) {
  
  ## In case X_sev/Y_sev are factors or character-like
  X_sev_num <- as.numeric(as.character(dat$X_sev))
  Y_sev_num <- as.numeric(as.character(dat$Y_sev))
  
  for (cp in cutpoints) {
    
    cp_lab <- safe_name(as.character(cp))
    
    x_name <- paste0("X_dich_", cp_lab)
    y_name <- paste0("Y_dich_", cp_lab)
    
    dat[[x_name]] <- as.integer(X_sev_num > cp)
    dat[[y_name]] <- as.integer(Y_sev_num > cp)
  }
  
  dat
}

## ------------------------------------------------------------
## Helper: build weighting subgroup grid
## ------------------------------------------------------------

make_subgroup_grid <- function(cutpoints = c(1.5, 2.5, 3.5, 4.5)) {
  
  cp_labs <- safe_name(as.character(cutpoints))
  
  L_vars <- c(
    "Y_sev",
    paste0("Y_dich_", cp_labs)
  )
  
  K_vars <- c(
    "X_sev",
    paste0("X_dich_", cp_labs)
  )
  
  L_cut <- c(NA, cutpoints)
  K_cut <- c(NA, cutpoints)
  
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
  
  subgroup_grid <- merge(L_info, K_info, all = TRUE)
  
  subgroup_grid$dich_label <- paste0(
    "L_", subgroup_grid$L_var,
    "__K_", subgroup_grid$K_var
  )
  
  subgroup_grid$dich_label <- safe_name(subgroup_grid$dich_label)
  
  subgroup_grid
}

## ------------------------------------------------------------
## Helper: choose association variables
## ------------------------------------------------------------

get_association_variables <- function(association_type) {
  
  if (association_type == "pearson") {
    
    ## Pearson uses the original continuous outcomes
    out <- list(
      Y_var = "Y",
      X_var = "X"
    )
    
  } else if (association_type == "spearman") {
    
    ## Spearman uses the ordinal severity outcomes
    out <- list(
      Y_var = "Y_sev",
      X_var = "X_sev"
    )
    
  } else if (association_type == "phi") {
    
    ## Phi uses binary versions of the original continuous outcomes
    out <- list(
      Y_var = "Y_bin",
      X_var = "X_bin"
    )
    
  } else {
    
    stop("association_type must be one of: pearson, spearman, phi")
  }
  
  out
}

## ------------------------------------------------------------
## Function to process one saved simulation chunk
## for one association type and one subgrouping choice
## ------------------------------------------------------------

process_one_sim_file <- function(row_id,
                                 association_type = "pearson",
                                 L_var,
                                 K_var,
                                 L_cut,
                                 K_cut,
                                 dich_label,
                                 res_dir,
                                 grid,
                                 cutpoints = c(1.5, 2.5, 3.5, 4.5)) {
  
  sim_file <- grid$out_file[row_id]
  
  cat("Processing row", row_id, "of", nrow(grid), "\n")
  cat("Reading:", sim_file, "\n")
  cat("Association type:", association_type, "\n")
  cat("Subgrouping:", dich_label, "\n")
  
  sim <- readRDS(sim_file)
  
  mdat_obs_list <- lapply(sim$data, function(dat) {
    
    dat <- dat[dat$R == 1, ]
    dat <- dat %>% dplyr::filter(size > 1)
    
    dat$cID <- as.integer(
      factor(paste0(dat$replicate, "_", dat$cluster))
    )
    
    ## Binary variables for Phi coefficient
    ## 1 = positive outcome, 0 = non-positive outcome
    dat$X_bin <- as.integer(dat$X > 0)
    dat$Y_bin <- as.integer(dat$Y > 0)
    
    ## Dichotomized severity subgrouping variables
    dat <- add_severity_dichotomies(
      dat = dat,
      cutpoints = cutpoints
    )
    
    dat
  })
  
  assoc_vars <- get_association_variables(association_type)
  
  Y_vars <- assoc_vars$Y_var
  X_vars <- assoc_vars$X_var
  
  L_vars <- L_var
  K_vars <- K_var
  
  weight_types <- c("no_weight", "CW", "PPW", "OPW", "MOPW")
  
  all_res <- do.call(
    rbind,
    lapply(seq_along(mdat_obs_list), function(rep_id) {
      
      mdat_obs <- mdat_obs_list[[rep_id]]
      
      assoc_results <- do.call(
        rbind,
        lapply(weight_types, function(wtype) {
          
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
        })
      )
      
      assoc_results$replicate <- sim$rep_start + rep_id - 1
      assoc_results$chunk_id  <- sim$chunk_id
      assoc_results$label     <- sim$label
      assoc_results$M         <- sim$M
      assoc_results$eta_x     <- sim$params$eta_x
      assoc_results$eta_y     <- sim$params$eta_y
      assoc_results$rho_xy    <- sim$params$rho_xy
      assoc_results$rho_uv    <- sim$params$rho_uv
      
      ## Explicit subgrouping labels
      assoc_results$association_type <- association_type
      assoc_results$dich_label       <- dich_label
      assoc_results$L_subgroup_var   <- L_var
      assoc_results$K_subgroup_var   <- K_var
      assoc_results$L_cut            <- L_cut
      assoc_results$K_cut            <- K_cut
      
      assoc_results
    })
  )
  
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
  
  saveRDS(all_res, out_file)
  
  cat("Saved:", out_file, "\n\n")
  
  rm(sim, mdat_obs_list, all_res)
  gc()
  
  invisible(out_file)
}

## ------------------------------------------------------------
## Function to run one association type and one subgrouping choice
## over the full simulation grid
## ------------------------------------------------------------

run_one_association_subgrouping <- function(association_type,
                                            subgroup_row,
                                            save_dir,
                                            grid,
                                            cores = 50,
                                            cutpoints = c(1.5, 2.5, 3.5, 4.5)) {
  
  L_var      <- subgroup_row$L_var
  K_var      <- subgroup_row$K_var
  L_cut      <- subgroup_row$L_cut
  K_cut      <- subgroup_row$K_cut
  dich_label <- subgroup_row$dich_label
  
  association_folder <- tools::toTitleCase(association_type)
  
  res_dir <- file.path(
    save_dir,
    "Association_Results",
    association_folder,
    dich_label
  )
  
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("\n")
  cat("============================================================\n")
  cat("Running association type:", association_type, "\n")
  cat("L subgroup variable:", L_var, "\n")
  cat("K subgroup variable:", K_var, "\n")
  cat("Saving to:", res_dir, "\n")
  cat("============================================================\n")
  
  res_files <- parallel::mclapply(
    seq_len(nrow(grid)),
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
    },
    mc.cores = cores
  )
  
  res_files <- unlist(res_files, use.names = FALSE)
  
  ## Save an index for this association type and subgrouping
  grid_out <- grid
  
  grid_out$association_type <- association_type
  grid_out$dich_label       <- dich_label
  grid_out$L_subgroup_var   <- L_var
  grid_out$K_subgroup_var   <- K_var
  grid_out$L_cut            <- L_cut
  grid_out$K_cut            <- K_cut
  grid_out$result_file      <- res_files
  
  index_file <- file.path(
    res_dir,
    paste0(
      association_type,
      "_",
      dich_label,
      "_results_index.rds"
    )
  )
  
  saveRDS(grid_out, index_file)
  
  cat("\n")
  cat("Saved index file:\n")
  cat(index_file, "\n")
  
  ## ----------------------------------------------------------
  ## Combine all result files for this association/subgrouping
  ## ----------------------------------------------------------
  
  cat("\n")
  cat("Combining all result files for:\n")
  cat("Association:", association_type, "\n")
  cat("Subgrouping:", dich_label, "\n")
  
  all_results <- do.call(
    rbind,
    lapply(res_files, readRDS)
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
  
  saveRDS(all_results, all_results_file)
  
  cat("Saved combined results:\n")
  cat(all_results_file, "\n")
  
  rm(all_results)
  gc()
  
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
}

## ------------------------------------------------------------
## Function to run one association type over all subgroupings
## ------------------------------------------------------------

run_simulation_association <- function(association_type,
                                       save_dir,
                                       grid,
                                       subgroup_grid,
                                       cores = 50,
                                       cutpoints = c(1.5, 2.5, 3.5, 4.5)) {
  
  if (!association_type %in% c("pearson", "spearman", "phi")) {
    stop("association_type must be one of: pearson, spearman, phi")
  }
  
  cat("\n")
  cat("############################################################\n")
  cat("Starting association type:", association_type, "\n")
  cat("Number of subgrouping choices:", nrow(subgroup_grid), "\n")
  cat("############################################################\n")
  
  subgroup_outputs <- lapply(
    seq_len(nrow(subgroup_grid)),
    function(sg_id) {
      
      run_one_association_subgrouping(
        association_type = association_type,
        subgroup_row = subgroup_grid[sg_id, ],
        save_dir = save_dir,
        grid = grid,
        cores = cores,
        cutpoints = cutpoints
      )
    }
  )
  
  names(subgroup_outputs) <- subgroup_grid$dich_label
  
  ## Save association-level index
  association_folder <- tools::toTitleCase(association_type)
  
  association_dir <- file.path(
    save_dir,
    "Association_Results",
    association_folder
  )
  
  association_index_file <- file.path(
    association_dir,
    paste0(association_type, "_all_subgroupings_index.rds")
  )
  
  saveRDS(subgroup_outputs, association_index_file)
  
  cat("\n")
  cat("Saved association-level subgrouping index:\n")
  cat(association_index_file, "\n")
  
  invisible(
    list(
      association_type = association_type,
      association_dir = association_dir,
      subgroup_outputs = subgroup_outputs,
      association_index_file = association_index_file
    )
  )
}

## ------------------------------------------------------------
## Run Pearson, Spearman, and Phi over all subgrouping choices
## ------------------------------------------------------------

cores <- 100

cutpoints <- c(1.5, 2.5, 3.5, 4.5)

subgroup_grid <- make_subgroup_grid(
  cutpoints = cutpoints
)

association_types <- c("pearson", "spearman", "phi")

association_outputs <- lapply(
  association_types,
  function(a_type) {
    run_simulation_association(
      association_type = a_type,
      save_dir = save_dir,
      grid = grid,
      subgroup_grid = subgroup_grid,
      cores = cores,
      cutpoints = cutpoints
    )
  }
)

names(association_outputs) <- association_types

## ------------------------------------------------------------
## Save a master index of all association result objects
## ------------------------------------------------------------

master_index_file <- file.path(
  save_dir,
  "Association_Results",
  "correlation_results_master_index.rds"
)

dir.create(
  dirname(master_index_file),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(association_outputs, master_index_file)

cat("\n")
cat("============================================================\n")
cat("Finished all correlation calculations.\n")
cat("Master index saved to:\n")
cat(master_index_file, "\n")
cat("============================================================\n")