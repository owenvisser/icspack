# Simulating Stouffer ISS tests over simulation grid
HiperGator <- TRUE

if (HiperGator) {
  custom_lib <- "/orange/somnath.datta/NHANES/Rpackages"
  
  .libPaths(c(custom_lib, .libPaths()))
  
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"))
  Sys.setenv(R_LIBS_USER = custom_lib)
  
  library(crspack)
  library(dplyr)
  library(parallel)
  
  save_dir <- "/orange/somnath.datta/NHANES/Simulation_RDS"
  grid <- readRDS(file.path(save_dir, "simulation_grid_index.rds"))
  
} else {
  
  library(crspack)
  library(dplyr)
  library(parallel)
  
  save_dir <- "C:/Users/owvis/OneDrive - University of Florida/nHANES/Simulation_RDS"
  grid <- readRDS(file.path(save_dir, "simulation_grid_index.rds"))
}


sub_ASD_pvals <- function(mdat,
                          m_sub = 20,
                          cores = 1,
                          seed = 1,
                          npctiles = 10,
                          K = 1000,
                          useZ = "percentiles",
                          Zvals = NULL,
                          HiperGator = FALSE) {
  
  clusters <- unique(mdat[, "cID"])
  
  cluster_groups <- split(
    clusters,
    ceiling(seq_along(clusters) / m_sub)
  )
  
  dat_groups <- lapply(cluster_groups, function(x) {
    mdat[mdat[, "cID"] %in% x, , drop = FALSE]
  })
  
  B <- length(dat_groups)
  
  set.seed(seed)
  seeds <- sample.int(1e8, B)
  
  task_list <- Map(
    function(dat, seed) {
      list(dat = dat, seed = seed)
    },
    dat_groups,
    seeds
  )
  
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  if (HiperGator) {
    parallel::clusterEvalQ(cl, {
      .libPaths(c("/orange/somnath.datta/NHANES/Rpackages", .libPaths()))
      library(crspack)
    })
  } else {
    parallel::clusterEvalQ(cl, {
      library(crspack)
    })
  }
  
  parallel::clusterExport(
    cl,
    c("useZ", "npctiles", "K", "Zvals"),
    envir = environment()
  )
  
  out <- parallel::parLapplyLB(cl, task_list, function(task) {
    
    set.seed(task$seed)
    
    x <- crspack::ASDpv(
      mdat = task$dat,
      useZ = useZ,
      Zvals = Zvals,
      npctiles = npctiles,
      K = K,
      parallel = FALSE
    )
    
    p <- pmin(pmax(x$pvalue, 1e-10), 1 - 1e-10)
    
    p
  })
  
  unlist(out)
}


stouffer_from_pvals <- function(pvals) {
  
  pvals <- pmin(pmax(pvals, 1e-10), 1 - 1e-10)
  
  zi <- qnorm(1 - pvals)
  
  z_stouffer <- sum(zi) / sqrt(length(zi))
  p_stouffer <- 1 - pnorm(z_stouffer)
  
  list(
    pvals = pvals,
    zi = zi,
    z_stouffer = z_stouffer,
    p_stouffer = p_stouffer,
    B = length(pvals)
  )
}


read_one_chunk <- function(out_file) {
  
  sim <- readRDS(out_file)
  
  mdat <- do.call(rbind, sim$data)
  
  mdat
}


prepare_chunk_obs <- function(mdat) {
  
  mdat_obs <- mdat[mdat$R == 1, ]
  mdat_obs <- mdat_obs %>% dplyr::filter(size > 1)
  
  mdat_obs$cID <- as.integer(
    factor(paste0(mdat_obs$replicate, "_", mdat_obs$cluster))
  )
  
  mdat_obs
}


make_test_dat <- function(mdat_obs,
                          direction = c("Z_equals_X", "Z_equals_Y")) {
  
  direction <- match.arg(direction)
  
  if (direction == "Z_equals_X") {
    
    test_dat <- cbind(
      cID = mdat_obs$cID,
      Y   = mdat_obs$Y,
      Z   = mdat_obs$X
    )
    
  } else if (direction == "Z_equals_Y") {
    
    test_dat <- cbind(
      cID = mdat_obs$cID,
      Y   = mdat_obs$X,
      Z   = mdat_obs$Y
    )
  }
  
  test_dat
}


run_one_chunk_direction <- function(out_file,
                                    direction = c("Z_equals_X", "Z_equals_Y"),
                                    cores = 100,
                                    HiperGator = TRUE) {
  
  direction <- match.arg(direction)
  
  mdat <- read_one_chunk(out_file)
  mdat_obs <- prepare_chunk_obs(mdat)
  
  test_dat <- make_test_dat(
    mdat_obs = mdat_obs,
    direction = direction
  )
  
  pvals <- sub_ASD_pvals(
    mdat = test_dat,
    m_sub = 10,
    K = 100,
    useZ = "percentiles",
    Zvals = NULL,
    npctiles = 5,
    cores = cores,
    HiperGator = HiperGator
  )
  
  n_clusters <- length(unique(test_dat[, "cID"]))
  
  rm(mdat, mdat_obs, test_dat)
  gc()
  
  list(
    pvals = pvals,
    n_clusters = n_clusters,
    B = length(pvals)
  )
}


## ------------------------------------------------------------
## Run only M = 20 scenarios
## ------------------------------------------------------------

grid_m20 <- grid %>%
  dplyr::filter(M == 20)

scenario_info <- grid_m20 %>%
  dplyr::group_by(scenario, eta_x, eta_y, rho_xy, rho_uv, M) %>%
  dplyr::summarise(
    out_files = list(out_file),
    n_chunks = dplyr::n(),
    .groups = "drop"
  )

cores <- 100

if (is.na(cores) || cores < 1) {
  cores <- 1
}

stouffer_results <- data.frame(
  scenario = rep(scenario_info$scenario, each = 2),
  eta_x    = rep(scenario_info$eta_x, each = 2),
  eta_y    = rep(scenario_info$eta_y, each = 2),
  rho_xy   = rep(scenario_info$rho_xy, each = 2),
  rho_uv   = rep(scenario_info$rho_uv, each = 2),
  M        = rep(scenario_info$M, each = 2),
  direction = rep(c("Z_equals_X", "Z_equals_Y"), times = nrow(scenario_info)),
  n_chunks = rep(scenario_info$n_chunks, each = 2),
  n_clusters = NA_integer_,
  B = NA_integer_,
  m_sub = NA_integer_,
  z_stouffer = NA_real_,
  p_stouffer = NA_real_
)

row_id <- 1

for (s in seq_len(nrow(scenario_info))) {
  
  cat("Running scenario", s, "of", nrow(scenario_info), "\n")
  cat("Original scenario ID:", scenario_info$scenario[s], "\n")
  cat("eta_x:", scenario_info$eta_x[s], "\n")
  cat("eta_y:", scenario_info$eta_y[s], "\n")
  cat("rho_xy:", scenario_info$rho_xy[s], "\n")
  cat("rho_uv:", scenario_info$rho_uv[s], "\n")
  cat("M:", scenario_info$M[s], "\n")
  cat("Chunks:", scenario_info$n_chunks[s], "\n")
  cat("Cores:", cores, "\n\n")
  
  out_files_now <- scenario_info$out_files[[s]]
  
  for (direction_now in c("Z_equals_X", "Z_equals_Y")) {
    
    cat("  Direction:", direction_now, "\n")
    
    all_pvals <- numeric(0)
    total_clusters <- 0L
    
    for (chunk_id in seq_along(out_files_now)) {
      
      cat("    Chunk", chunk_id, "of", length(out_files_now), "\n")
      
      chunk_res <- run_one_chunk_direction(
        out_file = out_files_now[chunk_id],
        direction = direction_now,
        cores = cores,
        HiperGator = HiperGator
      )
      
      all_pvals <- c(all_pvals, chunk_res$pvals)
      total_clusters <- total_clusters + chunk_res$n_clusters
      
      rm(chunk_res)
      gc()
      
      saveRDS(
        all_pvals,
        file.path(
          save_dir,
          paste0(
            "partial_pvals_scenario_",
            scenario_info$scenario[s],
            "_",
            direction_now,
            ".rds"
          )
        )
      )
    }
    
    res <- stouffer_from_pvals(all_pvals)
    
    stouffer_results$n_clusters[row_id] <- total_clusters
    stouffer_results$B[row_id] <- res$B
    stouffer_results$m_sub[row_id] <- 20
    stouffer_results$z_stouffer[row_id] <- res$z_stouffer
    stouffer_results$p_stouffer[row_id] <- res$p_stouffer
    
    cat("  Total ASD tests =", res$B, "\n")
    cat("  z_stouffer =", res$z_stouffer, "\n")
    cat("  p_stouffer =", res$p_stouffer, "\n\n")
    
    saveRDS(
      stouffer_results,
      file.path(save_dir, "simulation_stouffer_results_M20_both_directions_partial.rds")
    )
    
    write.csv(
      stouffer_results,
      file.path(save_dir, "simulation_stouffer_results_M20_both_directions_partial.csv"),
      row.names = FALSE
    )
    
    row_id <- row_id + 1
    
    rm(all_pvals, res)
    gc()
  }
}


saveRDS(
  stouffer_results,
  file.path(save_dir, "simulation_stouffer_results_M20_both_directions.rds")
)

write.csv(
  stouffer_results,
  file.path(save_dir, "simulation_stouffer_results_M20_both_directions.csv"),
  row.names = FALSE
)

stouffer_results