####################################################################################################
#
# 07_simulation_pvalues.R
#
# PURPOSE:
#   Run the ASD test on the simulated data created in Script 6.
#
#   For each simulation scenario:
#     1. Read each saved simulation chunk.
#     2. Keep observed observations (R == 1).
#     3. Remove clusters with size <= 1.
#     4. Create a unique cluster ID across simulation replicates.
#     5. Run ASDpv() on small groups of clusters.
#     6. Repeat in both directions:
#           Z = X, outcome = Y
#           Z = Y, outcome = X
#     7. Combine the resulting p-values using Stouffer's method.
#
# INPUT:
#   data/sim/simulation_grid_index.rds
#   data/sim/sim_reps_*.rds
#
# OUTPUT:
#   results/sim/simulation_stouffer_results.rds
#   results/sim/simulation_stouffer_results.csv
#
####################################################################################################


# ==================================================================================================
# 1. LOAD PACKAGES
# ==================================================================================================

library(crspack)
library(dplyr)
library(parallel)


# ==================================================================================================
# 2. FILE LOCATIONS
# ==================================================================================================

sim_dir <- file.path(
  "data",
  "sim"
)

results_dir <- file.path(
  "results",
  "sim"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Read the simulation index created in Script 6.
grid <- readRDS(
  file.path(
    sim_dir,
    "simulation_grid_index.rds"
  )
)


# ==================================================================================================
# 3. ASD SETTINGS
# ==================================================================================================

# Number of clusters placed into each ASD test.
m_sub <- 10

# Number of randomizations/permutations used by ASDpv().
K <- 100

# Z construction used by ASDpv().
useZ <- "percentiles"

# Number of percentile groups.
npctiles <- 5

# No prespecified Z values when useZ = "percentiles".
Zvals <- NULL

# Number of parallel R workers.
#
# We are deliberately keeping this modest.
# Script 6 showed that spawning a very large number of PSOCK workers
# locally can overwhelm the machine.
cores <- 4

if (is.na(cores) || cores < 1) {
  cores <- 1
}

# Base seed used to generate reproducible seeds for individual ASD tests.
base_seed <- 20260922


# ==================================================================================================
# 4. FUNCTION: RUN ASD TESTS ON GROUPS OF CLUSTERS
# ==================================================================================================

sub_ASD_pvals <- function(
    mdat,
    cl,
    m_sub = 10,
    seed = 1,
    npctiles = 5,
    K = 100,
    useZ = "percentiles",
    Zvals = NULL
) {

  # ------------------------------------------------------------------
  # Find all unique clusters.
  # ------------------------------------------------------------------

  clusters <- unique(
    mdat[, "cID"]
  )


  # ------------------------------------------------------------------
  # Break the clusters into groups of size m_sub.
  #
  # Example:
  #
  # clusters = 1,2,3,...,25
  # m_sub = 10
  #
  # group 1 = clusters  1-10
  # group 2 = clusters 11-20
  # group 3 = clusters 21-25
  # ------------------------------------------------------------------

  cluster_groups <- split(
    clusters,
    ceiling(
      seq_along(clusters) / m_sub
    )
  )


  # ------------------------------------------------------------------
  # Create one data set for each group of clusters.
  # ------------------------------------------------------------------

  dat_groups <- lapply(
    cluster_groups,
    function(cluster_ids) {

      mdat[
        mdat[, "cID"] %in% cluster_ids,
        ,
        drop = FALSE
      ]

    }
  )


  B <- length(dat_groups)


  if (B == 0) {
    return(numeric(0))
  }


  # ------------------------------------------------------------------
  # Give each ASD test its own reproducible random seed.
  # ------------------------------------------------------------------

  set.seed(seed)

  seeds <- sample.int(
    1e8,
    B
  )


  # ------------------------------------------------------------------
  # Pair each data set with its seed.
  # ------------------------------------------------------------------

  task_list <- Map(
    function(dat, task_seed) {

      list(
        dat = dat,
        seed = task_seed
      )

    },
    dat_groups,
    seeds
  )


  # ------------------------------------------------------------------
  # Run the ASD tests in parallel.
  #
  # The cluster "cl" was already created outside this function.
  # Therefore we do NOT repeatedly start and stop PSOCK workers.
  # ------------------------------------------------------------------

  out <- parallel::parLapplyLB(
    cl,
    task_list,
    function(
        task,
        useZ,
        Zvals,
        npctiles,
        K
    ) {

      set.seed(
        task$seed
      )


      fit <- crspack::ASDpv(
        mdat = task$dat,
        useZ = useZ,
        Zvals = Zvals,
        npctiles = npctiles,
        K = K,
        parallel = FALSE
      )


      p <- fit$pvalue


      # Make sure ASDpv returned one usable p-value.
      if (
        length(p) != 1 ||
        !is.finite(p)
      ) {

        stop(
          "ASDpv returned an invalid p-value."
        )

      }


      # Prevent exactly 0 or 1 because qnorm() would otherwise
      # become +/- infinity later when using Stouffer's method.
      p <- pmin(
        pmax(
          p,
          1e-10
        ),
        1 - 1e-10
      )


      p
    },
    useZ = useZ,
    Zvals = Zvals,
    npctiles = npctiles,
    K = K
  )


  unlist(
    out,
    use.names = FALSE
  )
}


# ==================================================================================================
# 5. FUNCTION: COMBINE P-VALUES USING STOUFFER'S METHOD
# ==================================================================================================

stouffer_from_pvals <- function(pvals) {

  if (length(pvals) == 0) {
    stop(
      "No p-values were supplied to Stouffer's method."
    )
  }


  # Prevent p = 0 or p = 1.
  pvals <- pmin(
    pmax(
      pvals,
      1e-10
    ),
    1 - 1e-10
  )


  # Convert each p-value into a standard-normal Z-score.
  zi <- qnorm(
    1 - pvals
  )


  # Stouffer combined Z statistic.
  z_stouffer <- sum(zi) /
    sqrt(
      length(zi)
    )


  # Convert the combined Z statistic back into a p-value.
  p_stouffer <- 1 -
    pnorm(
      z_stouffer
    )


  list(

    pvals = pvals,

    zi = zi,

    z_stouffer = z_stouffer,

    p_stouffer = p_stouffer,

    B = length(pvals)

  )
}


# ==================================================================================================
# 6. FUNCTION: READ ONE SIMULATION CHUNK
# ==================================================================================================

read_one_chunk <- function(out_file) {

  sim <- readRDS(
    out_file
  )


  # Each chunk contains a list called "data".
  #
  # Each element of sim$data is one simulated replicate.
  # Stack those replicates into one data frame.
  mdat <- do.call(
    rbind,
    sim$data
  )


  mdat
}


# ==================================================================================================
# 7. FUNCTION: PREPARE THE OBSERVED SIMULATION DATA
# ==================================================================================================

prepare_chunk_obs <- function(mdat) {

  # ------------------------------------------------------------------
  # Keep only observations that were retained:
  #
  # R = 1
  # ------------------------------------------------------------------

  mdat_obs <- mdat[
    mdat$R == 1,
  ]


  # ------------------------------------------------------------------
  # The simulator already defines "size" as the observed cluster size.
  #
  # Therefore size > 1 removes clusters containing only one observed
  # individual.
  # ------------------------------------------------------------------

  mdat_obs <- mdat_obs %>%
    dplyr::filter(
      size > 1
    )


  # ------------------------------------------------------------------
  # Cluster numbers repeat between simulation replicates.
  #
  # For example:
  #
  # replicate 1, cluster 1
  # replicate 2, cluster 1
  #
  # are two completely different clusters.
  #
  # Therefore create a new ID from:
  #
  # replicate + cluster
  # ------------------------------------------------------------------

  mdat_obs$cID <- as.integer(
    factor(
      paste0(
        mdat_obs$replicate,
        "_",
        mdat_obs$cluster
      )
    )
  )


  mdat_obs
}


# ==================================================================================================
# 8. FUNCTION: CREATE THE ASD INPUT DATA
# ==================================================================================================

make_test_dat <- function(
    mdat_obs,
    direction = c(
      "Z_equals_X",
      "Z_equals_Y"
    )
) {

  direction <- match.arg(
    direction
  )


  # ------------------------------------------------------------------
  # Direction 1:
  #
  # outcome = Y
  # Z       = X
  # ------------------------------------------------------------------

  if (direction == "Z_equals_X") {

    test_dat <- cbind(

      cID = mdat_obs$cID,

      Y = mdat_obs$Y,

      Z = mdat_obs$X

    )

  }


  # ------------------------------------------------------------------
  # Direction 2:
  #
  # outcome = X
  # Z       = Y
  # ------------------------------------------------------------------

  if (direction == "Z_equals_Y") {

    test_dat <- cbind(

      cID = mdat_obs$cID,

      Y = mdat_obs$X,

      Z = mdat_obs$Y

    )

  }


  test_dat
}


# ==================================================================================================
# 9. FUNCTION: RUN ONE CHUNK IN ONE DIRECTION
# ==================================================================================================

run_one_chunk_direction <- function(
    out_file,
    direction,
    scenario,
    chunk_id,
    cl,
    m_sub,
    K,
    useZ,
    Zvals,
    npctiles,
    base_seed
) {

  # Read the chunk.
  mdat <- read_one_chunk(
    out_file
  )


  # Keep observed observations and construct cID.
  mdat_obs <- prepare_chunk_obs(
    mdat
  )


  # Build the three-column data set expected by ASDpv().
  test_dat <- make_test_dat(
    mdat_obs = mdat_obs,
    direction = direction
  )


  # ------------------------------------------------------------------
  # Create a different reproducible seed for each:
  #
  # scenario
  # chunk
  # direction
  #
  # The old script always started with seed = 1 for every chunk.
  # Here we avoid reusing the same Monte Carlo seeds repeatedly.
  # ------------------------------------------------------------------

  direction_offset <- if (
    direction == "Z_equals_X"
  ) {

    0L

  } else {

    50000000L

  }


  seed_now <- base_seed +
    scenario * 1000L +
    chunk_id * 10L +
    direction_offset


  # Run the ASD tests.
  pvals <- sub_ASD_pvals(

    mdat = test_dat,

    cl = cl,

    m_sub = m_sub,

    seed = seed_now,

    K = K,

    useZ = useZ,

    Zvals = Zvals,

    npctiles = npctiles

  )


  # Number of clusters represented in this chunk.
  n_clusters <- length(
    unique(
      test_dat[, "cID"]
    )
  )


  # Remove large objects before moving on.
  rm(
    mdat,
    mdat_obs,
    test_dat
  )

  gc()


  list(

    pvals = pvals,

    n_clusters = n_clusters,

    B = length(pvals)

  )
}


# ==================================================================================================
# 10. ORGANIZE THE SIMULATION GRID BY SCENARIO
# ==================================================================================================

# Script 6 saved one row of the grid index for each CHUNK.
#
# We want Script 7 to operate scenario-by-scenario.
#
# Therefore group all chunk files belonging to the same scenario.

scenario_info <- grid %>%

  dplyr::arrange(
    scenario,
    chunk_id
  ) %>%

  dplyr::group_by(
    scenario,
    eta_x,
    eta_y,
    rho_xy,
    rho_uv,
    M
  ) %>%

  dplyr::summarise(

    out_files = list(
      out_file
    ),

    chunk_ids = list(
      chunk_id
    ),

    n_chunks = dplyr::n(),

    .groups = "drop"
  )


cat(
  "Number of simulation scenarios:",
  nrow(scenario_info),
  "\n"
)

cat(
  "Values of M:",
  paste(
    sort(
      unique(scenario_info$M)
    ),
    collapse = ", "
  ),
  "\n\n"
)


# ==================================================================================================
# 11. CREATE THE RESULTS DATA FRAME
# ==================================================================================================

# Each simulation scenario is tested in two directions.
#
# Therefore:
#
# number of result rows =
#
#     number of scenarios x 2

stouffer_results <- data.frame(

  scenario = rep(
    scenario_info$scenario,
    each = 2
  ),

  eta_x = rep(
    scenario_info$eta_x,
    each = 2
  ),

  eta_y = rep(
    scenario_info$eta_y,
    each = 2
  ),

  rho_xy = rep(
    scenario_info$rho_xy,
    each = 2
  ),

  rho_uv = rep(
    scenario_info$rho_uv,
    each = 2
  ),

  M = rep(
    scenario_info$M,
    each = 2
  ),

  direction = rep(
    c(
      "Z_equals_X",
      "Z_equals_Y"
    ),
    times = nrow(scenario_info)
  ),

  n_chunks = rep(
    scenario_info$n_chunks,
    each = 2
  ),

  n_clusters = NA_integer_,

  B = NA_integer_,

  m_sub = m_sub,

  K = K,

  npctiles = npctiles,

  z_stouffer = NA_real_,

  p_stouffer = NA_real_

)


# ==================================================================================================
# 12. CREATE ONE PARALLEL CLUSTER FOR THE ENTIRE SCRIPT
# ==================================================================================================

# Save the library paths from the main R session.
#
# This includes the renv library containing crspack.
current_libpaths <- .libPaths()


# ------------------------------------------------------------------
# IMPORTANT:
#
# "--vanilla" prevents the worker R sessions from reading the
# project's .Rprofile.
#
# That means the workers do NOT independently run:
#
#     source("renv/activate.R")
#
# which was the source of the huge startup problem we saw earlier.
# ------------------------------------------------------------------

cl <- parallel::makePSOCKcluster(

  cores,

  rscript_args = "--vanilla"

)


# Give the workers access to the same library paths as the main session.
parallel::clusterExport(

  cl,

  "current_libpaths",

  envir = environment()

)


parallel::clusterEvalQ(
  cl,
  {

    .libPaths(
      current_libpaths
    )

    library(crspack)

    NULL
  }
)


cat(
  "Parallel cluster created with",
  cores,
  "workers.\n\n"
)


# ==================================================================================================
# 13. RUN THE SIMULATION P-VALUE ANALYSIS
# ==================================================================================================

row_id <- 1


# tryCatch() guarantees that the PSOCK workers are shut down even
# if something fails during the simulation.

tryCatch(
  {


    for (
      s in seq_len(
        nrow(scenario_info)
      )
    ) {


      cat(
        "============================================================\n"
      )

      cat(
        "Running scenario",
        s,
        "of",
        nrow(scenario_info),
        "\n"
      )

      cat(
        "Scenario ID:",
        scenario_info$scenario[s],
        "\n"
      )

      cat(
        "eta_x:",
        scenario_info$eta_x[s],
        "\n"
      )

      cat(
        "eta_y:",
        scenario_info$eta_y[s],
        "\n"
      )

      cat(
        "rho_xy:",
        scenario_info$rho_xy[s],
        "\n"
      )

      cat(
        "rho_uv:",
        scenario_info$rho_uv[s],
        "\n"
      )

      cat(
        "M:",
        scenario_info$M[s],
        "\n"
      )

      cat(
        "Chunks:",
        scenario_info$n_chunks[s],
        "\n"
      )

      cat(
        "Cores:",
        cores,
        "\n\n"
      )


      # Files belonging to the current scenario.
      out_files_now <- scenario_info$out_files[[s]]

      chunk_ids_now <- scenario_info$chunk_ids[[s]]


      # ----------------------------------------------------------------
      # Run both directions.
      # ----------------------------------------------------------------

      for (
        direction_now in c(
          "Z_equals_X",
          "Z_equals_Y"
        )
      ) {


        cat(
          "  Direction:",
          direction_now,
          "\n"
        )


        # All ASD p-values from all chunks belonging to this scenario.
        all_pvals <- numeric(0)


        # Total number of clusters encountered across chunks.
        total_clusters <- 0L


        # --------------------------------------------------------------
        # Loop through the simulation chunks.
        # --------------------------------------------------------------

        for (
          chunk_index in seq_along(
            out_files_now
          )
        ) {


          chunk_id_now <- chunk_ids_now[
            chunk_index
          ]


          cat(
            "    Chunk",
            chunk_index,
            "of",
            length(out_files_now),
            "(chunk ID",
            chunk_id_now,
            ")\n"
          )


          chunk_res <- run_one_chunk_direction(

            out_file = out_files_now[
              chunk_index
            ],

            direction = direction_now,

            scenario = scenario_info$scenario[s],

            chunk_id = chunk_id_now,

            cl = cl,

            m_sub = m_sub,

            K = K,

            useZ = useZ,

            Zvals = Zvals,

            npctiles = npctiles,

            base_seed = base_seed

          )


          # Add this chunk's p-values to the scenario-level vector.
          all_pvals <- c(
            all_pvals,
            chunk_res$pvals
          )


          # Add the number of clusters.
          total_clusters <- total_clusters +
            chunk_res$n_clusters


          rm(
            chunk_res
          )

          gc()


          # ------------------------------------------------------------
          # SAVE PROGRESS AFTER EVERY CHUNK
          #
          # If the script crashes later, the p-values generated up to
          # this point are still available.
          # ------------------------------------------------------------

          saveRDS(

            list(

              scenario = scenario_info$scenario[s],

              direction = direction_now,

              last_chunk = chunk_id_now,

              pvals = all_pvals,

              n_clusters = total_clusters

            ),

            file.path(

              results_dir,

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


        # --------------------------------------------------------------
        # COMBINE ALL ASD P-VALUES FOR THIS SCENARIO
        # --------------------------------------------------------------

        res <- stouffer_from_pvals(
          all_pvals
        )


        # --------------------------------------------------------------
        # STORE THE RESULT
        # --------------------------------------------------------------

        stouffer_results$n_clusters[
          row_id
        ] <- total_clusters


        stouffer_results$B[
          row_id
        ] <- res$B


        stouffer_results$z_stouffer[
          row_id
        ] <- res$z_stouffer


        stouffer_results$p_stouffer[
          row_id
        ] <- res$p_stouffer


        cat(
          "\n"
        )

        cat(
          "  Total clusters =",
          total_clusters,
          "\n"
        )

        cat(
          "  Total ASD tests =",
          res$B,
          "\n"
        )

        cat(
          "  z_stouffer =",
          res$z_stouffer,
          "\n"
        )

        cat(
          "  p_stouffer =",
          res$p_stouffer,
          "\n\n"
        )


        # --------------------------------------------------------------
        # SAVE THE PARTIALLY COMPLETED RESULTS TABLE
        # --------------------------------------------------------------

        saveRDS(

          stouffer_results,

          file.path(
            results_dir,
            "simulation_stouffer_results_partial.rds"
          )

        )


        write.csv(

          stouffer_results,

          file.path(
            results_dir,
            "simulation_stouffer_results_partial.csv"
          ),

          row.names = FALSE
        )


        row_id <- row_id + 1


        rm(
          all_pvals,
          res
        )

        gc()

      }

    }

  },

  finally = {

    parallel::stopCluster(
      cl
    )

    cat(
      "\nParallel cluster stopped.\n"
    )

  }
)


# ==================================================================================================
# 14. SAVE FINAL RESULTS
# ==================================================================================================

saveRDS(

  stouffer_results,

  file.path(
    results_dir,
    "simulation_stouffer_results.rds"
  )

)


write.csv(

  stouffer_results,

  file.path(
    results_dir,
    "simulation_stouffer_results.csv"
  ),

  row.names = FALSE
)


# ==================================================================================================
# 15. DISPLAY RESULTS
# ==================================================================================================

stouffer_results