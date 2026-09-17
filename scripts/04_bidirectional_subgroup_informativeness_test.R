####################################################################################################
#
# NHANES SUBGROUP ASSOCIATION TESTING
#
# This script tests associations between the caries variables
#
#   DS, FS, DSI, FSI
#
# and the periodontal variables
#
#   maxCAL, maxPD
#
# using crspack::ASDpv().
#
# Subjects are divided into smaller groups of clusters, ASDpv() is
# run separately within each group, and the resulting p-values are
# combined using Stouffer's method.
#
####################################################################################################


# ==================================================================================================
# 1. LOAD PACKAGES
# ==================================================================================================

library(crspack)
library(dplyr)
library(parallel)


# ==================================================================================================
# 2. ANALYSIS SETTINGS
# ==================================================================================================

# Detect the number of physical CPU cores available.
detected_cores <- parallel::detectCores(logical = FALSE)

if (is.na(detected_cores)) {
  detected_cores <- 2
}


# Leave one physical core available for the operating system.
cores <- max(1, detected_cores - 1)


# Number of subjects included in each subgroup.
m_sub <- 10


# Number of permutations used by ASDpv().
#
# Start small when testing the script locally.
K <- 5

# For the full analysis:
# K <- 100


# Number of percentile cutoffs used by ASDpv().
npctiles <- 10


# Random seed used when assigning seeds to subgroup analyses.
seed <- 1


# ==================================================================================================
# 3. FILE PATHS
# ==================================================================================================

# Script 2 produces the tooth-level NHANES dataset needed here.
#
# The ASD analysis requires repeated tooth-level observations within
# each subject, so person_level_data.RData should NOT be used.

data_file <- file.path(
  "data",
  "Tooth_Level_data.RData"
)


# Store ASD-specific output inside the project's existing results folder.
results_dir <- file.path(
  "results",
  "ASD"
)

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Full R results.
save_file <- file.path(
  results_dir,
  "ASD_cross_set_results.rds"
)


# Simple summary table.
summary_file <- file.path(
  results_dir,
  "ASD_cross_set_results_summary.csv"
)


# ==================================================================================================
# 4. LOAD TOOTH-LEVEL DATA
# ==================================================================================================

load(data_file)


# Script 2 should have created a tooth-level data frame.
#
# This assumes the saved object is named df_tooth_level.
#
# If Script 2 saved it under a different object name, only this
# reference needs to be changed.

df_clean <- df_tooth_level


# ==================================================================================================
# 5. KEEP COMPLETE TOOTH-LEVEL OBSERVATIONS
# ==================================================================================================

# Keep only the variables required for this analysis before removing
# missing values.
#
# This is preferable to na.omit(df_tooth_level), because unrelated
# variables in the larger dataset should not cause a tooth to be removed.

df_clean <- df_clean %>%
  select(
    SEQN,
    DS,
    FS,
    DSI,
    FSI,
    maxCAL,
    maxPD
  ) %>%
  na.omit()


# Keep subjects having at least 10 usable tooth-level observations.
#
# This preserves the criterion used in the previous analysis.

df_clean <- df_clean %>%
  group_by(SEQN) %>%
  filter(n() > 9) %>%
  ungroup()


# ==================================================================================================
# 6. DEFINE VARIABLE SETS
# ==================================================================================================

# Caries variables.
set1 <- c(
  "DS",
  "FS",
  "DSI",
  "FSI"
)


# Periodontal variables.
set2 <- c(
  "maxCAL",
  "maxPD"
)


# ==================================================================================================
# 7. DEFINE PRESPECIFIED Z CUTPOINTS
# ==================================================================================================

# DS, FS, DSI, and FSI are binary variables.
#
# Therefore, when one of these variables is used as Z, the natural
# separating cutoff is 0.5.


# maxPD and maxCAL are discrete periodontal measurements.
#
# The cutoffs are placed halfway between consecutive observed values.
#
# For example:
#
#   observed values:
#
#       0, 1, 2, 3
#
#   resulting cutoffs:
#
#       0.5, 1.5, 2.5


pd_unique <- sort(
  unique(df_clean$maxPD)
)

pdvals <- pd_unique + 0.5

pdvals <- pdvals[
  -length(pdvals)
]


cal_unique <- sort(
  unique(df_clean$maxCAL)
)

calvals <- cal_unique + 0.5

calvals <- calvals[
  -length(calvals)
]


# Function returning the correct cutoffs according to which
# variable is being used as Z.

get_zvals <- function(zvar) {
  if (zvar %in% set1) {
    return(0.5)
  }

  if (zvar == "maxPD") {
    return(pdvals)
  }

  if (zvar == "maxCAL") {
    return(calvals)
  }

  stop(
    "No Zvals rule defined for: ",
    zvar
  )
}


# ==================================================================================================
# 8. DEFINE CROSS-SET TESTS
# ==================================================================================================

# Caries variable as Y and periodontal variable as Z.
pairs_1_to_2 <- expand.grid(
  Y = set1,
  Z = set2,
  stringsAsFactors = FALSE
)

# Periodontal variable as Y and caries variable as Z.
pairs_2_to_1 <- expand.grid(
  Y = set2,
  Z = set1,
  stringsAsFactors = FALSE
)

# Combine the two directions.
pairs_all <- rbind(
  pairs_1_to_2,
  pairs_2_to_1
)

# ==================================================================================================
# 9. DEFINE SUBGROUP ASD FUNCTION
# ==================================================================================================

sub_ASD <- function(
    mdat,
    cl = NULL,
    m_sub = 50,
    seed = 1,
    npctiles = 10,
    K = 50,
    useZ = "prespecified",
    Zvals = NULL) {


  # -----------------------------------------------------------------------------------------------
  # Find the unique clusters.
  #
  # Here, each cluster corresponds to one NHANES participant.
  # -----------------------------------------------------------------------------------------------

  clusters <- unique(
    mdat[, "cID"]
  )


  # -----------------------------------------------------------------------------------------------
  # Divide the subjects into groups containing approximately m_sub
  # subjects each.
  # -----------------------------------------------------------------------------------------------

  cluster_groups <- split(
    clusters,
    ceiling(
      seq_along(clusters) / m_sub
    )
  )


  # -----------------------------------------------------------------------------------------------
  # Construct one tooth-level dataset for each group of subjects.
  # -----------------------------------------------------------------------------------------------

  dat_groups <- lapply(
    cluster_groups,
    function(ids) {

      mdat[
        mdat[, "cID"] %in% ids,
        ,
        drop = FALSE
      ]

    }
  )


  # Number of subgroup analyses.
  B <- length(dat_groups)


  # -----------------------------------------------------------------------------------------------
  # Generate reproducible seeds for the individual subgroup analyses.
  # -----------------------------------------------------------------------------------------------

  set.seed(seed)

  seeds <- sample.int(
    1e8,
    B
  )


  # -----------------------------------------------------------------------------------------------
  # Sequential version.
  #
  # This is useful if cores = 1 or if parallelization causes problems.
  # -----------------------------------------------------------------------------------------------

  if (is.null(cl)) {

    out <- lapply(
      seq_len(B),
      function(b) {

        set.seed(
          seeds[b]
        )


        x <- crspack::ASDpv(
          mdat = dat_groups[[b]],
          useZ = useZ,
          Zvals = Zvals,
          npctiles = npctiles,
          K = K,
          parallel = FALSE
        )


        p <- pmin(
          pmax(
            x$pvalue,
            1e-10
          ),
          1 - 1e-10
        )


        return(p)

      }
    )


  } else {


    # ---------------------------------------------------------------------------------------------
    # Export everything needed by the workers.
    # ---------------------------------------------------------------------------------------------

    parallel::clusterExport(
      cl,
      varlist = c(
        "dat_groups",
        "useZ",
        "Zvals",
        "npctiles",
        "K",
        "seeds"
      ),
      envir = environment()
    )


    # ---------------------------------------------------------------------------------------------
    # Run subgroup analyses in parallel.
    #
    # ASDpv() itself uses parallel = FALSE because we are already
    # parallelizing at the subgroup level.
    # ---------------------------------------------------------------------------------------------

    out <- parallel::parLapply(
      cl,
      seq_len(B),
      function(b) {

        set.seed(
          seeds[b]
        )


        x <- crspack::ASDpv(
          mdat = dat_groups[[b]],
          useZ = useZ,
          Zvals = Zvals,
          npctiles = npctiles,
          K = K,
          parallel = FALSE
        )


        p <- pmin(
          pmax(
            x$pvalue,
            1e-10
          ),
          1 - 1e-10
        )


        return(p)

      }
    )

  }


  # -----------------------------------------------------------------------------------------------
  # Combine the subgroup p-values.
  # -----------------------------------------------------------------------------------------------

  pvals <- unlist(out)


  # Convert the one-sided p-values to Z statistics:
  #
  #   Z_b = Phi^(-1)(1 - p_b)

  zi <- qnorm(
    1 - pvals
  )


  # Stouffer's statistic:
  #
  #                 sum Z_b
  #   Z_Stouffer = -----------
  #                  sqrt(B)

  z_stouffer <- sum(zi) /
    sqrt(length(zi))


  # Corresponding one-sided p-value.

  p_stouffer <- 1 -
    pnorm(z_stouffer)


  return(
    list(
      pvals = pvals,
      zi = zi,
      z_stouffer = z_stouffer,
      p_stouffer = p_stouffer
    )
  )

}

# ==================================================================================================
# 10. CREATE PARALLEL CLUSTER
# ==================================================================================================

# Save the library paths from the current R session.
#
# The first path is the renv library for this project.
current_libpaths <- .libPaths()


# Start PSOCK workers without reading .Rprofile.
#
# This prevents every worker from independently trying to activate
# the same renv project at the same time.
cl <- parallel::makePSOCKcluster(
  cores,
  rscript_args = "--vanilla"
)


# Give each worker the same package-library paths as the main R session.
parallel::clusterCall(
  cl,
  function(paths) {
    .libPaths(paths)
  },
  current_libpaths
)


# Load crspack on every worker.
parallel::clusterEvalQ(
  cl,
  {
    library(crspack)
    NULL
  }
)

# ==================================================================================================
# 11. RUN ALL CROSS-SET TESTS
# ==================================================================================================

results <- list()


for (r in seq_len(nrow(pairs_all))) {


  # Variables for the current comparison.

  yvar <- pairs_all$Y[r]

  zvar <- pairs_all$Z[r]


  # Name used to store the result.

  result_name <- paste0(
    "Y_",
    yvar,
    "__Z_",
    zvar
  )


  # -----------------------------------------------------------------------------------------------
  # ASDpv() expects a matrix containing:
  #
  #   cID = cluster / subject ID
  #   Y   = response variable
  #   Z   = subgroup-defining covariate
  #
  # Each subject therefore appears multiple times, once for each
  # tooth contributing to the analysis.
  # -----------------------------------------------------------------------------------------------

  mdat <- cbind(
    cID = df_clean$SEQN,
    Y   = df_clean[[yvar]],
    Z   = df_clean[[zvar]]
  )


  # -----------------------------------------------------------------------------------------------
  # Run the subgroup ASD procedure.
  # -----------------------------------------------------------------------------------------------

  results[[result_name]] <- sub_ASD(
    mdat = mdat,
    cl = cl,
    m_sub = m_sub,
    seed = seed,
    npctiles = npctiles,
    K = K,
    useZ = "prespecified",
    Zvals = get_zvals(zvar)
  )


  # -----------------------------------------------------------------------------------------------
  # Save after every completed comparison.
  #
  # This prevents loss of completed work if a later comparison takes
  # a long time or the R session stops.
  # -----------------------------------------------------------------------------------------------

  saveRDS(
    list(
      results = results,
      pairs_all = pairs_all
    ),
    file = save_file
  )

}


# ==================================================================================================
# 12. STOP PARALLEL CLUSTER
# ==================================================================================================

parallel::stopCluster(
  cl
)


# ==================================================================================================
# 13. CREATE SUMMARY TABLE
# ==================================================================================================

summary_res <- data.frame(

  comparison = names(results),

  Y = pairs_all$Y,

  Z = pairs_all$Z,

  z_stouffer = sapply(
    results,
    function(x) {

      x$z_stouffer

    }
  ),

  p_stouffer = sapply(
    results,
    function(x) {

      x$p_stouffer

    }
  ),

  row.names = NULL

)


# ==================================================================================================
# 14. SAVE SUMMARY RESULTS
# ==================================================================================================

write.csv(
  summary_res,
  summary_file,
  row.names = FALSE
)


summary_res


# ==================================================================================================
# 15. CREATE AND SAVE RESULTS TABLE
# ==================================================================================================

library(gt)


results_table <- summary_res %>%
  select(
    Y,
    Z,
    z_stouffer,
    p_stouffer
  ) %>%
  gt() %>%
  cols_label(
    Y = "Y Variable",
    Z = "Z Variable",
    z_stouffer = "Stouffer Z",
    p_stouffer = "P-value"
  ) %>%
  fmt_number(
    columns = z_stouffer,
    decimals = 3
  ) %>%
  fmt_number(
    columns = p_stouffer,
    decimals = 4
  ) %>%
  tab_header(
    title = "NHANES Cross-Set Association Tests"
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 18,
    data_row.padding = px(6)
  )


# Display the table in RStudio.
results_table


# Save the table as an image.
gtsave(
  results_table,
  filename = file.path(
    results_dir,
    "ASD_cross_set_results_table.png"
  )
)
