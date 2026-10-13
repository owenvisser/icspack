####################################################################################################
#
# 04_bidirectional_subgroup_informativeness_test.R
#
# NHANES BIDIRECTIONAL DS / FS / PERIODONTAL SUBGROUP INFORMATIVENESS TESTING
#
# Tests DS, DSI, FS, and FSI against maxCAL and maxPD in both directions
# using crspack::ASDpv().
#
# When a binary caries measure defines the subgroup, the only cutpoint is 0.5.
# When maxCAL or maxPD defines the subgroup, every distinct observed split is tested.
# Because the periodontal measures are integer-valued, these cutpoints are formed by
# adding 0.5 to every observed value except the maximum.
#
# Subjects are divided into independent groups of 10 clusters. ASDpv() is run
# separately within each group and the resulting p-values are combined using
# Stouffer's method.
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
K <- 100

# For a quick local code check:
# K <- 5


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
  "Tooth_Level_Data.RData"
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

# Caries variables are binary, so zero versus nonzero is separated at 0.5.
#
# Periodontal variables use the same three clinically defined subgroup
# cutpoints used in the association analysis:
#
#   maxCAL: 2.5, 3.5, 4.5  -> CAL >= 3, >= 4, >= 5
#   maxPD:  3.5, 4.5, 5.5  -> PD  >= 4, >= 5, >= 6

get_zvals <- function(zvar) {
  if (zvar %in% set1) {
    return(0.5)
  }

  if (zvar == "maxCAL") {
    return(
      c(
        2.5,
        3.5,
        4.5
      )
    )
  }

  if (zvar == "maxPD") {
    return(
      c(
        3.5,
        4.5,
        5.5
      )
    )
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

pairs_all$direction <- c(
  rep(
    "Caries Y / Periodontal Z",
    nrow(pairs_1_to_2)
  ),
  rep(
    "Periodontal Y / Caries Z",
    nrow(pairs_2_to_1)
  )
)

# Record the cutpoints used for each subgroup-defining variable.
pairs_all$Z_cutpoints <- vapply(
  pairs_all$Z,
  function(zvar) {
    paste(
      get_zvals(zvar),
      collapse = ", "
    )
  },
  character(1)
)


# Create the comparison names once and verify that every requested test is unique.
comparison_names <- paste0(
  "Y_",
  pairs_all$Y,
  "__Z_",
  pairs_all$Z
)

if (anyDuplicated(comparison_names)) {
  stop("Duplicate comparison names were generated.")
}


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

current_libpaths <- .libPaths()

cl <- parallel::makePSOCKcluster(
  cores,
  rscript_args = "--vanilla"
)

parallel::clusterCall(
  cl,
  function(paths) {
    .libPaths(paths)
  },
  current_libpaths
)

parallel::clusterEvalQ(
  cl,
  {
    library(crspack)
    NULL
  }
)


# ==================================================================================================
# 11. RUN ALL BIDIRECTIONAL TESTS
# ==================================================================================================

# Preallocate one result slot for every requested comparison.
results <- setNames(
  vector("list", nrow(pairs_all)),
  comparison_names
)

# Always close the PSOCK workers, including after an error or user interrupt.
tryCatch(
  {
    for (r in seq_len(nrow(pairs_all))) {
      yvar <- pairs_all$Y[r]
      zvar <- pairs_all$Z[r]
      result_name <- comparison_names[r]

      cat(
        sprintf(
          "\n[%d/%d] Running %s\n",
          r,
          nrow(pairs_all),
          result_name
        )
      )
      flush.console()

      mdat <- cbind(
        cID = df_clean$SEQN,
        Y   = df_clean[[yvar]],
        Z   = df_clean[[zvar]]
      )

      results[[r]] <- sub_ASD(
        mdat = mdat,
        cl = cl,
        m_sub = m_sub,
        seed = seed,
        npctiles = npctiles,
        K = K,
        useZ = "prespecified",
        Zvals = get_zvals(zvar)
      )

      cat(
        sprintf(
          "[%d/%d] Completed %s\n",
          r,
          nrow(pairs_all),
          result_name
        )
      )
      flush.console()

      # Save a checkpoint after every completed comparison.
      saveRDS(
        list(
          results = results,
          pairs_all = pairs_all,
          comparison_names = comparison_names
        ),
        file = save_file
      )
    }
  },
  finally = {
    try(
      parallel::stopCluster(cl),
      silent = TRUE
    )
  }
)

# Do not construct a summary from an incomplete run. This prevents R from
# silently recycling a partial results vector across all rows of pairs_all.
completed <- !vapply(
  results,
  is.null,
  logical(1)
)

if (!all(completed)) {
  stop(
    sprintf(
      "Only %d of %d comparisons completed. The checkpoint RDS was saved, but no summary CSV was written.",
      sum(completed),
      length(completed)
    )
  )
}


# ==================================================================================================
# 13. CREATE SUMMARY TABLE
# ==================================================================================================

summary_res <- data.frame(
  comparison = comparison_names,
  direction = pairs_all$direction,
  Y = pairs_all$Y,
  Z = pairs_all$Z,
  Z_cutpoints = pairs_all$Z_cutpoints,
  z_stouffer = vapply(
    results,
    function(x) x$z_stouffer,
    numeric(1)
  ),
  p_stouffer = vapply(
    results,
    function(x) x$p_stouffer,
    numeric(1)
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
# 15. OPTIONAL RESULTS TABLE IMAGE
# ==================================================================================================

# The CSV above is the reproducible analysis output. Set this to TRUE only
# when a rendered PNG summary is also desired. Keeping it FALSE avoids
# possible webshot/chromote processes keeping the R session open.
save_table_png <- FALSE

if (save_table_png) {
  library(gt)

  results_table <- summary_res %>%
    select(
      direction,
      Y,
      Z,
      Z_cutpoints,
      z_stouffer,
      p_stouffer
    ) %>%
    gt() %>%
    cols_label(
      direction = "Direction",
      Y = "Y Variable",
      Z = "Z Variable",
      Z_cutpoints = "Z Cutpoints",
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
      title = "NHANES Bidirectional DS / FS / Periodontal Subgroup Informativeness Tests"
    ) %>%
    tab_options(
      table.font.size = 14,
      heading.title.font.size = 18,
      data_row.padding = px(6)
    )

  gtsave(
    results_table,
    filename = file.path(
      results_dir,
      "ASD_cross_set_results_table.png"
    )
  )
}
