####################################################################################################
#
# 04b_bidirectional_dfs_periodontal_informativeness_test.R
#
# NHANES BIDIRECTIONAL DFS / PERIODONTAL SUBGROUP INFORMATIVENESS TESTING
#
# This script tests associations between the caries variables
#
#   DFS, DFSI
#
# and the periodontal variables
#
#   maxCAL, maxPD
#
# in BOTH directions using crspack::ASDpv().
#
# When maxCAL or maxPD is the subgroup-defining variable Z, the prespecified
# cutpoints reproduce the original periodontal dichotomizations used in 05b/05d:
#
#   maxCAL >= 3, >= 4, >= 5  -> Z cutpoints 2.5, 3.5, 4.5
#   maxPD  >= 4, >= 5, >= 6  -> Z cutpoints 3.5, 4.5, 5.5
#
# The half-integer cutpoints are used because ASDpv() splits a discrete integer
# periodontal measure at Z <= z versus Z > z. Thus, for example, z = 2.5 gives
# exactly the same split as maxCAL < 3 versus maxCAL >= 3.
#
# When DFS or DFSI is the subgroup-defining variable Z, the natural binary
# separating cutoff is 0.5.
#
# Subjects are divided into smaller groups of clusters, ASDpv() is run separately
# within each group, and the resulting p-values are combined using Stouffer's method.
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

data_file <- file.path(
  "data",
  "Tooth_Level_Data.RData"
)


results_dir <- file.path(
  "results",
  "ASD"
)


dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


save_file <- file.path(
  results_dir,
  "ASD_DFS_periodontal_bidirectional_results.rds"
)


summary_file <- file.path(
  results_dir,
  "ASD_DFS_periodontal_bidirectional_results_summary.csv"
)


# ==================================================================================================
# 4. LOAD TOOTH-LEVEL DATA
# ==================================================================================================

load(data_file)


df_clean <- df_tooth_level


# ==================================================================================================
# 5. KEEP COMPLETE TOOTH-LEVEL OBSERVATIONS
# ==================================================================================================

# Retain only the variables required for this analysis before removing missing values.
df_clean <- df_clean %>%
  select(
    SEQN,
    DFS,
    DFSI,
    maxCAL,
    maxPD
  ) %>%
  na.omit()


# Keep subjects having at least 10 usable tooth-level observations,
# matching the criterion used in 04.
df_clean <- df_clean %>%
  group_by(SEQN) %>%
  filter(n() > 9) %>%
  ungroup()


# ==================================================================================================
# 6. DEFINE VARIABLE SETS
# ==================================================================================================

# Binary caries variables.
set1 <- c(
  "DFS",
  "DFSI"
)


# Original continuous tooth-level periodontal variables.
set2 <- c(
  "maxCAL",
  "maxPD"
)


# ==================================================================================================
# 7. DEFINE PRESPECIFIED Z CUTPOINTS
# ==================================================================================================

# DFS and DFSI are binary, so their separating cutoff is 0.5.
#
# For maxCAL and maxPD, use only the clinical dichotomizations used in the
# weighted association analyses rather than every possible observed split.
# Because the periodontal measures take integer values, the half-integer
# cutpoints reproduce the desired >= thresholds exactly.

get_zvals <- function(zvar) {

  if (zvar %in% set1) {
    return(0.5)
  }

  if (zvar == "maxCAL") {
    return(
      c(
        2.5,  # CAL >= 3
        3.5,  # CAL >= 4
        4.5   # CAL >= 5
      )
    )
  }

  if (zvar == "maxPD") {
    return(
      c(
        3.5,  # PD >= 4
        4.5,  # PD >= 5
        5.5   # PD >= 6
      )
    )
  }

  stop(
    "No Zvals rule defined for: ",
    zvar
  )
}


# ==================================================================================================
# 8. DEFINE BIDIRECTIONAL CROSS-SET TESTS
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


# Record the actual prespecified Z cutpoints for the output table.
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
  # Here, each cluster corresponds to one NHANES participant.
  # -----------------------------------------------------------------------------------------------

  clusters <- unique(
    mdat[, "cID"]
  )


  # -----------------------------------------------------------------------------------------------
  # Divide the subjects into groups containing approximately m_sub subjects each.
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
    # ASDpv() itself uses parallel = FALSE because parallelization already occurs here.
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
  # Combine the subgroup p-values using Stouffer's method.
  # -----------------------------------------------------------------------------------------------

  pvals <- unlist(out)


  zi <- qnorm(
    1 - pvals
  )


  z_stouffer <- sum(zi) /
    sqrt(length(zi))


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

results <- list()


for (r in seq_len(nrow(pairs_all))) {


  yvar <- pairs_all$Y[r]
  zvar <- pairs_all$Z[r]


  result_name <- paste0(
    "Y_",
    yvar,
    "__Z_",
    zvar
  )


  cat(
    "\nRunning:",
    result_name,
    "\n"
  )


  mdat <- cbind(
    cID = df_clean$SEQN,
    Y   = df_clean[[yvar]],
    Z   = df_clean[[zvar]]
  )


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


  # Save after every completed comparison so completed work is retained.
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

  direction = pairs_all$direction,

  Y = pairs_all$Y,

  Z = pairs_all$Z,

  Z_cutpoints = pairs_all$Z_cutpoints,

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
    title = "NHANES Bidirectional DFS / Periodontal Subgroup Informativeness Tests"
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 18,
    data_row.padding = px(6)
  )


results_table


gtsave(
  results_table,
  filename = file.path(
    results_dir,
    "ASD_DFS_periodontal_bidirectional_results_table.png"
  )
)
