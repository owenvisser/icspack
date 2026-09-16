# File for testing NHANES for subgroup correlation using Sams test.

sub_ASD <- function(mdat, m_sub = 50, cores = 1, seed = 1,
                    npctiles = 10, K = 50, useZ = "prespecified", Zvals = NULL) {
  
  clusters <- unique(mdat[, "cID"])
  
  cluster_groups <- split(
    clusters,
    ceiling(seq_along(clusters) / m_sub)
  )
  
  dat_groups <- lapply(cluster_groups, function(x) {
    mdat[mdat[, "cID"] %in% x, ]
  })
  
  B <- length(dat_groups)
  
  set.seed(seed)
  seeds <- sample.int(1e8, B)
  
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl))
  
  parallel::clusterEvalQ(cl, {
    .libPaths(c("/orange/somnath.datta/NHANES/Rpackages", .libPaths()))
    library(crspack)
  })
  
  parallel::clusterExport(
    cl,
    c("dat_groups", "useZ", "npctiles", "K", "Zvals", "seeds"),
    envir = environment()
  )
  
  out <- parallel::parLapply(cl, seq_len(B), function(b) {
    set.seed(seeds[b])
    
    x <- crspack::ASDpv(
      mdat = dat_groups[[b]],
      useZ = useZ,
      Zvals = Zvals,
      npctiles = npctiles,
      K = K,
      parallel = FALSE
    )
    
    p <- pmin(pmax(x$pvalue, 1e-10), 1 - 1e-10)
    p
  })
  
  pvals <- unlist(out)
  
  zi <- qnorm(1 - pvals)
  
  z_stouffer <- sum(zi) / sqrt(length(zi))
  p_stouffer <- 1 - pnorm(z_stouffer)
  
  list(
    pvals = pvals,
    zi = zi,
    z_stouffer = z_stouffer,
    p_stouffer = p_stouffer
  )
}


############################################################
# HiperGator setup
############################################################

custom_lib <- "/orange/somnath.datta/NHANES/Rpackages"

.libPaths(c(custom_lib, .libPaths()))

Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"))
Sys.setenv(R_LIBS_USER = custom_lib)

library(crspack)
library(dplyr)
library(parallel)

cores <- 50

base_dir <- "/orange/somnath.datta/NHANES"
data_file <- file.path(base_dir, "FullData", "Tooth_Level_Data.RData")

res_dir <- file.path(base_dir, "ASD_Results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

save_file <- file.path(res_dir, "ASD_cross_set_results_DISC_CUTS.RData")
summary_file <- file.path(res_dir, "ASD_cross_set_results_DISC_CUTS_summary.csv")

setwd(base_dir)

load(data_file)

############################################################
# Clean data
############################################################

df_clean <- na.omit(df_tooth_level)

df_clean <- df_clean %>%
  group_by(SEQN) %>%
  filter(n() > 9) %>%
  ungroup()

############################################################
# Define Z values
############################################################

pdvals <- sort(unique(df_clean$maxPD) + 0.5)[-length(unique(df_clean$maxPD))]
calvals <- sort(unique(df_clean$maxCAL) + 0.5)[-length(unique(df_clean$maxCAL))]

set1 <- c("DS", "FS", "DSI", "FSI")
set2 <- c("maxCAL", "maxPD")

get_zvals <- function(zvar) {
  if (zvar %in% set1) {
    return(c(0.5))
  }
  
  if (zvar == "maxPD") {
    return(pdvals)
  }
  
  if (zvar == "maxCAL") {
    return(calvals)
  }
  
  stop("No Zvals rule defined for: ", zvar)
}

############################################################
# Define cross-set tests
############################################################

pairs_1_to_2 <- expand.grid(
  Y = set1,
  Z = set2,
  stringsAsFactors = FALSE
)

pairs_2_to_1 <- expand.grid(
  Y = set2,
  Z = set1,
  stringsAsFactors = FALSE
)

pairs_all <- rbind(pairs_1_to_2, pairs_2_to_1)

############################################################
# Run tests
############################################################

results <- list()

for (r in seq_len(nrow(pairs_all))) {
  
  yvar <- pairs_all$Y[r]
  zvar <- pairs_all$Z[r]
  
  result_name <- paste0("Y_", yvar, "__Z_", zvar)
  
  cat("Running:", result_name, "\n")
  cat("Y =", yvar, "| Z =", zvar, "| cores =", cores, "\n")
  
  results[[result_name]] <- sub_ASD(
    mdat = cbind(
      cID = df_clean$SEQN,
      Y   = df_clean[[yvar]],
      Z   = df_clean[[zvar]]
    ),
    m_sub = 10,
    K = 100,
    useZ = "prespecified",
    Zvals = get_zvals(zvar),
    npctiles = 10,
    cores = cores
  )
  
  save(
    results,
    pairs_all,
    file = save_file
  )
}

############################################################
# Save final results
############################################################

save(
  results,
  pairs_all,
  file = save_file
)

############################################################
# Reload and summarize
############################################################

load(save_file)

summary_res <- data.frame(
  comparison = names(results),
  Y = pairs_all$Y,
  Z = pairs_all$Z,
  z_stouffer = sapply(results, function(x) x$z_stouffer),
  p_stouffer = sapply(results, function(x) x$p_stouffer),
  row.names = NULL
)

write.csv(summary_res, summary_file, row.names = FALSE)

print(summary_res)