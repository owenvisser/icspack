
# Generate the 9-setting simulation using clustPear's published methods only.

stopifnot(requireNamespace("here", quietly = TRUE))
stopifnot(requireNamespace("clustPear", quietly = TRUE))

source(here::here("R", "simulation.R"))

cfg <- simulation_config()
settings <- simulation_settings(cfg)

root <- here::here()
folder <- file.path(root, "results", "simulation", cfg$name, "chunks")

dir.create(folder, recursive = TRUE, showWarnings = FALSE)

nb <- ceiling(cfg$reps / cfg$block)

jobs <- expand.grid(
  setting = seq_len(nrow(settings)),
  block = seq_len(nb),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

needs <- !vapply(
  seq_len(nrow(jobs)),
  function(i) {
    file.exists(file.path(
      folder,
      sprintf(
        "setting_%02d_block_%03d.rds",
        jobs$setting[i],
        jobs$block[i]
      )
    ))
  },
  logical(1)
)

remaining <- jobs[needs, , drop = FALSE]

message(
  "Simulation: ", cfg$mode, "; ", nrow(settings), " settings x ",
  cfg$reps, " replicates x ", length(cfg$K), " K values"
)

message("Completed blocks: ", nrow(jobs) - nrow(remaining), "/", nrow(jobs))


run_remaining <- function() {
  if (!nrow(remaining)) {
    return(invisible(NULL))
  }

  workers <- min(as.integer(cfg$workers), nrow(remaining))

  if (workers > 1L) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      c("simulate_and_estimate", "run_simulation_block"),
      envir = environment(run_simulation_block)
    )

    parallel::clusterEvalQ(cl, library(clustPear))
  }

  waves <- split(
    seq_len(nrow(remaining)),
    ceiling(seq_len(nrow(remaining)) / workers)
  )

  for (j in seq_along(waves)) {
    block_jobs <- lapply(
      waves[[j]],
      function(k) as.list(remaining[k, ])
    )

    if (workers > 1L) {
      parallel::parLapplyLB(
        cl,
        block_jobs,
        run_simulation_block,
        cfg = cfg,
        settings = settings,
        folder = folder
      )
    } else {
      lapply(
        block_jobs,
        run_simulation_block,
        cfg = cfg,
        settings = settings,
        folder = folder
      )
    }

    message(sprintf(
      "Checkpoint blocks: %d/%d complete",
      nrow(jobs) - nrow(remaining) + sum(lengths(waves[seq_len(j)])),
      nrow(jobs)
    ))
  }
}

run_remaining()

message("All simulation blocks completed. Next: scripts/06_summarize_simulation.R")
