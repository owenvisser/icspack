
# clustPear::simulate_clustered() implements equations (dgp) and (retention).
# clustPear::cluster_pearson() provides every estimate and delta-method SE.

simulation_config <- function() {
  mode <- Sys.getenv("SIM_MODE", "pilot")

  if (!mode %in% c("pilot", "full")) {
    stop("SIM_MODE must be pilot or full")
  }

  full <- identical(mode, "full")

  list(
    mode = mode,
    name = if (full) "full" else "pilot",
    M = if (full) 100L else 20L,
    nmax = if (full) 100L else 30L,
    reps = if (full) 1000L else 3L,
    block = if (full) 100L else 3L,
    rho_uv = c(0, 0.25, 0.50),
    eta = c(0, 2, 4),
    sigma2 = 0.5,
    tau2 = 0.5,
    K = c(2L, 3L, 5L, 10L, 20L, 30L, 50L, 100L),
    er_max_pairs = if (full) 100000L else 4096L,
    seed = 20261103L,
    workers = if (full) 5L else 1L
  )
}

simulation_settings <- function(cfg) {
  expand.grid(
    rho_uv = cfg$rho_uv,
    eta = cfg$eta,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

# Run one generated dataset. All four methods share EXACTLY the same sample.

simulate_and_estimate <- function(s, r, cfg, settings) {
  rho_uv <- settings$rho_uv[s]
  eta <- settings$eta[s]
  seed <- as.integer(cfg$seed + s * 100000L + r)

  d <- clustPear::simulate_clustered(
    M = cfg$M,
    nmax = cfg$nmax,
    rho_uv = rho_uv,
    eta = eta,
    sigma2 = cfg$sigma2,
    tau2 = cfg$tau2,
    seed = seed
  )

  retained_fraction <- nrow(d) / (cfg$M * cfg$nmax)

  n <- table(d$cluster_id)

  d <- d[
    as.character(d$cluster_id) %in% names(n)[n >= 2L],
    ,
    drop = FALSE
  ]

  exclusion_fraction <- 1 - length(unique(d$cluster_id)) / cfg$M

  if (nrow(d) < 3L || length(unique(d$cluster_id)) < 2L) {
    stop("Too few retained observations (setting ", s, ", replicate ", r, ")")
  }

  ref_K <- 20L

  base <- clustPear::cluster_pearson(
    d,
    "cluster_id",
    "X",
    "Y",
    K = ref_K,
    er_max_pairs = cfg$er_max_pairs,
    method = "all"
  )

  # NW, CW and ER do not vary with K. EQ is recomputed for every other K.
  # No bootstrap, no jackknife: all SEs come from clustPear's delta method.

  data <- vector("list", length(cfg$K))

  for (j in seq_along(cfg$K)) {
    K <- cfg$K[j]

    eq <- if (K == ref_K) {
      base[base$method == "EQ", , drop = FALSE]
    } else {
      clustPear::cluster_pearson(
        d,
        "cluster_id",
        "X",
        "Y",
        method = "EQ",
        K = K,
        er_max_pairs = cfg$er_max_pairs
      )
    }

    others <- base[
      base$method %in% c("NW", "CW", "ER"),
      ,
      drop = FALSE
    ]

    fit <- rbind(others, eq)

    fit <- fit[
      match(c("NW", "CW", "EQ", "ER"), fit$method),
      ,
      drop = FALSE
    ]

    rho0 <- (rho_uv + cfg$sigma2 / 2) / (1 + cfg$sigma2)

    data[[j]] <- data.frame(
      rho_uv = rho_uv,
      eta = eta,
      rho0 = rho0,
      replicate = r,
      K = K,
      method = fit$method,
      estimate = fit$pearson,
      se = fit$se,
      n_clusters = fit$n_clusters,
      n_obs = fit$n_obs,
      mean_retained_n = nrow(d) / length(unique(d$cluster_id)),
      retained_fraction = retained_fraction,
      exclusion_fraction = exclusion_fraction,
      er_exact = fit$er_exact,
      abs_error = abs(fit$pearson - rho0),
      squared_error = (fit$pearson - rho0)^2,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, data)
}

run_simulation_block <- function(job, cfg, settings, folder) {
  s <- job$setting
  b <- job$block

  path <- file.path(
    folder,
    sprintf("setting_%02d_block_%03d.rds", s, b)
  )

  if (file.exists(path)) {
    return(basename(path))
  }

  first <- (b - 1L) * cfg$block + 1L
  last <- min(cfg$reps, b * cfg$block)

  rows <- lapply(seq.int(first, last), function(r) {
    simulate_and_estimate(s, r, cfg, settings)
  })

  obj <- list(
    setting = s,
    block = b,
    first = first,
    last = last,
    M = cfg$M,
    nmax = cfg$nmax,
    er_max_pairs = cfg$er_max_pairs,
    K = cfg$K,
    data = do.call(rbind, rows)
  )

  tmp <- paste0(path, ".tmp")

  saveRDS(obj, tmp)

  if (!file.rename(tmp, path)) {
    stop("Could not finalize checkpoint: ", path)
  }

  basename(path)
}

read_simulation_results <- function(cfg) {
  root <- here::here()
  folder <- file.path(root, "results", "simulation", cfg$name, "chunks")

  settings <- simulation_settings(cfg)
  nb <- ceiling(cfg$reps / cfg$block)

  paths <- unlist(lapply(seq_len(nrow(settings)), function(s) {
    file.path(
      folder,
      sprintf("setting_%02d_block_%03d.rds", s, seq_len(nb))
    )
  }))

  missing <- paths[!file.exists(paths)]

  if (length(missing)) {
    stop(length(missing), " simulation checkpoints missing; run 05 first")
  }

  parts <- lapply(paths, readRDS)

  for (z in parts) {
    if (
      z$M != cfg$M ||
      z$nmax != cfg$nmax ||
      z$er_max_pairs != cfg$er_max_pairs ||
      !identical(z$K, cfg$K)
    ) {
      stop("Checkpoint configuration mismatch. Use a new run name.")
    }
  }

  d <- do.call(rbind, lapply(parts, `[[`, "data"))

  rownames(d) <- NULL

  expect <- nrow(settings) * cfg$reps * length(cfg$K) * 4L

  if (nrow(d) != expect) {
    stop("Unexpected simulation row count: ", nrow(d))
  }

  d
}
