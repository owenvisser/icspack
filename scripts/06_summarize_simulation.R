
# Summaries for T1, F1, S1 and their numerical supplementary tables.

stopifnot(requireNamespace("here", quietly = TRUE))

source(here::here("R", "simulation.R"))

cfg <- simulation_config()
raw <- read_simulation_results(cfg)

root <- here::here()
out_dir <- file.path(root, "results", "simulation", cfg$name)

by_scenario <- split(
  raw,
  interaction(
    raw$rho_uv,
    raw$eta,
    raw$K,
    raw$method,
    drop = TRUE
  )
)

sum_one <- function(z) {
  data.frame(
    rho_uv = z$rho_uv[1],
    eta = z$eta[1],
    rho0 = z$rho0[1],
    K = z$K[1],
    method = z$method[1],
    n_replicates = nrow(z),
    mean_retained_n = mean(z$mean_retained_n),
    retention_fraction = mean(z$retained_fraction),
    exclusion_fraction = mean(z$exclusion_fraction),
    mean_estimate = mean(z$estimate),
    mean_se = mean(z$se),
    empirical_sd = stats::sd(z$estimate),
    mae = mean(z$abs_error),
    rmse = sqrt(mean(z$squared_error)),
    er_exact_fraction = if (z$method[1] == "ER") mean(z$er_exact) else NA_real_
  )
}

scenario <- do.call(rbind, lapply(by_scenario, sum_one))

scenario <- scenario[
  order(
    scenario$eta,
    scenario$rho_uv,
    scenario$K,
    match(scenario$method, c("NW", "CW", "EQ", "ER"))
  ),
]

rownames(scenario) <- NULL

utils::write.csv(
  scenario,
  file.path(out_dir, "scenario_summary.csv"),
  row.names = FALSE
)

by_pool <- split(
  scenario,
  interaction(scenario$K, scenario$method, drop = TRUE)
)

pooled <- do.call(rbind, lapply(by_pool, function(z) {
  data.frame(
    K = z$K[1],
    method = z$method[1],
    pooled_mae = mean(z$mae),
    pooled_rmse = mean(z$rmse)
  )
}))

pooled <- pooled[
  order(
    pooled$K,
    match(pooled$method, c("NW", "CW", "EQ", "ER"))
  ),
]

rownames(pooled) <- NULL

utils::write.csv(
  pooled,
  file.path(out_dir, "pooled_MAE_by_K.csv"),
  row.names = FALSE
)

# Paired within-replicate absolute EQ/ER differences; retained for paper text.

eq <- raw[
  raw$method == "EQ",
  c("rho_uv", "eta", "K", "replicate", "estimate")
]

er <- raw[
  raw$method == "ER",
  c("rho_uv", "eta", "K", "replicate", "estimate")
]

paired <- merge(
  eq,
  er,
  by = c("rho_uv", "eta", "K", "replicate"),
  suffixes = c("_EQ", "_ER")
)

paired$abs_difference <- abs(paired$estimate_EQ - paired$estimate_ER)

agreement <- aggregate(
  abs_difference ~ rho_uv + eta + K,
  paired,
  mean
)

names(agreement)[names(agreement) == "abs_difference"] <- "mean_abs_EQ_minus_ER"

utils::write.csv(
  agreement,
  file.path(out_dir, "eq_er_agreement.csv"),
  row.names = FALSE
)

message(
  "Saved simulation summaries for ", cfg$mode, ": ",
  nrow(scenario), " scenario rows"
)
