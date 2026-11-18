
# Manuscript NHANES correlations: clustPear only; delta-method sandwich SE.
# No PSU jackknife, bootstrap, or alternative correlation implementation.

stopifnot(requireNamespace("clustPear", quietly = TRUE))

if (utils::packageVersion("clustPear") < "0.1.2") {
  stop("clustPear 0.1.2 or newer is required for survey_weights")
}

source(here::here("R", "common.R"))

d <- load_teeth()
groups <- age_groups(d)

periodontal <- c("CAL", "PD")
caries <- c("DS", "DSI", "FS", "FSI", "DFS", "DFSI")

K <- 20L
max_pairs <- 100000L

out <- vector(
  "list",
  length(groups) * length(periodontal) * length(caries)
)

pos <- 0L

for (age in names(groups)) {
  g <- groups[[age]]

  for (x in periodontal) {
    for (y in caries) {
      pos <- pos + 1L

      message(sprintf(
        "[%d/%d] %s %s/%s",
        pos, length(out), age, x, y
      ))

      cols <- c("SEQN", x, y, "WTMEC4YR")

      z <- g[
        stats::complete.cases(g[cols]) &
          is.finite(as.numeric(g[[x]])) &
          is.finite(as.numeric(g[[y]])) &
          is.finite(g$WTMEC4YR) &
          g$WTMEC4YR > 0,
        ,
        drop = FALSE
      ]

      # Match the paper's eligible-cluster convention after pairwise deletion.
      n <- table(z$SEQN)

      z <- z[
        as.character(z$SEQN) %in% names(n)[n >= 2L],
        ,
        drop = FALSE
      ]

      if (!nrow(z)) {
        stop("No eligible complete pairs for ", age, " ", x, "/", y)
      }

      fit <- clustPear::cluster_pearson(
        data = z,
        cluster = "SEQN",
        x = x,
        y = y,
        survey_weights = "WTMEC4YR",
        K = K,
        er_max_pairs = max_pairs,
        method = "all"
      )

      # SE is clustPear's independent-patient sandwich delta-method SE;
      # this is not a design-based NHANES stratum/PSU variance estimate.
      out[[pos]] <- data.frame(
        age_group = age,
        periodontal = x,
        caries = y,
        method = fit$method,
        estimate = fit$pearson,
        standard_error = fit$se,
        n_clusters = fit$n_clusters,
        n_obs = fit$n_obs,
        er_exact = fit$er_exact,
        stringsAsFactors = FALSE
      )
    }
  }
}

result <- do.call(rbind, out)

if (nrow(result) != 192L) {
  stop("Expected 4 age groups x 12 pairs x 4 methods")
}

path <- file.path(
  repo_root(),
  "results",
  "nhanes",
  "correlations_delta.csv"
)

dir.create(
  dirname(path),
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.csv(result, path, row.names = FALSE)

if (any(result$method == "ER" & result$er_exact %in% FALSE)) {
  warning(
    "At least one ER calculation was capped; inspect er_exact before calling it exact"
  )
}

message("Saved only delta-method NHANES association results: ", path)
