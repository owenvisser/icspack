# Shared reproducible paths, labels, and strictly formatted LaTeX cells.
repo_root <- function() {
  if (!requireNamespace("here", quietly = TRUE)) stop("Install package: here")
  here::here()
}
load_teeth <- function() {
  path <- file.path(repo_root(), "data", "nhanes_teeth.rds")
  if (!file.exists(path)) stop("Run scripts/01, 02, and 03 first: ", path)
  d <- readRDS(path)
  needed <- c(
    "SEQN", "RIDAGEYR", "CAL", "PD", "DS", "DSI", "FS",
    "FSI", "DFS", "DFSI", "WTMEC4YR"
  )
  if (!all(needed %in% names(d))) {
    stop(
      "Missing processed tooth columns: ",
      paste(setdiff(needed, names(d)), collapse = ", ")
    )
  }
  d
}
age_groups <- function(d) {
  age <- as.numeric(d$RIDAGEYR)
  list(
    "All" = d,
    "30-49" = d[age >= 30 & age < 50, , drop = FALSE],
    "50-64" = d[age >= 50 & age < 65, , drop = FALSE],
    "65+" = d[age >= 65, , drop = FALSE]
  )
}
fmt_num <- function(x, digits = 3L) sprintf(paste0("%.", digits, "f"), x)
fmt_est_se <- function(est, se, digits = 3L) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), est, se)
}
write_tex <- function(lines, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, path, useBytes = TRUE)
  message("Wrote ", path)
}
