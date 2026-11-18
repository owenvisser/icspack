# package requirements
stopifnot(requireNamespace("here", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("foreign", quietly = TRUE))

root <- here::here()

read_cycle <- function(tag, folder) {
  components <- c(
    DEMO = "DEMO",
    OHXPER = "OHXPER",
    OHXDEN = "OHXDEN",
    SMQ = "SMQ",
    SMQRTU = "SMQRTU"
  )
  files <- file.path(
    root,
    "data-raw",
    folder,
    paste0(components, "_", tag, ".XPT")
  )

  missing <- files[!file.exists(files)]

  if (length(missing)) {
    stop("Missing raw NHANES files: ", paste(missing, collapse = ", "))
  }

  out <- lapply(files, foreign::read.xport)

  names(out) <- names(components)

  for (name in names(out)) {
    if (anyDuplicated(out[[name]]$SEQN)) {
      stop("Duplicate SEQN in ", name, "_", tag)
    }
  }
  out
}

a <- read_cycle("G", "2011_2012")
b <- read_cycle("H", "2013_2014")

peri <- dplyr::bind_rows(a$OHXPER, b$OHXPER)
dent <- dplyr::bind_rows(a$OHXDEN, b$OHXDEN)
shared_demo <- intersect(names(a$DEMO), names(b$DEMO))
demo <- dplyr::bind_rows(a$DEMO[shared_demo], b$DEMO[shared_demo])
smq <- dplyr::bind_rows(a$SMQ, b$SMQ)
smqrtu <- dplyr::bind_rows(a$SMQRTU, b$SMQRTU)

full <- dplyr::full_join(peri, dent, by = "SEQN") |>
  dplyr::full_join(demo, by = "SEQN") |>
  dplyr::full_join(smq, by = "SEQN") |>
  dplyr::full_join(smqrtu, by = "SEQN")

if (all(c("OHDEXSTS.x", "OHDEXSTS.y") %in% names(full))) {
  full$OHDEXSTS <- full$OHDEXSTS.x
  full[c("OHDEXSTS.x", "OHDEXSTS.y")] <- NULL
}

if (anyDuplicated(full$SEQN)) stop("Merged NHANES contains duplicated SEQN")
needed <- c("SEQN", "RIDAGEYR", "OHDEXSTS", "WTMEC2YR", "SDMVSTRA", "SDMVPSU")

if (!all(needed %in% names(full))) {
  stop("Missing merged demographic/design variables: ", paste(setdiff(needed, names(full)), collapse = ", "))
}
dir.create(file.path(root, "data"), showWarnings = FALSE, recursive = TRUE)

utils::write.csv(full, file.path(root, "data", "nhanes_full.csv"), row.names = FALSE)
utils::write.csv(peri, file.path(root, "data", "nhanes_perio.csv"), row.names = FALSE)


utils::write.csv(data.frame(Variable.Name = names(full)),
  file.path(root, "data", "variable_info.csv"),
  row.names = FALSE
)

message("Wrote merged raw NHANES data: ", nrow(full), " participants")
