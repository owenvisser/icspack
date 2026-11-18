stopifnot(requireNamespace("here", quietly = TRUE))
root <- here::here()

processed_path <- file.path(root, "data", "icspack_Tooth_Level_Data.RData")
raw_path <- file.path(root, "data", "nhanes_full.csv")
out_path <- file.path(root, "data", "nhanes_teeth.rds")

if (!file.exists(processed_path)) stop("Missing: ", processed_path)
if (!file.exists(raw_path)) stop("Missing: ", raw_path)
if (file.exists(out_path)) message("Replacing existing derived data: ", out_path)

# Source: 02_process_nhanes_dental_data.R saves object df_tooth_level.
loaded <- new.env(parent = baseenv())
load(processed_path, envir = loaded)
if (!exists("df_tooth_level", envir = loaded, inherits = FALSE) ||
  !is.data.frame(loaded$df_tooth_level)) {
  stop("icspack_Tooth_Level_Data.RData must contain the data.frame df_tooth_level")
}
tooth <- loaded$df_tooth_level
needed <- c(
  "SEQN", "ToothNumber", "Tooth_Count", "maxCAL", "maxPD",
  "DS_score", "DSI_score", "FS_score", "FSI_score",
  "WTMEC4YR", "SDMVSTRA", "SDMVPSU"
)
missing <- setdiff(needed, names(tooth))
if (length(missing)) stop("icspack tooth data missing: ", paste(missing, collapse = ", "))
if (anyDuplicated(tooth[c("SEQN", "ToothNumber")])) {
  stop("Duplicate SEQN/ToothNumber in icspack tooth data")
}

# Source: 01_build_nhanes_data.R saves merged demographics and tooth status.
raw <- utils::read.csv(raw_path, check.names = FALSE)
nums <- c(2:15, 18:31) # no third molars
status_cols <- paste0("OHX", sprintf("%02d", nums), "TC")
needed_raw <- c("SEQN", "RIDAGEYR", "OHDEXSTS", status_cols)
missing_raw <- setdiff(needed_raw, names(raw))
if (length(missing_raw)) {
  stop(
    "nhanes_full.csv missing: ",
    paste(missing_raw, collapse = ", ")
  )
}
if (anyDuplicated(raw$SEQN)) stop("Repeated SEQN in nhanes_full.csv")
raw <- raw[
  !is.na(raw$RIDAGEYR) & raw$RIDAGEYR >= 30 &
    !is.na(raw$OHDEXSTS) & raw$OHDEXSTS == 1,
  needed_raw,
  drop = FALSE
]

# Retain only teeth recorded as present in the original tooth-count rule
# (status 1, 2 or 3). Do not treat the 28 per-person candidate positions
# as 28 observed teeth; the icspack periodontal reshaping includes positions
# with missing teeth, which can otherwise distort patient eligibility.
present <- do.call(rbind, lapply(seq_along(nums), function(j) {
  idx <- !is.na(raw[[status_cols[j]]]) & raw[[status_cols[j]]] %in% c(1, 2, 3)
  data.frame(
    SEQN = raw$SEQN[idx], ToothNumber = nums[j],
    RIDAGEYR = raw$RIDAGEYR[idx]
  )
}))
if (!nrow(present)) stop("No eligible present teeth in nhanes_full.csv")
if (anyDuplicated(present[c("SEQN", "ToothNumber")])) {
  stop("Duplicate present tooth in raw data")
}
joined <- merge(tooth, present,
  by = c("SEQN", "ToothNumber"), all = FALSE,
  sort = FALSE
)
joined <- joined[!is.na(joined$Tooth_Count) & joined$Tooth_Count >= 10, ,
  drop = FALSE
]
if (!nrow(joined)) stop("No eligible teeth after merging and >=10 teeth rule")

# Use the continuous surface proportions (_score), NOT the binary DS/FS
# columns in icspack. DFS and DFSI follow the current manuscript's averaged
# normalized component convention rather than the original additive score.
num <- function(x) as.numeric(as.character(x))
tooth_level_data <- data.frame(
  SEQN = joined$SEQN,
  tooth_number = as.integer(joined$ToothNumber),
  RIDAGEYR = num(joined$RIDAGEYR),
  CAL = num(joined$maxCAL), PD = num(joined$maxPD),
  DS = num(joined$DS_score), DSI = num(joined$DSI_score),
  FS = num(joined$FS_score), FSI = num(joined$FSI_score),
  WTMEC4YR = num(joined$WTMEC4YR),
  SDMVSTRA = joined$SDMVSTRA, SDMVPSU = joined$SDMVPSU
)
tooth_level_data$DFS <- (tooth_level_data$DS + tooth_level_data$FS) / 2
tooth_level_data$DFSI <- (tooth_level_data$DSI + tooth_level_data$FSI) / 2
if (anyDuplicated(tooth_level_data[c("SEQN", "tooth_number")])) {
  stop("Duplicate tooth after conversion")
}
if (anyNA(tooth_level_data$RIDAGEYR) ||
  any(!is.finite(tooth_level_data$WTMEC4YR)) ||
  any(tooth_level_data$WTMEC4YR <= 0) ||
  anyNA(tooth_level_data$SDMVSTRA) ||
  anyNA(tooth_level_data$SDMVPSU)) {
  stop("Missing or invalid age, survey weights, strata, or PSUs")
}

# Retain eligible participants with >=10 present teeth, and check the
# actual data rather than assuming the source's prefilter is consistent.
counts <- table(tooth_level_data$SEQN)
if (any(counts < 10L)) stop("Some included participants have fewer than 10 teeth")
print(data.frame(
  patients = length(counts), teeth = nrow(tooth_level_data),
  complete_CAL_DFS = sum(stats::complete.cases(tooth_level_data[c("CAL", "DFS")])),
  complete_PD_DFS = sum(stats::complete.cases(tooth_level_data[c("PD", "DFS")]))
))


# The final analysis uses this exact simplified numeric data frame.
saveRDS(tooth_level_data, file = out_path, compress = "xz")
message("Saved processed NHANES tooth-level data: ", out_path)
