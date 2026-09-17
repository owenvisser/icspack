library(gtsummary)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(gt)
library(here)

# IMPORT DATA

var.info <- read.csv(here("data", "variable_info.csv"))
peri.uid <- read.csv(here("data", "nhanes_full.csv")) # only people with some periodontal data

# We want only the 2:15 & 18:31 teeth, ignoring wisdom teeth
ohx.nums <- paste0("OHX", sprintf("%02d", c(2:15, 18:31)))

# Periodontal measures: 
# CJ = free gingival margin to cementoenamel junction (CEJ),
# LA = loss of attachment, 
# PC = pocket depth;
# site letters are A = mesio-lingual, D = distal, L = mid-lingual,
#                  M = mid-facial, P = disto-lingual, S = mesio-facial
infovec <- c(
  "CJA", "CJD", "CJL", "CJM", "CJP", "CJS",
  "LAA", "LAD", "LAL", "LAM", "LAP", "LAS",
  "PCA", "PCD", "PCL", "PCM", "PCP", "PCS",
  "TC", "CSC", "CTC"
)

# A matrix of all tooth codes, dental caries surface codes/conditions, and periodontal measurements
teeth_list <- as.data.frame(
  t(
    sapply(ohx.nums, function(x) {
      sapply(infovec, function(y) paste0(x, y))
    })
  )
)

# FILTER

peri.uid <- peri.uid %>%
  # Retain adults age 30 and older
  filter(DMDHRAGE >= 30) %>%
  # Retain participants with a complete oral health examination
  filter(OHDEXSTS == 1) %>%
  rowwise() %>%
  mutate(
    # Tooth status:  1 = primary tooth, 2 = permanent tooth, 3 = implant,
    #                4 = tooth not present, 5 = permanent dental root fragment
    # Count teeth present as primary teeth, permanent teeth, or implants - Bruce/Levy
    Tooth_Count = sum(c_across(all_of(teeth_list$TC)) %in% c(1, 2, 3)),
    # Count teeth missing or represented only by a permanent dental root fragment - Bruce/Levy
    Teeth_Missing = sum(c_across(all_of(teeth_list$TC)) %in% c(4, 5))
  ) %>%
  ungroup() %>%
  # Retain participants with at least 9 teeth present - Requested by Bruce/Levy
  filter(Tooth_Count >= 9)

# additionally filter by age group
agedf <- peri.uid %>%
  mutate(
    # Recode exact age into 30-49, 50-64, and 65+ age groups
    Age = cut(
      DMDHRAGE,
      breaks = c(29, 49, 64, 80),
      labels = c("30-49", "50-64", "65+"),
      include.lowest = TRUE
    )
  ) %>%
  select(SEQN, Age)

# CATEGORIZE PERIODONTAL STATUS

# Calculate the mean after excluding missing values and NHANES code 99
avg_exclude <- function(x) {
  valid_values <- x[!is.na(x) & x != 99]

  if (length(valid_values) > 0) {
    mean(valid_values)
  } else {
    NA_real_
  }
}

# Find the maximum after excluding missing values and NHANES code 99
max_exclude <- function(x) {
  valid_values <- x[!is.na(x) & x != 99]

  if (length(valid_values) > 0) {
    max(valid_values)
  } else {
    NA_real_
  }
}

# Count how many valid observations equal the maximum observed value
max_exclude_count <- function(x) {
  valid_values <- x[!is.na(x) & x != 99]

  if (length(valid_values) > 0) {
    sum(valid_values == max(valid_values))
  } else {
    NA_real_
  }
}

# Indicate whether any valid measurement is greater than a specified threshold
get_over_n <- function(x, n) {
  valid_values <- x[!is.na(x) & x != 99]

  if (length(valid_values) > 0) {
    as.integer(any(valid_values > n))
  } else {
    NA_real_
  }
}


# CREATE LONG-FORM PERIODONTAL DATA

# Get the CAL, pocket depth, and FGM-to-CEJ variable names for all included teeth
perio_vars <- unlist(teeth_list %>% select(-TC, -CTC, -CSC), use.names = FALSE)

Perio_long <- peri.uid %>%
  # Retain participant ID and periodontal measurements for all included teeth
  select(
    SEQN,
    all_of(perio_vars)
  ) %>%
  # Pivot periodontal measurements from wide to long format
  pivot_longer(
    cols = all_of(perio_vars),
    names_to = "ToothCodename",
    values_to = "Measurement"
  ) %>%
  # Separate NHANES variable names into tooth number, measurement type, and measurement location
  mutate(
    Measurement = as.numeric(Measurement),
    ToothNumber = substr(ToothCodename, 4, 5),
    MeasurementCode = substr(ToothCodename, 6, 7),
    MeasurementLoc = substr(ToothCodename, 8, 8)
  ) %>%
  # Convert NHANES measurement codes into interpretable labels
  mutate(
    MeasurementType = case_when(
      MeasurementCode == "CJ" ~ "FGM to CEJ measurement (mm)",
      MeasurementCode == "LA" ~ "Clinical Loss of Attachment (mm)",
      MeasurementCode == "PC" ~ "Pocket Depth (mm)"
    ),
    MeasurementFace = case_when(
      MeasurementLoc == "A" ~ "mesio-lingual",
      MeasurementLoc == "D" ~ "distal",
      MeasurementLoc == "L" ~ "mid-lingual",
      MeasurementLoc == "M" ~ "mid-facial",
      MeasurementLoc == "P" ~ "distal-lingual",
      MeasurementLoc == "S" ~ "mesio-facial"
    )
  ) %>%
  # Remove FGM-to-CEJ measurements because CAL and pocket depth are used for periodontal classification
  filter(MeasurementType != "FGM to CEJ measurement (mm)") %>%
  # Retain variables needed for periodontal analyses
  select(SEQN, ToothNumber, MeasurementType, MeasurementFace, Measurement)


# CREATE TOOTH-REGION GROUPS

Perio_tooth <- Perio_long %>%
  mutate(
    # Group teeth into anatomical regions
    Tooth_type = case_when(
      ToothNumber %in% c("02", "03", "14", "15") ~ "Upper Molars",
      ToothNumber %in% c("04", "05", "12", "13") ~ "Upper Premolars",
      ToothNumber %in% c("06", "07", "08", "09", "10", "11") ~ "Upper Anterior",
      ToothNumber %in% c("18", "19", "30", "31") ~ "Lower Molars",
      ToothNumber %in% c("20", "21", "28", "29") ~ "Lower Premolars",
      ToothNumber %in% c("22", "23", "24", "25", "26", "27") ~ "Lower Anterior",
      TRUE ~ NA_character_
    )
  )

# Save long-format periodontal data for later analyses
save(
  Perio_tooth,
  file = here("data", "LongData_Filtered.RData")
)


# CREATE TOOTH-LEVEL CONTINUOUS MEASURES

Perio_tooth_level_allMeasures <- Perio_tooth %>%
  group_by(SEQN, ToothNumber, Tooth_type) %>%
  summarize(
    maxCAL = max_exclude(
      Measurement[MeasurementType == "Clinical Loss of Attachment (mm)"]
    ),
    avgCAL = avg_exclude(
      Measurement[MeasurementType == "Clinical Loss of Attachment (mm)"]
    ),
    maxPD = max_exclude(
      Measurement[MeasurementType == "Pocket Depth (mm)"]
    ),
    avgPD = avg_exclude(
      Measurement[MeasurementType == "Pocket Depth (mm)"]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ToothNumber = as.numeric(ToothNumber)
  )


# CREATE TOOTH-LEVEL BINARY THRESHOLD MEASURES

Perio_tooth_level <- Perio_tooth %>%
  group_by(SEQN, ToothNumber, Tooth_type) %>%
  summarize(
    `CAL >= 3mm` = get_over_n(
      Measurement[MeasurementType == "Clinical Loss of Attachment (mm)"],
      2
    ),
    `CAL >= 4mm` = get_over_n(
      Measurement[MeasurementType == "Clinical Loss of Attachment (mm)"],
      3
    ),
    `CAL >= 5mm` = get_over_n(
      Measurement[MeasurementType == "Clinical Loss of Attachment (mm)"],
      4
    ),
    `PD >= 4mm` = get_over_n(
      Measurement[MeasurementType == "Pocket Depth (mm)"],
      3
    ),
    `PD >= 5mm` = get_over_n(
      Measurement[MeasurementType == "Pocket Depth (mm)"],
      4
    ),
    `PD >= 6mm` = get_over_n(
      Measurement[MeasurementType == "Pocket Depth (mm)"],
      5
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ToothNumber = as.numeric(ToothNumber)
  )


# CREATE PERSON-LEVEL PERIODONTAL STATUS

Perio_long <- Perio_long %>%
  # Summarize each tooth separately for interproximal sites and all measured sites
  group_by(SEQN, ToothNumber, MeasurementType) %>%
  summarize(
    # Highest severity score among interproximal measurement locations
    ToothSeverity_Mesial_Distal = max_exclude(
      Measurement[
        MeasurementFace %in% c(
          "mesio-lingual",
          "distal",
          "distal-lingual",
          "mesio-facial"
        )
      ]
    ),

    # Number of interproximal sites equal to the maximum measurement on the tooth
    TSNum_Max_Mesial_Distal = max_exclude_count(
      Measurement[
        MeasurementFace %in% c(
          "mesio-lingual",
          "distal",
          "distal-lingual",
          "mesio-facial"
        )
      ]
    ),

    # Highest severity score across all measured sites on the tooth
    ToothSeverity_Total = max_exclude(Measurement),

    .groups = "drop_last"
  ) %>%
  group_by(SEQN) %>%
  summarize(
    # Severe periodontitis requires CAL >= 6 mm on at least two teeth and PD >= 5 mm on at least one tooth
    Severe_PD = ifelse(
      sum(
        ToothSeverity_Mesial_Distal[
          MeasurementType == "Clinical Loss of Attachment (mm)"
        ] %in% 6:12,
        na.rm = TRUE
      ) >= 2 &
        sum(
          ToothSeverity_Mesial_Distal[
            MeasurementType == "Pocket Depth (mm)"
          ] %in% 5:12,
          na.rm = TRUE
        ) >= 1,
      1,
      0
    ),

    # Moderate periodontitis requires CAL >= 4 mm on at least two teeth or PD >= 5 mm on at least two teeth
    Moderate_PD = ifelse(
      sum(
        ToothSeverity_Mesial_Distal[
          MeasurementType == "Clinical Loss of Attachment (mm)"
        ] %in% 4:12,
        na.rm = TRUE
      ) >= 2 |
        sum(
          ToothSeverity_Mesial_Distal[
            MeasurementType == "Pocket Depth (mm)"
          ] %in% 5:12,
          na.rm = TRUE
        ) >= 2,
      1,
      0
    ),

    # Mild periodontitis requires CAL >= 3 mm at two interproximal sites and PD meeting the mild disease threshold
    Mild_PD = ifelse(
      all(
        any(
          # CAL >= 3 mm on at least two different teeth
          sum(
            ToothSeverity_Mesial_Distal[
              MeasurementType == "Clinical Loss of Attachment (mm)"
            ] %in% 3:12,
            na.rm = TRUE
          ) >= 2,

          # Or two interproximal CAL measurements at the maximum value on the same tooth
          TSNum_Max_Mesial_Distal[
            ToothSeverity_Mesial_Distal[
              MeasurementType == "Clinical Loss of Attachment (mm)"
            ] %in% 3:12
          ] >= 2,

          na.rm = TRUE
        ),

        any(
          # PD >= 4 mm on at least two different teeth
          sum(
            ToothSeverity_Mesial_Distal[
              MeasurementType == "Pocket Depth (mm)"
            ] %in% 4:12,
            na.rm = TRUE
          ) >= 2,

          # Or PD >= 5 mm at any measured site
          sum(
            ToothSeverity_Total[
              MeasurementType == "Pocket Depth (mm)"
            ] %in% 5:12,
            na.rm = TRUE
          ) >= 1,

          na.rm = TRUE
        )
      ),
      1,
      0
    ),

    .groups = "drop"
  ) %>%
  # Assign each participant to the highest periodontal severity category only
  mutate(
    Moderate_PD = ifelse(Severe_PD == 1, 0, Moderate_PD),
    Mild_PD = ifelse(Severe_PD == 1 | Moderate_PD == 1, 0, Mild_PD),
    No_PD = ifelse(Mild_PD == 1 | Moderate_PD == 1 | Severe_PD == 1, 0, 1)
  )

# CATEGORIZE DENTAL CARIES

# Get the tooth-level caries variables for all included teeth
caries_vars <- unlist(teeth_list %>% select(CTC, CSC), use.names = FALSE)

Caries_long <- peri.uid %>%
  # Retain participant ID and caries variables for all included teeth
  select(
    SEQN,
    all_of(caries_vars)
  ) %>%
  # Pivot tooth-level caries variables from wide to long format
  pivot_longer(
    cols = all_of(caries_vars),
    names_to = "ToothCodename",
    values_to = "Measurement",
    values_transform = as.character
  ) %>%
  # Separate NHANES variable names into tooth number and measurement code
  mutate(
    ToothNumber = as.numeric(substr(ToothCodename, 4, 5)),
    MeasurementCode = substr(ToothCodename, 6, 8)
  ) %>%
  select(SEQN, ToothNumber, MeasurementCode, Measurement) %>%
  pivot_wider(
    names_from = MeasurementCode,
    values_from = Measurement
  ) %>%
  # Identify individual decayed and restored surfaces from the CSC surface code
  mutate(
    `0` = ifelse(CTC %in% c("Z", "K") & grepl("0", CSC), 1, 0),
    `1` = ifelse(CTC %in% c("Z", "K") & grepl("1", CSC), 1, 0),
    `2` = ifelse(CTC %in% c("Z", "K") & grepl("2", CSC), 1, 0),
    `3` = ifelse(CTC %in% c("Z", "K") & grepl("3", CSC), 1, 0),
    `4` = ifelse(CTC %in% c("Z", "K") & grepl("4", CSC), 1, 0),
    `5` = ifelse(CTC %in% c("Z", "K") & grepl("5", CSC), 1, 0),
    `6` = ifelse(CTC %in% c("Z", "K") & grepl("6", CSC), 1, 0),
    `7` = ifelse(CTC %in% c("Z", "K") & grepl("7", CSC), 1, 0),
    `8` = ifelse(CTC %in% c("Z", "K") & grepl("8", CSC), 1, 0),
    `9` = ifelse(CTC %in% c("Z", "K") & grepl("9", CSC), 1, 0)
  ) %>%
  select(-CSC) %>%
  # Calculate decayed, filled, and missing surface counts for each tooth
  mutate(
    DS_interprox = case_when(
      CTC %in% c("Z", "K") ~ `3` + `4`,
      CTC == "J" ~ 2,
      TRUE ~ 0
    ),
    DS = case_when(
      CTC %in% c("Z", "K") ~ `0` + `1` + `2` + `3` + `4`,
      CTC == "J" & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == "J" ~ 5,
      TRUE ~ 0
    ),
    FS_interprox = case_when(
      CTC %in% c("Z", "K") ~ `8` + `9`,
      CTC == "T" ~ 2,
      TRUE ~ 0
    ),
    FS = case_when(
      CTC %in% c("Z", "K") ~ `5` + `6` + `7` + `8` + `9`,
      CTC == "T" & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == "T" ~ 5,
      TRUE ~ 0
    ),
    MS_interprox = case_when(
      CTC %in% c("E", "R", "P") ~ 2,
      TRUE ~ 0
    ),
    MS = case_when(
      CTC %in% c("E", "R", "P") & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC %in% c("E", "R", "P") ~ 5,
      TRUE ~ 0
    )
  ) %>%
  # Create tooth-level surface counts for decay, restorations, and caries experience
  mutate(
    DS_Tooth = DS,
    DS_Tooth_Interprox = DS_interprox,
    FS_Tooth = FS,
    FS_Tooth_Interprox = FS_interprox,
    DFS_Tooth = DS + FS,
    DFS_Tooth_Interprox = DS_interprox + FS_interprox,
    DMFS_Tooth = DS + FS + MS,
    DMFS_Tooth_Interprox = DS_interprox + FS_interprox + MS_interprox
  )


# CREATE TOOTH-LEVEL CARIES DATA

Caries_tooth <- Caries_long %>%
  mutate(
    # Anterior teeth have four measured surfaces and posterior teeth have five
    Num_Surfaces = ifelse(
      ToothNumber %in% c(6:11, 22:27),
      4,
      5
    ),

    # Every tooth has two interproximal surfaces
    Num_Interprox_Surfaces = 2,

    # Continuous scores represent the proportion of tooth surfaces affected
    DS_score = DS_Tooth / Num_Surfaces,
    DSI_score = DS_Tooth_Interprox / Num_Interprox_Surfaces,
    FS_score = FS_Tooth / Num_Surfaces,
    FSI_score = FS_Tooth_Interprox / Num_Interprox_Surfaces,
    DFS_score = DFS_Tooth / Num_Surfaces,
    DFSI_score = DFS_Tooth_Interprox / Num_Interprox_Surfaces,
    DMFS_score = DMFS_Tooth / Num_Surfaces,
    DMFSI_score = DMFS_Tooth_Interprox / Num_Interprox_Surfaces,

    # Group teeth into anatomical regions
    Tooth_type = case_when(
      ToothNumber %in% c(2, 3, 14, 15) ~ "Upper Molars",
      ToothNumber %in% c(4, 5, 12, 13) ~ "Upper Premolars",
      ToothNumber %in% c(6:11) ~ "Upper Anterior",
      ToothNumber %in% c(18, 19, 30, 31) ~ "Lower Molars",
      ToothNumber %in% c(20, 21, 28, 29) ~ "Lower Premolars",
      ToothNumber %in% c(22:27) ~ "Lower Anterior",
      TRUE ~ NA_character_
    )
  )


# CREATE TOOTH-LEVEL BINARY CARIES INDICATORS

Caries_tooth_level <- Caries_tooth %>%
  mutate(
    # Binary indicators identify whether any affected surface is present on the tooth
    DS = ifelse(DS_score > 0, 1, 0),
    DSI = ifelse(DSI_score > 0, 1, 0),
    FS = ifelse(FS_score > 0, 1, 0),
    FSI = ifelse(FSI_score > 0, 1, 0),
    DFS = ifelse(DFS_score > 0, 1, 0),
    DFSI = ifelse(DFSI_score > 0, 1, 0),
    DMFS = ifelse(DMFS_score > 0, 1, 0),
    DMFSI = ifelse(DMFSI_score > 0, 1, 0)
  ) %>%
  select(
    SEQN,
    ToothNumber,
    Tooth_type,
    DS_score,
    DSI_score,
    FS_score,
    FSI_score,
    DFS_score,
    DFSI_score,
    DMFS_score,
    DMFSI_score,
    DS,
    DSI,
    FS,
    FSI,
    DFS,
    DFSI,
    DMFS,
    DMFSI
  )


# CREATE PERSON-LEVEL CARIES SCORES

Caries_long <- Caries_long %>%
  group_by(SEQN) %>%
  summarise(
    FS_score = sum(FS_Tooth, na.rm = TRUE),
    FSI_score = sum(FS_Tooth_Interprox, na.rm = TRUE),
    DS_score = sum(DS_Tooth, na.rm = TRUE),
    DSI_score = sum(DS_Tooth_Interprox, na.rm = TRUE),
    DFS_score = sum(DFS_Tooth, na.rm = TRUE),
    DFSI_score = sum(DFS_Tooth_Interprox, na.rm = TRUE),
    DMFS_score = sum(DMFS_Tooth, na.rm = TRUE),
    DMFSI_score = sum(DMFS_Tooth_Interprox, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Binary person-level indicators identify whether any affected surface is present
    `Surface Decay or Restoration (DFS)` = ifelse(DFS_score > 0, 1, 0),
    `Surface Decay or Restoration on Interproximal Surfaces (DFSI)` = ifelse(DFSI_score > 0, 1, 0),
    `Untreated Interproximal Dental Caries (DSI)` = ifelse(DSI_score > 0, 1, 0),
    `Untreated Dental Caries (DS)` = ifelse(DS_score > 0, 1, 0),
    `Treated Interproximal Dental Caries (FSI)` = ifelse(FSI_score > 0, 1, 0),
    `Treated Dental Caries (FS)` = ifelse(FS_score > 0, 1, 0),
    `Caries Experience (DMFS)` = ifelse(DMFS_score > 0, 1, 0),
    `Interproximal Caries Experience (DMFSI)` = ifelse(DMFSI_score > 0, 1, 0)
  )


# CREATE DEMOGRAPHIC COVARIATES

# NHANES prefixes: OHX = oral health, SMQ = smoking questionnaire, DEMO = demographics
df <- peri.uid %>%
  mutate(
    Gender = factor(
      case_when(
        RIAGENDR == 1 ~ "Male",
        RIAGENDR == 2 ~ "Female",
        TRUE ~ NA_character_
      )
    ),

    Race = factor(
      case_when(
        RIDRETH3 == 1 ~ "Mexican American",
        RIDRETH3 == 2 ~ "Other Hispanic",
        RIDRETH3 == 3 ~ "Non-Hispanic White",
        RIDRETH3 == 4 ~ "Non-Hispanic Black",
        RIDRETH3 == 6 ~ "Non-Hispanic Asian",
        RIDRETH3 == 7 ~ "Other Race - Including Multi-Racial",
        TRUE ~ NA_character_
      )
    ),

    Age = cut(
      DMDHRAGE,
      breaks = c(29, 49, 64, 80),
      labels = c("30-49", "50-64", "65+"),
      include.lowest = TRUE
    ),

    Income_Household = factor(
      case_when(
        INDHHIN2 == 1 ~ "$0 to $4,999",
        INDHHIN2 == 2 ~ "$5,000 to $9,999",
        INDHHIN2 == 3 ~ "$10,000 to $14,999",
        INDHHIN2 == 4 ~ "$15,000 to $19,999",
        INDHHIN2 == 5 ~ "$20,000 to $24,999",
        INDHHIN2 == 6 ~ "$25,000 to $34,999",
        INDHHIN2 == 7 ~ "$35,000 to $44,999",
        INDHHIN2 == 8 ~ "$45,000 to $54,999",
        INDHHIN2 == 9 ~ "$55,000 to $64,999",
        INDHHIN2 == 10 ~ "$65,000 to $74,999",
        INDHHIN2 == 12 ~ "$20,000 and Over",
        INDHHIN2 == 13 ~ "Under $20,000",
        INDHHIN2 == 14 ~ "$75,000 to $99,999",
        INDHHIN2 == 15 ~ "$100,000 and Over",
        INDHHIN2 == 77 ~ "Refused",
        INDHHIN2 == 99 ~ "Do not know",
        TRUE ~ NA_character_
      )
    ),

    Income_Family = factor(
      case_when(
        INDFMIN2 == 1 ~ "$0 to $4,999",
        INDFMIN2 == 2 ~ "$5,000 to $9,999",
        INDFMIN2 == 3 ~ "$10,000 to $14,999",
        INDFMIN2 == 4 ~ "$15,000 to $19,999",
        INDFMIN2 == 5 ~ "$20,000 to $24,999",
        INDFMIN2 == 6 ~ "$25,000 to $34,999",
        INDFMIN2 == 7 ~ "$35,000 to $44,999",
        INDFMIN2 == 8 ~ "$45,000 to $54,999",
        INDFMIN2 == 9 ~ "$55,000 to $64,999",
        INDFMIN2 == 10 ~ "$65,000 to $74,999",
        INDFMIN2 == 12 ~ "$20,000 and Over",
        INDFMIN2 == 13 ~ "Under $20,000",
        INDFMIN2 == 14 ~ "$75,000 to $99,999",
        INDFMIN2 == 15 ~ "$100,000 and Over",
        INDFMIN2 == 77 ~ "Refused",
        INDFMIN2 == 99 ~ "Do not know",
        TRUE ~ NA_character_
      )
    ),

    Smoking_Status = factor(
      case_when(
        SMQ040 == 1 ~ "Every day",
        SMQ040 == 2 ~ "Some days",
        SMQ040 == 3 ~ "Not at all",
        SMQ040 == 7 ~ "Refused",
        SMQ040 == 9 ~ "Do not know",
        TRUE ~ NA_character_
      )
    ),

    Missing_Teeth = factor(
      case_when(
        Teeth_Missing == 0 ~ "0",
        Teeth_Missing %in% 1:5 ~ "1-5",
        Teeth_Missing %in% 6:28 ~ "6-28",
        TRUE ~ NA_character_
      )
    )
  ) %>%
  select(
    SEQN,
    Age,
    Gender,
    Race,
    Income_Family,
    Income_Household,
    Smoking_Status,
    Missing_Teeth,
    Teeth_Missing,
    Tooth_Count
  )


# CREATE PERSON-LEVEL DATA

df_person_level <- df %>%
  left_join(Perio_long, by = "SEQN") %>%
  left_join(Caries_long, by = "SEQN")


# CREATE TOOTH-LEVEL DATA

df_tooth_level <- Perio_tooth_level_allMeasures %>%
  full_join(
    Perio_tooth_level,
    by = c("SEQN", "ToothNumber", "Tooth_type")
  ) %>%
  full_join(
    Caries_tooth_level,
    by = c("SEQN", "ToothNumber", "Tooth_type")
  ) %>%
  left_join(
    df %>% select(SEQN, Tooth_Count),
    by = "SEQN"
  )


# SAVE PROCESSED DATA

save(
  df_person_level,
  file = here("data", "Person_Level_Data.RData")
)

save(
  df_tooth_level,
  file = here("data", "Tooth_Level_Data.RData")
)

