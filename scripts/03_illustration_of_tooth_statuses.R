library(dplyr)
library(tidyr)
library(gt)
library(here)


# ============================================================
# 1. LOAD PROCESSED DATA
# ============================================================

load(here("data", "Person_Level_Data.RData"))
load(here("data", "Tooth_Level_Data.RData"))


# ============================================================
# 2. ADD AGE GROUP TO TOOTH-LEVEL DATA
# ============================================================

# Age is already contained in the person-level dataset.
# Create a one-row-per-person lookup table and join it onto
# the tooth-level dataset.

age_lookup <- df_person_level %>%
  select(SEQN, Age) %>%
  distinct()

df_tooth_level <- df_tooth_level %>%
  left_join(
    age_lookup,
    by = "SEQN"
  )


# ============================================================
# 3. HELPER FUNCTION FOR PREVALENCE
# ============================================================

# All prevalence variables used below are binary:
#
#   1 = condition is present
#   0 = condition is absent
#
# Therefore, the mean of the variable is the proportion with
# the condition. Multiplying by 100 converts this to percent.

get_percent <- function(x) {
  mean(x, na.rm = TRUE) * 100
}


# ============================================================
# 4. PERIODONTAL PREVALENCE - PERSON LEVEL
# ============================================================

# Each participant belongs to only one periodontal severity
# category because Severe_PD, Moderate_PD, and Mild_PD were
# made mutually exclusive in Script 2.
#
# Calculate prevalence within each age group.

perio_person_age <- df_person_level %>%
  group_by(Age) %>%
  summarise(
    Count = n(),
    Severe = get_percent(Severe_PD),
    Moderate = get_percent(Moderate_PD),
    Mild = get_percent(Mild_PD),
    .groups = "drop"
  )


# Calculate the same quantities for the entire sample.

perio_person_total <- df_person_level %>%
  summarise(
    Count = n(),
    Severe = get_percent(Severe_PD),
    Moderate = get_percent(Moderate_PD),
    Mild = get_percent(Mild_PD)
  ) %>%
  mutate(
    Age = "Total"
  )


# Combine age-specific and total results.

perio_person <- bind_rows(
  perio_person_age,
  perio_person_total
) %>%
  mutate(
    Age = as.character(Age)
  ) %>%
  select(
    Age,
    Count,
    Severe,
    Moderate,
    Mild
  )


# ============================================================
# 5. CREATE TOOTH-REGION VARIABLES
# ============================================================

# Tooth_type contains six categories:
#
#   Upper Anterior
#   Upper Premolars
#   Upper Molars
#   Lower Anterior
#   Lower Premolars
#   Lower Molars
#
# Separate this into:
#
#   Jaw    = Upper / Lower
#   Region = Anterior / Premolar / Molar
#
# This makes it easier to produce the U/L entries shown in
# the original table.

tooth_table_data <- df_tooth_level %>%
  mutate(
    Jaw = case_when(
      grepl("^Upper", Tooth_type) ~ "Upper",
      grepl("^Lower", Tooth_type) ~ "Lower",
      TRUE ~ NA_character_
    ),

    Region = case_when(
      grepl("Anterior$", Tooth_type) ~ "Anterior",
      grepl("Premolars$", Tooth_type) ~ "Premolar",
      grepl("Molars$", Tooth_type) ~ "Molar",
      TRUE ~ NA_character_
    )
  )


# ============================================================
# 6. PERIODONTAL PREVALENCE - TOOTH LEVEL
# ============================================================

# For the periodontal tooth-level table we want:
#
#   CAL >= 5 mm
#   PD  >= 4 mm
#
# separately for upper and lower teeth within each tooth region.
#
# For example:
#
#   Upper anterior CAL prevalence =
#
#       number of upper anterior teeth with CAL >= 5 mm
#       ------------------------------------------------ × 100
#       number of upper anterior teeth with observed CAL
#
#
# First calculate the prevalence separately for each
# Age × Region × Jaw combination.

perio_tooth_age <- tooth_table_data %>%
  filter(
    !is.na(Age),
    !is.na(Jaw),
    !is.na(Region)
  ) %>%
  group_by(
    Age,
    Region,
    Jaw
  ) %>%
  summarise(
    CAL = get_percent(`CAL >= 5mm`),
    PD = get_percent(`PD >= 4mm`),
    .groups = "drop"
  )


# Repeat for the full sample.

perio_tooth_total <- tooth_table_data %>%
  filter(
    !is.na(Jaw),
    !is.na(Region)
  ) %>%
  group_by(
    Region,
    Jaw
  ) %>%
  summarise(
    CAL = get_percent(`CAL >= 5mm`),
    PD = get_percent(`PD >= 4mm`),
    .groups = "drop"
  ) %>%
  mutate(
    Age = "Total"
  )


# Combine age-specific and total estimates.

perio_tooth <- bind_rows(
  perio_tooth_age,
  perio_tooth_total
) %>%
  mutate(
    Age = as.character(Age)
  )


# ============================================================
# 7. FORMAT PERIODONTAL TOOTH PREVALENCE AS UPPER/LOWER
# ============================================================

# The original table reports each cell as:
#
#       Upper % / Lower %
#
# For example:
#
#       3.8 / 5.7
#
# Pivot Upper and Lower into separate columns first.

perio_tooth_wide <- perio_tooth %>%
  pivot_wider(
    names_from = Jaw,
    values_from = c(CAL, PD)
  )


# Combine Upper and Lower percentages into one formatted value.

perio_tooth_wide <- perio_tooth_wide %>%
  mutate(
    CAL_UL = sprintf("%.1f/%.1f", CAL_Upper, CAL_Lower),
    PD_UL  = sprintf("%.1f/%.1f", PD_Upper, PD_Lower)
  ) %>%
  select(
    Age,
    Region,
    CAL_UL,
    PD_UL
  )


# Put Anterior, Premolar, and Molar into separate columns.

perio_CAL <- perio_tooth_wide %>%
  select(
    Age,
    Region,
    CAL_UL
  ) %>%
  pivot_wider(
    names_from = Region,
    values_from = CAL_UL
  )


perio_PD <- perio_tooth_wide %>%
  select(
    Age,
    Region,
    PD_UL
  ) %>%
  pivot_wider(
    names_from = Region,
    values_from = PD_UL
  )


# ============================================================
# 8. PERIODONTAL TABLE
# ============================================================

# Format person-level percentages to one decimal place.

perio_person_display <- perio_person %>%
  mutate(
    Severe = sprintf("%.1f", Severe),
    Moderate = sprintf("%.1f", Moderate),
    Mild = sprintf("%.1f", Mild)
  )


# Print person-level periodontal section.

perio_person_table <- perio_person_display %>%
  gt() %>%
  tab_header(
    title = "Periodontal Prevalence"
  ) %>%
  tab_spanner(
    label = "Periodontitis Prevalence (%)",
    columns = c(Severe, Moderate, Mild)
  ) %>%
  cols_label(
    Age = "Ages",
    Count = "Count",
    Severe = "Severe",
    Moderate = "Moderate",
    Mild = "Mild"
  )

perio_person_table


# Print tooth-level CAL section.

perio_CAL_table <- perio_CAL %>%
  gt() %>%
  tab_header(
    title = "Tooth-Level Periodontal Prevalence"
  ) %>%
  tab_spanner(
    label = "CAL >= 5 mm (Upper/Lower %)",
    columns = c(Anterior, Premolar, Molar)
  ) %>%
  cols_label(
    Age = "Ages",
    Anterior = "Anterior",
    Premolar = "Premolar",
    Molar = "Molar"
  )

perio_CAL_table


# Print tooth-level PD section.

perio_PD_table <- perio_PD %>%
  gt() %>%
  tab_spanner(
    label = "PD >= 4 mm (Upper/Lower %)",
    columns = c(Anterior, Premolar, Molar)
  ) %>%
  cols_label(
    Age = "Ages",
    Anterior = "Anterior",
    Premolar = "Premolar",
    Molar = "Molar"
  )

perio_PD_table


# ============================================================
# 9. CARIES PREVALENCE - PERSON LEVEL
# ============================================================

# Person-level caries prevalence represents ANY presence.
#
# For example:
#
#   DS = 1 if the person has at least one decayed surface
#
#   FS = 1 if the person has at least one filled surface
#
#   DFS = 1 if the person has at least one decayed OR
#         filled surface
#
#   DMFS = 1 if the person has any decayed, missing, or
#          filled surface
#
# These categories are NOT mutually exclusive.


caries_person_age <- df_person_level %>%
  group_by(Age) %>%
  summarise(
    Count = n(),

    DS = get_percent(
      `Untreated Dental Caries (DS)`
    ),

    DSI = get_percent(
      `Untreated Interproximal Dental Caries (DSI)`
    ),

    FS = get_percent(
      `Treated Dental Caries (FS)`
    ),

    FSI = get_percent(
      `Treated Interproximal Dental Caries (FSI)`
    ),

    DFS = get_percent(
      `Surface Decay or Restoration (DFS)`
    ),

    DFSI = get_percent(
      `Surface Decay or Restoration on Interproximal Surfaces (DFSI)`
    ),

    DMFS = get_percent(
      `Caries Experience (DMFS)`
    ),

    DMFSI = get_percent(
      `Interproximal Caries Experience (DMFSI)`
    ),

    .groups = "drop"
  )


# Overall prevalence.

caries_person_total <- df_person_level %>%
  summarise(
    Count = n(),

    DS = get_percent(
      `Untreated Dental Caries (DS)`
    ),

    DSI = get_percent(
      `Untreated Interproximal Dental Caries (DSI)`
    ),

    FS = get_percent(
      `Treated Dental Caries (FS)`
    ),

    FSI = get_percent(
      `Treated Interproximal Dental Caries (FSI)`
    ),

    DFS = get_percent(
      `Surface Decay or Restoration (DFS)`
    ),

    DFSI = get_percent(
      `Surface Decay or Restoration on Interproximal Surfaces (DFSI)`
    ),

    DMFS = get_percent(
      `Caries Experience (DMFS)`
    ),

    DMFSI = get_percent(
      `Interproximal Caries Experience (DMFSI)`
    )
  ) %>%
  mutate(
    Age = "Total"
  )


caries_person <- bind_rows(
  caries_person_age,
  caries_person_total
) %>%
  mutate(
    Age = as.character(Age)
  ) %>%
  select(
    Age,
    Count,
    DS,
    DSI,
    FS,
    FSI,
    DFS,
    DFSI,
    DMFS,
    DMFSI
  )


# ============================================================
# 10. CARIES PREVALENCE - TOOTH LEVEL
# ============================================================

# The original tooth-level table reports:
#
#   DS prevalence
#   DMFS prevalence
#
# by Anterior / Premolar / Molar and Upper / Lower.


caries_tooth_age <- tooth_table_data %>%
  filter(
    !is.na(Age),
    !is.na(Jaw),
    !is.na(Region)
  ) %>%
  group_by(
    Age,
    Region,
    Jaw
  ) %>%
  summarise(
    DS = get_percent(DS),
    DMFS = get_percent(DMFS),
    .groups = "drop"
  )


caries_tooth_total <- tooth_table_data %>%
  filter(
    !is.na(Jaw),
    !is.na(Region)
  ) %>%
  group_by(
    Region,
    Jaw
  ) %>%
  summarise(
    DS = get_percent(DS),
    DMFS = get_percent(DMFS),
    .groups = "drop"
  ) %>%
  mutate(
    Age = "Total"
  )


caries_tooth <- bind_rows(
  caries_tooth_age,
  caries_tooth_total
) %>%
  mutate(
    Age = as.character(Age)
  )


# ============================================================
# 11. FORMAT CARIES TOOTH PREVALENCE AS UPPER/LOWER
# ============================================================

caries_tooth_wide <- caries_tooth %>%
  pivot_wider(
    names_from = Jaw,
    values_from = c(DS, DMFS)
  ) %>%
  mutate(
    DS_UL = sprintf("%.1f/%.1f", DS_Upper, DS_Lower),
    DMFS_UL = sprintf("%.1f/%.1f", DMFS_Upper, DMFS_Lower)
  ) %>%
  select(
    Age,
    Region,
    DS_UL,
    DMFS_UL
  )


caries_DS <- caries_tooth_wide %>%
  select(
    Age,
    Region,
    DS_UL
  ) %>%
  pivot_wider(
    names_from = Region,
    values_from = DS_UL
  )


caries_DMFS <- caries_tooth_wide %>%
  select(
    Age,
    Region,
    DMFS_UL
  ) %>%
  pivot_wider(
    names_from = Region,
    values_from = DMFS_UL
  )


# ============================================================
# 12. CARIES TABLE
# ============================================================

caries_person_display <- caries_person %>%
  mutate(
    across(
      c(DS, DSI, FS, FSI, DFS, DFSI, DMFS, DMFSI),
      ~ sprintf("%.1f", .x)
    )
  )


caries_person_table <- caries_person_display %>%
  gt() %>%
  tab_header(
    title = "Caries Prevalence"
  ) %>%
  tab_spanner(
    label = "Caries Prevalence (%, I = Interproximal)",
    columns = c(
      DS,
      DSI,
      FS,
      FSI,
      DFS,
      DFSI,
      DMFS,
      DMFSI
    )
  ) %>%
  cols_label(
    Age = "Ages",
    Count = "Count"
  )

caries_person_table


# Tooth-level DS prevalence.

caries_DS_table <- caries_DS %>%
  gt() %>%
  tab_header(
    title = "Tooth-Level Caries Prevalence"
  ) %>%
  tab_spanner(
    label = "DS Prevalence (Upper/Lower %)",
    columns = c(Anterior, Premolar, Molar)
  ) %>%
  cols_label(
    Age = "Ages",
    Anterior = "Anterior",
    Premolar = "Premolar",
    Molar = "Molar"
  )

caries_DS_table


# Tooth-level DMFS prevalence.

caries_DMFS_table <- caries_DMFS %>%
  gt() %>%
  tab_spanner(
    label = "DMFS Prevalence (Upper/Lower %)",
    columns = c(Anterior, Premolar, Molar)
  ) %>%
  cols_label(
    Age = "Ages",
    Anterior = "Anterior",
    Premolar = "Premolar",
    Molar = "Molar"
  )

caries_DMFS_table


# ============================================================
# 13. SAVE TABLE IMAGES
# ============================================================

# Save all table images to the package figures folder.
# here() will resolve this relative to the project root.

gtsave(
  perio_person_table,
  here("figures", "perio_person_table.png"),
  vwidth = 900,
  vheight = 500
)

gtsave(
  perio_CAL_table,
  here("figures", "perio_CAL_table.png"),
  vwidth = 900,
  vheight = 500
)

gtsave(
  perio_PD_table,
  here("figures", "perio_PD_table.png"),
  vwidth = 900,
  vheight = 500
)

gtsave(
  caries_person_table,
  here("figures", "caries_person_table.png"),
  vwidth = 900,
  vheight = 500
)

gtsave(
  caries_DS_table,
  here("figures", "caries_DS_table.png"),
  vwidth = 900,
  vheight = 500
)

gtsave(
  caries_DMFS_table,
  here("figures", "caries_DMFS_table.png"),
  vwidth = 900,
  vheight = 500
)
