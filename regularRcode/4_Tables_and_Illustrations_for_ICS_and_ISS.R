library(gtsummary)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(gt)

####################################################################################################
#
# Import Data
#
####################################################################################################

#Set source as current file path - ON LAPTOP
setwd(r"(C:\Users\owvis\OneDrive - University of Florida\nHANES\Full Data)")
load(file = r"(C:\Users\owvis\OneDrive - University of Florida\nHANES\Full Data\Person_Level_Data.RData)")
load(file = r"(C:\Users\owvis\OneDrive - University of Florida\nHANES\Full Data\Tooth_Level_Data.RData)")


df_plot <- df_person_level %>%
  filter(Tooth_Count > 0) %>%
  mutate(
    PD_Status = case_when(
      No_PD == 1 ~ "No PD",
      Mild_PD == 1 ~ "Mild PD",
      Moderate_PD == 1 ~ "Moderate PD",
      Severe_PD == 1 ~ "Severe PD",
      TRUE ~ NA_character_
    ),
    PD_Status = factor(PD_Status,
                       levels = c("No PD", "Mild PD", "Moderate PD", "Severe PD"))
  )

boxplot(Tooth_Count ~ PD_Status, data = df_plot,
        xlab = "Periodontitis Severity", ylab = "Tooth Count")


df_plot2 <- df_person_level %>%
  filter(Tooth_Count > 0) %>%
  mutate(
    PD_Status = case_when(
      No_PD == 1 ~ "No PD",
      Mild_PD == 1 ~ "Mild PD",
      Moderate_PD == 1 ~ "Moderate PD",
      Severe_PD == 1 ~ "Severe PD",
      TRUE ~ NA_character_
    ),
    PD_Status = factor(PD_Status,
                       levels = c("No PD", "Mild PD", "Moderate PD", "Severe PD")),
    
    Tooth_Bin = cut(
      Tooth_Count,
      breaks = c(0, 6, 12, 18, 24, 28),
      labels = c("1-6", "7-12", "13-18", "19-24", "25-28"),
      include.lowest = TRUE
    )
  )

tab <- table(df_plot2$Tooth_Bin, df_plot2$PD_Status)


matplot(
  x = 1:nrow(tab),
  y = tab,
  type = "l",
  lty = 1,
  lwd = 2,
  xaxt = "n",
  xlab = "Tooth Count (Binned)",
  ylab = "Frequency",
  col = 1:4
)

axis(1, at = 1:nrow(tab), labels = rownames(tab))

legend("topleft",
       legend = colnames(tab),
       col = 1:4,
       lty = 1,
       lwd = 2)


df_plot3 <- df_person_level %>%
  mutate(
    DMFSI_bin = cut(
      DMFSI_score,
      breaks = seq(0, 128, by = 16),
      right = FALSE
    )
  )

boxplot(Tooth_Count ~ DMFSI_bin, data = df_plot3,
        xlab = "DMFSI Score",
        ylab = "Tooth Count")

#hence we have informative cluster size.



str(df_tooth_level)



# mild periodontitis defined as 2 or more interproximal sites with CAL 3 mm or greater -can be same tooth?-
# and 
# (
# 2 or more interproximal sites with PPD 4 mm or greater (not on the same tooth) 
# or 
# 1 or more sites with 5 mm or more.
# )

library(dplyr)
library(tidyr)
library(ggplot2)

df_cluster <- df_tooth_level %>%
  filter(!is.na(`CAL >= 3mm`) & !is.na(`CAL >= 4mm`) & !is.na(`CAL >= 5mm`)) %>%
  group_by(SEQN) %>%
  summarise(
    Tooth_Count = first(Tooth_Count),
    CAL_teeth_n = n(),
    CAL3 = sum(`CAL >= 3mm`, na.rm = TRUE),
    CAL4 = sum(`CAL >= 4mm`, na.rm = TRUE),
    CAL5 = sum(`CAL >= 5mm`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Tooth_Count > 0, CAL_teeth_n > 0) %>%
  mutate(
    Tooth_bin = cut(
      Tooth_Count,
      breaks = seq(0, 30, by = 5),
      include.lowest = TRUE,
      right = TRUE
    ),
    CAL3_prop = CAL3 / CAL_teeth_n,
    CAL4_prop = CAL4 / CAL_teeth_n,
    CAL5_prop = CAL5 / CAL_teeth_n
  ) %>%
  select(SEQN, Tooth_bin, CAL3_prop, CAL4_prop, CAL5_prop) %>%
  pivot_longer(
    cols = c(CAL3_prop, CAL4_prop, CAL5_prop),
    names_to = "Measure",
    values_to = "Prop"
  ) %>%
  mutate(
    Measure = factor(
      Measure,
      levels = c("CAL3_prop", "CAL4_prop", "CAL5_prop"),
      labels = c("CAL ≥ 3mm", "CAL ≥ 4mm", "CAL ≥ 5mm")
    ),
    Prop_bin = cut(
      Prop,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      right = TRUE
    )
  )

df_heat <- df_cluster %>%
  filter(!is.na(Tooth_bin), !is.na(Prop_bin)) %>%
  count(Measure, Tooth_bin, Prop_bin, name = "n") %>%
  group_by(Measure, Tooth_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_heat, aes(x = Prop_bin, y = Tooth_bin, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ Measure, ncol = 1) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(
    x = "Proportion of Teeth in Subgroup",
    y = "Total Tooth Count (Binned)",
    fill = "Row Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )






library(dplyr)
library(tidyr)
library(ggplot2)

df_cluster <- df_tooth_level %>%
  filter(!is.na(DS) & !is.na(FS) & !is.na(DMFS)) %>%
  group_by(SEQN) %>%
  summarise(
    Tooth_Count = first(Tooth_Count),
    Caries_teeth_n = n(),
    DS_n = sum(DS, na.rm = TRUE),
    FS_n = sum(FS, na.rm = TRUE),
    DMFS_n = sum(DMFS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Tooth_Count > 0, Caries_teeth_n > 0) %>%
  mutate(
    Tooth_bin = cut(
      Tooth_Count,
      breaks = seq(0, 30, by = 5),
      include.lowest = TRUE,
      right = TRUE
    ),
    DS_prop = DS_n / Caries_teeth_n,
    FS_prop = FS_n / Caries_teeth_n,
    DMFS_prop = DMFS_n / Caries_teeth_n
  ) %>%
  select(SEQN, Tooth_bin, DS_prop, FS_prop, DMFS_prop) %>%
  pivot_longer(
    cols = c(DS_prop, FS_prop, DMFS_prop),
    names_to = "Measure",
    values_to = "Prop"
  ) %>%
  mutate(
    Measure = factor(
      Measure,
      levels = c("DS_prop", "FS_prop", "DMFS_prop"),
      labels = c("DS", "FS", "DMFS")
    ),
    Prop_bin = cut(
      Prop,
      breaks = seq(0, 1, by = 0.1),
      include.lowest = TRUE,
      right = TRUE
    )
  )

df_heat <- df_cluster %>%
  filter(!is.na(Tooth_bin), !is.na(Prop_bin)) %>%
  count(Measure, Tooth_bin, Prop_bin, name = "n") %>%
  group_by(Measure, Tooth_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(df_heat, aes(x = Prop_bin, y = Tooth_bin, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ Measure, ncol = 1) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(
    x = "Proportion of Teeth in Subgroup",
    y = "Total Tooth Count (Binned)",
    fill = "Row Proportion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )



