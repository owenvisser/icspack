library(gtsummary)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(gt)

####################################################################################################
#
# IMPORT DATA - HELPER MATRIX
#
####################################################################################################

#Set source as current file path
setwd(r"(C:\Users\owvis\Desktop\nHANES\Full Data)")

var.info <- read.csv("variable.info.csv")[,-1]
peri.uid <- read.csv("All_UID.csv")[,-1] # only people with some periodontal data

# we want only the 2:15 & 18:31 teeth. ignoring wisdom teeth.
ohx.nums <-  paste0("OHX", sprintf("%02d", c(2:15, 18:31)))

infovec <- c('CJA', 'CJD', 'CJL', 'CJM', 'CJP', 'CJS',
             'LAA', 'LAD', 'LAL', 'LAM', 'LAP', 'LAS',  
             'PCA', 'PCD', 'PCL', 'PCM', 'PCP', 'PCS',
             'TC','CSC','CTC')

# A matrix of all tooth codes, dental caries surface codes/conditions, and periodontal measurements
teeth_list <- as.data.frame(t(sapply(ohx.nums, function(x){
  sapply(infovec, function(y) paste0(x, y))
})))


####################################################################################################
#
# FILTER
#
####################################################################################################

peri.uid <- peri.uid %>%
  #retain adults 30 and older
  filter(DMDHRAGE >= 30) %>%
  #retain complete status for oral exam.
  filter(OHDEXSTS == 1) %>% 
  rowwise() %>%
  mutate(
    # 1 = Primary tooth
    # 2 = Permanent tooth
    # 3 = Implant
    # 4 = Tooth not present
    # 5 = Permanent dental root fragment
    #Get the count of teeth that are present in a patient
    Tooth_Count = sum(c_across(all_of(teeth_list$TC)) %in% c(1, 2, 3)), 
    #get the number of teeth marked missing
    Teeth_Missing = sum(c_across(all_of(teeth_list$TC)) %in% c(4,5))) %>%
  ungroup() %>%
  #retain people with at least 9 observations
  filter(Tooth_Count >= 9)

agedf <- peri.uid %>% mutate(
  Age = as.factor(cut(as.numeric(peri.uid$DMDHRAGE),
                      breaks = c(29, 49, 64, 80),
                      labels = c("30-49", "50-64", "65+"),
                      include.lowest = TRUE))) %>%
  select(SEQN, Age)



####################################################################################################
#
# CATEGORIZE PERIO
#
####################################################################################################
# Define a custom function to find the mean value excluding NA and 99
avg_exclude <- function(x) {
  # Exclude NA and 99 from x
  valid_values <- x[!is.na(x) & x != 99]
  if (length(valid_values) > 0) {
    return(mean(valid_values))
  } else {
    return(NA)  # Return NA if no valid values are found
  }
}

# Define a custom function to find the maximum value excluding NA and 99
max_exclude <- function(x) {
  # Exclude NA and 99 from x
  valid_values <- x[!is.na(x) & x != 99]
  if (length(valid_values) > 0) {
    return(max(valid_values))
  } else {
    return(NA)  # Return NA if no valid values are found
  }
}

# Define a custom function to check the number of x = max(x) excluding NA and 99
max_exclude_count <- function(x) {
  # Exclude NA and 99 from x
  valid_values <- x[!is.na(x) & x != 99]
  if (length(valid_values) > 0) {
    return(sum(valid_values == max(valid_values)))
  } else {
    return(NA)  # Return NA if no valid values are found
  }
}


get_over_n <- function(x, n) {
  valid_values <- x[!is.na(x) & x != 99]
  if (length(valid_values) > 0) {
    return(ifelse(any(valid_values > n), 1, 0))
  } else {
    return(NA)  # Return NA if no valid values are found
  }
}
 

Perio_long <- peri.uid %>%
  # Want these columns
  select(c('SEQN', 
           unlist(as.vector(teeth_list %>% select(-c(TC, CTC, CSC))))  )) %>%
  # Set the Names
  setNames(c('SEQN',
             unlist(as.vector(teeth_list %>% select(-c(TC, CTC, CSC))))  )) %>% 
  # Pivot to a long matrix for categorization
  pivot_longer(cols = starts_with('OHX'),
               names_to = 'ToothCodename',
               values_to = 'Measurement',
               values_transform = as.character) %>%
  # Separate Tooth Number and Measurement Abbr.
  mutate(
    Measurement = as.numeric(Measurement),
    ToothNumber = substr(ToothCodename, 4, 5),
    MeasurementCode = substr(ToothCodename, 6, 7),
    MeasurementLoc = substr(ToothCodename, 8, 8)
  ) %>%
  mutate(
    # Measurement Type
    MeasurementType = case_when(
      MeasurementCode  == 'CJ' ~ 'FGM to CEJ measurement (mm)',
      MeasurementCode  == 'LA' ~ 'Clinical Loss of Attachment (mm)',
      MeasurementCode  == 'PC' ~ 'Pocket Depth (mm)'),
    # Face of Tooth Measurement is On
    MeasurementFace = case_when(
      MeasurementLoc == 'A' ~ 'mesio-lingual',
      MeasurementLoc == 'D' ~ 'distal',
      MeasurementLoc == 'L' ~ 'mid-lingual',
      MeasurementLoc == 'M' ~ 'mid-facial',
      MeasurementLoc == 'P' ~ 'distal-lingual',
      MeasurementLoc == 'S' ~ 'mesia-facial'
    )) %>%
  # Filter out Non-Important Measurements
  filter(MeasurementType != 'FGM to CEJ measurement (mm)') %>%
  # Filter out Non-Important Faces
  #filter(!(MeasurementFace %in% c('mid-lingual' , 'mid-facial'))) %>%
  # Clean Up Columns
  select(SEQN, ToothNumber, MeasurementType, MeasurementFace, Measurement)  

###
### Short section for tooth level
###
Perio_tooth <- Perio_long %>%
  mutate(
    Tooth_type = case_when(
      ToothNumber %in% c('02','03','14','15') ~ 'Upper Molars',
      ToothNumber %in% c('04','05','12','13') ~ 'Upper Premolars',
      ToothNumber %in% c('06', '11', '07','08','09','10') ~ 'Upper Anterior',
      ToothNumber %in% c('18','19','30','31') ~ 'Lower Molars',
      ToothNumber %in% c('20','21','28','29') ~ 'Lower Premolars',
      ToothNumber %in% c('22', '27', '23','24','25','26') ~ 'Lower Anterior',
      TRUE ~ NA
      )
    )

save(Perio_tooth, file = r"(C:\Users\owvis\Desktop\nHANES\Full Data\LongData_Filtered.RData)")

# Perio_site_level <- Perio_tooth %>%
#   group_by(SEQN, ToothNumber, MeasurementFace, Tooth_type) %>%
#   summarize(
#     `CAL >= 3mm` = get_over_n(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)'], 2),
#     `CAL >= 4mm` = get_over_n(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)'], 3),
#     `CAL >= 5mm` = get_over_n(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)'], 6),
#     `PD >= 4mm` = get_over_n(Measurement[MeasurementType == 'Pocket Depth (mm)'], 3),
#     `PD >= 5mm` = get_over_n(Measurement[MeasurementType == 'Pocket Depth (mm)'], 4),
#     `PD >= 6mm` = get_over_n(Measurement[MeasurementType == 'Pocket Depth (mm)'], 5)
#   ) %>%
#   ungroup() %>%
#   group_by(SEQN, ToothNumber, MeasurementFace) %>%
#   pivot_wider(names_from = Tooth_type,
#               values_from = c(`CAL >= 3mm`, `CAL >= 4mm`, `CAL >= 5mm`, `PD >= 4mm`, `PD >= 5mm`, `PD >= 6mm`))
# 
# Perio_site_level <- agedf %>% full_join(Perio_site_level) %>% select(-ToothNumber)


Perio_tooth_level_allMeasures <- Perio_tooth %>%
  group_by(SEQN, ToothNumber, Tooth_type) %>%
  summarize(.groups = "drop_last",
            maxCAL = max(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)']),
            avgCAL = mean(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)']),
            maxPD = max(Measurement[MeasurementType == 'Pocket Depth (mm)']),
            avgPD = mean(Measurement[MeasurementType == 'Pocket Depth (mm)'])
  ) %>%
  mutate(ToothNumber = as.numeric(ToothNumber),
         maxCAL = ifelse(maxCAL > 90, NA, maxCAL),
         avgCAL = ifelse(avgCAL > 90, NA, avgCAL),
         avgPD = ifelse(avgPD > 90, NA, avgPD),
         maxPD = ifelse(maxPD > 90, NA, maxPD),
         ) %>%
  ungroup()


Perio_tooth_level <- Perio_tooth %>%
  group_by(SEQN, ToothNumber, Tooth_type) %>%
  summarize(
    `CAL >= 3mm` = get_over_n(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)'], 2),
    `CAL >= 4mm` = get_over_n(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)'], 3),
    `CAL >= 5mm` = get_over_n(Measurement[MeasurementType == 'Clinical Loss of Attachment (mm)'], 4),
    `PD >= 4mm` = get_over_n(Measurement[MeasurementType == 'Pocket Depth (mm)'], 3),
    `PD >= 5mm` = get_over_n(Measurement[MeasurementType == 'Pocket Depth (mm)'], 4),
    `PD >= 6mm` = get_over_n(Measurement[MeasurementType == 'Pocket Depth (mm)'], 5)
  ) %>%
  mutate(ToothNumber = as.numeric(ToothNumber)) %>%
  ungroup()



###
### Finish person level
###


Perio_long <- Perio_long %>%
  # Keep the order
  group_by(SEQN, ToothNumber, MeasurementType) %>%
  summarize(
    #Highest severity score at mesial/distal locations
    ToothSeverity_Mesial_Distal = max_exclude(Measurement[MeasurementFace %in% c('mesio-lingual', 'distal', 'distal-lingual', 'mesio-facial')]),
    #Number of the highest score at mesial distal locations
    TSNum_Max_Mesial_Distal = max_exclude_count(Measurement[MeasurementFace %in% c('mesio-lingual', 'distal', 'distal-lingual', 'mesio-facial')]),
    #Highest severity score at any location
    ToothSeverity_Total = max_exclude(Measurement)
            ) %>%
  group_by(SEQN) %>% 
  summarize(
    Severe_PD = ifelse(sum(ToothSeverity_Mesial_Distal[MeasurementType == 'Clinical Loss of Attachment (mm)'] %in% 6:12, na.rm = TRUE) >= 2 &
                      sum(ToothSeverity_Mesial_Distal[MeasurementType == 'Pocket Depth (mm)'] %in% 5:12, na.rm = TRUE) >= 1,
                    1, 0),
    Moderate_PD = ifelse(sum(ToothSeverity_Mesial_Distal[MeasurementType == 'Clinical Loss of Attachment (mm)'] %in% 4:12, na.rm = TRUE) >= 2 |
                        sum(ToothSeverity_Mesial_Distal[MeasurementType == 'Pocket Depth (mm)'] %in% 5:12, na.rm = TRUE) >= 2,
                      1, 0),
    # mild periodontitis defined as 2 or more interproximal sites with CAL 3 mm or greater -can be same tooth?-
    # and 
    # (
    # 2 or more interproximal sites with PPD 4 mm or greater (not on the same tooth) 
    # or 
    # 1 or more sites with 5 mm or more.
    # )
    Mild_PD = ifelse(
      all(
        any(
          sum(ToothSeverity_Mesial_Distal[MeasurementType == 'Clinical Loss of Attachment (mm)'] %in% 3:12, na.rm = TRUE) >= 2, #two diff teeth
          TSNum_Max_Mesial_Distal[ToothSeverity_Mesial_Distal[MeasurementType == 'Clinical Loss of Attachment (mm)'] %in% 3:12] >= 2 #could be same tooth
        , na.rm = TRUE), 
        any(
          sum(ToothSeverity_Mesial_Distal[MeasurementType == 'Pocket Depth (mm)'] %in% 4:12, na.rm = TRUE) >= 2,
          sum(ToothSeverity_Total[MeasurementType == 'Pocket Depth (mm)'] %in% 5:12, na.rm = TRUE) >= 1
          , na.rm = TRUE)
        ), 1, 0)) %>%
  mutate(Moderate_PD = ifelse(Severe_PD == 1, 0, Moderate_PD),
         Mild_PD = ifelse(Severe_PD == 1 | Moderate_PD == 1 , 0, Mild_PD),
         No_PD = ifelse(Mild_PD == 1 | Moderate_PD == 1 | Severe_PD == 1, 0, 1))


####################################################################################################
#
# CATEGORIZE DENTAL CARIES
#
####################################################################################################

Caries_long <- peri.uid %>%
  # Want these columns
  select(c('SEQN',
           unlist(as.vector(teeth_list %>% select(CTC, CSC)))  )) %>%
  # Set the Names
  setNames(c('SEQN',
             unlist(as.vector(teeth_list %>% select(CTC, CSC)))  )) %>% 
  # Pivot to a long matrix for categorization
  pivot_longer(cols = starts_with('OHX'),
               names_to = 'ToothCodename',
               values_to = 'Measurement',
               values_transform = as.character) %>%
  # Separate Tooth Number and Measurement Abbr.
  mutate(
    ToothNumber = as.numeric(substr(ToothCodename, 4, 5)),
    MeasurementCode = substr(ToothCodename, 6, 8)
  ) %>%
  select(SEQN, ToothNumber, MeasurementCode, Measurement) %>%
  pivot_wider(
    names_from = MeasurementCode,
    values_from = Measurement
  ) %>%
  #assign surface codes to ctc == z cases.
  mutate(`0` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('0', CSC), 1, 0),
         `1` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('1', CSC), 1, 0),
         `2` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('2', CSC), 1, 0),
         `3` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('3', CSC), 1, 0),
         `4` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('4', CSC), 1, 0),
         `5` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('5', CSC), 1, 0),
         `6` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('6', CSC), 1, 0),
         `7` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('7', CSC), 1, 0),
         `8` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('8', CSC), 1, 0),
         `9` = ifelse((CTC == 'Z' |  CTC == 'K') & grepl('9', CSC), 1, 0))%>% 
  select(-CSC) %>%
  mutate(
    DS_interprox = case_when(
      CTC == 'Z' | CTC == 'K' ~ `3` + `4`,
      CTC == 'J' & ToothNumber %in% c(6:11, 22:27) ~ 2,
      CTC == 'J' & !(ToothNumber %in% c(6:11, 22:27)) ~ 2,
      TRUE ~ 0
    ),
    DS = case_when(
      CTC == 'Z' | CTC == 'K' ~ `0` + `1` + `2` + `3` + `4`,
      CTC == 'J' & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == 'J' & !(ToothNumber %in% c(6:11, 22:27)) ~ 5,
      TRUE ~ 0
      ),
    FS_interprox = case_when(
      CTC == 'Z' | CTC == 'K' ~ `8` + `9`,
      CTC == 'T' & ToothNumber %in% c(6:11, 22:27) ~ 2,
      CTC == 'T' & !(ToothNumber %in% c(6:11, 22:27)) ~ 2,
      TRUE ~ 0
    ),
    FS = case_when(
      CTC == 'Z' | CTC == 'K' ~ `5` + `6` + `7` + `8` + `9`,
      CTC == 'T' & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == 'T' & !(ToothNumber %in% c(6:11, 22:27)) ~ 5,
      TRUE ~ 0
    ),
    MS = case_when(
      CTC == 'E' & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == 'E' & !(ToothNumber %in% c(6:11, 22:27)) ~ 5,
      CTC == 'R' & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == 'R' & !(ToothNumber %in% c(6:11, 22:27)) ~ 5,
      CTC == 'P' & ToothNumber %in% c(6:11, 22:27) ~ 4,
      CTC == 'P' & !(ToothNumber %in% c(6:11, 22:27)) ~ 5,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    DS_Tooth_Interprox = DS_interprox,
    DS_Tooth = DS, 
    FS_Tooth = FS, 
    FS_Tooth_Interprox = FS_interprox,
    DFS_Tooth = DS + FS,
    DFS_Tooth_Interprox = DS_interprox + FS_interprox,
    DMFS_Tooth = DS + FS + MS,
    DMFS_Tooth_Interprox = DS_interprox + FS_interprox + MS
  )

Caries_tooth <- Caries_long


Caries_tooth <- Caries_tooth %>%
  mutate(
    Tooth_type = case_when(
      ToothNumber %in% c(2,3,14,15) ~ 'Upper Molars',
      ToothNumber %in% c(4,5,12,13) ~ 'Upper Premolars',
      ToothNumber %in% c(6:11) ~ 'Upper Anterior',
      ToothNumber %in% c(18,19,30,31) ~ 'Lower Molars',
      ToothNumber %in% c(20,21,28,29) ~ 'Lower Premolars',
      ToothNumber %in% c(22:27) ~ 'Lower Anterior',
      TRUE ~ NA
    ))

Caries_tooth_level <- Caries_tooth %>%
  group_by(SEQN, ToothNumber, Tooth_type) %>%
  summarise(
    FS_score = sum(FS_Tooth),
    FSI_score = sum(FS_Tooth_Interprox),
    DSI_score = sum(DS_Tooth_Interprox),
    DS_score = sum(DS_Tooth),
    DFS_score = sum(DFS_Tooth),
    DFSI_score = sum(DFS_Tooth_Interprox),
    DMFS_score = sum(DMFS_Tooth),
    DMFSI_score = sum(DMFS_Tooth_Interprox)
  ) %>%
  mutate(
    `DFS` = ifelse(DFS_score > 0, 1, 0),
    `DFSI` = ifelse(DFSI_score > 0, 1, 0),
    `DSI` = ifelse(DSI_score > 0, 1 , 0),
    `DS` = ifelse(DS_score > 0, 1 , 0),
    `FSI` = ifelse(FSI_score > 0, 1 , 0),
    `FS` = ifelse(FS_score > 0, 1 , 0),
    `DMFS` = ifelse(DMFS_score > 0, 1, 0),
    `DMFSI` = ifelse( DMFSI_score >0, 1, 0)) %>%
  select(SEQN, ToothNumber, Tooth_type, DS, DSI, FS, FSI, DFS, DFSI, DMFS, DMFSI) %>%
  ungroup()


Caries_long <- Caries_long %>%
  select(SEQN, FS_Tooth_Interprox, FS_Tooth, DS_Tooth_Interprox, DS_Tooth, DMFS_Tooth, DMFS_Tooth_Interprox, DFS_Tooth, DFS_Tooth_Interprox) %>%
  group_by(SEQN) %>%
  summarise(
    FS_score = sum(FS_Tooth),
    FSI_score = sum(FS_Tooth_Interprox),
    DSI_score = sum(DS_Tooth_Interprox),
    DS_score = sum(DS_Tooth),
    DFS_score = sum(DFS_Tooth),
    DFSI_score = sum(DFS_Tooth_Interprox),
    DMFS_score = sum(DMFS_Tooth),
    DMFSI_score = sum(DMFS_Tooth_Interprox)
  ) %>%
  mutate(
    `Surface Decay or Restoration (DFS)` = ifelse(DFS_score > 0, 1, 0),
    `Surface Decay or Restoration on Interproximal Surfaces (DFSI)` = ifelse(DFSI_score > 0, 1, 0),
    `Untreated Interproximal Dental Caries (DSI)` = ifelse(DSI_score > 0, 1 , 0),
    `Untreated Dental Caries (DS)` = ifelse(DS_score > 0, 1 , 0),
    `Treated Interproximal Dental Caries (FSI)` = ifelse(FSI_score > 0, 1 , 0),
    `Treated Dental Caries (FS)` = ifelse(FS_score > 0, 1 , 0),
    `Caries Experience (DMFS)` = ifelse(DMFS_score > 0, 1, 0),
    `Interproximal Caries Experience (DMFSI)` = ifelse( DMFSI_score >0, 1, 0)) %>%
  mutate(FSI_score = ifelse(FSI_score > 0, FSI_score, 0),
         FS_score = ifelse(FS_score > 0, FS_score, 0),
         DSI_score = ifelse(DSI_score > 0, DSI_score, 0),
         DS_score = ifelse(DS_score > 0, DS_score, 0),
         DMFS_score = ifelse(DMFS_score > 0, DMFS_score, 0),
         DFSI_score = ifelse(DFSI_score > 0, DFSI_score, 0),
         DFS_score = ifelse(DFS_score > 0, DFS_score, 0),
         DMFSI_score = ifelse(DMFSI_score > 0 , DMFSI_score, 0))
  

####################################################################################################
#
# FIND IMPORTANT COVARIATES & CREATING THE DATA FRAME
#
####################################################################################################

#TBX = uid
#OHX = teeth 
#SMQ = smoke 
#DEMO = demographic 
# For more DEMO variables see
# https://wwwn.cdc.gov/nchs/nhanes/2011-2012/DEMO_G.htm
# https://wwwn.cdc.gov/nchs/nhanes/2013-2014/DEMO_H.htm
# https://wwwn.cdc.gov/nchs/nhanes/2013-2014/SMQ_H.htm

#  Change DEMO to SMQ, OHX, TBX for more variable info.

df <- peri.uid %>% mutate(
  Gender = as.factor(ifelse(RIAGENDR == 1, 'Male', 'Female')),
  Race = as.factor(case_when(
    RIDRETH3 == 1 ~ 'Mexican American', 
    RIDRETH3 == 2 ~ 'Other Hispanic',
    RIDRETH3 == 3 ~ 'Non-Hispanic White',
    RIDRETH3 == 4 ~ 'Non-Hispanic Black',
    RIDRETH3 == 5 ~ 'Non-Hispanic Asian',
    RIDRETH3 == 6 ~ 'Other Race - Including Multi-Racial',
    TRUE ~ NA
  )),
  Age = as.factor(cut(as.numeric(peri.uid$DMDHRAGE),
            breaks = c(29, 49, 64, 80),
            labels = c("30-49", "50-64", "65+"),
            include.lowest = TRUE)),
  Income_Family = as.factor(case_when(
    INDHHIN2 == 1 ~ '$0 to $4,999',
    INDHHIN2 == 2 ~ '$ 5,000 to $ 9,999',
    INDHHIN2 == 3 ~ '$10,000 to $14,999',
    INDHHIN2 == 4 ~ '$15,000 to $19,999',
    INDHHIN2 == 5 ~ '$20,000 to $24,999',
    INDHHIN2 == 6 ~ '$25,000 to $34,999',
    INDHHIN2 == 7 ~ '$35,000 to $44,999',
    INDHHIN2 == 8 ~ '$45,000 to $54,999',
    INDHHIN2 == 9 ~ '$55,000 to $64,999',
    INDHHIN2 == 10 ~ '$65,000 to $74,999',
    INDHHIN2 == 12 ~ '$20,000 and Over',
    INDHHIN2 == 13 ~ 'Under $20,000',
    INDHHIN2 == 14 ~ '$75,000 to $99,999',
    INDHHIN2 == 15 ~ '$100,000 and Over',
    INDHHIN2 == 77 ~ 'Refused',
    INDHHIN2 == 99 ~ 'Do not know',
    TRUE ~ NA
  )),
  Income_Household = as.factor(case_when(
    INDFMIN2 == 1 ~ '$0 to $4,999',
    INDFMIN2 == 2 ~ '$ 5,000 to $ 9,999',
    INDFMIN2 == 3 ~ '$10,000 to $14,999',
    INDFMIN2 == 4 ~ '$15,000 to $19,999',
    INDFMIN2 == 5 ~ '$20,000 to $24,999',
    INDFMIN2 == 6 ~ '$25,000 to $34,999',
    INDFMIN2 == 7 ~ '$35,000 to $44,999',
    INDFMIN2 == 8 ~ '$45,000 to $54,999',
    INDFMIN2 == 9 ~ '$55,000 to $64,999',
    INDFMIN2 == 10 ~ '$65,000 to $74,999',
    INDFMIN2 == 12 ~ '$20,000 and Over',
    INDFMIN2 == 13 ~ 'Under $20,000',
    INDFMIN2 == 14 ~ '$75,000 to $99,999',
    INDFMIN2 == 15 ~ '$100,000 and Over',
    INDFMIN2 == 77 ~ 'Refused',
    INDFMIN2 == 99 ~ 'Do not know',
    TRUE ~ NA
  )),
  Smoking_Status = as.factor(case_when(
    SMQ040 == 1 ~ 'Every day',
    SMQ040 == 2 ~ 'Some days',
    SMQ040 == 3 ~ 'Not at all',
    SMQ040 == 4 ~ 'Refused',
    SMQ040 == 5 ~ 'Do not know',
    TRUE ~ NA
  )),
  Missing_Teeth = as.factor(case_when(
    Teeth_Missing %in% 0 ~ '0',
    Teeth_Missing %in% 1:5 ~ '1-5',
    Teeth_Missing %in% 6:28 ~ '6-27',
    TRUE ~ NA
  ))) %>%
  select(SEQN, Age, Gender, Race, Income_Family, Income_Household, Smoking_Status, Missing_Teeth, Teeth_Missing, Tooth_Count)



df_person_level <- df %>%
  left_join(Perio_long, by = 'SEQN') %>%
  left_join(Caries_long, by = 'SEQN') 


df_tooth_level <- Perio_tooth_level_allMeasures  %>%
  full_join(Caries_tooth_level, by = c('SEQN', "ToothNumber", "Tooth_type")) 

df_tooth_level <- df_tooth_level %>%
  left_join(
    df_person_level %>% select(SEQN, Tooth_Count),
    by = "SEQN"
  )

save(df_person_level,  file = r"(C:\Users\owvis\Desktop\nHANES\Full Data\Person_Level_Data.RData)")

save(df_tooth_level,  file = r"(C:\Users\owvis\Desktop\nHANES\Full Data\Tooth_Level_Data.RData)")





















# IDK some shit below i forgot what I used it for. 




library(gtsummary)
df1 %>%
  tbl_summary(by = Age,
              statistic = list(all_continuous() ~ "{mean} / {sd} ({p10}, {p25}, {p75}, {p90})", all_categorical() ~ "{n} ({p}%)"),
              # Handling missing data,
              digits = list(Severe_PD ~ c(0, 1),
                            Moderate_PD ~ c(0, 1),
                            Mild_PD ~ c(0, 1),
                            `Untreated Interproximal Dental Caries (DSI)` ~ c(0, 1),
                            `Untreated Dental Caries (DS)` ~ c(0, 1),
                            `Caries Experience (DMFS)` ~ c(0, 1),
                            `Treated Dental Caries (FS)` ~ c(0, 1), 
                            `Treated Interproximal Dental Caries (FSI)` ~ c(0, 1),
              `Surface Decay or Restoration (DFS)` ~ c(0, 1),
              `Surface Decay or Restoration on Interproximal Surfaces (DFSI)` ~ c(0, 1),
              `Interproximal Caries Experience (DMFSI)` ~ c(0, 1)),
              missing = "no",
              include = c(
                Severe_PD, 
                Mild_PD, 
                Moderate_PD,
                `Surface Decay or Restoration (DFS)`,
                `Surface Decay or Restoration on Interproximal Surfaces (DFSI)`,
                `Untreated Dental Caries (DS)`,
                `Untreated Interproximal Dental Caries (DSI)`,
                `Treated Dental Caries (FS)`, 
                `Treated Interproximal Dental Caries (FSI)`,
                `Caries Experience (DMFS)`,
                `Interproximal Caries Experience (DMFSI)`)) %>%
  add_overall()


#Tooth level stuff
tl <- Perio_tooth_level %>% tbl_summary(
  by = Age, 
  statistic = list(all_continuous() ~ "{mean} ", all_categorical() ~ "{n} ({p}%)"),
  digits = list(all_categorical() ~ c(0, 1)),
  missing = 'no',
  include = c(-SEQN)
) %>%
  add_overall()


#Tooth level stuff
tlc <- Caries_tooth_level %>% tbl_summary(
  by = Age, 
  statistic = list(all_continuous() ~ "{mean} ", all_categorical() ~ "{n} ({p}%)"),
  digits = list(all_categorical() ~ c(0, 1)),
  missing = 'no',
  include = c(-SEQN)
) %>%
  add_overall()

# 
# #Site level stuff
# sl <- Perio_site_level %>% tbl_summary(
#   by = Age, 
#   statistic = list(all_continuous() ~ "{mean} ", all_categorical() ~ "{n} ({p}%)"),
#   digits = list(all_categorical() ~ c(0, 1)),
#   missing = 'no',
#   include = c(-SEQN)
# ) %>%
#   add_overall()


latex_tl <- gt::as_latex(as_gt(tl))
latex_tlc <- gt::as_latex(as_gt(tlc))
# 
# latex_sl <- gt::as_latex(as_gt(sl))

cat(as.character(latex_tl), file = "tooth_level_table.tex")
cat(as.character(latex_tlc), file = "tooth_level_table_caries.tex")
# cat(as.character(latex_sl), file = "site_level_table.tex")

#See teeth_list definitions below:
# var.info %>%
#   filter(grepl('OHX02', Variable.Name)) %>%
#   filter(grepl('Loss of', Variable.Description)) %>%
#   select(Variable.Name, Variable.Description)


