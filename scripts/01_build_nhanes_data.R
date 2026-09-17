# 01_build_nhanes_data.R
#
# Purpose:
#   Read NHANES 2011-2012 and 2013-2014 source files, combine survey cycles, 
#   merge datasets by SEQN, collect NHANES variable descriptions, and save
#   the resulting datasets to the package data directory.
#
# Input: 
#    raw data files from NHANES survey database. These files are 
#    specific to the data values required for this study. 
#    See the *_Oral_Health_Recorders_Procedures_Manual.pdf
#    for more information.
#
# Output:
# data/variable_info.csv
#   containing all variable information
# data/nhanes_full.csv
#   containing all variable entries across all fields
# data/nhanes_perio.csv
#   containing only variable entries for periodontitis and caries related fields

# 1. LOAD REQUIRED PACKAGES

library(dplyr)
library(here)
library(rvest)
library(foreign)

# 2. DEFINE FILE PATHS

path_2011_2012 <- here("data-raw", "2011_2012")
path_2013_2014 <- here("data-raw", "2013_2014")
output_path <- here("data")

# 3. READ NHANES 2011-2012 DATA

demo_12 <- foreign::read.xport(file.path(path_2011_2012, "DEMO_G.XPT"))
peri_12 <- foreign::read.xport(file.path(path_2011_2012, "OHXPER_G.XPT"))
dent_12 <- foreign::read.xport(file.path(path_2011_2012, "OHXDEN_G.XPT"))
smok_cigUse_12 <- foreign::read.xport(file.path(path_2011_2012, "SMQ_G.XPT"))
smok_recent_12 <- foreign::read.xport(file.path(path_2011_2012, "SMQRTU_G.XPT"))

# 4. READ NHANES 2013-2014 DATA

demo_34 <- foreign::read.xport(file.path(path_2013_2014, "DEMO_H.XPT"))
peri_34 <- foreign::read.xport(file.path(path_2013_2014, "OHXPER_H.XPT"))
dent_34 <- foreign::read.xport(file.path(path_2013_2014, "OHXDEN_H.XPT"))
smok_cigUse_34 <- foreign::read.xport(file.path(path_2013_2014, "SMQ_H.XPT"))
smok_recent_34 <- foreign::read.xport(file.path(path_2013_2014, "SMQRTU_H.XPT"))

# 5. CHECK THAT SEQN IS UNIQUE WITHIN EACH DATASET

nhanes_files <- list(
  demo_12 = demo_12,
  peri_12 = peri_12,
  dent_12 = dent_12,
  smok_cigUse_12 = smok_cigUse_12,
  smok_recent_12 = smok_recent_12,
  demo_34 = demo_34,
  peri_34 = peri_34,
  dent_34 = dent_34,
  smok_cigUse_34 = smok_cigUse_34,
  smok_recent_34 = smok_recent_34
)

check_unique_seqn <- function(data) {
  length(data$SEQN) == length(unique(data$SEQN))
}

seqn_check <- vapply(
  nhanes_files,
  check_unique_seqn,
  logical(1)
)

print(seqn_check)

# 6. COMBINE NHANES SURVEY CYCLES

peri_all <- bind_rows(
  peri_12,
  peri_34
)

# Keep only demographic variables appearing in both cycles.
common_demo_variables <- intersect(
  names(demo_12),
  names(demo_34)
)

demo_all <- bind_rows(
  demo_12 %>% select(all_of(common_demo_variables)),
  demo_34 %>% select(all_of(common_demo_variables))
)

dent_all <- bind_rows(
  dent_12,
  dent_34
)

smok_cigUse_all <- bind_rows(
  smok_cigUse_12,
  smok_cigUse_34
)

smok_recent_all <- bind_rows(
  smok_recent_12,
  smok_recent_34
)

# 7. MERGE DATASETS BY PARTICIPANT ID

fulldata <- peri_all %>%
  full_join(dent_all, by = "SEQN") %>%
  full_join(demo_all, by = "SEQN") %>%
  full_join(smok_cigUse_all, by = "SEQN") %>%
  full_join(smok_recent_all, by = "SEQN")

# 8. DOWNLOAD NHANES VARIABLE INFORMATION

url12_exam <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Examination&Cycle=2011-2012"
url12_demo <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Demographics&Cycle=2011-2012"
url12_ques <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Questionnaire&Cycle=2011-2012"
url34_exam <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Examination&Cycle=2013-2014"
url34_demo <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Demographics&Cycle=2013-2014"
url34_ques <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Questionnaire&Cycle=2013-2014"

# 9. EXTRACT VARIABLE DESCRIPTION TABLES

get_variable_table <- function(url) {
  read_html(url) %>%
    html_elements("table") %>%
    html_table() %>%
    as.data.frame()
}

desc_exam <- bind_rows(
  get_variable_table(url12_exam),
  get_variable_table(url34_exam)
) %>%
  distinct(Variable.Name, .keep_all = TRUE)

desc_demo <- bind_rows(
  get_variable_table(url12_demo),
  get_variable_table(url34_demo)
) %>%
  distinct(Variable.Name, .keep_all = TRUE)

desc_ques <- bind_rows(
  get_variable_table(url12_ques),
  get_variable_table(url34_ques)
) %>%
  distinct(Variable.Name, .keep_all = TRUE)

# 10. COMBINE VARIABLE INFORMATION


var_info <- bind_rows(
  desc_exam,
  desc_demo,
  desc_ques
)

# 11. RESOLVE DUPLICATED OHDEXSTS VARIABLE

# OHDEXSTS appears in more than one merged NHANES component.
# Retain the value from the periodontal examination dataset.

if (
  "OHDEXSTS.x" %in% names(fulldata) &&
  "OHDEXSTS.y" %in% names(fulldata)
) {

  fulldata <- fulldata %>%
    mutate(
      OHDEXSTS = OHDEXSTS.x
    ) %>%
    select(
      -OHDEXSTS.x,
      -OHDEXSTS.y
    )
}

# 12. CHECK VARIABLE DOCUMENTATION

unknown_variables <- setdiff(
  names(fulldata),
  var_info$Variable.Name
)

if (length(unknown_variables) > 0) {
  message(
    "Variables without NHANES variable descriptions: ",
    paste(unknown_variables, collapse = ", ")
  )
}

# Keep documentation only for variables appearing in the dataset.

var_info <- var_info %>%
  filter(
    Variable.Name %in% names(fulldata)
  )

# 13. SAVE DATA

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  var_info,
  here("data", "variable_info.csv"),
  row.names = FALSE
)

write.csv(
  fulldata,
  here("data", "nhanes_full.csv"),
  row.names = FALSE
)

write.csv(
  peri_all,
  here("data", "nhanes_perio.csv"),
  row.names = FALSE
)