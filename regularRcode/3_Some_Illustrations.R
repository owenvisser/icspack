#Illustration
#Some summary statistics with the data we have.
setwd("C:/Users/owvis/OneDrive - University of Florida/nHANES/Full Data")
load("organized_data.RData")

setwd("C:/Users/owvis/OneDrive - University of Florida/nHANES/Full Data")
df <- read.csv("All_UID.csv")[,-1] # only people with some periodontal data


## Looking into Ages
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(tidyr)

anyDuplicated(names(df))

df <- df %>% filter(!if_all(everything(), is.na)) %>% mutate(across(everything(), as.factor)) 

#Patient ages
hist(as.numeric(df$Age),
     main = 'Frequency of Age for Dental Patients with Perio Data',
     xlab = 'Age')

age_groups <- cut(as.numeric(df$Age),
                  breaks = c(34, 39, 49, 59, 64, 69),
                  labels = c("35-39", "40-49", "50-59", "60-64", "65-69"),
                  include.lowest = TRUE)
age_groups <- na.omit(age_groups)
table(age_groups)
Full_age_groups <- cut(as.numeric(df$Age),
                      breaks = c(1, 15, 29, 39, 49, 59, 69, 79, 89),
                      labels = c( "2-15", "16-29", "30-39",
                                  "40-49", "50-59", "60-69", "70-79", "80-89"),
                      include.lowest = TRUE)
Full_age_groups <- na.omit(Full_age_groups)
table(Full_age_groups)
# Plot histogram of age groups
ggplot(data = data.frame(AgeGroup = age_groups), aes(x = AgeGroup)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Histogram of Age Groups", x = "Age Group", y = "Count") +
  theme_minimal()

# Create a table of counts
age_group_counts <- table(age_groups)

# Print the table
print(age_group_counts)

## Looking into Tooth Counts - ignoring wisdom teeth
tc <- paste0("OHX", sprintf("%02d", c(1:28)), "TC")
tc_df <- df %>% select(all_of(tc))

# Checking toothcount distribution
# labels = c('1 = Primary tooth', '2 = Permanent tooth', '3 = Implant', '4 = Tooth not present', '5 = Permanent dental root fragment'))
toothdist <- tc_df %>%
  mutate(Perm_Count = rowSums(across(everything(), ~ .x %in% c(1,2,3,5)), na.rm = T))
toothdist <- toothdist$Perm_Count
toothdist <- toothdist[!is.na(toothdist)]
hist(toothdist, breaks = 32)


tc_df1 <- lapply(tc_df, function(x) factor(x, levels = 1:5))
table_tc <- sapply(tc_df1, table)
table_tc <- as.data.frame(table_tc)
rownames(table_tc) <- 1:5

long_df <- table_tc %>%
  rownames_to_column(var = "Level") %>%
  pivot_longer(cols = -all_of('Level'), names_to = "Section", values_to = "Counts") %>%
  group_by(Level) %>%
  summarize(TotalCounts = sum(Counts))
# Plotting total counts per level

ggplot(long_df, aes(x = Level, y = TotalCounts, fill = Level)) +
  geom_bar(stat = "identity") +
  labs(x = "Level", y = "Total Counts", title = "Total Counts per Level Across All Teeth") +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set1"), 
                    labels = c('1 = Primary tooth', '2 = Permanent tooth', '3 = Implant', '4 = Tooth not present', '5 = Permanent dental root fragment'))+
  guides(fill = guide_legend(ncol = 3))

# 1 = Primary tooth
# 2 = Permanent tooth
# 3 = Implant
# 4 = Tooth not present
# 5 = Permanent dental root fragment

long_df2 <- pivot_longer(table_tc, cols = starts_with("OHX"), names_to = "Section", values_to = "Counts")
long_df2$Level <- as.factor(rep(1:5, each = 32))


ggplot(long_df2, aes(x = Section, y = Counts, fill = Level)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "free_y", ncol = 1) +  # One column, separate plot per level
  labs(x = "Tooth", y = "Counts", title = "Counts of Individual Tooth Level (free scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom',
        panel.spacing = unit(0.1, 'lines')) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set1"), 
                    labels = c('1 = Primary tooth', '2 = Permanent tooth', '3 = Implant', '4 = Tooth not present', '5 = Permanent dental root fragment'))+
  guides(fill = guide_legend(ncol = 3))

ggplot(long_df2, aes(x = Section, y = Counts, fill = Level)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed", ncol = 1) +  # One column, separate plot per level
  labs(x = "Tooth", y = "Counts", title = "Counts of Individual Tooth Level (fixed scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'bottom',
        panel.spacing = unit(0.1, 'lines')) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Set1"), 
                    labels = c('1 = Primary tooth', '2 = Permanent tooth', '3 = Implant', '4 = Tooth not present', '5 = Permanent dental root fragment'))+
  guides(fill = guide_legend(ncol = 3))





##################################### More stuff. 

loa_measures <- c('CJA', 'CJD', 'CJL', 'CJM', 'CJP', 'CJS')
tooth_counts <- c('TC')
coronal_info <- c('CSC','CTC')


summary(df_bytooth$OHX02, maxsum = 100)


#### UTILIZING GPT TO PRINT: CSC

supertab <- lapply(df_bytooth, function(x){
  summary(x %>% select(all_of(grep('CSC', colnames(x)))), maxsum = 100)
})

st2 <- do.call(cbind, supertab[1:28])

# Create the data frame
data <- data.frame(
  Tooth = c("OHX02CSC", "OHX03CSC", "OHX04CSC", "OHX05CSC", "OHX06CSC",
            "OHX07CSC", "OHX08CSC", "OHX09CSC", "OHX10CSC", "OHX11CSC",
            "OHX12CSC", "OHX13CSC", "OHX14CSC", "OHX15CSC", "OHX18CSC",
            "OHX19CSC", "OHX20CSC", "OHX21CSC", "OHX22CSC", "OHX23CSC",
            "OHX24CSC", "OHX25CSC", "OHX26CSC", "OHX27CSC", "OHX28CSC",
            "OHX29CSC", "OHX30CSC", "OHX31CSC"),
  `Decay only` = c(237, 215, 184, 206, 268,
                   304, 318, 318, 293, 255,
                   178, 173, 211, 296, 275,
                   191, 216, 184, 182, 110,
                   84, 110, 122, 190, 175,
                   219, 222, 294),
  `Filling only` = c(3112, 3508, 2325, 2049, 717,
                     963, 1014, 1007, 1050, 774,
                     2013, 2348, 3543, 3125, 3273,
                     3222, 2264, 1292, 342, 191,
                     141, 154, 166, 324, 1322,
                     2278, 3241, 3346),
  `Fillings and Decay` = c(25, 41, 19, 27, 14,
                           18, 24, 20, 20, 17,
                           24, 20, 43, 23, 57,
                           62, 18, 14, 7, 3,
                           2, 2, 3, 5, 20,
                           34, 65, 60),
  `NA's` = c(5660, 5270, 6506, 6752, 8035,
             7749, 7678, 7689, 7671, 7988,
             6819, 6493, 5237, 5590, 5429,
             5559, 6536, 7544, 8503, 8730,
             8807, 8768, 8743, 8515, 7517,
             6503, 5506, 5334)
)

# Convert to long format for easier plotting
data_long <- data %>%
  pivot_longer(cols = -Tooth, names_to = "Condition", values_to = "Count")

# Plot using ggplot2
ggplot(data_long, aes(x = Tooth, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Oral Health Data by Tooth",
       x = "Tooth",
       y = "Count",
       fill = "Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ONLY WHEN ABOVE IS "Z"
# For caries, the allowable codes are as follows:
# 0 = Lingual caries
# 1 = Occlusal caries
# 2 = Facial caries
# 3 = Mesial caries
# 4 = Distal caries
# For filled teeth or restorations, the allowable surface codes are as follows:
# 5 = Lingual restoration
# 6 = Occlusal restoration
# 7 = Facial restoration
# 8 = Mesial restoration
# 9 = Distal restoration
# C = Crown (short call for both primary and permanent teeth)


#### UTILIZING GPT TO PRINT: CTC

supertab <- lapply(df_bytooth, function(x){
  summary(x %>% select(all_of(grep('CTC', colnames(x)))), maxsum = 100)
})

do.call(cbind, supertab[1:28])

# Create a data frame containing the new data
data2 <- data.frame(
  Tooth = c("OHX02CTC", "OHX03CTC", "OHX04CTC", "OHX05CTC", "OHX06CTC", "OHX07CTC", "OHX08CTC", "OHX09CTC", "OHX10CTC", "OHX11CTC",
            "OHX12CTC", "OHX13CTC", "OHX14CTC", "OHX15CTC", "OHX18CTC", "OHX19CTC", "OHX20CTC", "OHX21CTC", "OHX22CTC", "OHX23CTC",
            "OHX24CTC", "OHX25CTC", "OHX26CTC", "OHX27CTC", "OHX28CTC", "OHX29CTC", "OHX30CTC", "OHX31CTC"),
  S = c(2073, 1664, 3250, 3397, 5474, 4799, 4623, 4653, 4772, 5457,
        3459, 3221, 1623, 1993, 1527, 1320, 3839, 5189, 6760, 6791,
        6744, 6696, 6770, 6739, 5159, 3782, 1309, 1515),
  Z = c(3374, 3764, 2526, 2282, 997, 1285, 1356, 1345, 1363, 1045,
        2215, 2540, 3796, 3444, 3605, 3475, 2490, 1490, 531, 304,
        227, 266, 291, 519, 1517, 2523, 3528, 3700),
  D = c(0, 0, 2, 0, 10, 1, 0, 0, 2, 1,
        0, 3, 0, 0, 0, 0, 11, 0, 2, 0,
        1, 1, 0, 2, 0, 11, 0, 0),
  K = c(0, 0, 2, 0, 2, 0, 0, 0, 0, 1,
        0, 1, 0, 0, 0, 0, 8, 0, 0, 0,
        0, 0, 0, 0, 0, 8, 0, 0),
  U = c(3, 0, 1, 3, 7, 16, 0, 0, 12, 8,
        3, 3, 0, 1, 0, 0, 4, 4, 0, 3,
        4, 2, 4, 1, 0, 3, 0, 3),
  E = c(1420, 1209, 822, 606, 243, 309, 290, 274, 316, 274,
        622, 820, 1215, 1428, 1931, 1973, 750, 394, 189, 233,
        284, 286, 254, 204, 388, 756, 1924, 1840),
  M = c(3, 5, 11, 398, 27, 43, 18, 16, 35, 22,
        370, 6, 4, 4, 3, 3, 15, 294, 4, 19,
        23, 30, 26, 9, 295, 9, 1, 5),
  R = c(36, 142, 179, 156, 44, 130, 149, 149, 128, 41,
        143, 189, 145, 21, 62, 270, 118, 33, 8, 51,
        76, 75, 48, 10, 43, 117, 284, 71),
  X = c(1, 2, 3, 4, 7, 33, 45, 38, 27, 3,
        3, 1, 1, 0, 1, 5, 2, 0, 1, 8,
        11, 14, 8, 1, 1, 3, 1, 0),
  P = c(1254, 1351, 1334, 1288, 1146, 1266, 1302, 1287, 1253, 1128,
        1300, 1326, 1371, 1261, 1034, 1099, 937, 815, 680, 805,
        835, 842, 805, 687, 815, 961, 1095, 1027),
  Q = c(34, 35, 37, 37, 36, 44, 54, 49, 42, 34,
        36, 34, 36, 35, 23, 22, 21, 21, 21, 25,
        29, 28, 26, 19, 20, 21, 21, 23),
  NAs = c(836, 862, 867, 863, 1041, 1108, 1197, 1223, 1084, 1020,
          883, 890, 843, 847, 848, 867, 839, 794, 838, 795,
          800, 794, 802, 843, 796, 840, 871, 850)
)

# Convert to long format for easier plotting
data2_long <- data2 %>%
  pivot_longer(cols = -Tooth, names_to = "Condition", values_to = "Count")

# Plot using ggplot2
ggplot(data2_long, aes(x = Tooth, y = Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Counts by Tooth and Condition",
       x = "Tooth",
       y = "Count",
       fill = "Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = 'bottom',
      panel.spacing = unit(0.1, 'lines')) +
  scale_fill_manual(values = brewer.pal(n = 12, name = "Paired"), 
                    labels = c(
                      'S = Sound permanent tooth (no decay or filling on any surface)',
                      'Z = Permanent tooth with surface condition *gives CSC plot',
                      'D = Sound primary (deciduous) tooth',
                      'K = Primary tooth with surface condition',
                      'U = Unerupted tooth',
                      'E = Missing due to dental disease (caries/periodontal disease)',
                      'M = Missing due to other causes (orthodontic/traumatic or other nondisease)',
                      'R = Missing due to dental disease but replaced by a fixed restoration',
                      'X = Missing due to other causes but replaced by a fixed restoration',
                      'P = Missing due to dental disease but replaced by a removable restoration',
                      'Q = Missing due to other causes but replaced by a removable restoration',
                      'NA = not available'
                    ))+
  guides(fill = guide_legend(ncol =3))

#Coronal Stuff:
# S = Sound permanent tooth (no decay or filling on any surface)
# Z = Permanent tooth with surface condition    *****This marks for the below section*****
# D = Sound primary (deciduous) tooth
# K = Primary tooth with surface condition
# U = Unerupted tooth
# E = Missing due to dental disease (caries/periodontal disease)
# M = Missing due to other causes (orthodontic/traumatic or other nondisease)
# R = Missing due to dental disease but replaced by a fixed restoration
# X = Missing due to other causes but replaced by a fixed restoration
# P = Missing due to dental disease but replaced by a removable restoration
# Q = Missing due to other causes but replaced by a removable restoration


#### UTILIZING GPT TO PRINT: CTC
loa_measures <- c('CJA', 'CJD', 'CJL', 'CJM', 'CJP', 'CJS')



supertab <- lapply(df_bytooth, function(x){
  presence_vector <- sapply(names(x), function(string) {
    any(sapply(loa_measures, function(loa_measures) grepl(loa_measures, string)))
  })
  data.frame(unclass(summary(x %>% select(all_of(which(presence_vector))), maxsum = 100)), check.names = FALSE, stringsAsFactors = FALSE)
})

data.frame(unclass(supertab), check.names = FALSE, stringsAsFactors = FALSE)
names(supertab) <- NULL
do.call(cbind, supertab[1:28])

supertab[[1]]$OHX02CJA


data_list <- list(
  OHX02CJA = c("-12 :   0  ", "-11 :   1  ", "-10 :   0  ", "-9  :   0  ", "-8  :   3  ",
               "-7  :   3  ", "-6  :   9  ", "-5  :  19  ", "-4  :  46  ", "-3  :  98  ",
               "-2  : 209  ", "-1  : 273  ", "0   :4107  ", "1   : 237  ", "2   :  26  ",
               "3   :   1  ", "4   :   0  ", "5   :   0  ", "6   :   0  ", "7   :   0  ",
               "8   :   0  ", "9   :   0  ", "10  :   0  ", "11  :   0  ", "12  :   0  ",
               "99  :2720  ", "NA's:1282  ")
)

# Extract Column Guides
get_guide <- function(column) {
  c <- sapply(column, function(x) strsplit(x, ":")[[1]][1])
  return(unlist(c))
}
# Extract Column Values
get_values <- function(column) {
  sapply(column, function(x) as.numeric(trimws(strsplit(trimws(x), ":")[[1]][2])))
}

# Create Data Frame by Extracting Guides
measure_summary <- do.call(cbind, lapply(supertab[1:28], function(y){
  values <- as.data.frame(sapply(y, get_values))
  rownames(values) <- get_guide(rownames(values))
  return(values)
}))

xy <- rowSums(measure_summary)
row_sums_df <- data.frame(Guide = names(xy), Count = as.numeric(xy))

# Load ggplot2
library(ggplot2)

# Plotting with ggplot2
ggplot(row_sums_df, aes(x = Guide, y = Count)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(
    title = "Sum of Measurements by Guide",
    x = "Guide",
    y = "Sum of Measurements"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))



#Allowable Entries for Loss of Attachment: 
# For each tooth site (distal, mid-facial, or mesio-facial) the following two spaces must be recorded:
# First Data Entry Space (FGM to CEJ measurement)
# -9 to 9 = Measurement in millimeters
# +A = Measurement is +10 millimeters
# +B = Measurement is +11 millimeters
# +C = Measurement is +12 millimeters
# 99 = Cannot be assessed
# Second Data Entry Space (FGM to sulcus base measurement)
# 0 – 9 = Measurement in millimeters
# A = Measurement is 10 millimeters
# B = Measurement is 11 millimeters
# C = Measurement is 12 millimeters
# 99 = Cannot be assessed


library(gtsummary)






# table(age_groups)

teeth_binary_df <- data.frame(
  apply(peri.uid %>% select(unlist(as.vector(teeth_list %>% select(TC)))), 2, function(x) {
    sapply(x, function(y){
      if (y %in% c(2, 3, 5)) {
        return(1)
      } else if(is.na(y)) {
        return(NA)
      } else {
        return(0)
      }})})) %>%
  setNames(unlist(as.vector(teeth_list %>% select(TC))))

teeth_df <- cbind(ages = df$Age, teeth_binary_df)

# Mapping of Universal Numbers to scientific names (partial example)
scientific_mapping <- c(
  "01" = "01 Maxillary Right Third Molar",
  "02" = "02 Maxillary Right Second Molar",
  "03" = "03 Maxillary Right First Molar",
  "04" = "04 Maxillary Right Second Premolar",
  "05" = "05 Maxillary Right First Premolar",
  "06" = "06 Maxillary Right Canine",
  "07" = "07 Maxillary Right Lateral Incisor",
  "08" = "08 Maxillary Right Central Incisor",
  "09" = "09 Maxillary Left Central Incisor",
  "10" = "10 Maxillary Left Lateral Incisor",
  "11" = "11 Maxillary Left Canine",
  "12" = "12 Maxillary Left First Premolar",
  "13" = "13 Maxillary Left Second Premolar",
  "14" = "14 Maxillary Left First Molar",
  "15" = "15 Maxillary Left Second Molar",
  "16" = "16 Maxillary Left Third Molar",
  "17" = "17 Mandibular Left Third Molar",
  "18" = "18 Mandibular Left Second Molar",
  "19" = "19 Mandibular Left First Molar",
  "20" = "20 Mandibular Left Second Premolar",
  "21" = "21 Mandibular Left First Premolar",
  "22" = "22 Mandibular Left Canine",
  "23" = "23 Mandibular Left Lateral Incisor",
  "24" = "24 Mandibular Left Central Incisor",
  "25" = "25 Mandibular Right Central Incisor",
  "26" = "26 Mandibular Right Lateral Incisor",
  "27" = "27 Mandibular Right Canine",
  "28" = "28 Mandibular Right First Premolar",
  "29" = "29 Mandibular Right Second Premolar",
  "30" = "30 Mandibular Right First Molar",
  "31" = "31 Mandibular Right Second Molar",
  "32" = "32 Mandibular Right Third Molar"
)

# Reshape the data to a long format
long_teeth_df <- teeth_df %>%
  pivot_longer(
    cols = starts_with("OHX"),  # Select columns that start with "OHX"
    names_to = "Tooth",
    values_to = "Presence"
  )

#Henry - denominator should not include the missing data. - just removed missing woo-hew!

# Extract the tooth number from the column name and apply scientific names
long_teeth_df <- long_teeth_df %>%
  mutate(Tooth = sub("OHX([0-9]+)TC", "\\1", Tooth)) %>%
  mutate(Tooth = scientific_mapping[Tooth])

# Calculate the percentage of presence for each Tooth-Age group combination
heatmap_data <- long_teeth_df %>%
  group_by(ages, Tooth) %>%
  summarize(Presence_Percentage = mean(Presence) * 100, .groups = "drop")

min_percentage <- min(heatmap_data$Presence_Percentage, na.rm = TRUE)
max_percentage <- max(heatmap_data$Presence_Percentage, na.rm = TRUE)

# Create the heatmap with ggplot2
heatmap_plot <- ggplot(heatmap_data, aes(x = ages, y = Tooth, fill = Presence_Percentage)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black",
                      name = paste0("Percentage\n(",
                                    round(min_percentage, 1),
                                    "%-", round(max_percentage, 1),
                                    "%)")) +
  labs(
    title = "Tooth Presence Percentage by Age Heatmap",
    x = "Ages",
    y = "Tooth"
  ) +
  theme_minimal()+
  theme(
    axis.text.y = element_text(hjust = 0)
  )

# Print the heatmap
print(heatmap_plot)


# Allowable Entries for Tooth Count:
# 1 = Primary tooth (deciduous)
# 2 = Permanent tooth
# 3 = Dental Implant
# 4 = Tooth not present
# 5 = Permanent dental root fragment


#scatterplot - heatmap
#
teeth_df2 <- teeth_df %>% mutate(tc = rowSums(select(., -all_of('ages')), na.rm = TRUE))
teeth_df2$Age_Groups <- age_groups

ggplot(teeth_df2, aes(x = as.factor(ages), y = as.factor(tc))) +
  geom_bin2d(binwidth = c(1, 1)) +  # Adjust the 'bins' parameter as needed
  scale_fill_gradient(low = "grey", high = "black") +  # Adjust colors
  theme_minimal() +
  labs(title = "Scattered Heatmap of Tooth Count Density by Individual Ages",
       x = bquote(bold('Age')~'(Years)'),
       y = bquote(bold('Count of Teeth') ~'(permanent, deciduous, dental implant)'),
       fill = 'Participant\nCount')+
  theme(
    axis.text.x = element_text(size = 7),  # Adjust x-axis text size
    axis.text.y = element_text(size = 7)   # Adjust y-axis text size
  )

ggplot(teeth_df2, aes(x = Age_Groups, y = as.factor(tc))) +
  geom_bin2d(binwidth = c(1, 1)) +  # Adjust the 'bins' parameter as needed
  scale_fill_gradient(low = "grey", high = "black") +  # Adjust colors
  theme_minimal() +
  labs(title = "Scattered Heatmap of Tooth Count Density by Age Groups",
       x = bquote(bold('Age')~'(Years)'),
       y = bquote(bold('Count of Teeth') ~'(permanent, deciduous, dental implant)'),
       fill = 'Participant\nCount')





# looking at clinical attachment loss average per tooth. 
CLA_df <- df %>% select(unlist(teeth_list %>% select(c('LAA', 'LAD', 'LAL', 'LAM', 'LAP', 'LAS')))) %>%
  mutate(across(everything(), ~ factor(.x, levels = c(c(-12:12), 99)))) %>%
  setNames(unlist(teeth_list %>% select(c('LAA', 'LAD', 'LAL', 'LAM', 'LAP', 'LAS'))))


ohx.nums <- sapply(c(2:15, 18:31), function(x){
  ifelse(x < 10,
         ret <- paste0('OHX0',x),
         ret <- paste0('OHX', x))
  return(ret)})

CLA_df <- cbind(ages = df$Age, # can replace with age_groups or df$Age
                as.data.frame(sapply(ohx.nums, function(x){
                  CLA_df %>% select(all_of(grep(pattern = x, colnames(CLA_df)))) %>%
                    apply(2, as.numeric) %>%
                    apply(1, function(y){
                      if(all(is.na(y))){ #if all NA, return NA
                        return(NA) 
                      } else if(!any(y == 99, na.rm = TRUE)) {#if there are no 99s give mean - ignore NA
                        return(mean(y, na.rm = TRUE))
                      } else if(all(y==99, na.rm = TRUE)) { #if there is only 99s give 99 - ignore NA
                        return(99)
                      } else {
                        return(mean(y[y!=99], na.rm = TRUE)) #if there is 99 and regular give mean - ignore NA and 99
                      }
                    }) %>% as.data.frame %>% setNames(NULL)})
                ))


# Reshape the data to a long format
long_cla_df <- CLA_df %>%
  pivot_longer(
    cols = starts_with("OHX"),  # Select columns that start with "OHX"
    names_to = "Tooth",
    values_to = "avgCLA"
  )

#Henry - denominator should not include the missing data. - just removed missing woo-hew!

# Extract the tooth number from the column name and apply scientific names
long_cla_df <- long_cla_df %>%
  mutate(Tooth = sub("OHX([0-9]+)", "\\1", Tooth)) %>%
  mutate(Tooth = scientific_mapping[Tooth])

# Calculate the percentage of presence for each Tooth-Age group combination
heatmap_data <- long_cla_df %>%
  group_by(ages, Tooth) %>% filter(avgCLA != 99) %>%
  summarize(super_avg_CLA = mean(avgCLA), .groups = "drop")

# Create the heatmap with ggplot2
heatmap_plot <- ggplot(heatmap_data, aes(x = ages, y = Tooth, fill = super_avg_CLA)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black",
                      name = "Distance(mm)") +
  labs(
    title = "Average Clinical Loss of Attachment per Tooth by Age Heatmap",
    x = "Ages",
    y = "Tooth"
  ) +
  theme_minimal()+
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(size = 10)
  )

# Print the heatmap
print(heatmap_plot)


