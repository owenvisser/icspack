### file for calculating the correlation values for Nhanes. 

# On laptop
load("C:/Users/owvis/OneDrive - University of Florida/nHANES/Full Data/Tooth_Level_Data.RData")
library(dplyr)
source("C:/Users/owvis/OneDrive - University of Florida/nHANES/Code/Weighting_Equations.R")
source("C:/Users/owvis/OneDrive - University of Florida/nHANES/Code/Correlation_Equations.R")


# # On HiperGator
# setwd("/orange/somnath.datta/NHANES/")
# load("FullData/Tooth_Level_Data.RData")
# library(dplyr)
# source("/blue/somnath.datta/NHANES/Weighting_Equations.R")
# source("/blue/somnath.datta/NHANES/Correlation_Equations.R")

#drop rows with NA values.

df_clean <- na.omit(df_tooth_level)

df_clean <- df_clean %>%
  group_by(SEQN) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(Tooth_Count > 9)

df_clean <- na.omit(df_clean)


df_clean <- df_clean %>%
  mutate(`Cal>=3` = ifelse(maxCAL >= 3, 1 , 0),
         `Cal>=4` = ifelse(maxCAL >= 4, 1 , 0),
         `Cal>=5` = ifelse(maxCAL >= 5, 1 , 0),
         `Pd>=4` = ifelse(maxPD >= 4, 1 , 0),
         `Pd>=5` = ifelse(maxPD >= 5, 1 , 0),
         `Pd>=6` = ifelse(maxPD >= 6, 1 , 0),
         )

# Setup of the grid of Y and X to use for association evaluation

# For spearman Coeff.
# Perio
Y_vars <- c("maxCAL", "maxCAL", "maxCAL", "maxPD", "maxPD", "maxPD" )
L_vars <- c("Cal>=3", "Cal>=4", "Cal>=5", "Pd>=4", "Pd>=5", "Pd>=6" )
# Caries
X_vars <- c("DS", "FS", "DSI", "FSI") #c("DS", "FS", "DFS", "DSI", "FSI", "DFSI")
K_vars <- X_vars

# For all weights
weight_types <- c("no_weight", "CW", "PPW", "OPW", "MOPW")

all_spearman_results <- do.call(
  rbind,
  lapply(weight_types, function(wtype) {
    
    run_association_grid(
      dat = df_clean,
      Y_vars = Y_vars,
      X_vars = X_vars,
      clusterID = "SEQN",
      clusterSize = "Tooth_Count",
      K_vars = K_vars,
      L_vars = L_vars,
      weight_type = wtype,
      association_type = "spearman"
    )
  })
)
all_spearman_results



# For phi Coeff.
# Perio
Y_vars <- c("Cal>=3", "Cal>=4", "Cal>=5", "Pd>=4", "Pd>=5", "Pd>=6" )
L_vars <- Y_vars
# Caries
X_vars <- c("DS", "FS", "DSI", "FSI") #c("DS", "FS", "DFS", "DSI", "FSI", "DFSI")
K_vars <- X_vars

# For all weights
weight_types <- c("no_weight", "CW", "PPW", "OPW", "MOPW")

all_phi_results <- do.call(
  rbind,
  lapply(weight_types, function(wtype) {
    
    run_association_grid(
      dat = df_clean,
      Y_vars = Y_vars,
      X_vars = X_vars,
      clusterID = "SEQN",
      clusterSize = "Tooth_Count",
      K_vars = K_vars,
      L_vars = L_vars,
      weight_type = wtype,
      association_type = "phi"
    )
  })
)

all_phi_results










library(dplyr)
library(tidyr)
library(stringr)

weight_order <- c("no_weight", "CW", "PPW", "OPW", "MOPW")
K_order <- c("FS", "FSI", "DS", "DSI")

wide_tab <- all_phi_results %>%
  select(L, K, weight_type, estimate, se) %>%
  rename(est = estimate) %>%
  mutate(
    weight_type = factor(weight_type, levels = weight_order),
    K = factor(K, levels = K_order),
    
    L_group = case_when(
      str_detect(L, "Cal") ~ "CAL",
      str_detect(L, "Pd")  ~ "PD",
      TRUE ~ NA_character_
    ),
    
    L_num = as.numeric(str_extract(L, "\\d+")),
    
    L_group = factor(L_group, levels = c("CAL", "PD"))
  ) %>%
  arrange(K, L_group, L_num) %>%
  pivot_wider(
    names_from = weight_type,
    values_from = c(est, se),
    names_glue = "{weight_type}_{.value}",
    names_vary = "slowest"
  ) %>%
  arrange(K, L_group, L_num)

print(wide_tab, n = 100)

fmt_est <- function(est, se) {
  ifelse(
    is.na(est) | is.na(se),
    "-",
    sprintf("$%.3f\\,(%.3f)$", est, se)
  )
}

latex_rows <- wide_tab %>%
  group_by(K) %>%
  mutate(
    K_row = row_number(),
    K_tex = ifelse(K_row == 1, as.character(K), "   ")
  ) %>%
  ungroup() %>%
  group_by(K, L_group) %>%
  mutate(
    within_group_row = row_number(),
    
    no_weight_tex = ifelse(
      within_group_row == 1,
      fmt_est(no_weight_est, no_weight_se),
      "-"
    ),
    
    CW_tex = ifelse(
      within_group_row == 1,
      fmt_est(CW_est, CW_se),
      "-"
    ),
    
    PPW_tex  = fmt_est(PPW_est, PPW_se),
    OPW_tex  = fmt_est(OPW_est, OPW_se),
    MOPW_tex = fmt_est(MOPW_est, MOPW_se)
  ) %>%
  ungroup() %>%
  mutate(
    L_tex = case_when(
      L_group == "CAL" ~ paste0("CAL$\\geq$", L_num),
      L_group == "PD"  ~ paste0("PD$\\geq$", L_num),
      TRUE ~ as.character(L)
    ),
    
    row_tex = sprintf(
      "%s & %s && %s & %s & %s & %s & %s \\\\",
      K_tex, L_tex,
      no_weight_tex,
      CW_tex,
      PPW_tex,
      OPW_tex,
      MOPW_tex
    )
  ) %>%
  group_by(K) %>%
  mutate(
    row_tex = ifelse(
      L_group == "PD" & row_number() == which(L_group == "PD")[1],
      paste0("   \\\\\n", row_tex),
      row_tex
    )
  ) %>%
  ungroup() %>%
  group_by(K) %>%
  summarise(
    block_tex = paste(row_tex, collapse = "\n"),
    .groups = "drop"
  ) %>%
  pull(block_tex)

cat(paste(latex_rows, collapse = "\n\\\\\n"))






# # For particular 
# # 
# pearson_results <- run_association_grid(
#   dat = dat,
#   Y_vars = Y_vars,
#   X_vars = X_vars,
#   clusterID = "clusterID",
#   clusterSize = "clusterSize",
#   K = "K",
#   L = "L",
#   weight_type = "CW",
#   association_type = "pearson"
# )
# 


# 
# all_pearson_results <- do.call(
#   rbind,
#   lapply(weight_types, function(wtype) {
#     
#     run_association_grid(
#       dat = dat,
#       Y_vars = Y_vars,
#       X_vars = X_vars,
#       clusterID = "clusterID",
#       clusterSize = "clusterSize",
#       K = "K",
#       L = "L",
#       weight_type = wtype,
#       association_type = "pearson"
#     )
#   })
# )
# 
# all_spearman_results <- do.call(
#   rbind,
#   lapply(weight_types, function(wtype) {
#     
#     run_association_grid(
#       dat = dat,
#       Y_vars = Y_vars,
#       X_vars = X_vars,
#       clusterID = "clusterID",
#       clusterSize = "clusterSize",
#       K = "K",
#       L = "L",
#       weight_type = wtype,
#       association_type = "spearman"
#     )
#   })
# )
# 
# all_phi_results <- do.call(
#   rbind,
#   lapply(weight_types, function(wtype) {
#     
#     run_association_grid(
#       dat = df_clean,
#       Y_vars = Y_vars,
#       X_vars = X_vars,
#       clusterID = "SEQN",
#       clusterSize = "Tooth_Count",
#       K_vars = K_vars,
#       L_vars = L_vars,
#       weight_type = wtype,
#       association_type = "phi"
#     )
#   })
# )








# 
# # Examples of making the dat frame for correlation.
# dat <- data.frame(
#   clusterID = df_clean$SEQN,
#   clusterSize = df_clean$Tooth_Count,
#   Y = df_clean$FS,
#   L = df_clean$FS,
#   X = df_clean$maxCAL,
#   K = df_clean$maxCAL
# )
# 
# w <- make_weights(
#   clusterID = dat$clusterID,
#   clusterSize = dat$clusterSize,
#   K = dat$K,
#   L = dat$L,
#   weight_type = "MOPW"
# )
# 
# pearson_association(
#   clusterID = dat$clusterID,
#   Y = dat$Y,
#   X = dat$X,
#   omega = w
# )
# 
# spearman_association(
#   clusterID = dat$clusterID,
#   Y = dat$Y,
#   X = dat$X,
#   omega = w
# )
# 
# df_clean2 <- df_clean %>%
#   mutate(`Cal>=3` = ifelse(maxCAL >= 3, 1 , 0),
#          `Cal>=4` = ifelse(maxCAL >= 4, 1 , 0),
#          `Cal>=5` = ifelse(maxCAL >= 5, 1 , 0),
#          `Pd>=4` = ifelse(maxPD >= 4, 1 , 0),
#          `Pd>=5` = ifelse(maxPD >= 5, 1 , 0),
#          `Pd>=6` = ifelse(maxPD >= 6, 1 , 0),
#          )
# 
# dat <- data.frame(
#   clusterID = df_clean2$SEQN,
#   clusterSize = df_clean2$Tooth_Count,
#   Y = df_clean2$FS,
#   L = df_clean2$FS,
#   X = df_clean2$`Cal>=3`,
#   K = df_clean2$`Cal>=3`
# )
# 
# phi_association(
#   clusterID = dat$clusterID,
#   Y = dat$Y,
#   X = dat$X,
#   omega = w
# )