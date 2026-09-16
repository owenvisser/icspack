#Packs
ll <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(message("package '",package, "' loaded."))
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}
ll("foreign")
ll("dplyr")
ll("tidyverse")
ll("rvest")
ll("xml2")

rm(list=ls())

#GetData
setwd(r"(C:\Users\owvis\Desktop\nHANES\2011_2012)")
list.files()
demo_12 <- read.xport("DEMO_G.XPT")
peri_12 <- read.xport("OHXPER_G.XPT")
dent_12 <- read.xport("OHXDEN_G.XPT")
smok_cigUse_12 <- read.xport("SMQ_G.XPT")
smok_recent_12 <- read.xport("SMQRTU_G.XPT")

setwd(r"(C:\Users\owvis\Desktop\nHANES\2013_2014)")
demo_34 <- read.xport("DEMO_H.XPT")
peri_34 <- read.xport("OHXPER_H.XPT")
dent_34 <- read.xport("OHXDEN_H.XPT")
smok_cigUse_34 <- read.xport("SMQ_H.XPT")
smok_recent_34 <- read.xport("SMQRTU_H.XPT")

#Check if all data have unique identifiers. (nHanes structure check)
files <- ls()[grep("_", ls())]
checkit <- function(d) length(d$SEQN) == length(unique(d$SEQN))
for(f in files) message(checkit(get(f)))
#Passes.

#MergeData
peri_all <- bind_rows(peri_12, peri_34)
demo_12 <- demo_12 %>% select(-setdiff(names(demo_12), names(demo_34)))
demo_all <- bind_rows(demo_12, demo_34)
dent_all <- bind_rows(dent_12, dent_34)
smok_cigUse_all <- bind_rows(smok_cigUse_12, smok_cigUse_34)
smok_recent_all <- bind_rows(smok_recent_12, smok_recent_34)

fulldata <- peri_all %>%
  full_join(dent_all,        by = "SEQN") %>%
  full_join(demo_all,        by = "SEQN") %>%
  full_join(smok_cigUse_all, by = "SEQN") %>%
  full_join(smok_recent_all, by = "SEQN")

#Clean ls - rm(list = ls()[which(!(ls() %in% ls()[grep("all", ls())]))])

#All survey data is recorded for periodontal data patients.
#Get Variable Info
url12_exam <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Examination&Cycle=2011-2012"
url12_demo <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Demographics&Cycle=2011-2012"
url12_ques <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Questionnaire&Cycle=2011-2012"
url34_exam <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Examination&Cycle=2013-2014"
url34_demo <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Demographics&Cycle=2013-2014"
url34_ques <- "https://wwwn.cdc.gov/nchs/nhanes/search/variablelist.aspx?Component=Questionnaire&Cycle=2013-2014"
gettable <- function(url) read_html(url) %>% html_nodes("table") %>% html_table() %>% as.data.frame()

desc_exam <- bind_rows(gettable(url12_exam), gettable(url34_exam)) %>% distinct(Variable.Name, .keep_all = TRUE)
desc_demo <- bind_rows(gettable(url12_demo), gettable(url34_demo)) %>% distinct(Variable.Name, .keep_all = TRUE)
desc_ques <- bind_rows(gettable(url12_ques), gettable(url34_ques)) %>% distinct(Variable.Name, .keep_all = TRUE)

#Clean ls
rm(list = c(ls()[grep("url", ls())], "gettable"))

#Join all Variable Info
var_info <- desc_exam %>% bind_rows(desc_demo) %>% bind_rows(desc_ques)

#Confirm all Variables Known
#
setdiff(colnames(fulldata), var_info$Variable.Name)
fulldata <- fulldata %>% mutate(OHDEXSTS = OHDEXSTS.x) %>% select(-OHDEXSTS.x, -OHDEXSTS.y)
var_info <- var_info %>% filter(Variable.Name %in% colnames(fulldata)) 

#Clean ls
rm(list = ls()[grep("desc", ls())])

#Write Data to .csv
setwd(r"(C:\Users\owvis\Desktop\nHANES\Full Data)")
write.csv(var_info, file = "variable.info.csv")
write.csv(fulldata, file = "All_UID.csv")
write.csv(peri_all, file = "Perio_UID.csv")

#Done. 