library(haven)
library(readxl)
library(dplyr)
library(ggplot2)
library(knitr)
library(stringr)
library(purrr)
library(lubridate)
library(readr)
library(tidyr)
library(haven)
library(lme4)
library(Matrix) # dependency for lme4
library(tableone)
library(kableExtra)
library(xtable)
library(gee)
library(tidyverse)
library(naniar)
library(fastDummies)


#Site names
site_names <- c("Ghana", "India-CMC", "India-SAS", "Kenya", "Pakistan", "Zambia")

###Import data
UploadDate = "2024-06-28"
path_to_data <- paste0('D:/Users/yipeng_wei/Documents/Stacked data/',UploadDate)# - for AWS data

uploadDate <- as.Date(UploadDate, format = "%Y-%m-%d")
setwd(path_to_data)

###Load data
########################################################################################
recreate_mnh_list <- FALSE #if I want to recreate the big data list then set this to TRUE
########################################################################################
tic <- Sys.time()
if(recreate_mnh_list == TRUE){
  mnh_names <- c()
  list_mnh <- dir(path = path_to_data, pattern = "*.csv", full.names = TRUE) #creates a list of all the csv files in the directory
  for (data_file in list_mnh[]) { #can test by just bringing in a small number (add 1:2 inside the bracket to do so)
    form_name <- substr(basename(data_file), 1,5) #substr pulls out the 1:5 spaces in a char (will pull out "mnh00" etc);
    #basename() pulls out just the name of the file from the entire directory/path.
    print(paste("Reading", form_name))
    assign(form_name, read.csv(data_file))
    mnh_names <- c(mnh_names, form_name) # it's buidling up a vector of obj names.
  }
  save(list=mnh_names, file = paste0("MNH_dfs_", basename(path_to_data), ".RDS"))
} else{
  load(file = paste0('MNH_dfs_',  basename(path_to_data), ".RDS"))
}

toc <- Sys.time()
print(toc - tic)

###Create a vector of women who are enrolled. 
###Get ReMAPP full cohort population
mnh02 <- mnh02 %>% 
  mutate(ENROLLED_FF = ifelse(M02_AGE_IEORRES==1 & M02_PC_IEORRES==1 & 
                                M02_CATCHMENT_IEORRES==1 & M02_CATCH_REMAIN_IEORRES==1 & 
                                M02_CONSENT_IEORRES==1, 1,0),
         remapp = case_when((SITE == "Ghana" & M02_SCRN_OBSSTDAT >= "2022-12-28") |
                              (SITE == "Kenya" & M02_SCRN_OBSSTDAT >= "2023-04-14") |
                              (SITE == "Zambia" & M02_SCRN_OBSSTDAT >= "2022-12-15") |
                              (SITE == "Pakistan" & M02_SCRN_OBSSTDAT >= "2022-09-22" & M02_SCRN_OBSSTDAT <= "2024-04-05") |
                              (SITE == "India-CMC" & M02_SCRN_OBSSTDAT >= "2023-06-20") |
                              (SITE == "India-SAS" & M02_SCRN_OBSSTDAT >= "2023-08-15") ~ 1,
                            TRUE ~ NA_real_)) %>% filter(ENROLLED_FF ==1 & remapp==1) # THIS DOES DROP SOME WOMEN!

ReMAPP.id<-mnh02 %>% select("SITE", "MOMID", "PREGID")

# Assuming ReMAPP.id already has unique MOMID and PREGID
ReMAPP.id <- distinct(ReMAPP.id, SITE, MOMID, PREGID)

data.nutr <- read_dta("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_NUTR.dta")
data.nutr1 <- read_dta("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_NUTR_RFA.dta")
data.risk <- read_dta("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_RISKS.dta")
data.anemia <- read_dta("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_ANEMIA.dta")
data.demographic<-read.csv("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_DEMOGRAPHIC.csv")
data.infection<-read.csv("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_INFECTION.csv")
data.supplement<-read.csv("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_SUPPLEMENT.csv")
data.rbc<-read.csv("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_RBC.csv")
data.fuel<-read.csv("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_FUEL.csv")

data.nutr<-left_join(ReMAPP.id, data.nutr, by = c("SITE", "MOMID","PREGID"))
data.nutr1<-left_join(ReMAPP.id, data.nutr1, by = c("SITE", "MOMID","PREGID"))
data.risk <-left_join(ReMAPP.id, data.risk, by = c("SITE", "MOMID", "PREGID"))
data.anemia<-left_join(ReMAPP.id, data.anemia, by = c("SITE", "MOMID","PREGID"))
data.demographic<-left_join(ReMAPP.id, data.demographic, by = c("SITE", "MOMID", "PREGID"))
data.infection<-left_join(ReMAPP.id, data.infection, by = c("SITE", "MOMID", "PREGID"))
data.infection <- distinct(data.infection, SITE, MOMID,.keep_all = TRUE)

data.aim3 <- data.anemia %>%
  left_join(data.nutr, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(data.risk, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(data.demographic, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(data.infection, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(data.supplement, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(data.rbc, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(data.fuel, by = c("SITE", "MOMID", "PREGID")) %>%
  mutate(
    ##Nutrition
    FERRITIN70_ANC20=case_when(
      FERRITIN70_ANC20==2 ~ 1,
      FERRITIN70_ANC20==1 ~ 0,
      TRUE ~ NA_real_),
    FERRITIN15_ANC20=case_when(
      FERRITIN15_ANC20==2 ~ 1,
      FERRITIN15_ANC20==1 ~ 0,
      TRUE ~ NA_real_),
    STFR_ANC20=case_when(
      STFR_ANC20==3 ~ "High",
      STFR_ANC20==2 ~ "Normal",
      STFR_ANC20==1 ~ "Low",
      TRUE ~ NA_character_),
    RBP4_ANC20=case_when(
      RBP4_ANC20==4 ~ "No_deficiency",
      RBP4_ANC20==3 ~ "Mild_deficiency",
      RBP4_ANC20==2 ~ "Moderate_deficiency",
      RBP4_ANC20==1 ~ "Severe_deficiency",
      TRUE ~ NA_character_),
    VITB12_COB_ANC20=case_when(
      VITB12_COB_ANC20==3 ~ "Deficient",
      VITB12_COB_ANC20==2 ~ "Insufficient",
      VITB12_COB_ANC20==1 ~ "Sufficient",
      TRUE ~ NA_character_),
    VITB12_HOL_ANC20=case_when(
      VITB12_HOL_ANC20==2 ~ 0,
      VITB12_HOL_ANC20==1 ~ 1,
      TRUE ~ NA_real_),
    FOL_SERUM_ANC20=case_when(
      FOL_SERUM_ANC20==4 ~ "Elevated",
      FOL_SERUM_ANC20==3 ~ "Normal",
      FOL_SERUM_ANC20==2 ~ "Possibly_deficient",
      FOL_SERUM_ANC20==1 ~ "Deficient",
      TRUE ~ NA_character_),
    HIGH_TG_44_ANC20=case_when(
      HIGH_TG_44_ANC20==2 ~ 1,
      HIGH_TG_44_ANC20==1 ~ 0,
      TRUE ~ NA_real_),
    ###Inflammation
    CRP_ANC20=case_when(
      CRP_ANC20==2 ~ 1,
      CRP_ANC20==1 ~ 0,
      TRUE ~ NA_real_),
    AGP_ANC20=case_when(
      AGP_ANC20==2 ~ 1,
      AGP_ANC20==1 ~ 0,
      TRUE ~ NA_real_),
    ##RBC
    MCV_ANC20=case_when(
      MCV_ANC20==3 ~ "Macrocytic",
      MCV_ANC20==2 ~ "Normal",
      MCV_ANC20==1 ~ "Microcytic",
      TRUE ~ NA_character_),
    RBC_SICKLE=case_when(
      M08_RBC_THALA_1==1|M08_RBC_THALA_3==1|M08_RBC_THALA_5==1|M08_RBC_THALA_15==1 ~ "Disease" ,
      M08_RBC_THALA_16==1 ~ "Trait" ,
      (M08_RBC_THALA_1==0|is.na(M08_RBC_THALA_1)) & (M08_RBC_THALA_3==0|is.na(M08_RBC_THALA_3)) & (M08_RBC_THALA_5==0|is.na(M08_RBC_THALA_5)) &(M08_RBC_THALA_15==0|is.na(M08_RBC_THALA_15)) ~ "Normal"
    ),
    RBC_Thalassemia=case_when(
      M08_RBC_THALA_6==1|M08_RBC_THALA_7==1|M08_RBC_THALA_8==1|M08_RBC_THALA_11==1|M08_RBC_THALA_12==1|M08_RBC_THALA_13==1|M08_RBC_THALA_14==1 ~ "Disease" ,
      (M08_RBC_THALA_6==0|is.na(M08_RBC_THALA_6)) & (M08_RBC_THALA_7==0|is.na(M08_RBC_THALA_7)) & (M08_RBC_THALA_8==0|is.na(M08_RBC_THALA_8)) & (M08_RBC_THALA_11==0|is.na(M08_RBC_THALA_11))& (M08_RBC_THALA_12==0|is.na(M08_RBC_THALA_12)) & (M08_RBC_THALA_13==0|is.na(M08_RBC_THALA_13)) & (M08_RBC_THALA_14==0|is.na(M08_RBC_THALA_14)) ~ "Normal" ,
    ),
    ##Supplement
    IRON_Supplement=case_when(
      M04_IRON_IV_CMOCCUR==1|M04_IRON_ORAL_CMOCCUR==1 ~ 1,
      (M04_IRON_IV_CMOCCUR==0|is.na(M04_IRON_IV_CMOCCUR))&(M04_IRON_ORAL_CMOCCUR==0|is.na(M04_IRON_ORAL_CMOCCUR)) ~ 0,
      TRUE ~ NA_real_,
    ),
    ##Demographic
    hh_smoke=case_when(
      hh_smoke==1 ~ 1,
      hh_smoke==0 ~ 0,
      TRUE ~ NA_real_),
    SCHOOL_MORE10=case_when(
      SCHOOL_MORE10==1 ~ 1,
      SCHOOL_MORE10==0 ~ 0,
      TRUE ~ NA_real_),
    water_improved = case_when(
      water_improved == 1 ~ 1,
      water_improved == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    toilet_improved = case_when(
      toilet_improved == 1 ~ 1,
      toilet_improved == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WEALTH_QUINT_5 = case_when(
      WEALTH_QUINT_5 == 1 ~ 1,
      WEALTH_QUINT_5 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WEALTH_QUINT_4 = case_when(
      WEALTH_QUINT_4 == 1 ~ 1,
      WEALTH_QUINT_4 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WEALTH_QUINT_3 = case_when(
      WEALTH_QUINT_3 == 1 ~ 1,
      WEALTH_QUINT_3 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WEALTH_QUINT_2 = case_when(
      WEALTH_QUINT_2 == 1 ~ 1,
      WEALTH_QUINT_2 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    WEALTH_QUINT_1 = case_when(
      WEALTH_QUINT_1 == 1 ~ 1,
      WEALTH_QUINT_1 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    ANEMIA_IND1 = case_when(
      ANEMIA_T1 %in% c(1, 2, 3) | ANEMIA_T2 %in% c(1, 2, 3) | ANEMIA_T3 %in% c(1, 2, 3) ~ 1,
      (ANEMIA_T1 ==55 |is.na(ANEMIA_T1) | ANEMIA_T1 == 0) & (ANEMIA_T2 ==55 | is.na(ANEMIA_T2) | ANEMIA_T2 == 0) & (ANEMIA_T3 ==55 | is.na(ANEMIA_T3) | ANEMIA_T3 == 0) ~ 0,
      (ANEMIA_T1 ==55 |is.na(ANEMIA_T1)) & (ANEMIA_T2 ==55 |is.na(ANEMIA_T2)) & (ANEMIA_T3 ==55 |is.na(ANEMIA_T3)) ~ NA_real_,
    ),
    ANEMIA_IND2 = case_when(
      ANEMIA_T1 %in% c(2, 3) | ANEMIA_T2 %in% c(2, 3) | ANEMIA_T3 %in% c(2, 3) ~ 1,
      (ANEMIA_T1 ==55 |is.na(ANEMIA_T1) | ANEMIA_T1 %in% c(0, 1)) & (ANEMIA_T2 ==55 | is.na(ANEMIA_T2) | ANEMIA_T2 %in% c(0, 1)) & (ANEMIA_T3 ==55 | is.na(ANEMIA_T3) | ANEMIA_T3 %in% c(0, 1)) ~ 0,
      (ANEMIA_T1 ==55 |is.na(ANEMIA_T1)) & (ANEMIA_T2 ==55 |is.na(ANEMIA_T2)) & (ANEMIA_T3 ==55 |is.na(ANEMIA_T3)) ~ NA_real_,
    )
  )

data.aim3$RBC_SICKLE<-ifelse(is.na(data.aim3$M08_RBC_THALA_1) & is.na(data.aim3$M08_RBC_THALA_3) & is.na(data.aim3$M08_RBC_THALA_5) & is.na(data.aim3$M08_RBC_THALA_15),NA,data.aim3$RBC_SICKLE)
data.aim3$RBC_Thalassemia<-ifelse(is.na(data.aim3$M08_RBC_THALA_6) & is.na(data.aim3$M08_RBC_THALA_7) & is.na(data.aim3$M08_RBC_THALA_8) & is.na(data.aim3$M08_RBC_THALA_11) & is.na(data.aim3$M08_RBC_THALA_12) & is.na(data.aim3$M08_RBC_THALA_13) & is.na(data.aim3$M08_RBC_THALA_14),NA,data.aim3$RBC_Thalassemia)
data.aim3$IRON_Supplement<-ifelse(is.na(data.aim3$M04_IRON_IV_CMOCCUR) & is.na(data.aim3$M04_IRON_ORAL_CMOCCUR),NA,data.aim3$IRON_Supplement)
data.aim3$ANEMIA_IND1<-ifelse((data.aim3$ANEMIA_T1 ==55 |is.na(data.aim3$ANEMIA_T1)) & (data.aim3$ANEMIA_T2 ==55 |is.na(data.aim3$ANEMIA_T2)) & (data.aim3$ANEMIA_T3 ==55 |is.na(data.aim3$ANEMIA_T3)),NA,data.aim3$ANEMIA_IND1)
data.aim3$ANEMIA_IND2<-ifelse((data.aim3$ANEMIA_T1 ==55 |is.na(data.aim3$ANEMIA_T1)) & (data.aim3$ANEMIA_T2 ==55 |is.na(data.aim3$ANEMIA_T2)) & (data.aim3$ANEMIA_T3 ==55 |is.na(data.aim3$ANEMIA_T3)),NA,data.aim3$ANEMIA_IND2)  

data.aim3<-data.aim3 %>% dummy_cols(select_columns = "STFR_ANC20")%>%
  dummy_cols(select_columns = "RBP4_ANC20")%>%
  dummy_cols(select_columns = "VITB12_COB_ANC20")%>%
  dummy_cols(select_columns = "FOL_SERUM_ANC20")%>%
  dummy_cols(select_columns = "MCV_ANC20")%>%
  dummy_cols(select_columns = "RBC_SICKLE")%>%
  dummy_cols(select_columns = "RBC_Thalassemia")

long_data.aim3 <- data.aim3 %>%
  pivot_longer(
    cols = c("ANEMIA_T1","ANEMIA_T2","ANEMIA_T3"),
    names_to = "ANEMIA_TRIMESTER",
    values_to = "ANEMIA_STATUS"
  )%>%
  mutate(
    ##Outcome
    ANEMIA_OUTCOME1=case_when(
      ANEMIA_STATUS %in% c(1,2,3) ~ 1,
      ANEMIA_STATUS == 0 ~ 0,
      TRUE ~ NA_real_),
    ANEMIA_OUTCOME2=case_when(
      ANEMIA_STATUS %in% c(2,3) ~ 1,
      ANEMIA_STATUS %in% c(0,1) ~ 0,
      TRUE ~ NA_real_))
