#****************************************************************************
#*Aim 1: Healthy cohort criterias
#*Author: Xiaoyan Hu
#****************************************************************************
rm(list = ls())

library(tidyverse)
library(lubridate)
library(naniar)

UploadDate = "2024-05-03"

#****************************************************************************
#0. Prepare original data
#****************************************************************************
#Load data with basic ID information
load(paste0("D:/Users/xyh/Documents/github/0-data-base/", UploadDate, "/df_mat.rda"))

#load mnh00 and keep necessary variables
mnh00 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh00_merged.csv")) %>% 
  select(SITE, SCRNID, M00_KNOWN_DOBYN_SCORRES, 
         M00_BRTHDAT, M00_ESTIMATED_AGE, M00_SCHOOL_YRS_SCORRES, M00_SCHOOL_SCORRES)

#load mnh03 and keep necessary variables
mnh03 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh03_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, M03_MARITAL_SCORRES,
         M03_SMOKE_OECOCCUR, M03_CHEW_BNUT_OECOCCUR, M03_CHEW_OECOCCUR, M03_DRINK_OECOCCUR)

#load mnh04 and keep necessary variables
mnh04 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh04_merged.csv")) %>% 
  filter(M04_TYPE_VISIT == 1) %>% 
  select(SITE, MOMID, PREGID, M04_PRETERM_RPORRES, M04_PH_PREV_RPORRES, M04_PH_PREVN_RPORRES, M04_PH_LIVE_RPORRES, 
         M04_MISCARRIAGE_RPORRES, M04_MISCARRIAGE_CT_RPORRES, M04_PH_OTH_RPORRES,M04_STILLBIRTH_RPORRES,
         M04_LOWBIRTHWT_RPORRES, M04_MALARIA_EVER_MHOCCUR, 
         M04_CANCER_EVER_MHOCCUR, M04_KIDNEY_EVER_MHOCCUR, M04_CARDIAC_EVER_MHOCCUR,
         M04_HIV_MHOCCUR, M04_HIV_EVER_MHOCCUR, M04_UNPL_CESARIAN_PROCCUR, M04_PREECLAMPSIA_RPORRES,
         M04_GEST_DIAB_RPORRES, M04_PREMATURE_RUPTURE_RPORRES,
         M04_MACROSOMIA_RPORRES, M04_OLIGOHYDRAMNIOS_RPORRES,
         M04_APH_RPORRES, M04_PPH_RPORRES)

#load mnh05 and keep necessary variables
mnh05 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh05_merged.csv")) %>% 
  filter(M05_TYPE_VISIT == 1) %>% 
  select(SITE, MOMID, PREGID, M05_ANT_PEDAT, M05_WEIGHT_PERES, M05_HEIGHT_PERES, M05_MUAC_PERES)

#load mnh06 and keep necessary variables
mnh06 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh06_merged.csv")) %>% 
  filter(M06_TYPE_VISIT == 1) %>% 
  select(SITE, MOMID, PREGID, M06_SINGLETON_PERES, 
         M06_BP_SYS_VSORRES_1, M06_BP_SYS_VSORRES_2, M06_BP_SYS_VSORRES_3,
         M06_BP_DIA_VSORRES_1, M06_BP_DIA_VSORRES_2, M06_BP_DIA_VSORRES_3,
         M06_MALARIA_POC_LBORRES, M06_MALARIA_POC_LBPERF, 
         M06_HBV_POC_LBORRES, M06_HBV_POC_LBPERF, M06_HCV_POC_LBORRES, M06_HCV_POC_LBPERF,
         M06_HIV_POC_LBORRES, M06_HIV_POC_LBPERF,
         num_range("M06_HB_POC_LBORRES_",1:12))

#load mnh08 and keep necessary variables
mnh08 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh08_merged.csv")) %>% 
  filter(M08_TYPE_VISIT == 1) %>% 
  select(SITE, MOMID, PREGID, 
         M08_MN_LBPERF_8, M08_FERRITIN_LBORRES, 
         M08_RBC_LBPERF_2, M08_RBC_THALA_LBORRES, M08_RBC_LBPERF_3, M08_RBC_GLUC6_LBORRES,
         M08_MN_LBPERF_12, M08_CRP_LBORRES, M08_MN_LBPERF_13, M08_AGP_LBORRES)

#load mnh09 and keep necessary variables
mnh09 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh09_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, num_range("M09_INFANTID_INF",1:4),
         num_range("M09_INFANTID_INF",1:4),
         num_range("M09_BIRTH_DSTERM_INF",1:4), 
         num_range("M09_DELIV_DSSTDAT_INF",1:4))

#merge maternal data after ReMAPP luanches
df_maternal <- df_mat %>%
  mutate(REMAPP_LAUNCH = ifelse((SITE == "Ghana" & M02_SCRN_OBSSTDAT >= "2022-12-28") |
                                  (SITE == "Kenya" & M02_SCRN_OBSSTDAT >= "2023-04-14") |
                                  (SITE == "Zambia" & M02_SCRN_OBSSTDAT >= "2022-12-15") |
                                  (SITE == "Pakistan" & M02_SCRN_OBSSTDAT >= "2022-09-22") |
                                  (SITE == "India-CMC" & M02_SCRN_OBSSTDAT >= "2023-06-20") |
                                  (SITE == "India-SAS" & M02_SCRN_OBSSTDAT >= "2023-08-15"), 1, 0)) %>% 
  filter(REMAPP_LAUNCH == 1) %>%
  left_join(mnh00, by = c("SITE", "SCRNID")) %>% 
  left_join(mnh03, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh04, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh05, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh06, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh08, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh09, by = c("SITE", "MOMID", "PREGID"))

#save data
save(df_maternal, file = "derived_data/df_maternal.rda")

#****************************************************************************
#1. define criteria
#****************************************************************************
#derive criteria
df_criteria <- df_maternal %>%
  mutate(
    # A. age at enrollment
    # Aged 18 to 34 years
    AGE_ENROLL = ifelse(M00_KNOWN_DOBYN_SCORRES == 1 &  M00_BRTHDAT != "1907-07-07", 
                        as.numeric(ymd(M02_SCRN_OBSSTDAT) - ymd(M00_BRTHDAT))/365,
                        ifelse(M00_KNOWN_DOBYN_SCORRES == 0 & M00_ESTIMATED_AGE != -7, M00_ESTIMATED_AGE, 99)),
    CRIT_AGE = ifelse((AGE_ENROLL > 0 & AGE_ENROLL < 18) | AGE_ENROLL > 34, 0,
                      ifelse(AGE_ENROLL >= 18 & AGE_ENROLL <= 34, 1, 55)
    ),
    
    # B. GA at enrollment
    # gestational age at enrollment - Gestational age <14 weeks 
    BASELINE_GA_WKS = floor(BOE_GA_DAYS_ENROLL/7),
    CRIT_GA = ifelse(BASELINE_GA_WKS > 0 & BASELINE_GA_WKS < 14, 1,
                     ifelse(BASELINE_GA_WKS >= 14 & BASELINE_GA_WKS <=26, 0,
                            ifelse(BASELINE_GA_WKS == -7 | is.na(BASELINE_GA_WKS), NA, 77))),
    
    # C. Pre-pregnancy or early pregnancy body mass index (BMI) of >18.5 and <30 kg/m2 AND mid-upper arm circumference (MUAC) > 23cm [45]
    # BMI
    BMI = M05_WEIGHT_PERES / M05_HEIGHT_PERES / M05_HEIGHT_PERES * 10000,
    
    TEMP_BMI = ifelse(BMI <= 18.5 | BMI >= 30, 0, 
                      ifelse(BMI > 18.5 & BMI < 30, 1, 55)
    ),
    # MUAC mid-upper arm circumference - MUAC
    TEMP_MUAC = ifelse(M05_MUAC_PERES <= 23, 0, 
                       ifelse(M05_MUAC_PERES > 23, 1, 55)
    ),
    CRIT_BMI_MUAC = case_when(
      TEMP_BMI == 1 & TEMP_MUAC == 1 ~ 1, 
      TEMP_BMI == 0 | TEMP_MUAC == 0 ~ 0, 
      TRUE ~ 55
    ),
    # D. Height ≥150 cm
    CRIT_HEIGHT = ifelse(M05_HEIGHT_PERES < 150, 0,
                         ifelse(M05_HEIGHT_PERES >= 150, 1, 55)
    ),
    # E. Singleton pregnancy
    CRIT_SINGLEPREG = ifelse(M06_SINGLETON_PERES == 0, 0,
                             ifelse(M06_SINGLETON_PERES == 1, 1, 55)
    ),
    # F. no iron deficiency (not iron deficient: serum ferritin > 15 mcg/L) data unit is ??g/dL couble check before use
    #convert unit from ug/dL to mcg/L
    FERRITIN_LBORRES = case_when(
      SITE %in% c("Ghana", "Kenya") ~ 10*M08_FERRITIN_LBORRES, #need add Pakistan and India conversion rate when we have data
      SITE == "Zambia" ~ M08_FERRITIN_LBORRES 
    ),
    CRIT_IRON = ifelse(FERRITIN_LBORRES > 15, 1,
                       ifelse(FERRITIN_LBORRES >0 & FERRITIN_LBORRES <= 15, 0,
                              ifelse(FERRITIN_LBORRES == 0, 0, 55))
    ),
    # G. no subclinical inflammation (CRP???5 and/or AGP???1) ??? check unit (mg/L for CRP and g/L for AGP in dd) double check the calculation before use
    CRIT_INFLAM = case_when(
      M08_CRP_LBORRES > 0 & M08_CRP_LBORRES <= 5 & M08_AGP_LBORRES >0 & M08_AGP_LBORRES <= 1 ~ 1,
      M08_CRP_LBORRES > 5 | M08_AGP_LBORRES > 1 ~ 0,
      M08_MN_LBPERF_12 == 0 | M08_MN_LBPERF_13 == 0 ~ 55,
      TRUE ~ 55
    )
  ) %>% 
  rowwise() %>% 
  mutate(
    # H.a. blood pressure
    M06_BP_SYS_1 = mean(c(M06_BP_SYS_VSORRES_1, M06_BP_SYS_VSORRES_2, M06_BP_SYS_VSORRES_3), na.rm = TRUE),
    M06_BP_DIA_1 = mean(c(M06_BP_DIA_VSORRES_1, M06_BP_DIA_VSORRES_2, M06_BP_DIA_VSORRES_3), na.rm = TRUE),
    
    CRIT_BP = ifelse(M06_BP_SYS_1 > 0 & M06_BP_SYS_1 < 140 & M06_BP_DIA_1 > 0 & M06_BP_DIA_1 < 90, 1,
                     ifelse(M06_BP_SYS_1 >= 140 | M06_BP_DIA_1 >= 90, 0, 55)
    )) %>% 
  ungroup() %>% 
  mutate(
    # H.b. no previous low birth weight delivery
    CRIT_LBW = ifelse(M04_LOWBIRTHWT_RPORRES == 1, 0,
                      ifelse(M04_PH_PREV_RPORRES == 0 | M04_LOWBIRTHWT_RPORRES == 0, 1,
                             ifelse(M04_LOWBIRTHWT_RPORRES == 99, 0, 55))
    ),
    # H.c. No previous reported stillbirth
    CRIT_STILLBIRTH = ifelse(M04_STILLBIRTH_RPORRES == 1, 0, #stillbirth,
                             ifelse(M04_PH_PREV_RPORRES == 0 | 
                                      M04_PH_OTH_RPORRES == 0 | #no fetal loss
                                      M04_STILLBIRTH_RPORRES == 0, 1,
                                    ifelse(M04_STILLBIRTH_RPORRES == 99, 0, 55))  
    ),
    # H.d. No previous reported unplanned cesarean delivery
    CRIT_UNPL_CESARIAN = case_when(
      M04_UNPL_CESARIAN_PROCCUR == 1 ~ 0, 
      M04_PH_PREV_RPORRES == 0 | M04_UNPL_CESARIAN_PROCCUR == 0 ~ 1,
      M04_UNPL_CESARIAN_PROCCUR == 99 ~ 0,
      TRUE ~ 55 
    ),
    # I. No hemoglobinopathies: SS, SC, SE, EE, CC, SD-Punjab, Sβthal, Eβthal, Cβthal, CD-Punjab, ED-Punjab, D-D-Punjab, 
    # D-Punjabβthal, Thalassemia major, Thalassemia intermedia, glucose-6-phosphate dehydrogenase deficiency, or Alpha thalassemia
    CRIT_HEMOGLOBINOPATHIES = ifelse(M08_RBC_THALA_LBORRES == 0 & M08_RBC_GLUC6_LBORRES == 0, 1,
                                     ifelse(M08_RBC_THALA_LBORRES == 1 | M08_RBC_GLUC6_LBORRES == 1, 0,
                                            ifelse(M08_RBC_LBPERF_2 == 0 | M08_RBC_LBPERF_3 == 0, 55, 55))
    ),
    #J. No reported cigarette smoking, tobacco chewing, or betel nut use during pregnancy
    CRIT_SMOKE = case_when(
      SITE == "Zambia" & (M03_SMOKE_OECOCCUR == 1 | M03_CHEW_OECOCCUR == 1) ~ 0,
      SITE == "Zambia" & (M03_SMOKE_OECOCCUR == 0 & M03_CHEW_OECOCCUR == 0) ~ 1,
      M03_SMOKE_OECOCCUR == 1 | M03_CHEW_BNUT_OECOCCUR == 1 | M03_CHEW_OECOCCUR == 1 ~ 0,
      M03_SMOKE_OECOCCUR == 0 & M03_CHEW_BNUT_OECOCCUR == 0 & M03_CHEW_OECOCCUR == 0 ~ 1,
      TRUE ~ 55
    ),
    #K. No reported alcohol consumption during pregnancy
    CRIT_DRINK = ifelse(SITE == "Pakistan", 666,
                        ifelse(M03_DRINK_OECOCCUR == 1, 0,
                               ifelse(M03_DRINK_OECOCCUR == 0, 1,
                                      ifelse(M03_DRINK_OECOCCUR == 66, 0,
                                             ifelse(M03_DRINK_OECOCCUR == 77, 0, 55)))) #temporary code for Kenya, check for other country
    ), 
    #L. No known history or current chronic disease including cancer, kidney disease, and cardiac conditions
    CRIT_CHRONIC = ifelse(M04_CANCER_EVER_MHOCCUR == 1 | M04_KIDNEY_EVER_MHOCCUR == 1 | 
                            M04_CARDIAC_EVER_MHOCCUR == 1, 0,
                          ifelse(M04_CANCER_EVER_MHOCCUR == 0 & M04_KIDNEY_EVER_MHOCCUR == 0 & 
                                   M04_CARDIAC_EVER_MHOCCUR == 0, 1,
                                 ifelse(M04_CANCER_EVER_MHOCCUR == 99 | M04_KIDNEY_EVER_MHOCCUR == 99 | 
                                          M04_CARDIAC_EVER_MHOCCUR == 99, 0, 55))
    ),
    #M. No known history or current HIV
    # if "Record HIV results" = positive, then CRIT_HIV=0 (ineligible) [M06_HIV_POC_LBORRES]
    CRIT_HIV = ifelse(M06_HIV_POC_LBORRES == 1, 0, 
   # if "Record HIV results" = negative, then CRIT_HIV=1 (eligible) [M06_HIV_POC_LBORRES]
    ifelse(M06_HIV_POC_LBORRES == 0, 1,
   # if "Record HIV results" = 55, then CRIT_HIV=55 (pending) [M06_HIV_POC_LBORRES]
    ifelse(M06_HIV_POC_LBORRES == 55, 55,  
    # if "Have you ever been diagnosed with HIV?" = yes OR 
    # if "Have you had any of the following issues since becoming pregnant with the current pregnancy, HIV" = yes, then CRIT_HIV=0 (ineligible)
    ifelse(M04_HIV_EVER_MHOCCUR == 1 |  M04_HIV_MHOCCUR == 1, 0,
    # if "Have you ever been diagnosed with HIV?" = no AND 
    # if "Have you had any of the following issues since becoming pregnant with the current pregnancy, HIV" = no, then CRIT_HIV=1 (eligible)
    ifelse(M04_HIV_EVER_MHOCCUR == 0 & M04_HIV_MHOCCUR == 0, 1,
    # if "Have you ever been diagnosed with HIV?" = 55, then CRIT_HIV=55 (pending) OR
    # if "Have you had any of the following issues since becoming pregnant with the current pregnancy, HIV" = 55, then CRIT_HIV=55 (pending)
    ifelse(M04_HIV_EVER_MHOCCUR == 55 | M04_HIV_MHOCCUR == 55, 55,  
    # if "Have you ever been diagnosed with HIV?" = 77/99 AND 
    # if "Have you had any of the following issues since becoming pregnant with the current pregnancy, HIV" = 77/99, then CRIT_HIV=55 (pending)
    ifelse(M04_HIV_EVER_MHOCCUR %in% c(0,99) & M04_HIV_MHOCCUR %in% c(0,99), 0, 55)))))) 
    ),
    #N. No current malaria infection (per rapid diagnostic test)
    CRIT_MALARIA = case_when(
      M06_MALARIA_POC_LBORRES == 1 ~ 0,
      M06_MALARIA_POC_LBORRES == 0 ~ 1,
      M06_MALARIA_POC_LBPERF == 0 ~ 0,
      TRUE ~ 55
    ),
    #O. No current Hepatitis B virus infection (per rapid diagnostic test)
    CRIT_HEPATITISB = ifelse(M06_HBV_POC_LBORRES == 1, 0,
                             ifelse(M06_HBV_POC_LBORRES == 0, 1,
                                    ifelse(M06_HBV_POC_LBPERF == 0, 55, 55))
    ),
    #P. No current Hepatitis C virus infection (per rapid diagnostic test)
    CRIT_HEPATITISC = ifelse(M06_HCV_POC_LBORRES == 1, 0,
                             ifelse(M06_HCV_POC_LBORRES == 0, 1,
                                    ifelse(M06_HCV_POC_LBPERF == 0, 55, 55)))
  ) 
#After enrollment, participants will be excluded from the final analysis if any of the following occur: 
#Multiple pregnancies not identified at recruitment
#Severe conditions not evident at recruitment including cancer, HIV, TB, or Malaria
#Severe pregnancy-related conditions requiring hospital admission including eclampsia or severe pre-eclampsia
save(df_criteria, file = "derived_data/df_criteria.rda")

#**************************************************************************************
#*2. check eligibility and save df_healthy.rda
#**************************************************************************************
#code 666 for any not applicable by site
healthyOutcome <- df_criteria %>% 
  rowwise() %>%
  mutate(HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE)) %>% 
  mutate(
    HEALTHY_ELIGIBLE = case_when(
      if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
      if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
      HEALTHY_CHECK < 19 ~ 3), #19 criterias
    #!!!!!! temp code for healthy_eligible due to small eligible sample
    HEALTHY_ELIGIBLE = case_when(
      CRIT_AGE == 1 &
        # CRIT_GA == 1 &
        CRIT_BMI_MUAC == 1 &
        CRIT_HEIGHT == 1 &
        CRIT_SINGLEPREG == 1 &
        # CRIT_IRON == 1 
        # CRIT_INFLAM == 1 
        CRIT_BP == 1 &
        CRIT_LBW == 1 &
        CRIT_STILLBIRTH == 1 &
        CRIT_UNPL_CESARIAN == 1 &
        # CRIT_HEMOGLOBINOPATHIES == 1 & 
        CRIT_SMOKE == 1 & 
        CRIT_DRINK %in% c(1,666) &
        CRIT_CHRONIC == 1 &
        CRIT_HIV == 1 &
        CRIT_MALARIA == 1 & 
        CRIT_HEPATITISB == 1 & 
        CRIT_HEPATITISC == 1
      ~ 1, 
      TRUE ~ 0
    ) ) %>%
  ungroup() 

df_healthy <- healthyOutcome %>% 
  filter(HEALTHY_ELIGIBLE == 1)
table(df_healthy$SITE)

#save data
save(healthyOutcome, file = "derived_data/healthyOutcome.rda")
save(df_healthy, file = "derived_data/df_healthy.rda")
#*****************************************************************************
#*3. HB data --> long and wide
#*****************************************************************************
#*prepare data
prep_hb1 <- df_healthy %>%
  dplyr:: select("MOMID", "PREGID", "SITE",
                 EST_CONCEP_DATE,
                 M03_SMOKE_OECOCCUR) %>%
  #replace 7s and 5s with NA
  replace_with_na_all(condition = ~.< 0) %>%
  replace_with_na_all(condition = ~. %in% c("1907-07-07", "1907-05-05")) %>%
  distinct()

#mnh08_hb
mnh08_hb <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh08_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, 
         M08_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES) %>% 
  #remove duplicates if there is any
  group_by(SITE, MOMID, PREGID, M08_TYPE_VISIT) %>% 
  mutate(n = n()) %>% 
  filter(n == 1, M08_TYPE_VISIT %in% c(1:12)) %>% 
  replace_with_na_all(condition = ~.< 0) %>%
  replace_with_na_all(condition = ~. %in% c("1907-07-07", "1907-05-05")) 


#long hb data
df_hb_long1_all <- mnh08_hb %>% 
  left_join(prep_hb1, by = c("SITE", "MOMID", "PREGID")) %>% 
  #derive variables
  mutate(
    ga_wks = case_when(
      M08_TYPE_VISIT >= 6 ~ NA_real_,
      M08_TYPE_VISIT < 6 ~ as.numeric(ymd(M08_LBSTDAT) - EST_CONCEP_DATE)/7
    ),
    trimester = case_when(
      M08_TYPE_VISIT >= 6 ~ NA_real_,
      ga_wks > 0 & ga_wks < 14 ~ 1,
      ga_wks >= 14 & ga_wks < 27 ~ 2,
      ga_wks >= 27 & ga_wks <= 40 ~ 3,
      TRUE ~ NA_real_
    ), 
    hb_alti = case_when(
      SITE == "Kenya" ~ M08_CBC_HB_LBORRES - 0.2,
      SITE == "Zambia" ~ M08_CBC_HB_LBORRES - 0.5,
      !is.na(M08_CBC_HB_LBORRES) ~ M08_CBC_HB_LBORRES,
      TRUE ~ NA_real_
    ), 
    #adjust for both smoke and altitude
    hb = case_when(
      M03_SMOKE_OECOCCUR == 1 ~ hb_alti - 0.3,
      M03_SMOKE_OECOCCUR == 0 ~ hb_alti,
      TRUE ~ NA_real_
    )) %>% 
  #Remove rows that have gestational weeks less than 10 or greater than 50 (42-week pregnancy + 8-week postpartum).
  #Remove rows that have hemoglobin less than 5 or greater than 18.
  #Remove rows that have missing values on gestational weeks.
  filter(!is.na(EST_CONCEP_DATE) & !is.na(M08_LBSTDAT) & !is.na(M08_CBC_HB_LBORRES)) 

#wide hb data for all healthy cohort with ga_wks info and hb values
df_hb_wide1_all <- df_hb_long1_all %>%
  pivot_wider(
    names_from = M08_TYPE_VISIT,
    values_from = -c("MOMID","PREGID","SITE", EST_CONCEP_DATE, M03_SMOKE_OECOCCUR)
  ) %>%
  rowwise() %>%
  ungroup()

 #remove outliers
df_hb_long1 <- df_hb_long1_all %>% 
  filter(ga_wks >= 6) %>%
  filter(hb >= 5 & hb <= 18) 

#wide hb data 
df_hb_wide1 <- df_hb_long1 %>%
  pivot_wider(
    names_from = M08_TYPE_VISIT,
    values_from = -c("MOMID","PREGID","SITE", EST_CONCEP_DATE, M03_SMOKE_OECOCCUR)
  ) %>%
  rowwise() %>%
  ungroup()

#save data
save(df_hb_long1_all, file = "derived_data/df_hb_long1_all.rda")
save(df_hb_wide1_all, file = "derived_data/df_hb_wide1_all.rda")
save(df_hb_long1, file = "derived_data/df_hb_long1.rda")
save(df_hb_wide1, file = "derived_data/df_hb_wide1.rda")

#*****************************************************************************
#*4. sensitivity analysis data
#*****************************************************************************
#*prepare data
prep_sensitive_long <- df_hb_long1 %>% 
  left_join(df_healthy %>% 
              dplyr::select(MOMID, PREGID, SITE,
                            num_range("M09_INFANTID_INF",1:4), 
                            num_range("M09_DELIV_DSSTDAT_INF",1:4), 
                            num_range("M09_BIRTH_DSTERM_INF",1:4)), 
            by = c("MOMID", "PREGID", "SITE")) %>% 
  replace_with_na_all(condition = ~.== "1907-07-07") %>%
  replace_with_na_all(condition = ~.x %in% c("n/a")) %>%
  replace_with_na_at(.vars = c("M09_BIRTH_DSTERM_INF1", "M09_BIRTH_DSTERM_INF2",
                               "M09_BIRTH_DSTERM_INF3", "M09_BIRTH_DSTERM_INF4"),
                     condition = ~.x == 77) %>%
  distinct() 

#long data
df_sensitive_long <- prep_sensitive_long %>% 
  mutate(DOB_mom = 
           pmin(M09_DELIV_DSSTDAT_INF1, M09_DELIV_DSSTDAT_INF2, 
                M09_DELIV_DSSTDAT_INF3, M09_DELIV_DSSTDAT_INF4, na.rm = TRUE)) %>%
  mutate(GESTAGEBIRTH_BOE_DAYS = as.numeric(ymd(DOB_mom) - EST_CONCEP_DATE), 
         GESTAGEBIRTH_BOE = floor(GESTAGEBIRTH_BOE_DAYS/7)) %>%
  mutate(PRETERMBIRTH_LT37 = ifelse(GESTAGEBIRTH_BOE >= 20 & GESTAGEBIRTH_BOE < 37, 1, 0)) 

#save data
save(df_sensitive_long, file = "derived_data/df_sensitive_long.rda")

#*****************************************************************************
#*5. data for eligiblity and missingness
#*****************************************************************************
df_eli <- df_criteria %>%
  dplyr::select(MOMID, PREGID, SITE, starts_with("CRIT_")) %>%
  mutate(
    C1 = ifelse(!is.na(CRIT_AGE),CRIT_AGE,55), 
    C2 = ifelse(!is.na(CRIT_GA),CRIT_GA,55), 
    C3 = ifelse(!is.na(CRIT_BMI_MUAC),CRIT_BMI_MUAC,55),
    C4 = ifelse(!is.na(CRIT_HEIGHT),CRIT_HEIGHT,55), 
    C5 = ifelse(!is.na(CRIT_SINGLEPREG),CRIT_SINGLEPREG,55), 
    C6 = ifelse(!is.na(CRIT_LBW),CRIT_LBW,55),
    C7 = ifelse(!is.na(CRIT_STILLBIRTH),CRIT_STILLBIRTH,55), 
    C8 = ifelse(!is.na(CRIT_UNPL_CESARIAN),CRIT_UNPL_CESARIAN,55), 
    C9 = ifelse(!is.na(CRIT_SMOKE),CRIT_SMOKE,55), 
    C10 = case_when(
      CRIT_DRINK == 1|CRIT_DRINK == 0 ~ CRIT_DRINK,
      CRIT_DRINK == 666 ~ 1,
      TRUE~55), 
    C11 = ifelse(!is.na(CRIT_CHRONIC),CRIT_CHRONIC,55),
    C12 = ifelse(!is.na(CRIT_HIV),CRIT_HIV,55), 
    C13 = ifelse(!is.na(CRIT_BP),CRIT_BP,55), 
    C14 = ifelse(!is.na(CRIT_MALARIA),CRIT_MALARIA,55), 
    C15 = ifelse(!is.na(CRIT_HEPATITISB),CRIT_HEPATITISB,55),
    C16 = ifelse(!is.na(CRIT_HEPATITISC),CRIT_HEPATITISC,55),
    C17 = ifelse(!is.na(CRIT_IRON),CRIT_IRON,55),
    C18 = ifelse(!is.na(CRIT_INFLAM),CRIT_INFLAM,55), 
    C19 = ifelse(!is.na(CRIT_HEMOGLOBINOPATHIES),CRIT_HEMOGLOBINOPATHIES,55)
  ) %>% 
  mutate(across(matches("C\\d+"),
                function(x) 
                  factor(x, 
                         levels = c(1,0,55),
                         labels = c("Eligible", "Ineligible", "Pending")
                  )))


#make a long dataset
df_eli_long <- df_eli %>% 
  dplyr::select(-starts_with("CRIT_")) %>% 
  pivot_longer(cols = -c(MOMID, PREGID, SITE), 
               names_to = "Variable", 
               values_to = "Value")

#order variable value
df_eli_long$Variable <- factor(
  df_eli_long$Variable,
  levels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", 
             "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19"))

#save data
save(df_eli_long, file = "derived_data/df_eli_long.rda")


