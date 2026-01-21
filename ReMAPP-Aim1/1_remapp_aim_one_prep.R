#****************************************************************************
#*Aim 1: Healthy cohort criteria data preparation
#*Author: Xiaoyan Hu
#*Email: xyh@gwu.edu
#****************************************************************************
rm(list = ls())

library(tidyverse)
library(lubridate)
library(naniar)
library(BRINDA)
library(dplyr)
library(readxl)

UploadDate = "2025-10-31"

# Define base directory  
base_dir <- "D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim1"


# Create the full path for today's directory  
new_dir <- file.path(base_dir, UploadDate)

# Create the directory if it doesn't already exist  
if (!dir.exists(new_dir)) {
  dir.create(new_dir, recursive = TRUE)
}

# Set the new directory as the working directory  
setwd(new_dir)

# Print confirmation  
print(paste("Working directory set to:", getwd()))

derived_data_dir <- file.path(new_dir, "derived_data")

# Create Derived_Data directory if it doesn't exist  
if (!dir.exists(derived_data_dir)) {
  dir.create(derived_data_dir)
}

#****************************************************************************
#1. LOAD AND MERGE DATASET ----
#****************************************************************************
##load maternal enrolment ----
if (file.exists("derived_data/MAT_ENROLL.rda")) {
  load("derived_data/MAT_ENROLL.rda")   # loads MAT_ENROLL
} else {
  MAT_ENROLL <- read.csv(paste0("Z:/Outcome Data/", UploadDate, "/MAT_ENROLL.csv"))
  save(MAT_ENROLL, file = "derived_data/MAT_ENROLL.rda")
}

##load inf outcomes (composite, preterm37, lbw2500, sga10)----
if (file.exists("derived_data/INF_OUTCOMES.rda")) {
  load("derived_data/INF_OUTCOMES.rda")   # loads INF_OUTCOMES
} else {
  INF_OUTCOMES <- read.csv(paste0("Z:/Outcome Data/", UploadDate, "/INF_OUTCOMES.csv"))
  save(INF_OUTCOMES, file = "derived_data/INF_OUTCOMES.rda")
}

INF_OUTCOMES <- INF_OUTCOMES %>% 
  filter(!is.na(INFANTID)) %>% 
  select(SITE, MOMID, PREGID, INFANTID, PRETERMBIRTH_LT37) %>% 
  group_by(SITE, MOMID, PREGID) %>% 
  summarise(ifpreterm37mom = any(PRETERMBIRTH_LT37 == 1)) %>% 
  mutate(preterm37_mom = ifelse(ifpreterm37mom, 1, 0)) %>%
  ungroup()

##load gwg dataset----
if (file.exists("derived_data/MAT_GWG.rda")) {
  load("derived_data/MAT_GWG.rda")   
  
} else {
  MAT_GWG <- read.csv(paste0("Z:/Outcome Data/", UploadDate, "/GWG_OUTCOME_LONG_", UploadDate, ".csv"))
  save(MAT_GWG, file = "derived_data/MAT_GWG.rda")  
}

bmi_df <- MAT_GWG %>%
  select(
    SITE, MOMID, PREGID,
    WEIGHT_IMPUTED_FLAG,
    BMI_IMPUTTED = BMI,
    WEIGHT_ENROLL,
    HEIGHT_ENROLL
  ) %>%
  group_by(SITE, MOMID, PREGID, BMI_IMPUTTED) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    # Convert height to meters
    HEIGHT_M = HEIGHT_ENROLL / 100,
    
    # Manual BMI calculation
    BMI_MANUAL = WEIGHT_ENROLL / (HEIGHT_M ^ 2),
    
    # Rounded versions
    BMI_MANUAL_R = round(BMI_MANUAL, 1),
    BMI_IMPUTTED_R = round(BMI_IMPUTTED, 1),
    
    # Difference between rounded BMIs
    BMI_DIFF = BMI_IMPUTTED_R - BMI_MANUAL_R
  ) %>% 
  filter (!is.na(WEIGHT_IMPUTED_FLAG))


##load mat flowchart dataset----
if (file.exists("derived_data/MAT_FLOWCHART.rda")) {
  load("derived_data/MAT_FLOWCHART.rda")   # loads INF_OUTCOMES
} else {
  MAT_FLOWCHART <- read.csv(paste0("Z:/Outcome Data/", UploadDate, "/MAT_FLOWCHART.csv"))
  save(MAT_FLOWCHART, file = "derived_data/MAT_FLOWCHART.rda")
}

  
###save screening remapp data ----
screen_df_clean <- MAT_FLOWCHART %>%
  filter(REMAPP_PRESCRN == 1 | REMAPP_SCRN == 1)  %>%
  filter (!(REMAPP_PRESCRN == 1 & REMAPP_SCRN != 1 & ELIGIBLE == 1))  %>%
  mutate(REMAPP_PRESCRN = case_when(REMAPP_SCRN == 1 ~ 1, 
                                     TRUE ~ REMAPP_PRESCRN)) %>%
  distinct(SITE, SCRNID, MOMID, PREGID, .keep_all = TRUE)

save(screen_df_clean, file = "derived_data/screen_df_clean.rda")

##load mnh00 and keep necessary variables----
if (file.exists("derived_data/mnh00.rda")) {
  load("derived_data/mnh00.rda")   # loads mnh00
} else {
  mnh00 <- read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh00_merged.csv")) 
  save(mnh00, file = "derived_data/mnh00.rda")
}

mnh00_df <- mnh00 %>% 
  select(SITE, SCRNID,  M00_BRTHDAT, M00_ESTIMATED_AGE) %>%
  mutate(
    M00_BRTHDAT     = ymd(parse_date_time(
    M00_BRTHDAT,
      c("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d", "%d-%b-%y")
    )),
    M00_BRTHDAT = ymd(M00_BRTHDAT)
  ) %>%
  mutate(
    M00_ESTIMATED_AGE     = as.numeric(M00_ESTIMATED_AGE),
    M00_ESTIMATED_AGE= case_when(M00_ESTIMATED_AGE <= 0 ~ NA,
                                 TRUE ~ M00_ESTIMATED_AGE)
  )

##load mnh03 and keep necessary variables----
if (file.exists("derived_data/mnh03.rda")) {
  load("derived_data/mnh03.rda")   # loads mnh03
} else {
  mnh03 <- read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh03_merged.csv")) %>%
    dplyr::select(
      SITE, MOMID, PREGID, M03_MARITAL_SCORRES,
      M03_SMOKE_OECOCCUR, M03_CHEW_BNUT_OECOCCUR,
      M03_CHEW_OECOCCUR, M03_DRINK_OECOCCUR
    )
  save(mnh03, file = "derived_data/mnh03.rda")
}

## load mnh04 and keep necessary variables----
if (file.exists("derived_data/mnh04_raw.rda")) {
  load("derived_data/mnh04_raw.rda")   # loads mnh04_raw
} else {
  mnh04_raw <- read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh04_merged.csv"))
  save(mnh04_raw, file = "derived_data/mnh04_raw.rda")
}

mnh04<-  mnh04_raw %>% 
  filter(M04_TYPE_VISIT == 1) %>% 
  select(SITE, MOMID, PREGID, M04_PRETERM_RPORRES, M04_PH_PREV_RPORRES, M04_PH_PREVN_RPORRES, M04_PH_LIVE_RPORRES, 
         M04_MISCARRIAGE_RPORRES, M04_MISCARRIAGE_CT_RPORRES, M04_PH_OTH_RPORRES,M04_STILLBIRTH_RPORRES,
         M04_LOWBIRTHWT_RPORRES, M04_MALARIA_EVER_MHOCCUR, 
         M04_CANCER_EVER_MHOCCUR, M04_KIDNEY_EVER_MHOCCUR, M04_CARDIAC_EVER_MHOCCUR,
         M04_HIV_MHOCCUR, M04_HIV_EVER_MHOCCUR, M04_UNPL_CESARIAN_PROCCUR, M04_PREECLAMPSIA_RPORRES,
         M04_GEST_DIAB_RPORRES, M04_PREMATURE_RUPTURE_RPORRES,
         M04_MACROSOMIA_RPORRES, M04_OLIGOHYDRAMNIOS_RPORRES,
         M04_APH_RPORRES, M04_PPH_RPORRES)

#test for duplicates
test_dup_04 <- mnh04 %>%
  group_by(SITE, MOMID, PREGID) %>%
  filter(n() > 1) %>%
  ungroup()


##load mnh05 and keep necessary variables----
if (file.exists("derived_data/mnh05.rda")) {
  load("derived_data/mnh05.rda")   # loads mnh05
} else {
  mnh05 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh05_merged.csv")) %>% 
    filter(M05_TYPE_VISIT == 1) %>% 
    select(SITE, MOMID, PREGID, M05_ANT_PEDAT, M05_WEIGHT_PERES, M05_HEIGHT_PERES, M05_MUAC_PERES)
    save(mnh05, file = "derived_data/mnh05.rda")
}

test_dup_05 <- mnh05 %>%
  group_by(SITE, MOMID, PREGID) %>%
  filter(n() > 1) %>%
  ungroup()

##load mnh06 and keep necessary variables----
if (file.exists("derived_data/mnh06.rda")) {
  load("derived_data/mnh06.rda")   # loads mnh06
} else {
  mnh06 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh06_merged.csv")) %>% 
    filter(M06_TYPE_VISIT == 1) %>% 
    select(SITE, MOMID, PREGID, M06_SINGLETON_PERES, 
           M06_BP_SYS_VSORRES_1, M06_BP_SYS_VSORRES_2, M06_BP_SYS_VSORRES_3,
           M06_BP_DIA_VSORRES_1, M06_BP_DIA_VSORRES_2, M06_BP_DIA_VSORRES_3,
           M06_MALARIA_POC_LBORRES, M06_MALARIA_POC_LBPERF, 
           M06_HBV_POC_LBORRES, M06_HBV_POC_LBPERF, M06_HCV_POC_LBORRES, M06_HCV_POC_LBPERF,
           M06_HIV_POC_LBORRES, M06_HIV_POC_LBPERF,
           num_range("M06_HB_POC_LBORRES_",1:12)) #update to include all visits
  save(mnh06, file = "derived_data/mnh06.rda")
}

test_dup_06 <- mnh06 %>%
  group_by(SITE, MOMID, PREGID) %>%
  filter(n() > 1) %>%
  ungroup()

##load mnh08 and keep necessary variables----
if (file.exists("derived_data/mnh08_raw.rda")) {
  load("derived_data/mnh08_raw.rda")   # loads mnh08_raw
} else {
  mnh08_raw <- read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh08_merged.csv"))
  save(mnh08_raw, file = "derived_data/mnh08_raw.rda")
}


# For MNH08 quanys variables - phasing out the use of TYPE_VISIT, therefore, 
## we are taking the first non missing lab value measured between 0 and 139 days

mnh08 <- mnh08_raw %>%
  left_join(MAT_ENROLL %>% select(SITE, MOMID, PREGID, PREG_START_DATE),
            by = c("SITE", "MOMID", "PREGID")) %>%
  select(SITE, MOMID, PREGID, M08_LBSTDAT, M08_FERRITIN_LBORRES, 
          M08_CRP_LBORRES, M08_AGP_LBORRES, PREG_START_DATE) %>%
  mutate(M08_LBSTDAT = ymd(parse_date_time(M08_LBSTDAT, 
                                           c("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d", "%d-%b-%y"))),
         PREG_START_DATE = ymd(PREG_START_DATE))  %>%
  mutate(gest_age = as.numeric(M08_LBSTDAT - PREG_START_DATE)) %>%
  filter (gest_age >= 0)

# helper: first non-missing value in a vector
first_non_na <- function(x) {
  i <- which(!is.na(x))[1]
  if (is.na(i)) NA_real_ else x[i]
}

mnh08 <- mnh08_raw %>%
  left_join(
    MAT_ENROLL %>% select(SITE, MOMID, PREGID, PREG_START_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  select(
    SITE, MOMID, PREGID, M08_LBSTDAT,
    M08_FERRITIN_LBORRES, M08_CRP_LBORRES, M08_AGP_LBORRES,
    PREG_START_DATE
  ) %>%
  mutate(
    M08_LBSTDAT     = ymd(parse_date_time(
      M08_LBSTDAT,
      c("%d/%m/%Y", "%d-%m-%Y", "%Y-%m-%d", "%d-%b-%y")
    )),
    PREG_START_DATE = ymd(PREG_START_DATE)
  ) %>%
  mutate(
    gest_age = as.numeric(M08_LBSTDAT - PREG_START_DATE)
  ) %>%
  # keep only 0–139 days
  filter(dplyr::between(gest_age, 0, 139)) %>%
  # clean lab values: -5, -7, or < 0 --> NA
  mutate(
    across(
      c(M08_FERRITIN_LBORRES, M08_CRP_LBORRES, M08_AGP_LBORRES),
      ~ ifelse(. %in% c(-5, -7) | . < 0, NA_real_, .)
    )
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    # earliest GA in window for reference
    gest_age = min(gest_age, na.rm = TRUE),
    M08_FERRITIN_LBORRES = first_non_na(M08_FERRITIN_LBORRES),
    M08_CRP_LBORRES      = first_non_na(M08_CRP_LBORRES),
    M08_AGP_LBORRES      = first_non_na(M08_AGP_LBORRES),
    .groups = "drop"
  ) %>%
  # drop moms with no usable lab data at all
  dplyr::filter(
    if_any(
      c(M08_FERRITIN_LBORRES, M08_CRP_LBORRES, M08_AGP_LBORRES),
      ~ !is.na(.)
    )
  )


### g6pd Criteria----
# Step 1: Create G6PD criteria
g6pd_criteria_raw <- mnh08_raw %>% 
  select(SITE, MOMID, PREGID, M08_TYPE_VISIT, M08_RBC_LBPERF_3, M08_RBC_G6PD_LBORRES, M08_LBSTDAT) %>% 
  #filter(M08_RBC_LBPERF_3 == 1 & M08_RBC_G6PD_LBORRES > 0) %>%
  filter(M08_RBC_G6PD_LBORRES > 0) %>% 
  mutate(M08_LBSTDAT =  ymd(parse_date_time(M08_LBSTDAT, c("%d/%m/%Y","%d-%m-%Y","%Y-%m-%d", "%d-%b-%y"))),
    CRIT_G6PD = case_when(
    M08_RBC_G6PD_LBORRES == 77 ~ 55,  # Fix wrong value
    M08_RBC_G6PD_LBORRES >= 6.1 ~ 1,
    M08_RBC_G6PD_LBORRES >= 0 & M08_RBC_G6PD_LBORRES < 6.1 ~ 0,
    TRUE ~ 55
  ))

# Step 2: Select the EARLIEST G6PD result by date
g6pd_criteria <- g6pd_criteria_raw %>%
  group_by(SITE, MOMID, PREGID) %>%
  arrange(M08_LBSTDAT, .by_group = TRUE) %>%  # earliest first
  slice(1) %>%
  ungroup() %>% 
  select(SITE, MOMID, PREGID, M08_RBC_G6PD_LBORRES, M08_LBSTDAT, CRIT_G6PD)

# Query - Step 3: Flag if participant had conflicting G6PD criteria across dates
g6pd_conflict_flag <- g6pd_criteria_raw %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(n_unique = n_distinct(CRIT_G6PD), .groups = "drop") %>%
  mutate(conflict_flag = ifelse(n_unique > 1, TRUE, FALSE)) %>% 
  filter (conflict_flag == TRUE) %>% 
  left_join(g6pd_criteria_raw, by= c("SITE", "MOMID", "PREGID"))

dup_counts <- g6pd_criteria_raw %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

# See how many participants have duplicates
nrow(dup_counts)
### rbc morphology----
rbc_morph_raw  <- mnh08_raw  %>% 
  select(SITE, MOMID, PREGID, M08_RBC_LBPERF_1, M08_RBC_LBPERF_2, M08_RBC_THALA_LBORRES, starts_with("M08_RBC_THALA"), 
         M08_RBC_SICKLE_LBORRES, M08_RBC_SPFY_THALA, M08_LBSTDAT, M08_TYPE_VISIT) %>%
  filter (M08_RBC_LBPERF_1 == 1 | M08_RBC_LBPERF_2 == 1) %>%
  mutate (
    # J. no hemoglobinopathies: SS, SC, SE, EE, CC, SD-Punjab, Sβthal, Eβthal, 
    #Cβthal, CD-Punjab, ED-Punjab, D-D-Punjab, D-Punjabβthal, Thalassemia major, Thalassemia intermedia, or Alpha thalassemia
    CRIT_HEMOGLOBINOPATHIES = case_when(
      # Case 1: Any of the M08_RBC_THALA_x variables is 1 OR grepl() condition is met
      (M08_RBC_THALA_1 == 1 | M08_RBC_THALA_2 == 1 | M08_RBC_THALA_3 == 1 | M08_RBC_THALA_4 == 1 |
         M08_RBC_THALA_5 == 1 | M08_RBC_THALA_6 == 1 | M08_RBC_THALA_7 == 1 | M08_RBC_THALA_8 == 1 |
         M08_RBC_THALA_9 == 1 | M08_RBC_THALA_10 == 1 | M08_RBC_THALA_11 == 1 | M08_RBC_THALA_12 == 1 |
         M08_RBC_THALA_13 == 1 | M08_RBC_THALA_14 == 1 ) ~ 0,
      
      # Case 2: grepl condition with M08_RBC_THALA_19
      (grepl("Interme|Diseas|Major|HbD Punjab", M08_RBC_SPFY_THALA, ignore.case = TRUE) & M08_RBC_THALA_19 == 1) ~ 0,
      
      # Case 3: If M08_RBC_SICKLE_LBORRES is 1, assign 0
      M08_RBC_SICKLE_LBORRES == 1 ~ 0,
      
      # Case 4: All M08_RBC_THALA_x are 0, but M08_RBC_THALA_16, 17, or 18 is 1 OR thala test results are 0
      ((M08_RBC_THALA_1 %in% c(0,77) & M08_RBC_THALA_2 %in% c(0,77) & M08_RBC_THALA_3 %in% c(0,77) & M08_RBC_THALA_4 %in% c(0,77) &
          M08_RBC_THALA_5 %in% c(0,77) & M08_RBC_THALA_6 %in% c(0,77) & M08_RBC_THALA_7 %in% c(0,77) & M08_RBC_THALA_8 %in% c(0,77) &
          M08_RBC_THALA_9 %in% c(0,77) & M08_RBC_THALA_10 %in% c(0,77) & M08_RBC_THALA_11 %in% c(0,77) & M08_RBC_THALA_12 %in% c(0,77) &
          M08_RBC_THALA_13 %in% c(0,77) & M08_RBC_THALA_14 %in% c(0,77)) & 
         (M08_RBC_THALA_15 == 1 | M08_RBC_THALA_16 == 1 | M08_RBC_THALA_17 == 1 | M08_RBC_THALA_18 == 1)) |
        M08_RBC_THALA_LBORRES == 0 ~ 1,
      
      # Case 5: grepl condition with M08_RBC_THALA_19 if it has trait/any regular hemoglobanopathy without disease
      (grepl("TRAIT|AF|FC|AE|AS|HbG|Normal|HbJ|HbF|", M08_RBC_SPFY_THALA, ignore.case = TRUE) & M08_RBC_THALA_19 == 1) ~ 1,
      
      
      # Default case
      TRUE ~ 55
    ))
# Step 2: Remove duplicates (keep one row per participant)
# (Assume: no specific date available ??? pick the record where CRIT_HEMOGLOBINOPATHIES is not 55 first if possible)
rbc_morph_criteria <- rbc_morph_raw %>%
  group_by(SITE, MOMID, PREGID) %>%
  arrange(CRIT_HEMOGLOBINOPATHIES) %>%  # Prioritize 0/1 over 55
  slice(1) %>%
  ungroup()

# Step 3: Flag participants with discrepancies (differing CRIT_HEMOGLOBINOPATHIES across records)
rbc_morph_conflict_flag <- rbc_morph_raw %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(n_unique = n_distinct(CRIT_HEMOGLOBINOPATHIES), .groups = "drop") %>%
  mutate(conflict_flag = ifelse(n_unique > 1, TRUE, FALSE))

##load mnh09 and keep necessary variables---
if (file.exists("derived_data/mnh09.rda")) {
  load("derived_data/mnh09.rda")   # loads mnh09
} else {
  mnh09 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh09_merged.csv")) %>% 
    select(SITE, MOMID, PREGID, num_range("M09_INFANTID_INF",1:4),
           num_range("M09_INFANTID_INF",1:4),
           num_range("M09_BIRTH_DSTERM_INF",1:4), 
           num_range("M09_DELIV_DSSTDAT_INF",1:4))
  save(mnh09, file = "derived_data/mnh09.rda")
}


##load in site reported health cohort ----
# Identify sheets C1–C20
library(readxl)
library(dplyr)
library(purrr)
library(stringr)

healthycohort_file <- "Z:/ReMAPP Healthy Cohort/ReMAPP_Healthy_Cohort_IDs.xlsx"

# Sheets C1–C20 only
target_sheets <- excel_sheets(healthycohort_file)
target_sheets <- target_sheets[target_sheets %in% paste0("C", 1:20)]

# Read each sheet, keep SITE/MOMID/PREGID + its SITE_REPORTED_C*
sheet_list <- map(
  target_sheets,
  ~ {
    sheet_name <- .x
    
    df <- read_excel(healthycohort_file, sheet = sheet_name) %>%
      # drop any "Comments" column
      select(-matches("(?i)^comments$"))
    
    # find the SITE_REPORTED column in this sheet
    site_col <- names(df)[str_detect(names(df), "^SITE_REPORTED")]
    
    # standardize that column name to SITE_REPORTED_C#
    desired_name <- paste0("SITE_REPORTED_", sheet_name)
    
    if (length(site_col) == 0) {
      df[[desired_name]] <- NA
    } else {
      df <- df %>%
        rename(!!desired_name := all_of(site_col[1]))
    }
    
    df %>%
      select(SITE, MOMID, PREGID, all_of(desired_name))
  }
)

# Merge all sheets wide by SITE/MOMID/PREGID
site_reported_hc <- reduce(sheet_list, full_join, by = c("SITE", "MOMID", "PREGID"))




##merge maternal data after remapp luanches----
df_maternal <- screen_df_clean %>%
  filter (ELIGIBLE == 1 & ENROLL_NO_ISSUES == 1) %>% 
  left_join(MAT_ENROLL, by = c("SITE", "MOMID", "PREGID", "SCRNID")) %>% 
  left_join(mnh00_df, by = c("SITE", "SCRNID")) %>% 
  left_join(mnh03, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh04, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh05, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh06, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh08, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(g6pd_criteria, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(rbc_morph_criteria, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh09, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(site_reported_hc, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(bmi_df, by = c("SITE", "MOMID", "PREGID"))

#save data
save(df_maternal, file = "derived_data/df_maternal.rda")

#****************************************************************************
#2. DEFINE HEALTHY CRITERIA ----
#****************************************************************************
##derive criteria intro
df_maternal <- df_maternal %>%
  mutate(FERRITIN_LBORRES  = M08_FERRITIN_LBORRES)

### adjusting ferritin units----
p50_ferritin_gh  <- median(df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Ghana"], na.rm = TRUE)
p50_ferritin_cmc  <- median(df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "India-CMC"], na.rm = TRUE)
p50_ferritin_sas  <- median(df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "India-SAS"], na.rm = TRUE)
p50_ferritin_ky  <- median(df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Kenya"], na.rm = TRUE)
p50_ferritin_pk  <- median(df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Pakistan"], na.rm = TRUE)
p50_ferritin_zm  <- median(df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Zambia"], na.rm = TRUE)

if (p50_ferritin_gh < 10 ) {
  df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Ghana"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Ghana"] * 10
} else  { df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Ghana"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Ghana"]
}


if (p50_ferritin_cmc < 10 ) {
  df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "India-CMC"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "India-CMC"] * 10
} else  { df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "India-CMC"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "India-CMC"]
}


if (p50_ferritin_sas < 10 ) {
  df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "India-SAS"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "India-SAS"] * 10
} else  { df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "India-SAS"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "India-SAS"]
}


if (p50_ferritin_ky < 15 ) {
  df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Kenya"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Kenya"] * 10
} else  { df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Kenya"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Kenya"]
}

if (p50_ferritin_pk < 10 ) {
  df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Pakistan"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Pakistan"] * 10
} else  { df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Pakistan"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Pakistan"]
}


if (p50_ferritin_zm < 10 ) {
  df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Zambia"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Zambia"] * 10
} else  { df_maternal$FERRITIN_LBORRES[df_maternal$SITE == "Zambia"] <- df_maternal$M08_FERRITIN_LBORRES[df_maternal$SITE == "Zambia"]
}

##prep criteria df ----
prep_criteria <- df_maternal %>%
  mutate(
    ## A. age at enrollment----
    # Aged 18 to 34 years
    # Clean ENROLL and BIRTH dates
    ENROLL_SCRN_DATE = ymd(ENROLL_SCRN_DATE),
    M00_BRTHDAT = ymd(M00_BRTHDAT),
    
    # Compute age in years when dates are valid
    AGE_ENROLL = case_when(
      ENROLL_SCRN_DATE %in% c(ymd("1907-07-07"), ymd("1905-05-05")) ~ NA_real_,
      M00_BRTHDAT %in% c(ymd("1907-07-07"), ymd("1905-05-05")) ~ NA_real_,
      TRUE ~ as.numeric(difftime(ENROLL_SCRN_DATE, M00_BRTHDAT, units = "days")) / 365
    ),
    
    # Use estimated age if AGE_ENROLL is NA and estimate is positive
    AGE_ENROLL = if_else(is.na(AGE_ENROLL) & M00_ESTIMATED_AGE > 0, M00_ESTIMATED_AGE, AGE_ENROLL),
    
    # Criteria for age eligibility
    CRIT_AGE = case_when(
      AGE_ENROLL >= 18 & AGE_ENROLL <= 34 ~ 1,
      AGE_ENROLL > 0 & AGE_ENROLL < 18 ~ 0,
      AGE_ENROLL > 34 ~ 0,
      TRUE ~ 55  # missing or unusable
    ),
    ## B. ga at enrollment----
    # gestational age at enrollment - Gestational age <14 weeks 
    BASELINE_GA_WKS = floor(BOE_GA_DAYS_ENROLL/7),
    CRIT_GA = ifelse(BASELINE_GA_WKS > 0 & BASELINE_GA_WKS < 14, 1,
                     ifelse(BASELINE_GA_WKS >= 14 & BASELINE_GA_WKS <=26, 0, 55)),
    
    ## C. pre-pregnancy or early pregnancy body mass index (BMI) of >18.5 and <30 kg/m2 AND mid-upper arm circumference (MUAC) > 23cm [45]----
    # BMI
    BMI = case_when(
      M05_WEIGHT_PERES > 0 & M05_HEIGHT_PERES > 0 ~  M05_WEIGHT_PERES / M05_HEIGHT_PERES / M05_HEIGHT_PERES * 10000, 
      TRUE ~ 55
    ),
    
    BMI_IMPUT = case_when (BMI_IMPUTTED_R <= 18.5 | BMI_IMPUTTED_R >= 30 ~ 0, 
                           BMI_IMPUTTED_R > 18.5 & BMI_IMPUTTED_R < 30 ~ 1, 
                           TRUE ~ 55),
    
    TEMP_BMI = case_when (BMI <= 18.5 | BMI >= 30 ~ 0, 
                          BMI > 18.5 & BMI < 30 ~ 1, 
                          TRUE ~ 55),
    
    # TEMP_BMI_INPUT = ifelse(BMI_IMPUT <= 18.5 | BMI_IMPUT >= 30, 0, 
    #                   ifelse(BMI_IMPUT > 18.5 & BMI_IMPUT < 30, 1, 55)
    # ),
    
    # MUAC mid-upper arm circumference - MUAC
    TEMP_MUAC = case_when (M05_MUAC_PERES > 0 & M05_MUAC_PERES <= 23 ~ 0, 
                           M05_MUAC_PERES > 23 ~ 1, 
                           TRUE ~ 55),
    
    CRIT_BMI_MUAC = case_when(
      #Using the imputted BMI first
      BMI_IMPUT == 1 & TEMP_MUAC == 1 ~ 1, 
      BMI_IMPUT == 0 | TEMP_MUAC == 0 ~ 0, 
      #Then afterward use raw BMI for missing imputted BMI 
      ##This is a temp fix, to follow up with Lili
      TEMP_BMI == 1 & TEMP_MUAC == 1 ~ 1, 
      TEMP_BMI == 0 | TEMP_MUAC == 0 ~ 0, 
      TRUE ~ 55
    ),
    
    # CRIT_BMI_IMPUTTED_MUAC = case_when(
    #   BMI_IMPUT == 1 & TEMP_MUAC == 1 ~ 1, 
    #   BMI_IMPUT == 0 | TEMP_MUAC == 0 ~ 0, 
    #   TRUE ~ 55
    # ),
    ## D. height 150 cm----
    CRIT_HEIGHT = ifelse(M05_HEIGHT_PERES > 0 & M05_HEIGHT_PERES < 150, 0,
                         ifelse(M05_HEIGHT_PERES >= 150, 1, 55)
    ),
    ## E. singleton pregnancy----
    CRIT_SINGLEPREG = ifelse(M06_SINGLETON_PERES == 0, 0,
                             ifelse(M06_SINGLETON_PERES == 1, 1, 55)
    ),
    ## G. no subclinical inflammation (CRP<=5 and/or AGP<=1) check unit (mg/L for CRP and g/L for AGP in dd) double check the calculation before use----
    CRIT_INFLAM = case_when(
      M08_CRP_LBORRES >= 0 & M08_CRP_LBORRES <= 5 & M08_AGP_LBORRES > 0 & M08_AGP_LBORRES <= 1 ~ 1,
      (M08_CRP_LBORRES>5 & !is.na(M08_CRP_LBORRES)) | (M08_AGP_LBORRES>1 & !is.na(M08_AGP_LBORRES)) ~ 0,
      is.na(M08_CRP_LBORRES) | is.na(M08_AGP_LBORRES)  ~ 55,
      TRUE ~ 55
    )
  ) %>% 
  rowwise() %>% 
  mutate_at(vars(starts_with("M06_BP_")), ~ ifelse(. < 0, NA, .)) %>% 
  mutate(

    ### H.a. blood pressure----
    M06_BP_SYS_1 = mean(c(M06_BP_SYS_VSORRES_1, M06_BP_SYS_VSORRES_2, M06_BP_SYS_VSORRES_3), na.rm = TRUE),
    M06_BP_DIA_1 = mean(c(M06_BP_DIA_VSORRES_1, M06_BP_DIA_VSORRES_2, M06_BP_DIA_VSORRES_3), na.rm = TRUE),
    
    CRIT_BP = case_when(
      # If systolic < 140 and diastolic < 90, and both > 0 → eligible
      M06_BP_SYS_1 > 0 & M06_BP_SYS_1 < 140 & 
        M06_BP_DIA_1 > 0 & M06_BP_DIA_1 < 90 ~ 1,
      
      # If systolic ≥ 140 or diastolic ≥ 90 → ineligible
      M06_BP_SYS_1 >= 140 | M06_BP_DIA_1 >= 90 ~ 0,
      
      # All other cases (missing, invalid, zero, etc.) → pending
      TRUE ~ 55)) %>%
  ungroup() %>% 
  mutate(
    ### H.b. no previous low birth weight delivery----
    CRIT_LBW = case_when(
      M04_LOWBIRTHWT_RPORRES == 1 ~ 0,  # had LBW → ineligible
      M04_PH_PREV_RPORRES == 0 | M04_LOWBIRTHWT_RPORRES == 0 ~ 1,  # no previous pregnancy or no LBW → eligible
      M04_LOWBIRTHWT_RPORRES == 99 ~ 0,  # don't know → assume ineligible
      TRUE ~ 55  # all other cases → pending
    ),
    ### H.c. no previous reported stillbirth----
    CRIT_STILLBIRTH = case_when(
      M04_STILLBIRTH_RPORRES == 1 ~ 0,  # stillbirth reported → ineligible
      M04_PH_PREV_RPORRES == 0 | 
        M04_PH_OTH_RPORRES == 0 | 
        M04_STILLBIRTH_RPORRES == 0 ~ 1,  # no prior pregnancy, no other history, or no stillbirth → eligible
      TRUE ~ 55  # all other cases → pending
    ),
    ### H.d. no previous reported unplanned cesarean delivery----
    CRIT_UNPL_CESARIAN = case_when(
      M04_UNPL_CESARIAN_PROCCUR == 1 ~ 0, 
      M04_PH_PREV_RPORRES == 0 | M04_UNPL_CESARIAN_PROCCUR == 0 ~ 1,
      M04_UNPL_CESARIAN_PROCCUR == 99 ~ 0,
      SITE_REPORTED_C8 %in% c(1,0) ~ SITE_REPORTED_C8,
      TRUE ~ 55 
    ),
    ## I. normal glucose-6-phosphate dehydrogenase (≥6.1 U/g Hb)----
    CRIT_G6PD = case_when(
      !is.na(CRIT_G6PD) ~ CRIT_G6PD,
      SITE_REPORTED_C20 %in% c(0, 1)   ~ SITE_REPORTED_C20,
      TRUE                             ~ 55
    ),
    ##J. hemoglobinopathies (made previously)----
    CRIT_HEMOGLOBINOPATHIES = case_when(
      !is.na(CRIT_HEMOGLOBINOPATHIES) ~ CRIT_HEMOGLOBINOPATHIES,
      SITE_REPORTED_C19 %in% c(0, 1)   ~ SITE_REPORTED_C19,
      TRUE                             ~ 55
    ),
    
    ##K. no reported cigarette smoking, tobacco chewing, or betel nut use during pregnancy----
    CRIT_SMOKE = case_when(
      SITE == "Zambia" & (M03_SMOKE_OECOCCUR == 1 | M03_CHEW_OECOCCUR == 1) ~ 0,
      SITE == "Zambia" & (M03_SMOKE_OECOCCUR == 0 & M03_CHEW_OECOCCUR == 0) ~ 1,
      M03_SMOKE_OECOCCUR == 1 | M03_CHEW_BNUT_OECOCCUR == 1 | M03_CHEW_OECOCCUR == 1 ~ 0,
      M03_SMOKE_OECOCCUR == 0 & M03_CHEW_BNUT_OECOCCUR == 0 & M03_CHEW_OECOCCUR == 0 ~ 1,
      SITE_REPORTED_C9 %in% c(1,0) ~ SITE_REPORTED_C9,
      TRUE ~ 55
    ),
    ##L. no reported alcohol consumption during pregnancy----
    CRIT_DRINK = case_when(
      SITE == "Pakistan" ~ 666,  # Not applicable in Pakistan
      M03_DRINK_OECOCCUR == 1 ~ 0,  # Reported drinking → ineligible
      M03_DRINK_OECOCCUR == 0 ~ 1,  # No drinking → eligible
      M03_DRINK_OECOCCUR == 66 ~ 0,  # "Don't ask" → ineligible
      M03_DRINK_OECOCCUR == 77 ~ 0,  # Refused → ineligible
      SITE_REPORTED_C9 %in% c(1,0) ~ SITE_REPORTED_C9,
      TRUE ~ 55  # All other/missing → pending
    ), 
    ## M. no known history or current chronic disease ----
    ##including cancer, kidney disease, and cardiac conditions
    CRIT_CHRONIC = case_when(
      M04_CANCER_EVER_MHOCCUR == 1 | 
        M04_KIDNEY_EVER_MHOCCUR == 1 | 
        M04_CARDIAC_EVER_MHOCCUR == 1 ~ 0,  # Any chronic illness → ineligible
      
      M04_CANCER_EVER_MHOCCUR == 0 & 
        M04_KIDNEY_EVER_MHOCCUR == 0 & 
        M04_CARDIAC_EVER_MHOCCUR == 0 ~ 1,  # no history of chronic illness → eligible
      
      M04_CANCER_EVER_MHOCCUR == 99 | 
        M04_KIDNEY_EVER_MHOCCUR == 99 | 
        M04_CARDIAC_EVER_MHOCCUR == 99 ~ 0,  # Unknown history → ineligible
      
      TRUE ~ 55  # All other cases → pending
    ),
    ##N. no known history or current HIV----
    # if "Record HIV results" = positive, then CRIT_HIV=0 (ineligible) [M06_HIV_POC_LBORRES]
    CRIT_HIV = case_when(
      # if "Record HIV results" = positive (1), then CRIT_HIV=0 (ineligible)
      M06_HIV_POC_LBORRES == 1 ~ 0,
      
      # if "Record HIV results" = negative (0), then CRIT_HIV=1 (eligible)
      M06_HIV_POC_LBORRES == 0 ~ 1,
      
      # if "Have you ever been diagnosed with HIV?" = yes (1)
      # OR if "Have you had any of the following issues since becoming pregnant, HIV" = yes (1),
      # then CRIT_HIV=0 (ineligible)
      M04_HIV_EVER_MHOCCUR == 1 | M04_HIV_MHOCCUR == 1 ~ 0,
      
      # if "Have you ever been diagnosed with HIV?" = no (0)
      # AND if "Have you had any of the following issues since becoming pregnant, HIV" = no (0,77) skip patter allowed
      #PJW comment: Allowing 99 - assuming they have never been diagnosed, I don't think they will receive treatment
      # then CRIT_HIV=1 (eligible)
      M04_HIV_EVER_MHOCCUR == 0 & M04_HIV_MHOCCUR %in% c(0, 77, 99) ~ 1,
      
      # if "Have you ever been diagnosed with HIV?" = 77/99
      # AND if "Have you had any of the following issues since becoming pregnant, HIV" = 99,
      # then CRIT_HIV=0 (assume ineligible)
      M04_HIV_EVER_MHOCCUR %in% c(77, 99) & M04_HIV_MHOCCUR == 99 ~ 0,
      
      # if "Record HIV results" = 55 OR
      # if "Have you ever been diagnosed with HIV?" = 55 OR
      # if "Have you had any of the following issues since becoming pregnant, HIV" = 55,
      # then CRIT_HIV=55 (pending)
      M06_HIV_POC_LBORRES == 55 |
        M04_HIV_EVER_MHOCCUR == 55 |
        M04_HIV_MHOCCUR == 55 ~ 55,
      
      # all other cases → pending
      TRUE ~ 55
    ),
    ##O. no current malaria infection ----
    #(per rapid diagnostic test)
    CRIT_MALARIA = case_when(
      M06_MALARIA_POC_LBORRES == 1  ~ 0,
      M06_MALARIA_POC_LBORRES == 0 | (SITE %in% c("India-CMC", "India-SAS", "Pakistan")) ~ 1,
      TRUE ~ 55
    ),
    ## P. no current Hepatitis B virus infection ----
    #(per rapid diagnostic test)
    CRIT_HEPATITISB = case_when(
      M06_HBV_POC_LBORRES == 1 ~ 0,  # Positive result → ineligible
      M06_HBV_POC_LBORRES == 0 ~ 1,  # Negative result → eligible
      TRUE ~ 55  # All other or missing → pending
    ),
    ## Q. no current Hepatitis C virus infection ----
    #(per rapid diagnostic test)
    CRIT_HEPATITISC = case_when(
      M06_HCV_POC_LBORRES == 1 ~ 0,  # Positive result → ineligible
      M06_HCV_POC_LBORRES == 0 ~ 1,  # Negative result → eligible
      TRUE ~ 55  # All other or missing → pending
    )
  ) 
   
df_criteria <- prep_criteria %>% 
  mutate(
    ## F. no iron deficiency ----
    #(not iron deficient: serum ferritin > 15 mcg/L(Ug/L)) 
    CRIT_IRON = case_when(
      (CRIT_INFLAM ==0 & FERRITIN_LBORRES <70) | (CRIT_INFLAM == 1 & FERRITIN_LBORRES<15) ~ 0, ## 0, ineligible if high inflammation & ferritin <70 OR normal inflammation & ferritin <15
      (CRIT_INFLAM ==0 & FERRITIN_LBORRES >=70) | (CRIT_INFLAM == 1 & FERRITIN_LBORRES>=15) ~ 1, ## 1, eligible
      CRIT_INFLAM == 55 ~  55, ## if inflammation status is pending, 55, pending
      TRUE ~ NA_real_
    ))

save(df_criteria, file = "derived_data/df_criteria.rda")

#**************************************************************************************
#3. ELIGIBILITY DATSETS CREATED  ----
#**************************************************************************************
#code 666 for any not applicable by site
healthyOutcome <- df_criteria %>% 
  rowwise() %>%
  mutate(
    HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE),
  ) %>% 
  # mutate(HEALTHY_ELIGIBLE = case_when(
  #   if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
  #   if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
  #   HEALTHY_CHECK < 20 ~ 3 #20 criteria
  # ) ) %>%
  mutate(
    #!!! temp code for healthy_eligible to allow exclude some criteria
    HEALTHY_ELIGIBLE = case_when(
      CRIT_AGE == 1 &
        # CRIT_GA == 1 & #with GA excluded
        CRIT_BMI_MUAC == 1 &
        CRIT_HEIGHT == 1 &
        CRIT_SINGLEPREG == 1 &
        CRIT_LBW == 1 &
        CRIT_STILLBIRTH == 1 &
        CRIT_UNPL_CESARIAN == 1 &
        CRIT_SMOKE == 1 & 
        CRIT_DRINK %in% c(1,666) &
        CRIT_CHRONIC == 1 &
        CRIT_HIV == 1 &
        CRIT_BP == 1 &
        CRIT_MALARIA == 1 & 
        CRIT_HEPATITISB == 1 & 
        CRIT_HEPATITISC == 1 &
        (CRIT_IRON == 1 | CRIT_IRON == 55) &
        # CRIT_IRON == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_INFLAM == 1 | CRIT_INFLAM == 55) & 
        # CRIT_INFLAM == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_HEMOGLOBINOPATHIES == 1 | CRIT_HEMOGLOBINOPATHIES == 55) &  
        # CRIT_HEMOGLOBINOPATHIES == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_G6PD == 1 | CRIT_G6PD == 55) 
      # CRIT_G6PD == 1 #update to this part once data collection completes (no 55 allowed)
      ~ 1, 
      TRUE ~ 0
    ),
    HEALTHY_ELIGIBLE_GA = case_when(
      CRIT_AGE == 1 &
      CRIT_GA == 1 & #with GA included as a criteria
        CRIT_BMI_MUAC == 1 &
        CRIT_HEIGHT == 1 &
        CRIT_SINGLEPREG == 1 &
        CRIT_LBW == 1 &
        CRIT_STILLBIRTH == 1 &
        CRIT_UNPL_CESARIAN == 1 &
        CRIT_SMOKE == 1 & 
        CRIT_DRINK %in% c(1,666) &
        CRIT_CHRONIC == 1 &
        CRIT_HIV == 1 &
        CRIT_BP == 1 &
        CRIT_MALARIA == 1 & 
        CRIT_HEPATITISB == 1 & 
        CRIT_HEPATITISC == 1 &
       (CRIT_IRON == 1 | CRIT_IRON == 55) &
        # CRIT_IRON == 1 & #update to this part once data collection completes (no 55 allowed)
      (CRIT_INFLAM == 1 | CRIT_INFLAM == 55) & 
        # CRIT_INFLAM == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_HEMOGLOBINOPATHIES == 1 | CRIT_HEMOGLOBINOPATHIES == 55) &  
        # CRIT_HEMOGLOBINOPATHIES == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_G6PD == 1 | CRIT_G6PD == 55) 
      # CRIT_G6PD == 1 #update to this part once data collection completes (no 55 allowed)
      ~ 1, 
      TRUE ~ 0
    ),
    HEALTHY_ELIGIBLE_G6PD = case_when(
      CRIT_AGE == 1 &
        # CRIT_GA == 1 & #with GA excluded
        CRIT_BMI_MUAC == 1 &
        CRIT_HEIGHT == 1 &
        CRIT_SINGLEPREG == 1 &
        CRIT_LBW == 1 &
        CRIT_STILLBIRTH == 1 &
        CRIT_UNPL_CESARIAN == 1 &
        CRIT_SMOKE == 1 & 
        CRIT_DRINK %in% c(1,666) &
        CRIT_CHRONIC == 1 &
        CRIT_HIV == 1 &
        CRIT_BP == 1 &
        CRIT_MALARIA == 1 & 
        CRIT_HEPATITISB == 1 & 
        CRIT_HEPATITISC == 1 &
        (CRIT_IRON == 1 | CRIT_IRON == 55) &
        # CRIT_IRON == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_INFLAM == 1 | CRIT_INFLAM == 55) & 
        # CRIT_INFLAM == 1 & #update to this part once data collection completes (no 55 allowed)
        (CRIT_HEMOGLOBINOPATHIES == 1 | CRIT_HEMOGLOBINOPATHIES == 55)
        # CRIT_HEMOGLOBINOPATHIES == 1 & #update to this part once data collection completes (no 55 allowed)
        #(CRIT_G6PD == 1 | CRIT_G6PD == 55) #REMOVING THIS CRIT
      ~ 1, 
      TRUE ~ 0
    )) %>% 
      mutate(
        HEALTHY_COHORT_COMP = case_when(
          # If any "CRIT_" column equals 55, assign 55
          any(c_across(starts_with("CRIT_")) == 55, na.rm = TRUE) ~ 55,
          
          # If all "CRIT_" columns are 1 or 666, assign 1
          all(c_across(starts_with("CRIT_")) %in% c(1, 666), na.rm = TRUE) ~ 1,
          
          # If any "CRIT_" column is 0, assign 0
          any(c_across(starts_with("CRIT_")) == 0, na.rm = TRUE) ~ 0,
          
          # Default case (optional)
          TRUE ~ NA_real_
        )
      ) %>%
  ungroup()  # Remove row-wise grouping for further operations

df_healthy <- healthyOutcome %>% 
  filter(HEALTHY_ELIGIBLE == 1)

df_healthy_comp <- healthyOutcome %>% 
  filter(HEALTHY_COHORT_COMP == 1)

df_healthy_ga <- healthyOutcome %>% 
  filter(HEALTHY_ELIGIBLE_GA == 1)

df_healthy_g6pd <- healthyOutcome %>% 
  filter(HEALTHY_ELIGIBLE_G6PD == 1)

healthyOutcome_clean <- healthyOutcome %>% 
  select (SCRNID, MOMID, PREGID, SITE, HEALTHY_ELIGIBLE, HEALTHY_ELIGIBLE_GA, 
          HEALTHY_COHORT_COMP,HEALTHY_ELIGIBLE_G6PD,starts_with("CRIT", ignore.case= TRUE))

#for some reason there are Ids (4) which are not in the MAT_ENROLL file. 
#So, to keep things harmonized, I am modifying this
#Using the flowchart dataset, removing the MAT_ENROLL check
all_screen_raw <- screen_df_clean %>% 
  left_join (healthyOutcome_clean, by = c("SCRNID", "MOMID", "PREGID", "SITE")) %>%
  distinct(SCRNID, SITE, .keep_all = TRUE)

#Also removing duplicate eligiblity with the same pregid
# Split rows by presence of keys
with_keys <- all_screen_raw %>% filter(!is.na(SITE) & !is.na(PREGID) & PREGID != "")
without_keys <- all_screen_raw %>% filter(is.na(PREGID) | PREGID == "")

# Deduplicate by SITE + PREGID + ELIGIBLE (drop duplicate eligibilities)
with_keys_dedup <- with_keys %>%
  distinct(SITE, PREGID, ELIGIBLE, .keep_all = TRUE)

# Recombine: keep untouched rows (with missing keys) + deduped rows
all_screen <- bind_rows(without_keys, with_keys_dedup)


table(df_healthy$SITE)
table(df_healthy_comp$SITE)
table(df_healthy_ga$SITE)
table(df_healthy_g6pd$SITE)

save(healthyOutcome, file = "derived_data/healthyOutcome.rda")
save(df_healthy, file = "derived_data/df_healthy.rda")
save(all_screen, file = "derived_data/all_screen.rda")

#*****************************************************************************
#4. HB DATA --> long and wide ----
#*****************************************************************************
#*prepare data
prep_hb1 <- mnh08_raw %>% 
  select(SITE, MOMID, PREGID, 
         M08_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES) %>% 
  #remove duplicates if there is any (PJW: changing to duplicates by date, to account for multiple unscheduled visits)
  group_by(SITE, MOMID, PREGID, M08_LBSTDAT) %>% 
  # group_by(SITE, MOMID, PREGID, M08_TYPE_VISIT) %>% 
  mutate(n = n()) %>% 
  filter(n == 1, M08_TYPE_VISIT %in% c(1:5, 13)) %>% #update to include all visits
  replace_with_na_all(condition = ~.< 0) %>%
  replace_with_na_all(condition = ~. %in% c("1907-07-07", "1905-05-05")) 


#long hb data for all data available
df_hb_long1_all <- prep_hb1 %>%
  right_join(
    healthyOutcome %>% 
      select(
        SITE, MOMID, PREGID, PREG_START_DATE, M03_SMOKE_OECOCCUR,
        BOE_GA_DAYS_ENROLL, HEALTHY_ELIGIBLE, HEALTHY_ELIGIBLE_G6PD,
        HEALTHY_COHORT_COMP, HEALTHY_CHECK
      ),
    by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  #derive variables
  mutate(
    ga_wks = case_when(
      M08_TYPE_VISIT >= 6 ~ NA_real_,
      M08_TYPE_VISIT < 6 ~ as.numeric(ymd(M08_LBSTDAT) - ymd(PREG_START_DATE))/7
    ),
    ga_days = case_when(
      M08_TYPE_VISIT >= 6 ~ NA_real_,
      M08_TYPE_VISIT < 6 ~ as.numeric(ymd(M08_LBSTDAT) - ymd(PREG_START_DATE))
    ),
    trimester = case_when(
      M08_TYPE_VISIT >= 6 ~ NA_real_,
      ga_wks > 0 & ga_wks < 14 ~ 1,
      ga_wks >= 14 & ga_wks < 28 ~ 2,
      ga_wks >= 28 & ga_wks <= 40 ~ 3,
      TRUE ~ NA_real_
    ), 
    EXPECTED_TYPE_VISIT = case_when(
      ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL <= 125 ~ 2, # Only participants who are ≤ 17 weeks at enrollment will have this visit
      ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
      ga_days >= 182 & ga_days <= 216 ~ 3,
      ga_days >= 217 & ga_days <= 237 ~ 4, # Might be updating the 237 number
      ga_days >= 238 & ga_days <= 300 ~ 5, 
      TRUE ~  M08_TYPE_VISIT
    ),
    hb_alti = case_when(
      #new adjustment
      SITE %in% c("Kenya", "Zambia") ~ M08_CBC_HB_LBORRES - 0.8,
      !is.na(M08_CBC_HB_LBORRES) ~ M08_CBC_HB_LBORRES,
      TRUE ~ NA_real_
    ), 
    #adjust for both smoke and altitude
    hb = case_when(
      M03_SMOKE_OECOCCUR == 1 ~ hb_alti - 0.3,
      M03_SMOKE_OECOCCUR == 0 ~ hb_alti,
      TRUE ~ NA_real_
    )) %>% 
  filter(!is.na(PREG_START_DATE) & !is.na(M08_LBSTDAT) & !is.na(M08_CBC_HB_LBORRES))  %>% 
  #to make the TYPE_VISIT more accurate and to remove duplicates
  #we want to take duplicates by actual EXPECTED_VISIT_TYPE with the worst Hb values
  group_by(SITE, MOMID, PREGID, EXPECTED_TYPE_VISIT) %>% 
  slice_max(order_by = hb, n = 1, with_ties = FALSE) %>%  # Keep the row with the highest hb value
  ungroup() %>%
  select(-n)


#wide hb data for all data available 
df_hb_wide1_all <- df_hb_long1_all %>%
  pivot_wider(
    names_from = EXPECTED_TYPE_VISIT,
    values_from = -c("MOMID","PREGID","SITE", PREG_START_DATE, M03_SMOKE_OECOCCUR,
                     HEALTHY_ELIGIBLE,HEALTHY_ELIGIBLE_G6PD, HEALTHY_COHORT_COMP,
                     HEALTHY_CHECK)) %>%
  rowwise() %>%
  ungroup()

#Now keep all data available by if they are healthy cohort (so we do _hc)
df_hb_wide1_hc <- df_hb_wide1_all %>%
  filter (HEALTHY_ELIGIBLE == 1)

df_hb_long1_hc <- df_hb_long1_all %>%
  filter (HEALTHY_ELIGIBLE == 1)

#Healthy cohort hb with G6PD criteria removed
df_hb_wide1_g6pd <- df_hb_wide1_all %>%
  filter (HEALTHY_ELIGIBLE_G6PD == 1)

df_hb_long1_g6pd <- df_hb_long1_all %>%
  filter (HEALTHY_ELIGIBLE_G6PD == 1)

#remove outliers
df_hb_long1 <- df_hb_long1_hc %>% 
  filter(ga_wks >= 5 & ga_wks < 45) %>%
  filter(hb >= 5 & hb <= 18) 

#wide hb data 
df_hb_wide1 <- df_hb_long1 %>%
  pivot_wider(
    names_from = EXPECTED_TYPE_VISIT,
    values_from = -c("MOMID","PREGID","SITE", PREG_START_DATE, M03_SMOKE_OECOCCUR)
  ) %>%
  rowwise() %>%
  ungroup()


test <- df_hb_wide1_all %>% filter (!(PREGID %in% df_hb_wide1$PREGID))

#save data
save(df_hb_long1_all, file = "derived_data/df_hb_long1_all.rda")
save(df_hb_wide1_all, file = "derived_data/df_hb_wide1_all.rda")

save(df_hb_long1_hc, file = "derived_data/df_hb_long1_hc.rda")
save(df_hb_wide1_hc, file = "derived_data/df_hb_wide1_hc.rda")

save(df_hb_long1_g6pd, file = "derived_data/df_hb_long1_g6pd.rda")
save(df_hb_wide1_g6pd, file = "derived_data/df_hb_wide1_g6pd.rda")

save(df_hb_long1, file = "derived_data/df_hb_long1.rda")
save(df_hb_wide1, file = "derived_data/df_hb_wide1.rda")

#*****************************************************************************
#5. SENSITIVITY ANALYISIS DATA ----
#*****************************************************************************
#*prepare data
df_sensitive_long <- df_hb_long1 %>% 
  left_join(INF_OUTCOMES) 

#save data
save(df_sensitive_long, file = "derived_data/df_sensitive_long.rda")

#*****************************************************************************
#6. FINAL ELIGIBILITY DATAFRAMES----
#*****************************************************************************
## Define c1 to c20----
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
      CRIT_DRINK == 1| CRIT_DRINK == 0 ~ CRIT_DRINK,
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
    C19 = ifelse(!is.na(CRIT_HEMOGLOBINOPATHIES),CRIT_HEMOGLOBINOPATHIES,55),
    C20 = ifelse(!is.na(CRIT_G6PD),CRIT_G6PD,55)
  ) %>% 
  mutate(across(matches("C\\d+"),
                function(x) 
                  factor(x, 
                         levels = c(1,0,55),
                         labels = c("Eligible", "Ineligible", "Pending")
                  )))

# df_eli_wide <- df_eli
# save(df_eli_wide, file = "derived_data/df_eli_wide.rda")
# library(openxlsx)
# write.xlsx(df_eli_wide, file = "derived_data/df_eli_wide.xlsx", overwrite = TRUE)

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
             "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20"))

save(df_eli_long, file = "derived_data/df_eli_long.rda")

# test <- df_hb_long1_all %>% filter (SITE == "Kenya")
# table (test$EXPECTED_TYPE_VISIT)

#*****************************************************************************
#7. Demographic Table Prep----
#*************************************************s****************************

if (file.exists("derived_data/inf_perinatal_df.rda")) {
  load("derived_data/inf_perinatal_df.rda") 
  print ("inf_perinatal_df dataset saved already")
} else {
  MAT_INF_OUTCOMES <- read.csv(paste0("Z:/Outcome Data/",UploadDate,"/INF_OUTCOMES.csv"))
  MAT_INF_OUTCOMES <- MAT_INF_OUTCOMES[order(MAT_INF_OUTCOMES$MOMID, MAT_INF_OUTCOMES$PREGID, MAT_INF_OUTCOMES$SITE), ]
  inf_perinatal_df <- MAT_INF_OUTCOMES[!duplicated(MAT_INF_OUTCOMES$PREGID), ]
  
  inf_perinatal_df <- inf_perinatal_df %>% 
    select(MOMID, PREGID, SITE, LIVEBIRTH, PRETERMBIRTH_LT37,
           LBW2500_ANY, BWEIGHT_ANY, SGA_CENTILE, SEX, INF_DTH)
  
  inf_perinatal_df[inf_perinatal_df == -5] <- NA
  inf_perinatal_df[inf_perinatal_df == -7] <- NA
  save(inf_perinatal_df, file = "derived_data/inf_perinatal_df.rda")
}


if (file.exists("derived_data/MAT_DEMOGRAPHIC.rda")) {
  load("derived_data/MAT_DEMOGRAPHIC.rda")   # loads MAT_DEMOGRAPHIC
  print ("MAT_DEMOGRAPHIC dataset saved already")
} else {
  MAT_DEMOGRAPHIC <- read.csv(paste0("Z:/Outcome Data/",UploadDate,"/MAT_DEMOGRAPHIC.csv"))
  save(MAT_DEMOGRAPHIC, file = "derived_data/MAT_DEMOGRAPHIC.rda")
}


