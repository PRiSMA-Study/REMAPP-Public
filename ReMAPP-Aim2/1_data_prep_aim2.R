#*****************************************************************************
#ReMAPP Aim2 Analysis
#Author: Xiaoyan Hu
#Email; xyh@gwu.edu
#Modified: Precious Williams
#Email: williams_pj@gwu.edu
#*****************************************************************************
rm(list = ls())

library(tidyverse)
library(lubridate)
library(growthstandards) ## INTERGROWTH PACKAGE
library(haven)
library(dplyr)
library(readxl)

UploadDate = "2025-10-31"

# Define base directory  
base_dir <- "D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2"

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
iso_data_dir <- file.path(new_dir, "iso_results")

# Create Derived_Data directory if it doesn't exist  
if (!dir.exists(derived_data_dir)) {
  dir.create(derived_data_dir)
}

if (!dir.exists(iso_data_dir)) {
  dir.create(iso_data_dir)s
}
#*****************************************************************************
#1. Load and process all required datasets ----
#*****************************************************************************
# helper: load existing RDA or build+save
load_or_build <- function(obj_name, rda_path, build_fun) {
  if (file.exists(rda_path)) {
    # load returns the names it loaded; we then return the requested object
    load(rda_path, envir = environment())
    return(get(obj_name, envir = environment()))
  } else {
    obj <- build_fun()
    assign(obj_name, obj, envir = environment())
    save(list = obj_name, file = rda_path)
    return(obj)
  }
}


## load MAT_ENROLL ----
MAT_ENROLL <- load_or_build(
  "MAT_ENROLL",
  "derived_data/MAT_ENROLL.rda",
  function() read_dta(paste0("Z:/Outcome Data/", UploadDate, "/MAT_ENROLL.dta"))
)

## load MAT_HDP (kept as RAW for traceability)  ----
MAT_HDP_RAW <- load_or_build(
  "MAT_HDP_RAW",
  "derived_data/MAT_HDP_RAW.rda",
  function() read_dta(paste0("Z:/Outcome Data/", UploadDate, "/MAT_HDP.dta"))
)

# processed MAT_HDP (derived from RAW; no separate RDA unless you want one)
MAT_HDP <- MAT_HDP_RAW %>%
  select(SITE, MOMID, PREGID, PREECLAMPSIA, PREECLAMPSIA_POSTPARTUM, 
         PREECLAMPSIA_GA, PREECLAMPSIA_GA_WK, PREECLAMPSIA_DATE)

## load MAT_NEAR_MISS  ----
MAT_NEAR_MISS <- load_or_build(
  "MAT_NEAR_MISS",
  "derived_data/MAT_NEAR_MISS.rda",
  function() read_dta(paste0("Z:/Outcome Data/", UploadDate, "/MAT_NEAR_MISS.dta"))
)

## load INF_OUTCOMES (composite, preterm37, lbw2500, sga10)  ----
INF_OUTCOMES <- load_or_build(
  "INF_OUTCOMES",
  "derived_data/INF_OUTCOMES.rda",
  function() read.csv(paste0("Z:/Outcome Data/", UploadDate, "/INF_OUTCOMES.csv"))
)

## load INF_COD  ----
INF_COD <- load_or_build(
  "INF_COD",
  "derived_data/INF_COD.rda",
  function() read.csv(paste0("Z:/Outcome Data/", UploadDate, "/INF_COD.csv"))
)

## load MAT_COD  ----
MAT_COD <- load_or_build(
  "MAT_COD",
  "derived_data/MAT_COD.rda",
  function() read.csv(paste0("Z:/Outcome Data/", UploadDate, "/MAT_COD.csv"))
)

## load MAT_ENDPOINTS  ----
MAT_ENDPOINTS <- load_or_build(
  "MAT_ENDPOINTS",
  "derived_data/MAT_ENDPOINTS_raw.rda",   # keep a raw cache, then select below
  function() read_dta(paste0("Z:/Outcome Data/", UploadDate, "/MAT_ENDPOINTS.dta"))
) %>%
  select(SITE, MOMID, PREGID, PREG_END_DATE, PREG_END)


## load MNH03 for smoking adjustment of hb  ----
mnh03 <- load_or_build(
  "mnh03",
  "derived_data/mnh03.rda",
  function() read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh03_merged.csv")) %>%
    select(SITE, MOMID, PREGID, M03_SMOKE_OECOCCUR)
)

## load mnh06 and process ----
# raw_mnh06 for smoking adjustment of hb
raw_mnh06 <- load_or_build(
  "raw_mnh06",
  "derived_data/raw_mnh06.rda",
  function()  read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh06_merged.csv")))


#Basically we want to make sure the right type_visit is being applied accordingly
raw_mnh06_2 <- raw_mnh06 %>% 
  select(SITE, MOMID, PREGID, 
         M06_TYPE_VISIT, M06_HB_POC_LBORRES, M06_SPHB_LBORRES, M06_DIAG_VSDAT) %>% 
  left_join(
    MAT_ENROLL %>% select(SITE, MOMID, PREGID, BOE_GA_DAYS_ENROLL, PREG_START_DATE, ENROLL),
    by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  left_join(
    MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  mutate_all(~ if_else(. < 0, NA, .)) %>%
  mutate(
    date_diff = as.numeric(
      ymd(M06_DIAG_VSDAT) -
        if_else(M06_TYPE_VISIT %in% c(1:5, 13), ymd(PREG_START_DATE), ymd(PREG_END_DATE))
    ),
    ga_days = if_else(M06_TYPE_VISIT %in% c(1:5, 13), date_diff, NA_real_),
    age_days = if_else(M06_TYPE_VISIT %in% c(6:12, 14), date_diff, NA_real_))

prep_mnh06  <- raw_mnh06_2 %>%
  filter(
    ENROLL == 1 &
      (is.na(M06_DIAG_VSDAT) | ymd(M06_DIAG_VSDAT) != ymd("1907-07-07")) &
      (# Pregnancy visits: If PREG_START_DATE is missing, skip this condition
        (!is.na(PREG_START_DATE) &
           ymd(M06_DIAG_VSDAT) >= ymd(PREG_START_DATE) &
           (is.na(PREG_END_DATE) | ymd(M06_DIAG_VSDAT) <= ymd(PREG_END_DATE)) &
           M06_TYPE_VISIT %in% c(1:5, 13)) |
          
          # Post-pregnancy visits: Must be after PREG_END_DATE
          (ymd(M06_DIAG_VSDAT) >= ymd(PREG_END_DATE) &
             M06_TYPE_VISIT %in% c(6:12, 14)))) %>%
  mutate(
    EXPECTED_TYPE_VISIT = case_when(
      # Gestational visit logic
      # Only participants who are 17 weeks at enrollment will have this visit
      ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL < ga_days ~ 2,
      ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
      # ga_days >= 13 & ga_days <= 181 ~ 2,  # Commented out logic
      ga_days >= 182 & ga_days <= 216 ~ 3,
      ga_days >= 217 & ga_days <= 237 ~ 4,  # Might be updating the 237 number
      ga_days >= 238 & ga_days <= 300 ~ 5,
      
      # Postnatal visit types
      M06_TYPE_VISIT == 6 | age_days < 3 ~ 6,
      age_days >= 3 & age_days <= 5 ~ 7,
      age_days >= 7 & age_days <= 14 ~ 8,
      age_days >= 28 & age_days <= 35 ~ 9,
      age_days >= 42 & age_days <= 104 ~ 10,
      age_days >= 104 & age_days <= 279 ~ 11,
      age_days >= 279 & age_days <= 454 ~ 12,
      
      # .5 visits between age-based windows
      age_days > 5 & age_days < 7 ~ 7.5,
      age_days > 14 & age_days < 28 ~ 8.5,
      age_days > 35 & age_days < 42 ~ 9.5,
      
      TRUE ~ 55 ),
    
    EXPECTED_TYPE_VISIT = if_else( EXPECTED_TYPE_VISIT %% 1 == 0,
                                   as.integer(EXPECTED_TYPE_VISIT),EXPECTED_TYPE_VISIT)) %>% 
  filter(ga_days < 300 | age_days < 455) %>% 
  filter (EXPECTED_TYPE_VISIT %in% c (1,2,3,4,5,6,7,8,9,10,11,12)) #till we decide what to do with extra data i.e 7.5,8.5,9.5

#store all non-missing hb data from all sites 
all_poc_hb <- prep_mnh06 %>% filter (!(is.na(M06_HB_POC_LBORRES) & is.na(M06_SPHB_LBORRES))) %>% 
  select(SITE, MOMID, PREGID, -M06_TYPE_VISIT, EXPECTED_TYPE_VISIT, 
         M06_DIAG_VSDAT, M06_HB_POC_LBORRES, M06_SPHB_LBORRES) %>% 
  rename (M06_TYPE_VISIT = EXPECTED_TYPE_VISIT)

### function: resolve_mnh06_hb ----
# Purpose: For each MOMID-PREGID-EXPECTED_TYPE_VISIT group in the mnh08 dataset,
#          this function selects a single Hb value (M06_HB_POC_LBORRES/SPHB) by:
#          - Calculating the median Hb within the group
#          - Measuring the distance of each Hb from that median
#          - Selecting the value farthest from the median (most "extreme")
#          - Randomly choosing one if there are ties
# Notes:
#   - Groups with all missing Hb values are excluded from processing
#   - This helps resolve duplicates by keeping one meaningful, extreme Hb value


resolve_mnh06_hb <- function(df, hb_var, keep_vars = c("SITE", "MOMID", "PREGID", "EXPECTED_TYPE_VISIT", "M06_DIAG_VSDAT")) {
  hb_var <- rlang::ensym(hb_var)
  
  set.seed(100)  # For reproducibility
  df %>%
    group_by(MOMID, PREGID, EXPECTED_TYPE_VISIT) %>%
    filter(any(!is.na(!!hb_var))) %>%
    mutate(
      visit_median = median(!!hb_var, na.rm = TRUE),
      hb_dist = abs(!!hb_var - visit_median),
      max_dist = max(hb_dist, na.rm = TRUE)
    ) %>%
    filter(hb_dist == max_dist) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(all_of(keep_vars), !!hb_var)
}

# Apply to POC Hb values
# For POC
resolved_poc <- prep_mnh06 %>%
  filter(!is.na(M06_HB_POC_LBORRES)) %>%
  resolve_mnh06_hb("M06_HB_POC_LBORRES")

# For SpHb
resolved_sphb <- prep_mnh06 %>%
  filter(!is.na(M06_SPHB_LBORRES)) %>%
  resolve_mnh06_hb("M06_SPHB_LBORRES")


resolved_poc_hb <- full_join(resolved_poc, resolved_sphb,
                               by = c("SITE", "MOMID", "PREGID", "EXPECTED_TYPE_VISIT")) %>%
  mutate(
    M06_DIAG_VSDAT = coalesce(M06_DIAG_VSDAT.x, M06_DIAG_VSDAT.y)
  ) %>%
  select(SITE, MOMID, PREGID, EXPECTED_TYPE_VISIT, M06_DIAG_VSDAT,
         M06_HB_POC_LBORRES, M06_SPHB_LBORRES)


mnh06 <- resolved_poc_hb %>% 
  select(SITE, MOMID, PREGID, EXPECTED_TYPE_VISIT, M06_DIAG_VSDAT, M06_SPHB_LBORRES, M06_HB_POC_LBORRES) %>% 
  rename(M06_TYPE_VISIT = EXPECTED_TYPE_VISIT) %>% 
  pivot_wider(
    names_from = M06_TYPE_VISIT,
    values_from = c(M06_DIAG_VSDAT, M06_SPHB_LBORRES, M06_HB_POC_LBORRES)
  )

## load mnh08 and process ----

#mnh08 for hemoglobin 
raw_mnh08 <- load_or_build(
  "raw_mnh08",
  "derived_data/raw_mnh08.rda",
  function()  read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh08_merged.csv")))


raw_mnh08_2  <- raw_mnh08 %>%
  left_join(
    MAT_ENROLL %>% select(SITE, MOMID, PREGID, BOE_GA_DAYS_ENROLL, PREG_START_DATE, ENROLL),
    by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  left_join(
    MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  mutate_all(~ if_else(. < 0, NA, .)) %>%
  mutate(
    date_diff = as.numeric(
      ymd(M08_LBSTDAT) -
        if_else(M08_TYPE_VISIT %in% c(1:5, 13), ymd(PREG_START_DATE), ymd(PREG_END_DATE))
    ),
    ga_days = if_else(M08_TYPE_VISIT %in% c(1:5, 13), date_diff, NA_real_),
    age_days = if_else(M08_TYPE_VISIT %in% c(6:12, 14), date_diff, NA_real_))
  
prep_mnh08  <- raw_mnh08_2 %>%
  filter(
    ENROLL == 1 &
      (is.na(M08_LBSTDAT) | ymd(M08_LBSTDAT) != ymd("1907-07-07")) &
      (# Pregnancy visits: If PREG_START_DATE is missing, skip this condition
        (!is.na(PREG_START_DATE) &
           ymd(M08_LBSTDAT) >= ymd(PREG_START_DATE) &
           (is.na(PREG_END_DATE) | ymd(M08_LBSTDAT) <= ymd(PREG_END_DATE)) &
           M08_TYPE_VISIT %in% c(1:5, 13)) |
          
          # Post-pregnancy visits: Must be after PREG_END_DATE
          (ymd(M08_LBSTDAT) >= ymd(PREG_END_DATE) &
             M08_TYPE_VISIT %in% c(6:12, 14)))) %>%
      mutate(
        EXPECTED_TYPE_VISIT = case_when(
          # Gestational visit logic
          # Only participants who are 17 weeks at enrollment will have this visit
          ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL < ga_days ~ 2,
          ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
          # ga_days >= 13 & ga_days <= 181 ~ 2,  # Commented out logic
          ga_days >= 182 & ga_days <= 216 ~ 3,
          ga_days >= 217 & ga_days <= 237 ~ 4,  # Might be updating the 237 number
          ga_days >= 238 & ga_days <= 300 ~ 5,
          
          # Postnatal visit types
          M08_TYPE_VISIT == 6 | age_days < 3 ~ 6,
          age_days >= 3 & age_days <= 5 ~ 7,
          age_days >= 7 & age_days <= 14 ~ 8,
          age_days >= 28 & age_days <= 35 ~ 9,
          age_days >= 42 & age_days <= 104 ~ 10,
          age_days >= 104 & age_days <= 279 ~ 11,
          age_days >= 279 & age_days <= 454 ~ 12,
          
          # .5 visits between age-based windows
          age_days > 5 & age_days < 7 ~ 7.5,
          age_days > 14 & age_days < 28 ~ 8.5,
          age_days > 35 & age_days < 42 ~ 9.5,
          
          TRUE ~ 55 ),
        
        EXPECTED_TYPE_VISIT = if_else( EXPECTED_TYPE_VISIT %% 1 == 0,
                              as.integer(EXPECTED_TYPE_VISIT),EXPECTED_TYPE_VISIT)) %>% 
        filter(ga_days < 300 | age_days < 455) %>% 
        filter (EXPECTED_TYPE_VISIT %in% c (1,2,3,4,5,6,7,8,9,10,11,12)) #till we decide what to do with extra data i.e 7.5,8.5,9.5

#store all non-missing hb data from all sites 
all_cbc_hb <- prep_mnh08 %>% filter (!is.na(M08_CBC_HB_LBORRES)) %>% 
  select(SITE, MOMID, PREGID, -M08_TYPE_VISIT, EXPECTED_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES) %>% 
  rename (M08_TYPE_VISIT = EXPECTED_TYPE_VISIT)


### function: resolve_mnh08_hb ----
# Purpose: For each MOMID-PREGID-EXPECTED_TYPE_VISIT group in the mnh08 dataset,
#          this function selects a single Hb value (M08_CBC_HB_LBORRES) by:
#          - Calculating the median Hb within the group
#          - Measuring the distance of each Hb from that median
#          - Selecting the value farthest from the median (most "extreme")
#          - Randomly choosing one if there are ties
# Notes:
#   - Groups with all missing Hb values are excluded from processing
#   - This helps resolve duplicates by keeping one meaningful, extreme Hb value


resolve_mnh08_hb <- function(df) {
  set.seed(100)  # For reproducibility
  df %>%
    group_by(MOMID, PREGID, EXPECTED_TYPE_VISIT) %>%
    # Proceed only for groups with at least one non-NA Hb value
    filter(any(!is.na(M08_CBC_HB_LBORRES))) %>%
    mutate(
      visit_median = median(M08_CBC_HB_LBORRES, na.rm = TRUE),
      hb_dist = abs(M08_CBC_HB_LBORRES - visit_median),
      max_dist = max(hb_dist, na.rm = TRUE)
    ) %>%
    # Filter rows where the distance is equal to the max
    filter(hb_dist == max_dist) %>%
    # If multiple tied, keep one randomly
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-visit_median, -hb_dist, -max_dist)
}

#applying the function to remove the duplicates
mnh08_nodup <- resolve_mnh08_hb(prep_mnh08)

# Check for remaining duplicates
# any_duplicated <- mnh08_nodup %>%
#   add_count(MOMID, PREGID, EXPECTED_TYPE_VISIT) %>%
#   filter(n > 1)
# print(nrow(any_duplicated))  # Should be 0 if duplicates were removed

### mnh08 wide dataset creation ----

mnh08 <- mnh08_nodup %>% 
select(SITE, MOMID, PREGID, -M08_TYPE_VISIT, EXPECTED_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES) %>% 
rename(M08_TYPE_VISIT = EXPECTED_TYPE_VISIT) %>% 
pivot_wider(
  names_from = M08_TYPE_VISIT,
  values_from = c(M08_LBSTDAT, M08_CBC_HB_LBORRES)
)


## mnh25 processing ----
#mnh25 for depression 
#PW changed (restructured based on Savannah's code)
raw_mnh25 <- load_or_build(
  "raw_mnh25",
  "derived_data/raw_mnh25.rda",
  function()  read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh25_merged.csv")))

prep_mnh25 <- raw_mnh25  %>% 
  ungroup() %>% 
  filter(M25_MAT_VISIT_MNH25 %in% c (1,2)) %>% 
  left_join(MAT_ENROLL %>% select(SITE, MOMID, PREGID, PREG_START_DATE, BOE_GA_DAYS_ENROLL, ENROLL)) %>% 
  left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE)) %>% 
  mutate_all(~ if_else(. < 0, NA, .)) %>% 
  ungroup() 

#Prep MNH25 to remove weird date occurrence 
restruc_mnh25 <- prep_mnh25 %>%
  filter(
    ENROLL == 1 &
      (is.na(M25_OBSSTDAT) | ymd(M25_OBSSTDAT) != ymd("1907-07-07")) &
      (
        # Pregnancy visits: If PREG_START_DATE is missing, skip this condition
        (!is.na(PREG_START_DATE) &
           ymd(M25_OBSSTDAT) >= ymd(PREG_START_DATE) &
           (is.na(PREG_END_DATE) | ymd(M25_OBSSTDAT) <= ymd(PREG_END_DATE)) &
           M25_TYPE_VISIT %in% c(1:5, 13)) |
          
          # Post-pregnancy visits: Must be after PREG_END_DATE
          (ymd(M25_OBSSTDAT) >= ymd(PREG_END_DATE) &
             M25_TYPE_VISIT %in% c(6:12, 14))
      )
  ) %>%
  mutate(
    date_diff = as.numeric(
      ymd(M25_OBSSTDAT) -
        if_else(M25_TYPE_VISIT %in% c(1:5, 13), ymd(PREG_START_DATE), ymd(PREG_END_DATE))
    ),
    ga_days = if_else(M25_TYPE_VISIT %in% c(1:5, 13), date_diff, NA_real_),
    age_days = if_else(M25_TYPE_VISIT %in% c(6:12, 14), date_diff, NA_real_),
    
    EXPECTED_TYPE_VISIT = case_when(
      # Gestational visit logic
      # Only participants who are 17 weeks at enrollment will have this visit
      ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL < ga_days ~ 2,
      ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
      # ga_days >= 13 & ga_days <= 181 ~ 2,  # Commented out logic
      ga_days >= 182 & ga_days <= 216 ~ 3,
      ga_days >= 217 & ga_days <= 237 ~ 4,  # Might be updating the 237 number
      ga_days >= 238 & ga_days <= 300 ~ 5,
      
      # Postnatal visit types
      M25_TYPE_VISIT == 6 | age_days < 3 ~ 6,
      age_days >= 3 & age_days <= 5 ~ 7,
      age_days >= 7 & age_days <= 14 ~ 8,
      age_days >= 28 & age_days <= 35 ~ 9,
      age_days >= 42 & age_days <= 104 ~ 10,
      age_days >= 104 & age_days <= 279 ~ 11,
      age_days >= 279 & age_days <= 454 ~ 12,
      
      # .5 visits between age-based windows
      age_days > 5 & age_days < 7 ~ 7.5,
      age_days > 14 & age_days < 28 ~ 8.5,
      age_days > 35 & age_days < 42 ~ 9.5,
      
      TRUE ~ 55 ),
    
    EXPECTED_TYPE_VISIT = if_else( EXPECTED_TYPE_VISIT %% 1 == 0,
                          as.integer(EXPECTED_TYPE_VISIT),EXPECTED_TYPE_VISIT)) %>% 
  
    filter(ga_days < 301 | age_days < 455)

#Step 1: Copy each original variable to a new "_R" column
epds_vars <- paste0("M25_EPDS01", sprintf("%02d", 1:10))  # Q1-Q10

for (var in epds_vars) {
  restruc_mnh25[[paste0(var, "_R")]] <- restruc_mnh25[[var]]
}

# Step 2: Recode based on wording and country
# positive-worded items (Q1, Q2, Q4)
# These questions indicate less depression when answered affirmatively

pos_vars <- c("M25_EPDS0101_R", "M25_EPDS0102_R", "M25_EPDS0104_R")

neg_vars <- c(
  "M25_EPDS0103_R", "M25_EPDS0105_R", "M25_EPDS0106_R",
  "M25_EPDS0107_R", "M25_EPDS0108_R", "M25_EPDS0109_R", "M25_EPDS0110_R"
)

restruc_mnh25 <- restruc_mnh25 %>%
  mutate(across(all_of(pos_vars), ~ case_when(
    SITE %in% c("Ghana", "Kenya") ~ .,
    . == 1 ~ 0,
    . == 2 ~ 1,
    . == 3 ~ 2,
    . == 4 ~ 3,
    . %in% c(55, 77) ~ NA_real_,
    TRUE ~ .
  )))

# Negative-worded items (Q3, Q5-Q10)
# These indicate more depression when answered affirmatively.
restruc_mnh25 <- restruc_mnh25 %>%
  mutate(across(all_of(neg_vars), ~ case_when(
    SITE %in% c("Ghana", "Kenya") ~ .,
    . == 1 ~ 3,
    . == 2 ~ 2,
    . == 3 ~ 1,
    . == 4 ~ 0,
    . %in% c(55, 77) ~ NA_real_,
    TRUE ~ .
  )))

# Step 3: Kenya-specific rule
# In Kenya, all questions are scored the same way:
# 1 = 0, 2 = 1, 3 = 2, 4 = 3, and invalid values become NA.

kenya_vars <- paste0("M25_EPDS01", sprintf("%02d", 1:10), "_R")

restruc_mnh25 <- restruc_mnh25 %>%
  mutate(across(all_of(kenya_vars), ~ case_when(
    SITE == "Kenya" & . == 1 ~ 0,
    SITE == "Kenya" & . == 2 ~ 1,
    SITE == "Kenya" & . == 3 ~ 2,
    SITE == "Kenya" & . == 4 ~ 3,
    SITE == "Kenya" & . %in% c(55, 77, 777) ~ NA_real_,
    TRUE ~ .
  )))

# ???????? Step 4: Ghana-specific rule
# In Ghana, you have separate variables like "M25_EPDS0101_Y" and "M25_EPDS0101_N" instead of one coded response.
# 
# Q1 & Q2 are positive, so:

for (var in c("M25_EPDS0101", "M25_EPDS0102")) {
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana"] <- NA
  
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_N")]] == 2] <- 3
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_N")]] == 1] <- 2
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_Y")]] == 1] <- 1
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_Y")]] == 2] <- 0
}

#Q3-Q10 are negative, so:

for (var in paste0("M25_EPDS01", sprintf("%02d", 3:10))) {
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana"] <- NA
  
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_Y")]] == 2] <- 3
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_Y")]] == 1] <- 2
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_N")]] == 1] <- 1
  restruc_mnh25[[paste0(var, "_R")]][restruc_mnh25$SITE == "Ghana" & restruc_mnh25[[paste0(var, "_N")]] == 2] <- 0
}


# List of recoded variables (Q1-Q10)
epds_recoded <- paste0("M25_EPDS01", sprintf("%02d", 1:10), "_R")

# Calculate total score and number of answered items
restruc_mnh25 <- restruc_mnh25 %>%
  rowwise() %>%
  mutate(
    dep = sum(c_across(all_of(epds_recoded)), na.rm = TRUE),
    q_answered = sum(!is.na(c_across(all_of(epds_recoded)))),
    epds_score = ifelse(q_answered > 0, round((dep / q_answered) * 10, 1), NA_real_),
    depress = case_when(
      SITE == "Ghana" & epds_score >= 11 ~ 1, 
      SITE == "India-CMC" & epds_score >= 11 ~ 1,
      SITE == "India-SAS" & epds_score >= 11 ~ 1,
      SITE == "Kenya" & epds_score >= 11 ~ 1,
      SITE == "Pakistan" & epds_score >= 11 ~ 1,
      SITE == "Zambia" & epds_score >= 11 ~ 1,
      epds_score >= 0 ~ 0,
      TRUE ~ NA_real_)
  ) %>%
  ungroup()

#prep_mnh25 the distribution by site to see if the recoded variables works (seems to)
library(ggplot2)

ggplot(restruc_mnh25, aes(x = epds_score, fill = SITE)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~SITE, scales = "free_y") +
  labs(
    title = "Distribution of EPDS Scores by Site",
    x = "EPDS Score (Scaled 0-10)",
    y = "Count"
  ) +
  theme_minimal()

mnh25_long <- restruc_mnh25 %>%
  mutate(
    ga_wks_25 = case_when(
      EXPECTED_TYPE_VISIT >= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT < 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), ymd(PREG_START_DATE), units = "days")) / 7)
    ), 
    pst_wks_25 = case_when(
      EXPECTED_TYPE_VISIT <= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT > 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), PREG_END_DATE, units = "days")) / 7)
    ),
    M25_TYPE_VISIT = EXPECTED_TYPE_VISIT
  ) %>%
  
  select(SITE, MOMID, PREGID, M25_TYPE_VISIT,  ga_wks_25, pst_wks_25, epds_score, 
         q_answered, depress)

mnh25 <- mnh25_long %>% 
  group_by(MOMID, PREGID, M25_TYPE_VISIT) %>% 
  slice_max(order_by = epds_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pivot_wider(names_from = M25_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID))


depr_mnh25_anc <- restruc_mnh25 %>%
  mutate(
    ga_wks_25 = case_when(
      EXPECTED_TYPE_VISIT >= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT < 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), ymd(PREG_START_DATE), units = "days")) / 7)
    ),
    pst_wks_25 = case_when(
      EXPECTED_TYPE_VISIT <= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT > 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), PREG_END_DATE, units = "days")) / 7)
    ),
    M25_TYPE_VISIT = EXPECTED_TYPE_VISIT
  ) %>%
  select(SITE, MOMID, PREGID, M25_TYPE_VISIT, ga_wks_25, pst_wks_25, epds_score, 
         q_answered, depress) %>%
  filter(!is.na(ga_wks_25) & ga_wks_25 > 28 & !is.na(epds_score)) %>%
  group_by(MOMID, PREGID, SITE) %>%
  slice_max(order_by = ga_wks_25, n = 1, with_ties = FALSE) %>%
  ungroup() 

depr_mnh25_pnc <- restruc_mnh25 %>%
  mutate(
    ga_wks_25 = case_when(
      EXPECTED_TYPE_VISIT >= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT < 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), ymd(PREG_START_DATE), units = "days")) / 7)
    ),
    pst_wks_25 = case_when(
      EXPECTED_TYPE_VISIT <= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT > 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), PREG_END_DATE, units = "days")) / 7)
    ),
    M25_TYPE_VISIT = EXPECTED_TYPE_VISIT
  ) %>%
  select(SITE, MOMID, PREGID, M25_TYPE_VISIT, ga_wks_25, pst_wks_25, epds_score, 
         q_answered, depress) %>%
  filter(!is.na(pst_wks_25) & !is.na(epds_score)) %>%
  group_by(MOMID, PREGID, SITE) %>%
  slice_max(order_by = pst_wks_25, n = 1, with_ties = FALSE) %>%
  ungroup() 

depr_mnh25_pnc6 <- restruc_mnh25 %>%
  mutate(
    ga_wks_25 = case_when(
      EXPECTED_TYPE_VISIT >= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT < 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), ymd(PREG_START_DATE), units = "days")) / 7)
    ),
    pst_wks_25 = case_when(
      EXPECTED_TYPE_VISIT <= 6 ~ NA_real_,
      EXPECTED_TYPE_VISIT > 6 ~ floor(as.numeric(difftime(ymd(M25_OBSSTDAT), PREG_END_DATE, units = "days")) / 7)
    ),
    M25_TYPE_VISIT = EXPECTED_TYPE_VISIT
  ) %>%
  select(SITE, MOMID, PREGID, M25_TYPE_VISIT, ga_wks_25, pst_wks_25, epds_score, 
         q_answered, depress) %>%
  filter(!is.na(pst_wks_25) & !is.na(epds_score)) %>%
  filter(pst_wks_25 >= 6 & pst_wks_25 <= 12) %>%
  group_by(MOMID, PREGID, SITE) %>%
  slice_max(order_by = pst_wks_25, n = 1, with_ties = FALSE) %>%
  ungroup() 


# library(writexl)
# file_path <- paste0("Z:/Outcome Data/", UploadDate, "/Depress_04_04_df.xlsx")
# write_xlsx(mnh25, path = file_path)

## mnh26 processing Prep----
#mnh26 for fatigue
raw_mnh26 <- load_or_build(
  "raw_mnh26",
  "derived_data/raw_mnh26.rda",
  function()  read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh26_merged.csv")))

prep_mnh26 <-  raw_mnh26 %>%
  filter(M26_MAT_VISIT_MNH26 %in% c (1,2)) %>%
  left_join(MAT_ENROLL %>% select(SITE, MOMID, PREGID, PREG_START_DATE, BOE_GA_DAYS_ENROLL, ENROLL)) %>% 
  left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE)) %>% 
  ungroup() %>% 
  mutate_all(~ if_else(. < 0, NA, .)) 

#Prep MNH26 to remove weird date occurrence 
restruc_mnh26 <- prep_mnh26 %>%
  filter(
    ENROLL == 1 &
      (is.na(M26_FTGE_OBSTDAT) | ymd(M26_FTGE_OBSTDAT) != ymd("1907-07-07")) &
      (
        # Pregnancy visits: If PREG_START_DATE is missing, skip this condition
        (!is.na(PREG_START_DATE) &
           ymd(M26_FTGE_OBSTDAT) >= ymd(PREG_START_DATE) &
           (is.na(PREG_END_DATE) | ymd(M26_FTGE_OBSTDAT) <= ymd(PREG_END_DATE)) &
           M26_TYPE_VISIT %in% c(1:5, 13)) |
          
          # Post-pregnancy visits: Must be after PREG_END_DATE
          (ymd(M26_FTGE_OBSTDAT) >= ymd(PREG_END_DATE) &
             M26_TYPE_VISIT %in% c(6:12, 14))
      )
  ) %>%
  mutate(
    date_diff = as.numeric(
      ymd(M26_FTGE_OBSTDAT) -
        if_else(M26_TYPE_VISIT %in% c(1:5, 13), ymd(PREG_START_DATE), ymd(PREG_END_DATE))
    ),
    ga_days = if_else(M26_TYPE_VISIT %in% c(1:5, 13), date_diff, NA_real_),
    age_days = if_else(M26_TYPE_VISIT %in% c(6:12, 14), date_diff, NA_real_),
    
    EXPECTED_TYPE_VISIT = case_when(
      # Gestational visit logic
      # Only participants who are 17 weeks at enrollment will have this visit
      ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL < ga_days ~ 2,
      ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
      # ga_days >= 13 & ga_days <= 181 ~ 2,  # Commented out logic
      ga_days >= 182 & ga_days <= 216 ~ 3,
      ga_days >= 217 & ga_days <= 237 ~ 4,  # Might be updating the 237 number
      ga_days >= 238 & ga_days <= 300 ~ 5,
      
      # Postnatal visit types
      M26_TYPE_VISIT == 6 | age_days < 3 ~ 6,
      age_days >= 3 & age_days <= 5 ~ 7,
      age_days >= 7 & age_days <= 14 ~ 8,
      age_days >= 28 & age_days <= 35 ~ 9,
      age_days >= 42 & age_days <= 104 ~ 10,
      age_days >= 104 & age_days <= 279 ~ 11,
      age_days >= 279 & age_days <= 454 ~ 12,
      
      # .5 visits between age-based windows
      age_days > 5 & age_days < 7 ~ 7.5,
      age_days > 14 & age_days < 28 ~ 8.5,
      age_days > 35 & age_days < 42 ~ 9.5,
      
      TRUE ~ 55 ),
    
    EXPECTED_TYPE_VISIT = if_else( EXPECTED_TYPE_VISIT %% 1 == 0,
                                   as.integer(EXPECTED_TYPE_VISIT),EXPECTED_TYPE_VISIT)) %>% 
  
  filter(ga_days < 301 | age_days < 455)

restruc_mnh26 <- restruc_mnh26 %>% 
  mutate(
    M26_TYPE_VISIT = EXPECTED_TYPE_VISIT,
    ga_wks_26 = case_when(
      M26_TYPE_VISIT >= 6 ~ NA_real_,
      M26_TYPE_VISIT < 6 ~ floor(as.numeric(difftime(ymd(M26_FTGE_OBSTDAT), 
                                ymd(PREG_START_DATE), units = "days")) / 7)), 
    pst_wks_26 = case_when(
      M26_TYPE_VISIT <= 6 ~ NA_real_,
      M26_TYPE_VISIT > 6 ~ floor(as.numeric(difftime(ymd(M26_FTGE_OBSTDAT), 
                                 PREG_END_DATE, units = "days")) / 7))) %>% 
  select(-M26_FTGE_ASSIST, M26_FTGE_OBSTDAT) 


# Assuming your dataframe is called restruc_mnh26
mnh26_df <- restruc_mnh26

# Reverse-coded items don't include AN5 and AN7
reverse_code_vars <- c("M26_FTGE_AN2", "M26_FTGE_HI7", "M26_FTGE_HI12", "M26_FTGE_AN1V",
                       "M26_FTGE_AN3", "M26_FTGE_AN4", "M26_FTGE_AN8", "M26_FTGE_AN12",
                       "M26_FTGE_AN14", "M26_FTGE_AN15", "M26_FTGE_AN16")

# Items not needing reverse code
no_reverse_needed <- c("M26_FTGE_AN5", "M26_FTGE_AN7")

# Destring numeric fields and clean out-of-range values
mnh26_df <- mnh26_df %>%
  mutate(across(all_of(c(reverse_code_vars, no_reverse_needed)), ~ as.numeric(.))) %>%
  mutate(across(all_of(c(no_reverse_needed)), ~ ifelse(!between(., 0, 4), NA, .))) %>%
  mutate(across(all_of(reverse_code_vars), ~ ifelse(!between(., 0, 4), NA, .))) 

# Create reverse-coded variables
for (v in reverse_code_vars) {
  new_var <- paste0(v, "_R")
  mnh26_df[[new_var]] <- 4 - mnh26_df[[v]]
}

# Define final item list (with reversed and original items)
fatigue_items <- c(paste0(reverse_code_vars, "_R"), "M26_FTGE_AN5", "M26_FTGE_AN7")

# Calculate number of items answered, total score, and scaled score
mnh26_df <- mnh26_df %>%
  rowwise() %>%
  mutate(
    questionsanswered = sum(!is.na(c_across(all_of(fatigue_items)))),
    fatigueraw = sum(c_across(all_of(fatigue_items)), na.rm = TRUE),
    fatigue_score = ifelse(questionsanswered > 0, floor(fatigueraw * 13 / questionsanswered), NA_real_),
    fatigued = case_when(fatigue_score <=30 ~ 1, fatigue_score > 30 ~ 0, 
                         TRUE ~ NA_real_)) %>% ungroup()

mnh26 <- mnh26_df %>%
  group_by(MOMID, PREGID, M26_TYPE_VISIT) %>% 
  slice_max(order_by = fatigue_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, M26_TYPE_VISIT, ga_wks_26, pst_wks_26, fatigue_score) %>% 
  pivot_wider(names_from = M26_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID))


fat_mnh26_anc <- mnh26_df %>%
  select(SITE, MOMID, PREGID, M26_TYPE_VISIT, ga_wks_26, pst_wks_26, fatigue_score, 
         fatigued) %>%
  filter(!is.na(ga_wks_26) & ga_wks_26 > 28 & !is.na(fatigue_score)) %>%
  group_by(MOMID, PREGID, SITE) %>%
  slice_max(order_by = ga_wks_26, n = 1, with_ties = FALSE) %>%
  ungroup() 

fat_mnh26_pnc <- mnh26_df %>%
  select(SITE, MOMID, PREGID, M26_TYPE_VISIT, ga_wks_26, pst_wks_26, fatigue_score, 
         fatigued) %>%
  filter(!is.na(pst_wks_26) & !is.na(fatigue_score)) %>%
  group_by(MOMID, PREGID, SITE) %>%
  slice_max(order_by = pst_wks_26, n = 1, with_ties = FALSE) %>%
  ungroup() 

fat_mnh26_pnc6 <- mnh26_df %>%
  select(SITE, MOMID, PREGID, M26_TYPE_VISIT, ga_wks_26, pst_wks_26, fatigue_score, 
         fatigued) %>%
  filter(!is.na(pst_wks_26) & !is.na(fatigue_score)) %>%
  filter(pst_wks_26 >= 6 & pst_wks_26 <= 12) %>%
  group_by(MOMID, PREGID, SITE) %>%
  slice_max(order_by = pst_wks_26, n = 1, with_ties = FALSE) %>%
  ungroup() 

# library(writexl)
# file_path <- paste0("Z:/Outcome Data/", UploadDate, "/Fatigue_04_04_df.xlsx")
# write_xlsx(mnh26, path = file_path)


#prep_mnh26 the distribution by site to see if the recoded variables works (seems to)
library(ggplot2)
ggplot(mnh26_df, aes(x = fatigue_score, fill = SITE)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~SITE, scales = "free_y") +
  labs(
    title = "Distribution of EPDS Scores by Site",
    x = "Fatigue Score",
    y = "Count"
  ) +
  theme_minimal()


## load MAT_HEMORRHAGE (HEM_PPH) ----
MAT_HEMORRHAGE <- load_or_build(
  "MAT_HEMORRHAGE",
  "derived_data/MAT_HEMORRHAGE.rda",
  function() read.csv(paste0("Z:/Outcome Data/",UploadDate,"/MAT_HEMORRHAGE.csv")) %>%
    select(SITE, MOMID, PREGID, HEM_PPH))


## load MAT_PRETERM (PPROM_OCCUR) ----
# MAT_PRETERM <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_PRETERM.dta")) 
MAT_PRETERM <- load_or_build(
  "MAT_PRETERM",
  "derived_data/MAT_PRETERM.rda",
  function() MAT_PRETERM <- 
    read_dta(paste0("Z:/Outcome Data/","2025-04-18","/MAT_PRETERM.dta")))


## combined df_maternal data set ----
df_maternal <- MAT_ENROLL %>%
  mutate(remapp = case_when(
    SITE == "Ghana"      & ENROLL_SCRN_DATE >= ymd("2022-12-28") & ENROLL_SCRN_DATE <= ymd("2024-10-29") ~ 1,
    SITE == "Kenya"      & ENROLL_SCRN_DATE >= ymd("2023-04-03") & ENROLL_SCRN_DATE <= ymd("2025-03-11") ~ 1,
    SITE == "Zambia"     & ENROLL_SCRN_DATE >= ymd("2022-12-15") & ENROLL_SCRN_DATE <= ymd("2025-03-20") ~ 1,
    SITE == "Pakistan"   & ENROLL_SCRN_DATE >= ymd("2022-09-22") & ENROLL_SCRN_DATE <= ymd("2024-01-17") ~ 1,
    SITE == "India-CMC"  & ENROLL_SCRN_DATE >= ymd("2023-06-20") & ENROLL_SCRN_DATE <= ymd("2025-07-01") ~ 1,
    SITE == "India-SAS"  & ENROLL_SCRN_DATE >= ymd("2023-08-15") & ENROLL_SCRN_DATE <= ymd("2025-03-06") ~ 1,
    TRUE ~ NA_real_)) %>%
  filter(remapp == 1) %>%
  left_join(MAT_ENDPOINTS, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh03, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh08, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh06, by = c("SITE", "MOMID", "PREGID"))

## creating a masterlist for all inclusion criteria ----
# enrolment check
### check for Hb missingness (1–12) ----
# Define column groups
cols_all  <- paste0("M08_CBC_HB_LBORRES_", 1:12)
cols_1_5  <- paste0("M08_CBC_HB_LBORRES_", 1:5)
cols_6_12 <- paste0("M08_CBC_HB_LBORRES_", 6:12)

mnh08_track <- mnh08

# safety: ensure columns exist
missing_cols <- setdiff(c(cols_all, "SITE", "MOMID", "PREGID"), names(mnh08_track))
if (length(missing_cols) > 0) {
  stop("These columns are missing from mnh08: ", paste(missing_cols, collapse = ", "))
}

# Add checks (use mnh08_track, not df; and keep 1–12 naming consistent)
mnh08_track$all_missing_1_12 <- rowSums(!is.na(mnh08_track[cols_all]))  == 0
mnh08_track$all_missing_1_5  <- rowSums(!is.na(mnh08_track[cols_1_5]))  == 0
mnh08_track$all_missing_6_12 <- rowSums(!is.na(mnh08_track[cols_6_12])) == 0
mnh08_track$count_real_1_5   <- apply(mnh08_track[cols_1_5], 1, function(x) sum(!is.na(x)))
mnh08_track$count_real_6_12  <- apply(mnh08_track[cols_6_12], 1, function(x) sum(!is.na(x)))

# Keep IDs and check results (select from mnh08_track and use the 1–12 column names)
hb_data_count <- mnh08_track[, c("SITE", "MOMID", "PREGID","all_missing_1_5","all_missing_6_12",
                                 "all_missing_1_12", "count_real_1_5", "count_real_6_12")]

# enrolment check
mat_enrolled_hb <- df_maternal %>% 
  select (MOMID, PREGID, SITE, REMAPP_ENROLL, PREG_END) %>% 
  left_join(hb_data_count, by = c("SITE", "MOMID", "PREGID")) %>% 
  mutate (all_missing_1_12 = case_when (is.na (count_real_6_12) & 
                                          is.na (count_real_1_5) ~ TRUE,
                                        TRUE ~ all_missing_1_12),
          all_missing_1_5 = case_when (is.na (count_real_1_5) ~ TRUE,
                                        TRUE ~ all_missing_1_5),
          all_missing_6_12 = case_when (is.na (count_real_6_12) ~ TRUE,
                                        TRUE ~ all_missing_6_12),
          
          count_real_1_5 = case_when(is.na(count_real_1_5) ~ 0,
                                     TRUE ~ count_real_1_5),
          count_real_6_12 = case_when(is.na(count_real_6_12) ~ 0,
                                      TRUE ~ count_real_6_12)) 
## maternal remapp masterlist ----
df_mat_mlist <-  mat_enrolled_hb   %>% 
  mutate (MISS_HB_ANC = case_when(all_missing_1_5 == TRUE ~ 1,
                                  all_missing_1_5 == FALSE ~ 0,
                                  TRUE ~ 55),
          MISS_HB_PNC = case_when(all_missing_6_12 == TRUE ~ 1,
                                  all_missing_6_12 == FALSE ~ 0,
                                  TRUE ~ 55),
          MISS_HB_ALL = case_when(all_missing_1_12 == TRUE ~ 1,
                                  all_missing_1_12 == FALSE ~ 0,
                                  TRUE ~ 55),
          HB_COUNT_ANC = count_real_1_5,
          HB_COUNT_PNC = count_real_6_12) %>% 
  select (MOMID, PREGID, SITE, REMAPP_ENROLL, PREG_END, MISS_HB_ANC, 
          MISS_HB_PNC, MISS_HB_ALL, HB_COUNT_ANC, HB_COUNT_PNC) %>% 
  distinct(MOMID, PREGID, .keep_all = TRUE)

##maternal analysis masterlist ----
df_mat_analysis <- df_mat_mlist %>% filter (REMAPP_ENROLL == 1 & MISS_HB_ANC == 0) %>% 
  distinct(MOMID, PREGID, .keep_all = TRUE)
  
##infant remapp masterlist ----
inf_enrolled_hb <- INF_OUTCOMES %>% 
  select (SITE, MOMID, PREGID, INFANTID, LIVEBIRTH, 
          BIRTH_OUTCOME_REPORTED, STILLBIRTH_20WK, STILLBIRTH_DENOM) %>% 
  filter(PREGID %in% df_maternal$PREGID) %>%
  left_join(hb_data_count, by = c("SITE", "MOMID", "PREGID"))%>% 
  mutate (all_missing_1_12 = case_when (is.na (count_real_6_12) & 
                                          is.na (count_real_1_5) ~ TRUE,
                                        TRUE ~ all_missing_1_12),
          all_missing_1_5 = case_when (is.na (count_real_1_5) ~ TRUE,
                                       TRUE ~ all_missing_1_5),
          all_missing_6_12 = case_when (is.na (count_real_6_12) ~ TRUE,
                                        TRUE ~ all_missing_6_12),
          count_real_1_5 = case_when(is.na(count_real_1_5) ~ 0,
                                     TRUE ~ count_real_1_5),
          count_real_6_12 = case_when(is.na(count_real_6_12) ~ 0,
                                      TRUE ~ count_real_6_12)) 

df_inf_mlist <-  inf_enrolled_hb   %>% 
  mutate (MISS_HB_ANC = case_when(all_missing_1_5 == TRUE ~ 1,
                                  all_missing_1_5 == FALSE ~ 0,
                                  TRUE ~ 55),
          MISS_HB_PNC = case_when(all_missing_6_12 == TRUE ~ 1,
                                  all_missing_6_12 == FALSE ~ 0,
                                  TRUE ~ 55),
          MISS_HB_ALL = case_when(all_missing_1_12 == TRUE ~ 1,
                                  all_missing_1_12 == FALSE ~ 0,
                                  TRUE ~ 55),
          HB_COUNT_ANC = count_real_1_5,
          HB_COUNT_PNC = count_real_6_12,
          REMAPP_ENROLL = 1) %>% 
  select (SITE, MOMID, PREGID, INFANTID, REMAPP_ENROLL, LIVEBIRTH,
          BIRTH_OUTCOME_REPORTED, STILLBIRTH_DENOM, STILLBIRTH_20WK,
          MISS_HB_ANC, MISS_HB_PNC, MISS_HB_ALL, HB_COUNT_ANC, 
          HB_COUNT_PNC) 


##infant analysis masterlist ----
df_inf_analysis <- df_inf_mlist %>% 
  filter ((STILLBIRTH_20WK == 1 | LIVEBIRTH == 1) & 
            MISS_HB_ANC == 0 & !is.na(INFANTID))

#*****************************************************************************
#2. Hemoglobin data processing ----
#*****************************************************************************
#*prepare data
prep_hb2 <- df_maternal %>% 
  dplyr:: select("SCRNID", "MOMID", "PREGID", "SITE",
                 PREG_START_DATE, PREG_END_DATE,
                 M03_SMOKE_OECOCCUR,
                 num_range("M08_TYPE_VISIT_",1:12), #these include all visits except PNC gap days
                 num_range("M08_CBC_HB_LBORRES_",1:12), #these include all visits except PNC gap days
                 num_range("M08_LBSTDAT_",1:12), #these include all visits except PNC gap days
                 num_range("M06_TYPE_VISIT_",1:12), #these include all visits except PNC gap days
                 num_range("M06_HB_POC_LBORRES_",1:12), #these include all visits except PNC gap days
                 num_range("M06_SPHB_LBORRES_",1:12), #these include all visits except PNC gap days
                 num_range("M06_DIAG_VSDAT_",1:12), #these include all visits except PNC gap days
  ) %>% 
  mutate(across(where(is.integer), ~ as.integer(.))) %>%
  #replace 7s and 5s with NA hb can't be 0
  mutate(across(where(is.numeric), ~ ifelse(. < 0, NA, .))) %>%
  mutate(across(where(is.character), ~ case_when(
    . %in% c("1907-07-07", "1905-05-05") ~ NA_character_,
    TRUE ~ .
  ))) 

##long hb data - basic hb data (no filter)----
df_hb_long2 <- prep_hb2 %>% 
  #to long format
  pivot_longer(
    -c("SCRNID","MOMID","PREGID","SITE", PREG_START_DATE, PREG_END_DATE, M03_SMOKE_OECOCCUR, 
    ),
    names_to = c(".value", "TYPE_VISIT"), 
    names_pattern = "^M\\d{2}_(.+)_(\\d+)"
  ) %>% 
  mutate(
    visit_type = as.numeric(TYPE_VISIT),
    #ga_wks is calculated from MNH08 for general use of hemoglobin (cbc hb)
    ga_wks = case_when(
      visit_type >= 6 ~ NA_real_,
      visit_type < 6 ~ (as.numeric(ymd(LBSTDAT) - ymd(PREG_START_DATE)) / 7)
    ),
    age_wks = case_when(
      visit_type >= 6 ~ floor(as.numeric(ymd(LBSTDAT) - ymd(PREG_END_DATE)) / 7),
      visit_type < 6 ~ NA_real_
    ),
    #ga_wks_06 is calculated from MNH06 for pochb and sphb
    ga_wks_06 = case_when(
      visit_type >= 6 ~ NA_real_,
      visit_type < 6 ~ (as.numeric(ymd(DIAG_VSDAT) - ymd(PREG_START_DATE))/7)
    ),
    #ga_wks_86 is calculated for hb-level or anemia
    ga_wks_86 = case_when(
      !is.na(CBC_HB_LBORRES) ~ ga_wks,
      !is.na(HB_POC_LBORRES) ~ ga_wks_06, 
      TRUE ~ NA_real_
    ),
    trimester = case_when(
      visit_type >= 6 ~ NA_real_,
      ga_wks > 0 & ga_wks < 14 ~ 1,
      ga_wks >= 14 & ga_wks < 28 ~ 2,
      ga_wks >= 28 & ga_wks <= 40 ~ 3,
      TRUE ~ NA_real_
    ), 
    trimester_06 = case_when(
      visit_type >= 6 ~ NA_real_,
      ga_wks_06 > 0 & ga_wks_06 < 14 ~ 1,
      ga_wks_06 >= 14 & ga_wks_06 < 28 ~ 2,
      ga_wks_06 >= 28 & ga_wks_06 <= 40 ~ 3,
      TRUE ~ NA_real_
    ), 
    trimester_86 = case_when(
      !is.na(CBC_HB_LBORRES) ~ trimester,
      !is.na(HB_POC_LBORRES) ~ trimester_06, 
      TRUE ~ NA_real_
    ), 
    # altitude adjustments
    hb_alti     = case_when(
      SITE %in% c("Kenya","Zambia") ~ CBC_HB_LBORRES - 0.8,
      !is.na(CBC_HB_LBORRES)        ~ CBC_HB_LBORRES,
      TRUE ~ NA_real_
    ),
    poc_hb_alti = case_when(
      SITE %in% c("Kenya","Zambia") ~ HB_POC_LBORRES - 0.8,
      !is.na(HB_POC_LBORRES)        ~ HB_POC_LBORRES,
      TRUE ~ NA_real_
    ),
    sphb_alti   = case_when(
      SITE %in% c("Kenya","Zambia") ~ SPHB_LBORRES - 0.8,
      !is.na(SPHB_LBORRES)          ~ SPHB_LBORRES,
      TRUE ~ NA_real_
    ),
    
    # single adjustment value that safely handles missing smoking status
    hb_adj = if_else(M03_SMOKE_OECOCCUR == 1, 0.3, 0, missing = 0),
    
    # apply adjustment
    hb     = hb_alti     - hb_adj,
    poc_hb = poc_hb_alti - hb_adj,
    sphb   = sphb_alti   - hb_adj,
    
    hb_cbc_poc = coalesce(hb, poc_hb),   # simpler/clearer
    hb_source  = case_when(
      !is.na(hb)     ~ "CBC",
      !is.na(poc_hb) ~ "POC",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    #will only calculate anemia < 42 days PNC (visit10 = pnc6)
    hb_level = case_when(
      trimester %in% c(1,3) & hb_cbc_poc > 0 & hb_cbc_poc < 7 ~ 1,
      trimester %in% c(1,3) & hb_cbc_poc >= 7 & hb_cbc_poc < 10 ~ 2,
      trimester %in% c(1,3) & hb_cbc_poc >= 10 & hb_cbc_poc < 11  ~ 3,
      trimester %in% c(1,3) & hb_cbc_poc >= 11 & hb_cbc_poc <= 13 ~ 4, 
      trimester %in% c(1,2,3) & hb_cbc_poc > 13 & hb_cbc_poc < 15 ~ 5, 
      trimester %in% c(1,2,3) & hb_cbc_poc >= 15 ~ 6, 
      trimester == 2 & hb_cbc_poc > 0 & hb_cbc_poc < 7 ~ 1,
      trimester == 2 & hb_cbc_poc >= 7 & hb_cbc_poc < 9.5 ~ 2,
      trimester == 2 & hb_cbc_poc >= 9.5 & hb_cbc_poc < 10.5  ~ 3,
      trimester == 2 & hb_cbc_poc >= 10.5 & hb_cbc_poc <= 13 ~ 4, 
      visit_type <= 9 & hb_cbc_poc > 0 & hb_cbc_poc < 7 ~ 1, 
      visit_type <= 9 & hb_cbc_poc >= 7 & hb_cbc_poc < 10 ~ 2, 
      visit_type <= 9 & hb_cbc_poc >= 10 & hb_cbc_poc < 11 ~ 3, 
      visit_type <= 9 & hb_cbc_poc >= 11 & hb_cbc_poc <= 13 ~ 4, 
      visit_type <= 9 & hb_cbc_poc > 13 & hb_cbc_poc < 15 ~ 5, 
      visit_type <= 9 & hb_cbc_poc >= 15 ~ 6, 
      TRUE ~ NA_real_
    )
  ) %>%
  select(-c(hb_alti, poc_hb_alti, sphb_alti, M03_SMOKE_OECOCCUR, visit_type)) 

df_hb_long = df_hb_long2 

df_hb_long$hb_level <- factor(
  df_hb_long$hb_level, 
  levels = c(1,2,3,4,5,6),
  labels = c("Severe", "Moderate", "Mild", "Normal", "High hb 13-<15g/dl", "High hb >=15g/dl")
)

df_hb_long2 <- df_hb_long2 %>% 
  filter(hb > 0)

#wide hb data 
df_hb_wide <- df_hb_long %>%
  group_by(SITE, MOMID, PREGID, TYPE_VISIT) %>%
  summarise(across(everything(), ~ first(na.omit(.))), .groups = "drop") %>%
  pivot_wider(
    names_from = TYPE_VISIT,
    values_from = -c(SCRNID, MOMID, PREGID, SITE, PREG_START_DATE, PREG_END_DATE),
    names_glue = "{.value}_{TYPE_VISIT}"
  )

df_hb_wide2 <- df_hb_wide %>%
select(-c(starts_with("hb_poc", ignore.case = TRUE), starts_with("hb_cbc", ignore.case = TRUE),
          starts_with("sphb", ignore.case = TRUE)))

## median Hb Value for ANC Visits ----
#************calculate median of the hb value for ANC visits per SITE***********************************
med_hb_anc_site <- df_hb_long2 %>% 
  filter(TYPE_VISIT < 6) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_exm_anc with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_exm_anc <- df_hb_long2 %>% 
  filter(TYPE_VISIT %in% 1:5) %>%
  left_join(med_hb_anc_site, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_anc, file = "derived_data/df_hb_exm_anc.rda")

table(df_hb_exm_anc$trimester)

 #Test if the randomization is working welll and the right trimeter is being chosen
test_trimester_match <- df_hb_exm_anc %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, trimester, ga_wks, TYPE_VISIT) %>%
  left_join(
    df_hb_long2 %>%
      filter(TYPE_VISIT %in% 1:5) %>%
      select(SITE, MOMID, PREGID, hb_orig = hb, trimester, visit_type = TYPE_VISIT ),
    by = c("SITE", "MOMID", "PREGID","trimester")
  ) %>% filter(visit_type == TYPE_VISIT) %>%
mutate (hb_diff = hb_orig - hb)

hist(test_trimester_match$hb_diff)


##median Hb Value for Trimester 1 & 2 ANC Visits ----
#************calculate median of the hb value for T1 & T2 ANC visits per SITE***********************************
med_hb_t12_site <- df_hb_long2 %>% 
  filter(trimester %in% c(1,2)) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_exm_anc with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_exm_t12 <- df_hb_long2 %>% 
  filter(trimester %in% c(1,2)) %>%
  left_join(med_hb_t12_site, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_t12, file = "derived_data/df_hb_exm_t12.rda")

table(df_hb_exm_t12$trimester)

#Test if the randomization is working welll and the right trimeter is being chosen
test_trimester_match_t12 <- df_hb_exm_t12 %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, trimester, ga_wks, TYPE_VISIT) %>%
  left_join(
    df_hb_long2 %>%
      filter(TYPE_VISIT %in% 1:5) %>%
      select(SITE, MOMID, PREGID, hb_orig = hb, trimester, visit_type = TYPE_VISIT ),
    by = c("SITE", "MOMID", "PREGID","trimester")
  ) %>% filter(visit_type == TYPE_VISIT) %>%
  mutate (hb_diff = hb_orig - hb)

hist(test_trimester_match_t12$hb_diff)



##median Hb Value for Trim 1 Visits ----
#************calculate median of the hb value for trimester 1***********************************
med_hb_trim1 <- df_hb_long2 %>% 
  filter(trimester == 1) %>% 
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#********************df_hb_exm_trim1 with one hb_exm value per mom********************
df_hb_exm_trim1 <- df_hb_long2 %>%
  filter(trimester == 1) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  left_join(med_hb_trim1, by = "SITE") %>%
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min)  %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_trim1, file = "derived_data/df_hb_exm_trim1.rda")


##median Hb Value for Trim 2 Visits ----
#************calculate median of the hb value for trimester2***********************************
med_hb_trim2 <- df_hb_long2 %>% 
  filter(trimester == 2) %>% 
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#********************df_hb_exm_trim2 with one hb_exm value per mom********************
df_hb_exm_trim2 <- df_hb_long2 %>% 
  #use ANC hb value only for current task
  filter(trimester == 2) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  left_join(med_hb_trim2, by = "SITE") %>%
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min) %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_trim2, file = "derived_data/df_hb_exm_trim2.rda")
##median Hb Value for Trim 3 Visits ----
#************calculate median of the hb value for trimester 3***********************************
med_hb_trim3 <- df_hb_long2 %>% 
  filter(trimester == 3) %>% 
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#********************df_hb_exm_trim3 with one hb_exm value per mom********************
df_hb_exm_trim3 <- df_hb_long2 %>% 
  #use ANC hb value only for current task
  filter(trimester == 3) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  left_join(med_hb_trim3, by = "SITE") %>%
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min)  %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_trim3, file = "derived_data/df_hb_exm_trim3.rda")

##median Hb Value for 6weeks visits ----
#************calculate median of the hb value for pnc 6weeks***********************************
med_hb_pnc6 <- df_hb_long2 %>% 
  filter(age_wks >= 6 & age_wks <=12) %>% 
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#********************df_hb_exm_trim2 with one hb_exm value per mom********************
df_hb_exm_pnc6 <- df_hb_long2 %>% 
  #use PNC hb value only for current task
  filter(age_wks >= 6 & age_wks <=12) %>% 
  select(SITE, MOMID, PREGID, hb, age_wks) %>% 
  left_join(med_hb_trim2, by = "SITE") %>%
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, age_wks) %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_pnc6, file = "derived_data/df_hb_exm_pnc6.rda")


##maternal anemia dataset creation ----
##New: We want to make sure Anemia severity is being defined as the same across the analysis
#Basically making sure that 
hb_trim1 <- df_hb_exm_trim1 %>% 
  mutate (trimester = 1,
          suffix = "T1") %>% 
  select (SITE, MOMID, PREGID, trimester, suffix, hb)

hb_trim2 <- df_hb_exm_trim2 %>% 
  mutate (trimester = 2,
          suffix = "T2") %>% 
  select (SITE, MOMID, PREGID, trimester, suffix, hb)

hb_trim3 <- df_hb_exm_trim3 %>% 
  mutate (trimester = 3,
          suffix = "T3") %>% 
  select (SITE, MOMID, PREGID, trimester, suffix, hb)

anc_hb_df <- bind_rows(hb_trim1, hb_trim2, hb_trim3) %>%
  mutate ( HB_CBC_PC = as.numeric(hb),
           HB_SOURCE = case_when(is.na(hb) ~ "CBC", TRUE  ~ NA_character_),
           ANEMIA = case_when( trimester %in% c(1,3) & HB_CBC_PC > 0 & HB_CBC_PC < 7 ~ 3, #severe anemia
           trimester %in% c(1,3) & HB_CBC_PC >= 7 & HB_CBC_PC < 10 ~ 2, #moderate anemia
           trimester %in% c(1,3) & HB_CBC_PC >= 10 & HB_CBC_PC < 11  ~ 1, #mild anemia
           trimester %in% c(1,3) & HB_CBC_PC >= 11  ~ 0, #no anemia

           trimester == 2 & HB_CBC_PC > 0 & HB_CBC_PC < 7 ~ 3,
           trimester == 2 & HB_CBC_PC >= 7 & HB_CBC_PC < 9.5 ~ 2,
           trimester == 2 & HB_CBC_PC >= 9.5 & HB_CBC_PC < 10.5  ~ 1,
           trimester == 2 & HB_CBC_PC >= 10.5  ~ 0))

anc_anemia_long <- anc_hb_df %>% 
  mutate (hb_level =  factor(
    ANEMIA, 
    levels = c(3,2,1,0),
    labels = c("Severe", "Moderate", "Mild", "Normal")
  ))
anc_hb_wide <- anc_hb_df %>%
  pivot_wider(
    id_cols = c(SITE, MOMID, PREGID),
    names_from = suffix,
    values_from = c(hb, HB_CBC_PC, HB_SOURCE, ANEMIA),
    names_glue = "{.value}_{suffix}"
  )

anc_hb_wide <- anc_hb_wide %>%
  mutate(
    ANEMIA_ANC = pmax(ANEMIA_T1, ANEMIA_T2, ANEMIA_T3, na.rm = TRUE)
  )

pnc_hb_df <- df_hb_long %>% 
  mutate (TYPE_VISIT = as.numeric(TYPE_VISIT)) %>%
  filter(TYPE_VISIT %in% c(9, 10, 11, 12)) %>% 
  select(SITE, MOMID, PREGID, TYPE_VISIT, hb_cbc_poc, hb_source, hb_level) %>%
  mutate(
    HB_CBC_PC = as.numeric(hb_cbc_poc),
    HB_SOURCE = hb_source,
    ANEMIA = case_when(
      HB_CBC_PC > 0 & HB_CBC_PC < 7 ~ 3, 
      HB_CBC_PC >= 7 & HB_CBC_PC < 10 ~ 2, 
      HB_CBC_PC >= 10 & HB_CBC_PC < 11 ~ 1, 
      HB_CBC_PC >= 11 ~ 0),
    suffix = case_when(
      TYPE_VISIT == 9 ~ "PNC4",
      TYPE_VISIT == 10 ~ "PNC6",
      TYPE_VISIT == 11 ~ "PNC26",
      TYPE_VISIT == 12 ~ "PNC52",  
      TRUE ~ NA_character_
    )
  )%>% 
  select(SITE, MOMID, PREGID, suffix, HB_CBC_PC, HB_SOURCE, ANEMIA, hb_cbc_poc) 


pnc_hb_wide <- pnc_hb_df %>%
  pivot_wider(
    id_cols = c(SITE, MOMID, PREGID),
    names_from = suffix,
    values_from = c(HB_CBC_PC, HB_SOURCE, ANEMIA),
    names_glue = "{.value}_{suffix}"
  )

MAT_ANEMIA <- df_maternal %>% select(SITE, SCRNID, MOMID, PREGID) %>% 
  left_join(anc_hb_wide, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(pnc_hb_wide, by = c("SITE", "MOMID", "PREGID")) 


##function for hb distribution 
plot_hb_by_outcome <- function(data, outcome_var) {
  outcome_sym <- rlang::sym(outcome_var)
  
  # Filter for non-missing outcome and available Hb
  df <- data %>%
    filter(!is.na(!!outcome_sym), !is.na(hb)) %>%
    mutate(outcome = as.factor(!!outcome_sym))
  
  # Summary stats
  summary_stats <- df %>%
    group_by(outcome) %>%
    summarise(
      count = n(),
      mean_hb = mean(hb, na.rm = TRUE),
      median_hb = median(hb, na.rm = TRUE),
      sd_hb = sd(hb, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_stats)
  
  # Plot
  ggplot(df, aes(x = outcome, y = hb, fill = outcome)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    labs(title = paste("Hb Distribution by", outcome_var),
         x = outcome_var,
         y = "Extreme Hb Value (hb)") +
    theme_minimal()
}

#*****************************************************************************
#3. Maternal outcome dataset processing ----
#*****************************************************************************
##maternal composite outcome----

#Potential pregnancy cause of death
# One regex to rule them all (case-insensitive, handles hyphens/spaces/UK-US spellings)
preg_direct_pattern <- paste0(
  "(?i)\\b(",
  paste(c(
    # Hemorrhage
    "obstetric\\s*ha?emorrh?age",
    "post\\s*partum\\s*ha?emorrh?age|\\bPPH\\b",
    "ante\\s*partum\\s*ha?emorrh?age|\\bAPH\\b",
    "placenta\\s*previa",
    "placental?\\s*abruption|abruptio\\s*placentae",
    "retained\\s*placenta",
    # Hypertensive disorders
    "pre[-\\s]*eclampsia|preeclampsia|eclampsia",
    "gestational\\s*hypertension",
    "pregnancy[-\\s]*induced\\s*hypertension|\\bPIH\\b",
    # Sepsis
    "puerperal\\s*sepsis|post\\s*partum\\s*sepsis|septic(a?e|e)?mia|chorioamnionitis",
    # Abortion-related
    "abortion|unsafe\\s*abortion|septic\\s*abortion|incomplete\\s*abortion|termination|miscarriage",
    # Ectopic pregnancy
    "ectopic\\s*pregnan(cy)?",
    # Embolism (pulmonary, amniotic, thromboembolism)
    "pulmonary\\s*embol(ism)?|amniotic\\s*fluid\\s*embol(ism)?|\\bAFE\\b",
    "thromboembol(ism)?|venous\\s*thromboembol(ism)?|\\bVTE\\b|\\bDVT\\b",
    # Other direct obstetric complications
    "uterine\\s*rupture",
    "obstruct(ed)?\\s*labou?r",
    "an?esthesia\\s*complication"
  ), collapse = "|"),
  ")\\b"
)


library(lubridate)
# Robust date cleaner
clean_date <- function(x) {
  # Already Date?
  if (inherits(x, "Date")) return(x)
  
  # Excel serials (numeric)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  
  # Coerce to character, trim, set blanks to NA
  x_chr <- as.character(x)
  x_chr <- str_trim(x_chr)
  x_chr[x_chr == ""] <- NA_character_
  
  # Try multiple common orders (with/without time)
  parsed <- suppressWarnings(parse_date_time(
    x_chr,
    orders = c(
      "Ymd","Y-m-d","Y.m.d",
      "mdY","m/d/Y","m.d.Y","m-B-Y","m-b-Y",
      "dmY","d/m/Y","d.m.Y","d-B-Y","d-b-Y",
      "Ymd HMS","mdy HMS","dmy HMS"
    ),
    exact = FALSE, tz = "UTC"
  ))
  
  # Fallback: 5-digit strings that look like Excel serials
  need_num <- is.na(parsed) & str_detect(x_chr, "^[0-9]{5}$")
  if (any(need_num)) {
    parsed[need_num] <- as.POSIXct(as.Date(as.numeric(x_chr[need_num]), origin = "1899-12-30"))
  }
  
  as.Date(parsed)
}

mat_dth <- MAT_COD %>%
  filter(DEATH_DATE_MISS == 0, COD_MISS == 0) %>%
  select(SITE, MOMID, PREGID, DEATH_DATE, COD, COD_TEXT) %>%
  mutate(
    PREG_CAUSE = case_when(
      COD %in% c(2, 3, 4, 12, 13) ~ 1L,
      COD == 17 & grepl(preg_direct_pattern, COD_TEXT, perl = TRUE) ~ 1L,
      COD == 17 & grepl("(?i)pregnan|maternal|post\\s*partum|puerper(iu|)m", COD_TEXT, perl = TRUE) ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  left_join(MAT_ENDPOINTS, by = c("MOMID","PREGID","SITE")) %>%
  mutate(
    DEATH_DATE    = clean_date(DEATH_DATE),
    PREG_END_DATE = clean_date(PREG_END_DATE),
    DAYS_PP = as.numeric(DEATH_DATE - PREG_END_DATE),
    DEATH_PP42 = case_when(
      # died during pregnancy or within 42 days postpartum
      !is.na(DEATH_DATE) & !is.na(PREG_END_DATE) &
        (DEATH_DATE <= PREG_END_DATE | (DAYS_PP >= 0 & DAYS_PP <= 42)) ~ 1L,
      TRUE ~ 0L
    )
  ) %>% 
  select (SITE, MOMID, PREGID, PREG_CAUSE, DEATH_PP42)

prep_mat_compo <- df_hb_exm_anc %>% 
  left_join(MAT_NEAR_MISS, by = c("SITE", "MOMID","PREGID")) %>% 
  left_join(mat_dth %>% select(SITE, MOMID, PREGID, PREG_CAUSE),
            by = c("SITE", "MOMID","PREGID")) %>% 
  left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE, PREG_END),
            by = c("SITE", "MOMID","PREGID")) %>%
  mutate (  UPLOADDATE    = as.Date(UploadDate),
            DAYS_PP42 = as.numeric(UPLOADDATE - PREG_END_DATE)) %>%
  filter (PREG_END == 1 & DAYS_PP42 >= 42 & is.na(NEARMISS_MISS_FORMS)) %>%
  mutate (  mat_out_compo = case_when (
            #mat_death from preg causes
            PREG_CAUSE == 1 | 
            #severe complications: missing sepsis
            HEM_PPH_SEV==1 | PREECLAMPSIA_SEV == 1 | UTERINE_RUP == 1 |
            #critical intervention
            MAT_ICU == 1 | HYSTERECTOMY == 1 | LAPAROTOMY == 1| TRANSFUSION_NEARMISS == 1|
            #Life threatning conditions: Organ failure/dysfunction
            ORG_FAIL == 1 | ORG_FAIL_HRT == 1 | ORG_FAIL_RESP == 1 | 
            ORG_FAIL_RENAL ==1 |ORG_FAIL_OTHR ==1| ORG_FAIL_LIVER == 1 | 
            ORG_FAIL_NEUR == 1 | ORG_FAIL_UTER == 1 | ORG_FAIL_HEM == 1 ~ 1L,
            TRUE ~ 0L))  %>%
  filter(mat_out_compo %in% c(0, 1))

df_mat_compo  <- prep_mat_compo  %>% 
  select(SITE, MOMID, PREGID, mat_out_compo, starts_with("hb"))
  
df_mat_compo_trim1 <- df_hb_exm_trim1 %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, mat_out_compo)) %>% 
  filter(mat_out_compo %in% c(0, 1))

df_mat_compo_trim2 <- df_hb_exm_trim2 %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, mat_out_compo)) %>% 
  filter(mat_out_compo %in% c(0, 1))

df_mat_compo_trim3 <- df_hb_exm_trim3 %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, mat_out_compo)) %>% 
  filter(mat_out_compo %in% c(0, 1))


df_mat_analysis <- df_mat_analysis %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, MAT_COMPO=mat_out_compo),
            by = c("SITE", "MOMID","PREGID"))

##postpartum hemorrhage----
df_mat_pph <- df_hb_exm_anc %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH),
                                by = c("SITE", "MOMID","PREGID")) %>% 
  filter(HEM_PPH %in% c(0, 1))

df_mat_pph_trim1 <- df_hb_exm_trim1 %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH)) %>% 
  filter(HEM_PPH %in% c(0, 1))

df_mat_pph_trim2 <- df_hb_exm_trim2 %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH)) %>% 
  filter(HEM_PPH %in% c(0, 1))

df_mat_pph_trim3 <- df_hb_exm_trim3 %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH)) %>% 
  filter(HEM_PPH %in% c(0, 1))


df_mat_analysis <- df_mat_analysis %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH),
            by = c("SITE", "MOMID","PREGID"))

##maternal postpartum anemia at PNC6----
df_mat_anemia <- MAT_ANEMIA %>% 
  select(SITE, MOMID, PREGID, ANEMIA_PNC6, ANEMIA_PNC26) %>%
  mutate(
    ppa_pnc6 = case_when(
      ANEMIA_PNC6 %in% c(1:3) ~ 1,
      ANEMIA_PNC6 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    ppa_pnc26 = case_when(
      ANEMIA_PNC26 %in% c(1:3) ~ 1,
      ANEMIA_PNC26 == 0 ~ 0,
      TRUE ~ NA_real_
    )) 

df_mat_ppa_pnc6 <- df_hb_exm_anc %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

df_mat_ppa_pnc6_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

df_mat_ppa_pnc6_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

df_mat_ppa_pnc6_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc6 %in% c(0, 1))


df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_anemia %>% select(SITE, MOMID, PREGID, 
                                     PPA_PNC6 = ppa_pnc6, PPA_PNC26 = ppa_pnc26),
            by = c("SITE", "MOMID","PREGID"))

##maternal postpartum anemia at PNC26----
df_mat_ppa_pnc26 <- df_hb_exm_anc %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_pnc6 <- df_hb_exm_pnc6 %>% 
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

##preterm premature rupture of membranes----
df_mat_pprom <- df_hb_exm_anc %>%
  left_join(MAT_PRETERM %>% select(SITE, MOMID, PREGID, PPROM_OCCUR)) %>%
  mutate(
    pprom = case_when(
      PPROM_OCCUR %in% c(1,0) ~ PPROM_OCCUR,
      TRUE ~ NA_real_
    )) %>%
  filter(pprom >= 0) 

df_mat_pprom_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_pprom %>% select(SITE, MOMID, PREGID, pprom), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(pprom %in% c(0, 1))

df_mat_pprom_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_pprom %>% select(SITE, MOMID, PREGID, pprom), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(pprom %in% c(0, 1))

df_mat_pprom_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_pprom %>% select(SITE, MOMID, PREGID, pprom), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(pprom %in% c(0, 1))

df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_pprom %>% select(SITE, MOMID, PREGID, PPROM = pprom),
            by = c("SITE", "MOMID","PREGID"))

##preeclampsia----
df_mat_preclamp <- df_hb_exm_anc %>%
  left_join(MAT_HDP, by = c("SITE", "MOMID","PREGID")) %>%
  mutate(
    preclamp = case_when(
      PREECLAMPSIA == 1 & 
      PREECLAMPSIA_POSTPARTUM == 0 & 
      PREECLAMPSIA_GA_WK > ga_wks ~ 1,
      PREECLAMPSIA == 0 ~ 0,
      TRUE ~ NA)) %>%
  filter(preclamp >= 0) 

df_mat_preclamp_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_preclamp %>% select(SITE, MOMID, PREGID, preclamp), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(preclamp %in% c(0, 1))

df_mat_preclamp_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_preclamp %>% select(SITE, MOMID, PREGID, preclamp), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(preclamp %in% c(0, 1))

df_mat_preclamp_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_preclamp %>% select(SITE, MOMID, PREGID, preclamp), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(preclamp %in% c(0, 1))

df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_preclamp %>% select(SITE, MOMID, PREGID, PREECLAMPSIA = preclamp),
            by = c("SITE", "MOMID","PREGID"))

##depression----
# ANC Hb with conservative GA weeks
hb_anc_long <- df_hb_long2 %>%
  filter(TYPE_VISIT %in% 1:5, !is.na(hb)) %>%
  select(SITE, MOMID, PREGID, TYPE_VISIT, LBSTDAT, ga_wks, hb, trimester) %>%
  mutate(
    TYPE_VISIT = as.numeric(TYPE_VISIT),
    ga_wks     = floor(as.numeric(ga_wks))   # conservative to weeks (round down)
  )

#Basically we want to choose the dpr score measurement closest to the hb measurement
df_hb_anc_dpr <- hb_anc_long %>%
  full_join(
    restruc_mnh25 %>%
      select(
        SITE, MOMID, PREGID, depress, epds_score,
        TYPE_VISIT = EXPECTED_TYPE_VISIT,
        M25_OBSSTDAT, ga_days_25 = ga_days
      ),
    by = c("SITE", "MOMID", "PREGID", "TYPE_VISIT")
  ) %>%
  filter(hb > 0, depress >= 0) %>%
  mutate(
    ga_wks_25 = floor(as.numeric(ga_days_25) / 7)  # conservative: days -> weeks, round down
  ) %>%
  filter(
    !is.na(ga_wks), !is.na(ga_wks_25),
    abs(ga_wks - ga_wks_25) <= 1                    # within one week of each other
  )

#then or the analysis we want one hb per mom
med_hb_anc_dpr <- df_hb_anc_dpr %>%
  filter(is.finite(hb)) %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")


df_mat_anc_dpr <- df_hb_anc_dpr %>%
  inner_join(med_hb_anc_dpr, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  ungroup() %>%
  select(SITE, MOMID, PREGID, epds_score, depress ,starts_with("hb", ignore.case = TRUE)) %>% 
  filter (PREGID %in% df_maternal$PREGID)


# Make a numeric trimester column once (no function)
df_hb_anc_dpr_trim <- df_hb_anc_dpr %>%
  mutate(
    trim_num = if (is.numeric(trimester)) as.integer(trimester)
    else suppressWarnings(as.integer(str_extract(trimester, "[123]")))
  )  %>%
  filter(hb > 0, epds_score >= 0)

### trimester 1----
df_hb_anc_dpr_trim1 <- df_hb_anc_dpr_trim %>%
  filter(trim_num == 1)

med_hb_anc_dpr_trim1 <- df_hb_anc_dpr_trim1 %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_anc_dpr_trim1 <- df_hb_anc_dpr_trim1 %>%
  inner_join(med_hb_anc_dpr_trim1, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, epds_score, depress ,starts_with("hb", ignore.case = TRUE))


###trimester 2 ----
df_hb_anc_dpr_trim2 <- df_hb_anc_dpr_trim %>%
  filter(trim_num == 2)

med_hb_anc_dpr_trim2 <- df_hb_anc_dpr_trim2 %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_anc_dpr_trim2 <- df_hb_anc_dpr_trim2 %>%
  inner_join(med_hb_anc_dpr_trim2, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, epds_score, depress ,starts_with("hb", ignore.case = TRUE))


###trimester 3 ----
df_hb_anc_dpr_trim3 <- df_hb_anc_dpr_trim %>%
  filter(trim_num == 3)

med_hb_anc_dpr_trim3 <- df_hb_anc_dpr_trim3 %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_anc_dpr_trim3 <- df_hb_anc_dpr_trim3 %>%
  inner_join(med_hb_anc_dpr_trim3, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, epds_score, depress ,starts_with("hb", ignore.case = TRUE))

#Optional: quick duplicate checks
# dupes1 <- df_mat_anc_dpr_trim1 %>% count(SITE, MOMID, PREGID) %>% filter(n > 1)
# dupes2 <- df_mat_anc_dpr_trim2 %>% count(SITE, MOMID, PREGID) %>% filter(n > 1)
# dupes3 <- df_mat_anc_dpr_trim3 %>% count(SITE, MOMID, PREGID) %>% filter(n > 1)

### pnc6 ----
hb_pnc_long <- df_hb_long2 %>%
  filter(TYPE_VISIT == 10 , !is.na(hb)) %>%
  select(SITE, MOMID, PREGID, TYPE_VISIT, LBSTDAT, age_wks, hb) %>%
  mutate(
    TYPE_VISIT = as.numeric(TYPE_VISIT),
    age_wks     = (as.numeric(age_wks))   # conservative to weeks (round down)
  )

df_dpr_pnc <- restruc_mnh25 %>%
  select(
    SITE, MOMID, PREGID, depress, epds_score,
    M25_TYPE_VISIT, EXPECTED_TYPE_VISIT,
    M25_OBSSTDAT, age_days_25 = age_days
  ) %>% mutate(
    age_wks_25 = floor(as.numeric(age_days_25) / 7)  # conservative: days -> weeks, round down
  ) %>%
  filter(EXPECTED_TYPE_VISIT > 6, depress >= 0)
  
#Basically we want to choose the dpr score measurement closest to the hb measurement
df_hb_pnc_dpr <- hb_pnc_long %>%
  left_join(
    df_dpr_pnc,
    by = c("SITE", "MOMID", "PREGID")) %>%
  filter(hb > 0, depress >= 0, TYPE_VISIT > 6 ) %>%
  filter(
    !is.na(age_wks), !is.na(age_wks_25),
    (abs(age_wks - age_wks_25) <= 2)) 

#then or the analysis we want one hb per mom
med_hb_pnc_dpr <- df_hb_pnc_dpr %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_pnc_dpr_pnc6 <- df_hb_pnc_dpr %>%
  inner_join(med_hb_pnc_dpr, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  ungroup() %>%
  select(SITE, MOMID, PREGID, epds_score, depress ,starts_with("hb", ignore.case = TRUE))


df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_anc_dpr %>% select(SITE, MOMID, PREGID, DEPRESS_ANC = depress),
            by = c("SITE", "MOMID","PREGID"))

df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_pnc_dpr_pnc6 %>% select(SITE, MOMID, PREGID, DEPRESS_PNC = depress),
            by = c("SITE", "MOMID","PREGID"))


##fatigue----

#Basically we want to choose the fat score measurement closest to the hb measurement
df_hb_anc_fat <- hb_anc_long %>%
  full_join(
    mnh26_df %>%
      select(
        SITE, MOMID, PREGID, fatigued, fatigue_score,
        TYPE_VISIT = EXPECTED_TYPE_VISIT,
        M26_FTGE_OBSTDAT, ga_days_26 = ga_days
      ),
    by = c("SITE", "MOMID", "PREGID", "TYPE_VISIT")
  ) %>%
  filter(hb > 0, fatigued >= 0) %>%
  mutate(
    ga_wks_26 = floor(as.numeric(ga_days_26) / 7)  # conservative: days -> weeks, round down
  ) %>%
  filter(
    !is.na(ga_wks), !is.na(ga_wks_26),
    abs(ga_wks - ga_wks_26) <= 1                    # within one week of each other
  ) %>%
  filter(hb > 0, fatigue_score >= 0)

#then or the analysis we want one hb per mom
med_hb_anc_fat <- df_hb_anc_fat %>%
  filter(is.finite(hb)) %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")


df_mat_anc_fat <- df_hb_anc_fat  %>%
  inner_join(med_hb_anc_fat, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  ungroup() %>%
  select(SITE, MOMID, PREGID, fatigue_score, fatigued ,starts_with("hb", ignore.case = TRUE))


# Make a numeric trimester column once (no function)
df_hb_anc_fat_trim <- df_hb_anc_fat %>%
  mutate(
    trim_num = as.integer(trimester)
  )  %>%
  filter(hb > 0, fatigue_score >= 0)

### trimester 1 -----
df_hb_anc_fat_trim1 <- df_hb_anc_fat_trim %>%
  filter(trim_num == 1)

med_hb_anc_fat_trim1 <- df_hb_anc_fat_trim1 %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_anc_fat_trim1 <- df_hb_anc_fat_trim1 %>%
  inner_join(med_hb_anc_fat_trim1, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, fatigue_score, fatigued ,starts_with("hb", ignore.case = TRUE))

### trimester 2 --------
df_hb_anc_fat_trim2 <- df_hb_anc_fat_trim %>%
  filter(trim_num == 2)

med_hb_anc_fat_trim2 <- df_hb_anc_fat_trim2 %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_anc_fat_trim2 <- df_hb_anc_fat_trim2 %>%
  inner_join(med_hb_anc_fat_trim2, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, fatigue_score, fatigued ,starts_with("hb", ignore.case = TRUE))


### trimester 3 --------
df_hb_anc_fat_trim3 <- df_hb_anc_fat_trim %>%
  filter(trim_num == 3)

med_hb_anc_fat_trim3 <- df_hb_anc_fat_trim3 %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_anc_fat_trim3 <- df_hb_anc_fat_trim3 %>%
  inner_join(med_hb_anc_fat_trim3, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SITE, MOMID, PREGID, fatigue_score, fatigued ,starts_with("hb", ignore.case = TRUE))


#Optional: quick duplicate checks
# dupes1 <- df_mat_anc_fat_trim1 %>% count(SITE, MOMID, PREGID) %>% filter(n > 1)
# dupes2 <- df_mat_anc_fat_trim2 %>% count(SITE, MOMID, PREGID) %>% filter(n > 1)
# dupes3 <- df_mat_anc_fat_trim3 %>% count(SITE, MOMID, PREGID) %>% filter(n > 1)

### pnc6 ----
df_fat_pnc <- mnh26_df %>%
  select(
    SITE, MOMID, PREGID, fatigued, fatigue_score,
    M26_TYPE_VISIT, EXPECTED_TYPE_VISIT,
    M26_FTGE_OBSTDAT, age_days_26 = age_days
  ) %>% mutate(
    age_wks_26 = floor(as.numeric(age_days_26) / 7)  # conservative: days -> weeks, round down
  ) %>%
  filter(EXPECTED_TYPE_VISIT > 6, fatigue_score >= 0)

#Basically we want to choose the fat score measurement closest to the hb measurement
df_hb_pnc_fat <- hb_pnc_long %>%
  left_join(
    df_fat_pnc,
    by = c("SITE", "MOMID", "PREGID")) %>%
  filter(hb > 0, fatigued >= 0, TYPE_VISIT > 6 ) %>%
  filter(
    !is.na(age_wks), !is.na(age_wks_26),
    (abs(age_wks - age_wks_26) <= 2)) 

#then or the analysis we want one hb per mom
med_hb_pnc_fat <- df_hb_pnc_fat %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_pnc_fat_pnc6 <- df_hb_pnc_fat %>%
  # inner_join(med_hb_pnc_fat, by = "SITE") %>%
  # mutate(
  #   signed_dev = hb - site_hb_median,
  #   abs_dev    = abs(signed_dev)
  # ) %>%
  # group_by(SITE, MOMID, PREGID) %>%
  # slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  # ungroup() %>%
  select(SITE, MOMID, PREGID, fatigue_score, fatigued ,starts_with("hb", ignore.case = TRUE))


df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_anc_fat %>% select(SITE, MOMID, PREGID, FATIGUE_ANC = fatigued),
            by = c("SITE", "MOMID","PREGID"))

df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_pnc_fat_pnc6 %>% select(SITE, MOMID, PREGID, FATIGUE_PNC = fatigued),
            by = c("SITE", "MOMID","PREGID"))


ggplot(df_mat_anc_fat_trim1, aes(x = hb, y = fatigue_score)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(x = "Hemoglobin (g/dL)", 
       y = "Fatigue Score",
       title = "Hemoglobin vs Fatigue Score") +
  theme_minimal()

ggplot(df_mat_anc_fat_trim2, aes(x = hb, y = fatigue_score)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(x = "Hemoglobin (g/dL)", 
       y = "Fatigue Score",
       title = "Hemoglobin vs Fatigue Score") +
  theme_minimal()


ggplot(df_mat_anc_fat_trim3, aes(x = hb, y = fatigue_score)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(x = "Hemoglobin (g/dL)", 
       y = "Fatigue Score",
       title = "Hemoglobin vs Fatigue Score") +
  theme_minimal()

ggplot(df_mat_pnc_fat_pnc6, aes(x = hb, y = fatigue_score)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(x = "Hemoglobin (g/dL)", 
       y = "Fatigue Score",
       title = "Hemoglobin vs Fatigue Score") +
  theme_minimal()
#*****************************************************************************
#4. Infant outcome dataset processing  ----
#*****************************************************************************
df_infant <- INF_OUTCOMES %>% 
  filter(!is.na(INFANTID)) %>% 
  select(SITE, MOMID, PREGID, INFANTID, LIVEBIRTH, BIRTH_OUTCOME_REPORTED,
         LBW2500_PRISMA, LBW1500_PRISMA, LBW_CAT_ANY, BWEIGHT_ANY,
         PRETERMBIRTH_CAT, PRETERMBIRTH_LT37, PRETERMBIRTH_LT34, 
         SGA_CAT, SGA_CENTILE, 
         INF_ASPH,
         INF_PSBI_IPC, INF_PSBI_PNC0, INF_PSBI_PNC1, INF_PSBI_PNC4, INF_PSBI_DENOM, 
         STILLBIRTH_20WK,  STILLBIRTH_22WK, STILLBIRTH_24WK, STILLBIRTH_28WK,
         INF_HYPERBILI_TCB15_24HR, INF_HYPERBILI_TCB15_5DAY, INF_HYPERBILI_TCB15_14DAY,
         INF_HYPERBILI_AAP_24HR, INF_HYPERBILI_AAP_5DAY, INF_HYPERBILI_AAP_14DAY
         ) %>% 
  semi_join(df_maternal, by = "PREGID") %>%      # keep pregnancies present in df_maternal
  filter(LIVEBIRTH == 1) %>%
  inner_join(                                      # be explicit: we only keep rows with an hb match
    df_hb_exm_anc, by = c("SITE", "MOMID", "PREGID")
  ) 

##compo data (preterm37, lbw2500, sga10)----
prep_compo <- df_infant %>% 
  mutate(
    preterm37 = case_when(
      PRETERMBIRTH_CAT %in% c(12,13,14,15) ~ 1,
      PRETERMBIRTH_CAT %in% c(10,11) ~ 0, 
      TRUE ~ NA_real_
    ), 
    lbw2500 = case_when(
      LBW_CAT_ANY %in% c(11,12) ~ 1,
      LBW_CAT_ANY %in% c(13,14) ~ 0, 
      TRUE ~ NA_real_
    ), 
    sga10 = case_when(
      SGA_CAT %in% c(11,12) ~ 1,
      SGA_CAT %in% c(13,14) ~ 0, 
      TRUE ~ NA_real_
    ),
    compo_pre_lbw_sga = case_when(
      preterm37 == 1 | lbw2500 == 1 | sga10 == 1 ~ 1,
      preterm37 == 0 & lbw2500 == 0 & sga10 == 0 ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>% 
  filter(hb > 0)

df_inf_compo <- prep_compo %>% 
  filter(compo_pre_lbw_sga >= 0)

df_inf_compo_trim1 <- df_inf_compo %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_compo_trim2 <- df_inf_compo %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_compo_trim3 <- df_inf_compo %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

plot_hb_by_outcome(df_inf_compo, "compo_pre_lbw_sga")

df_inf_analysis <- df_inf_analysis %>% 
  #select (-COMPOSITE, -LBW2500, -SGA10, -PRETERM37 ) %>%
  left_join(df_inf_compo %>% select(SITE, MOMID, PREGID, INFANTID, 
                                    COMPOSITE = compo_pre_lbw_sga,
                                    LBW2500 = lbw2500,
                                    SGA10 = sga10,
                                    PRETERM37 = preterm37),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##preterm37----
df_inf_preterm37 <- prep_compo %>% 
  filter(preterm37 >= 0 & hb > 0)

df_inf_preterm37_trim1 <- df_inf_preterm37 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_preterm37_trim2 <- df_inf_preterm37 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_preterm37_trim3 <- df_inf_preterm37 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

#plot_hb_by_outcome(df_inf_preterm37, "preterm37")

##lbw2500----

df_inf_lbw2500 <- prep_compo %>% 
  filter(lbw2500 >= 0)

df_inf_lbw2500_trim1 <- df_inf_lbw2500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_lbw2500_trim2 <- df_inf_lbw2500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_lbw2500_trim3 <- df_inf_lbw2500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

##sga10----
df_inf_sga10 <- prep_compo %>% 
  filter(sga10 >= 0)

df_inf_sga10_trim1 <- df_inf_sga10 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_sga10_trim2 <- df_inf_sga10 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_sga10_trim3 <- df_inf_sga10 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)


##preterm34-----

df_inf_preterm34 <- df_infant %>% 
  mutate(preterm34 = case_when(
    PRETERMBIRTH_CAT %in% c(13,14,15) ~ 1,
    PRETERMBIRTH_CAT %in% c(10,11,12) ~ 0, 
    TRUE ~ NA_real_
  )) %>% 
  filter(preterm34 >= 0 & hb > 0)

df_inf_preterm34_trim1 <- df_inf_preterm34 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_preterm34_trim2 <- df_inf_preterm34 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_preterm34_trim3 <- df_inf_preterm34 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-PRETERM34) %>%
  left_join(df_inf_preterm34 %>% select(SITE, MOMID, PREGID, INFANTID, 
                                    PRETERM34 = preterm34),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##lbw1500----

df_inf_lbw1500 <- df_infant %>% 
  mutate(lbw1500 = case_when(
    LBW_CAT_ANY %in% c(11) ~ 1,
    LBW_CAT_ANY %in% c(12,13,14) ~ 0, 
    TRUE ~ NA_real_
  )) %>% 
  filter(lbw1500 >= 0 & hb > 0)

df_inf_lbw1500_trim1 <- df_inf_lbw1500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_lbw1500_trim2 <- df_inf_lbw1500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_lbw1500_trim3 <- df_inf_lbw1500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-LBW1500) %>%
  left_join(df_inf_lbw1500 %>% select(SITE, MOMID, PREGID, INFANTID, 
                                    LBW1500 = lbw1500),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##sga3----
df_inf_sga3 <- df_infant %>% 
  mutate(
    sga3 = case_when(
      SGA_CAT %in% c(11) ~ 1,
      SGA_CAT %in% c(12,13,14) ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>% 
  filter(sga3 >= 0 & hb > 0)

df_inf_sga3_trim1 <- df_inf_sga3 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_sga3_trim2 <- df_inf_sga3 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_sga3_trim3 <- df_inf_sga3 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-SGA3) %>%
  left_join(df_inf_sga3 %>% select(SITE, MOMID, PREGID, INFANTID, 
                                      SGA3 = sga3),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##neonatal psbi (defined at IPC and PNC0 visits)----
df_inf_psbi <- df_infant %>% 
  mutate(
    inf_psbi = case_when(
      INF_PSBI_IPC == 1 | INF_PSBI_PNC0 == 1 | INF_PSBI_PNC1 == 1 | INF_PSBI_PNC4 == 1 ~ 1,
      INF_PSBI_IPC == 0 & INF_PSBI_PNC0 == 0 & INF_PSBI_PNC1 == 0 & INF_PSBI_PNC4 == 0 ~ 0, 
      TRUE ~ NA_real_
    ) 
  ) %>% 
  filter(hb > 0 & inf_psbi >= 0 & INF_PSBI_DENOM == 1)

df_inf_psbi_trim1 <- df_inf_psbi %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_psbi_trim2 <- df_inf_psbi %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_psbi_trim3 <- df_inf_psbi %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-PSBI) %>%
  left_join(df_inf_psbi %>% select(SITE, MOMID, PREGID, INFANTID, 
                                   PSBI = inf_psbi),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##neonatal asphyxia (defined at IPC and PNC0 visits)----
df_inf_asph <- df_infant %>% 
  mutate(inf_asph = ifelse(INF_ASPH %in% c(0,1), INF_ASPH, NA_real_)) %>% 
  filter(hb > 0 & inf_asph >= 0)

df_inf_asph_trim1 <- df_inf_asph %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_asph_trim2 <- df_inf_asph %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_asph_trim3 <- df_inf_asph %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-ASPH) %>%
  left_join(df_inf_asph %>% select(SITE, MOMID, PREGID, INFANTID, 
                                   ASPH = inf_asph),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##stillbirth20wks: death before delivery at GA >=20----

###median Hb Value for Stillbirth 20weeks Visits
#************calculate median of the hb value for ANC visits per SITE***********************************
med_hb_stbirth20_site <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks <= 20) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth20 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth20 <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks <= 20) %>%
  left_join(med_hb_stbirth20_site, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm) 

save(df_hb_stbirth20, file = "derived_data/df_hb_stbirth20.rda")



###median Hb Value for Stillbirth 20weeks Trimester 1 Visits ----
#*****calculate median hb value for Stillbirth Trimester 1 visits per SITE******
med_hb_stbirth20_trim1 <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks < 14) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth20_trim1 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth20_trim1 <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks < 14) %>%
  left_join(med_hb_stbirth20_trim1, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth20_trim1, file = "derived_data/df_hb_stbirth20_trim1.rda")

table(df_hb_stbirth20_trim1$trimester)

###median Hb Value for Stillbirth 20weeks Trimester 2 Visits ----
#*****calculate median hb value for Stillbirth Trimester 2 visits per SITE******
med_hb_stbirth20_trim2 <- df_hb_long2 %>% 
  filter(ga_wks >= 14 & ga_wks < 20) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth20_trim1 with one hb_exm value per mom ************
set.seed(100)
df_hb_stbirth20_trim2 <- df_hb_long2 %>% 
  filter(ga_wks >= 14 & ga_wks < 20) %>%
  left_join(med_hb_stbirth20_trim2, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth20_trim2, file = "derived_data/df_hb_stbirth20_trim2.rda")

table(df_hb_stbirth20_trim2$trimester)


infant_allbirths <- INF_OUTCOMES %>% 
  filter (STILLBIRTH_20WK == 1 | LIVEBIRTH == 1) %>% 
  filter (PREGID %in% df_maternal$PREGID)

#Stillbirth 20weeks outcome dataset generation
df_inf_stillbirth20 <- infant_allbirths %>% 
  left_join(       # be explicit: we only keep rows with an hb match
    df_hb_stbirth20, by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  mutate(inf_stillbirth20 = case_when(
    STILLBIRTH_20WK == 1 ~ 1,
    STILLBIRTH_20WK == 0 ~ 0, 
    TRUE ~ NA_real_
  )) %>% 
  filter(inf_stillbirth20 %in% c(0,1) & !is.na(hb)) %>% 
  select (SITE, MOMID, PREGID, INFANTID, inf_stillbirth20, starts_with("hb"))


test_still <- df_inf_stillbirth20 %>% 
  filter ((inf_stillbirth20 %in% c(0,1) & is.na(hb)))

df_inf_analysis <- df_inf_analysis %>% 
  #select (-INF_STILLBIRTH20) %>%
  left_join(df_inf_stillbirth20 %>% select(SITE, MOMID, PREGID, INFANTID, 
                                           INF_STILLBIRTH20 = inf_stillbirth20),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

df_inf_stillbirth20_trim1 <- df_inf_stillbirth20 %>% 
  select(-starts_with("hb")) %>% 
  inner_join(df_hb_stbirth20_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_stillbirth20_trim2 <- df_inf_stillbirth20 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_stbirth20_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)


##stillbirth28 weeks: death before delivery at GA >=28----
#Median Hb Value for Stillbirth 28weeks Visits
#************calculate median of the hb value for ANC visits per SITE***********************************
med_hb_stbirth28_site <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks <= 28) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth28 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth28 <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks <= 28) %>%
  left_join(med_hb_stbirth28_site, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth28, file = "derived_data/df_hb_stbirth28.rda")

table(df_hb_stbirth28$trimester)

###median Hb Value for Stillbirth 28weeks Trimester 1 Visits ----
#*calculate median hb value for Stillbirth Trimester 1 visits per SITE
med_hb_stbirth28_trim1 <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks < 14) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth28_trim1 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth28_trim1 <- df_hb_long2 %>% 
  filter(ga_wks > 0 & ga_wks < 14) %>%
  left_join(med_hb_stbirth28_trim1, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth28_trim1, file = "derived_data/df_hb_stbirth28_trim1.rda")

table(df_hb_stbirth28_trim1$trimester)

###median Hb Value for Stillbirth 28weeks Trimester 2 Visits ----
#*****calculate median hb value for Stillbirth Trimester 2 visits per SITE******
med_hb_stbirth28_trim2 <- df_hb_long2 %>% 
  filter(ga_wks >= 14 & ga_wks < 28) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth28_trim1 with one hb_exm value per mom ************
set.seed(100)
df_hb_stbirth28_trim2 <- df_hb_long2 %>% 
  filter(ga_wks >= 14 & ga_wks < 28) %>%
  left_join(med_hb_stbirth28_trim2, by = "SITE") %>%  # Join median hb early
  mutate(hb_dis = abs(hb - hb_med)) %>%       # Calculate deviation
  
  # Get max and min hb values per pregnancy
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>%
  
  # Find record with maximum deviation
  filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
  slice_sample(n = 1) %>%  # Randomly select if ties
  
  # Final formatting
  mutate(hb_exm = round(hb, 1)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, 
         trimester, ga_wks, TYPE_VISIT) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth28_trim2, file = "derived_data/df_hb_stbirth28_trim2.rda")

table(df_hb_stbirth28_trim2$trimester)


#Stillbirth 28weeks outcome dataset generation
df_inf_stillbirth28 <- infant_allbirths %>% 
  left_join(                                      # be explicit: we only keep rows with an hb match
    df_hb_stbirth28, by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  mutate(inf_stillbirth28 = case_when(
    STILLBIRTH_28WK == 1 ~ 1,
    STILLBIRTH_28WK == 0 ~ 0, 
    TRUE ~ NA_real_
  )) %>% 
  filter(hb > 0 & inf_stillbirth28 >= 0 ) %>% 
  select (SITE, MOMID, PREGID, INFANTID, inf_stillbirth28, starts_with("hb"))


df_inf_analysis <- df_inf_analysis %>% 
  #select (-INF_STILLBIRTH28) %>%
  left_join(df_inf_stillbirth28 %>% select(SITE, MOMID, PREGID, INFANTID, 
                                   INF_STILLBIRTH28 = inf_stillbirth28),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))


df_inf_stillbirth28_trim1 <- df_inf_stillbirth28 %>% 
  select(-starts_with("hb")) %>% 
  inner_join(df_hb_stbirth28_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_stillbirth28_trim2 <- df_inf_stillbirth28 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_stbirth28_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

##neonatal hyperbilirubinemia----
df_inf_hyperbili <- df_infant %>% 
  mutate(hyperbili = case_when(
    INF_HYPERBILI_AAP_24HR == 1 | INF_HYPERBILI_AAP_5DAY == 1 | INF_HYPERBILI_AAP_14DAY == 1 ~ 1,
    INF_HYPERBILI_AAP_24HR == 0 & INF_HYPERBILI_AAP_5DAY == 0 & INF_HYPERBILI_AAP_14DAY == 0 ~ 0,
    TRUE ~ NA_real_
  )) %>% 
  filter(hb > 0 & hyperbili >= 0)

df_inf_hyperbili_trim1 <- df_inf_hyperbili %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_hyperbili_trim2 <- df_inf_hyperbili %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_hyperbili_trim3 <- df_inf_hyperbili %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-HYPERBILI) %>%
  left_join(df_inf_hyperbili %>% select(SITE, MOMID, PREGID, INFANTID, 
                                           HYPERBILI = hyperbili),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))

##neonatal mortality----
inf_mortality <- INF_COD %>% 
  filter (AGE_DEATH >= 0 & AGE_DEATH <=28) %>% 
  mutate (NEONATAL_DTH = 1)

df_inf_mortality <- df_infant %>%
  select(SITE, MOMID, PREGID, INFANTID, LIVEBIRTH, hb, hb_exm, hb_max, hb_min, trimester, ga_wks) %>%
  left_join(
    dplyr::select(inf_mortality, SITE, MOMID, PREGID, INFANTID, NEONATAL_DTH),
    by = c("SITE", "MOMID", "PREGID", "INFANTID")
  ) %>%
  filter(LIVEBIRTH == 1) %>%
  mutate(neo_mortality = case_when(NEONATAL_DTH == 1 ~ 1, TRUE ~ 0)) %>%
  filter(hb > 0)


df_inf_mortality_trim1 <- df_inf_mortality %>% 
  select(-starts_with("hb"), -trimester, -ga_wks) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_mortality_trim2 <- df_inf_mortality %>% 
  select(-starts_with("hb"), -trimester, -ga_wks) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_mortality_trim3 <- df_inf_mortality %>% 
  select(-starts_with("hb"), -trimester, -ga_wks) %>% 
  left_join(df_hb_exm_trim3, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_analysis <- df_inf_analysis %>% 
  #select (-NEO_MORTALITY) %>%
  left_join(df_inf_mortality %>% select(SITE, MOMID, PREGID, INFANTID, 
                                        NEO_MORTALITY = neo_mortality),
            by = c("SITE", "MOMID","PREGID", "INFANTID"))
#*****************************************************************************
#5. Other Hemoglobin dataset ----
#*****************************************************************************
##long data including hb_type 
df_hb_3type <- df_hb_long %>%
  select(-hb_cbc_poc) %>%
  rename(
    CBC_HB = hb,
    POC_HB = poc_hb,
    SPHB = sphb
  ) %>%
  pivot_longer(c("CBC_HB", "POC_HB", "SPHB"),
               names_to = c("hb_type"),
               values_to = "hb_value") %>%
  #redefine ga_wks by using different ga_wks for different hb type
  mutate(ga_wks = case_when(
    hb_type == "CBC_HB" ~ ga_wks,
    hb_type %in% c("POC_HB", "SPHB") ~ ga_wks_06,
    TRUE ~ NA_real_
  )) %>% 
  select(MOMID, PREGID, SITE, hb_type, hb_value, visit_type = TYPE_VISIT, ga_wks)


##anemia status in each trimeseter ----
df_anemia_trim <- MAT_ANEMIA %>%
  select(SITE, MOMID, PREGID, ANEMIA_T1, ANEMIA_T2, ANEMIA_T3,
         ANEMIA_PNC6, ANEMIA_PNC26)   %>%
  filter(MOMID %in% df_maternal$MOMID) %>% 
  select(-c(ANEMIA_PNC6, ANEMIA_PNC26)) %>% 
  pivot_longer(c(ANEMIA_T1, ANEMIA_T2, ANEMIA_T3), 
               names_to = "trimester",
               values_to = "anemia_status") %>% 
  mutate(
    trimester = paste0("Trimester ", gsub("ANEMIA_T", "", trimester)), 
    anemia = case_when(
      anemia_status %in% c(1:3) ~ 1,
      anemia_status == 0 ~ 0,
      TRUE ~ NA_real_
    )) %>% 
  filter(!is.na(anemia)) 


##anemia for moms who complete pregnancy ----

df_anemia_complet_preg <- MAT_ANEMIA %>%
  select(SITE, MOMID, PREGID, ANEMIA_T1, ANEMIA_T2, ANEMIA_T3,
         ANEMIA_ANC, ANEMIA_PNC6, ANEMIA_PNC26) %>%
  right_join(df_maternal, by = c ("SITE", "MOMID", "PREGID")) %>%
  left_join(MAT_ENDPOINTS) %>% 
  filter(PREG_END == 1) 

test <- df_anemia_complet_preg %>% filter(is.na(ANEMIA_PNC6))

table (df_anemia_complet_preg$SITE)
#*****************************************************************************
#6. Save data ----
#*****************************************************************************
##maternal data(ReMAPP) ----
save(df_maternal, file = "derived_data/df_maternal.rda")
save(all_cbc_hb, file = "derived_data/all_cbc_hb.rda")
save(mat_dth, file = "derived_data/mat_dth.rda")

##hb data----
save(df_hb_long2, file = "derived_data/df_hb_long2.rda")
save(df_hb_wide2, file = "derived_data/df_hb_wide2.rda")
save(df_hb_exm_anc, file = "derived_data/df_hb_exm_anc.rda")
save(df_hb_long, file = "derived_data/df_hb_long.rda")
save(MAT_ANEMIA, file = "derived_data/MAT_ANEMIA.rda")
save(df_anemia_complet_preg, file = "derived_data/df_anemia_complet_preg.rda")
save(df_anemia_trim, file = "derived_data/df_anemia_trim.rda")
save(df_hb_3type, file = "derived_data/df_hb_3type.rda")
save(anc_anemia_long, file = "derived_data/anc_anemia_long.rda")


##matoutcome data----
# save(df_mat_dpr, file = "derived_data/df_mat_dpr.rda")
# save(df_mat_ftg, file = "derived_data/df_mat_ftg.rda")
save(df_mat_pph, file = "derived_data/df_mat_pph.rda")
save(df_mat_ppa_pnc6, file = "derived_data/df_mat_ppa_pnc6.rda")
save(df_mat_ppa_pnc26, file = "derived_data/df_mat_ppa_pnc26.rda")
save(df_mat_pprom, file = "derived_data/df_mat_pprom.rda")
save(df_mat_preclamp, file = "derived_data/df_mat_preclamp.rda")
save(df_mat_compo, file = "derived_data/df_mat_compo.rda")
save(df_mat_analysis, file = "derived_data/df_mat_analysis.rda")


###save fatigue  ----
save(df_mat_anc_fat,  file = "derived_data/df_mat_anc_fat.rda")

###3 Save ANC fatigue data by trimester
save(df_mat_anc_fat_trim1, file = "derived_data/df_mat_anc_fat_trim1.rda")
save(df_mat_anc_fat_trim2, file = "derived_data/df_mat_anc_fat_trim2.rda")
save(df_mat_anc_fat_trim3,file = "derived_data/df_mat_anc_fat_trim3.rda")

#### Save PNC fatigue data (all timepoints)
save(df_mat_pnc_fat_pnc6,file = "derived_data/df_mat_pnc_fat_pnc6.rda")



###save depression ----
save(df_mat_anc_dpr,  file = "derived_data/df_mat_anc_dpr.rda")
# Save PNC depression data by trimester
save(df_mat_anc_dpr_trim1, file = "derived_data/df_mat_anc_dpr_trim1.rda")
save(df_mat_anc_dpr_trim2, file = "derived_data/df_mat_anc_dpr_trim2.rda")
save(df_mat_anc_dpr_trim3,file = "derived_data/df_mat_anc_dpr_trim3.rda")

#### Save PNC depression data (all timepoints)
save(df_mat_pnc_dpr_pnc6,file = "derived_data/df_mat_pnc_dpr_pnc6.rda")

###trim1----
save(df_mat_pph_trim1, file = "derived_data/df_mat_pph_trim1.rda")

save(df_mat_ppa_pnc6_trim1, file = "derived_data/df_mat_ppa_pnc6_trim1.rda")
save(df_mat_ppa_pnc26_trim1, file = "derived_data/df_mat_ppa_pnc26_trim1.rda")

save(df_mat_pprom_trim1, file = "derived_data/df_mat_pprom_trim1.rda")
save(df_mat_preclamp_trim1, file = "derived_data/df_mat_preclamp_trim1.rda")
save(df_mat_compo_trim1, file = "derived_data/df_mat_compo_trim1.rda")

###trim2----
save(df_mat_pph_trim2, file = "derived_data/df_mat_pph_trim2.rda")
save(df_mat_ppa_pnc6_trim2, file = "derived_data/df_mat_ppa_pnc6_trim2.rda")
save(df_mat_ppa_pnc26_trim2, file = "derived_data/df_mat_ppa_pnc26_trim2.rda")
save(df_mat_pprom_trim2, file = "derived_data/df_mat_pprom_trim2.rda")
save(df_mat_preclamp_trim2, file = "derived_data/df_mat_preclamp_trim2.rda")
save(df_mat_compo_trim2, file = "derived_data/df_mat_compo_trim2.rda")

###trim3----
save(df_mat_pph_trim3, file = "derived_data/df_mat_pph_trim3.rda")
save(df_mat_ppa_pnc6_trim3, file = "derived_data/df_mat_ppa_pnc6_trim3.rda")
save(df_mat_ppa_pnc26_trim3, file = "derived_data/df_mat_ppa_pnc26_trim3.rda")
save(df_mat_pprom_trim3, file = "derived_data/df_mat_pprom_trim3.rda")
save(df_mat_preclamp_trim3, file = "derived_data/df_mat_preclamp_trim3.rda")
save(df_mat_compo_trim3, file = "derived_data/df_mat_compo_trim3.rda")
save(df_mat_ppa_pnc26_pnc6, file = "derived_data/df_mat_ppa_pnc26_pnc6.rda")

##infant data----
save(df_infant, file = "derived_data/df_infant.rda")

###infoutcome data----
save(df_inf_compo, file = "derived_data/df_inf_compo.rda")
save(df_inf_preterm37, file = "derived_data/df_inf_preterm37.rda")
save(df_inf_lbw2500, file = "derived_data/df_inf_lbw2500.rda")
save(df_inf_sga10, file = "derived_data/df_inf_sga10.rda")

save(df_inf_preterm34, file = "derived_data/df_inf_preterm34.rda")
save(df_inf_lbw1500, file = "derived_data/df_inf_lbw1500.rda")
save(df_inf_sga3, file = "derived_data/df_inf_sga3.rda")

save(df_inf_psbi, file = "derived_data/df_inf_psbi.rda")

save(df_inf_asph, file = "derived_data/df_inf_asph.rda")

save(df_inf_stillbirth20, file = "derived_data/df_inf_stillbirth20.rda")

save(df_inf_stillbirth28, file = "derived_data/df_inf_stillbirth28.rda")

save(df_inf_hyperbili, file = "derived_data/df_inf_hyperbili.rda")

save(df_inf_mortality, file = "derived_data/df_inf_mortality.rda")

save(df_inf_analysis, file = "derived_data/df_inf_analysis.rda")

###trim1----
save(df_inf_compo_trim1, file = "derived_data/df_inf_compo_trim1.rda")
save(df_inf_preterm37_trim1, file = "derived_data/df_inf_preterm37_trim1.rda")
save(df_inf_lbw2500_trim1, file = "derived_data/df_inf_lbw2500_trim1.rda")
save(df_inf_sga10_trim1, file = "derived_data/df_inf_sga10_trim1.rda")

save(df_inf_preterm34_trim1, file = "derived_data/df_inf_preterm34_trim1.rda")
save(df_inf_lbw1500_trim1, file = "derived_data/df_inf_lbw1500_trim1.rda")
save(df_inf_sga3_trim1, file = "derived_data/df_inf_sga3_trim1.rda")

save(df_inf_psbi_trim1, file = "derived_data/df_inf_psbi_trim1.rda")

save(df_inf_asph_trim1, file = "derived_data/df_inf_asph_trim1.rda")

save(df_inf_stillbirth20_trim1, file = "derived_data/df_inf_stillbirth20_trim1.rda")

save(df_inf_stillbirth28_trim1, file = "derived_data/df_inf_stillbirth28_trim1.rda")

save(df_inf_hyperbili_trim1, file = "derived_data/df_inf_hyperbili_trim1.rda")

save(df_inf_mortality_trim1, file = "derived_data/df_inf_mortality_trim1.rda")

###trim2----
save(df_inf_compo_trim2, file = "derived_data/df_inf_compo_trim2.rda")
save(df_inf_preterm37_trim2, file = "derived_data/df_inf_preterm37_trim2.rda")
save(df_inf_lbw2500_trim2, file = "derived_data/df_inf_lbw2500_trim2.rda")
save(df_inf_sga10_trim2, file = "derived_data/df_inf_sga10_trim2.rda")

save(df_inf_preterm34_trim2, file = "derived_data/df_inf_preterm34_trim2.rda")
save(df_inf_lbw1500_trim2, file = "derived_data/df_inf_lbw1500_trim2.rda")
save(df_inf_sga3_trim2, file = "derived_data/df_inf_sga3_trim2.rda")

save(df_inf_psbi_trim2, file = "derived_data/df_inf_psbi_trim2.rda")

save(df_inf_asph_trim2, file = "derived_data/df_inf_asph_trim2.rda")

save(df_inf_stillbirth20_trim2, file = "derived_data/df_inf_stillbirth20_trim2.rda")

save(df_inf_stillbirth28_trim2, file = "derived_data/df_inf_stillbirth28_trim2.rda")

save(df_inf_hyperbili_trim2, file = "derived_data/df_inf_hyperbili_trim2.rda")

save(df_inf_mortality_trim2, file = "derived_data/df_inf_mortality_trim2.rda")

###trim3----
save(df_inf_compo_trim3, file = "derived_data/df_inf_compo_trim3.rda")
save(df_inf_preterm37_trim3, file = "derived_data/df_inf_preterm37_trim3.rda")
save(df_inf_lbw2500_trim3, file = "derived_data/df_inf_lbw2500_trim3.rda")
save(df_inf_sga10_trim3, file = "derived_data/df_inf_sga10_trim3.rda")

save(df_inf_preterm34_trim3, file = "derived_data/df_inf_preterm34_trim3.rda")
save(df_inf_lbw1500_trim3, file = "derived_data/df_inf_lbw1500_trim3.rda")
save(df_inf_sga3_trim3, file = "derived_data/df_inf_sga3_trim3.rda")
save(df_inf_psbi_trim3, file = "derived_data/df_inf_psbi_trim3.rda")
save(df_inf_asph_trim3, file = "derived_data/df_inf_asph_trim3.rda")
save(df_inf_hyperbili_trim3, file = "derived_data/df_inf_hyperbili_trim3.rda")
save(df_inf_mortality_trim3, file = "derived_data/df_inf_mortality_trim3.rda")
#*****************************************************************************
#7. heatmap data - added 2024-07-29 ----
#*****************************************************************************
##Heat map data for infant outcome
df_heat_inf <- df_inf_lbw2500 %>%
  mutate(hb = round(hb)) %>% 
  group_by(hb) %>% 
  reframe(outcome = "Low birth weight (<2500g)", 
          count = sum(lbw2500 == 1), 
          hb = first(hb)) %>% 
  ungroup() %>% 
  bind_rows(df_inf_lbw1500 %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Low birth weight (<1500g)",
                      count = sum(lbw1500 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_inf_preterm37 %>%
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Preterm (<37 weeks)",
                      count = sum(preterm37 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_inf_preterm34 %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Preterm (<34 weeks)",
                      count = sum(preterm34 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_inf_sga10 %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "SGA (<10th)",
                      count = sum(sga10 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_inf_sga3 %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "SGA (<3th)",
                      count = sum(sga3 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_inf_psbi %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Neonatal PSBI",
                      count = sum(inf_psbi == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_inf_asph %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Birth asphyxia",
                      count = sum(inf_asph == 1), 
                      hb = first(hb)) %>% 
              ungroup()) 

outcome_order <- c("Preterm (<37 weeks)", "Low birth weight (<2500g)", "SGA (<10th)",
                   "Preterm (<34 weeks)","Low birth weight (<1500g)", "SGA (<3th)",
                   "Neonatal PSBI", "Birth asphyxia")

df_heat_inf$outcome <- factor(df_heat_inf$outcome,levels = rev(outcome_order))

save(df_heat_inf, file = "derived_data/df_heat_inf.rda")

##Heat map data for maternal outcome ----
df_heat_mat <- df_mat_pph %>% 
  mutate(hb = round(hb)) %>% 
  group_by(hb) %>% 
  reframe(outcome = "Postpartum Hemorrhage", 
          count = sum(HEM_PPH == 1), 
          hb = first(hb),
  ) %>% 
  ungroup() %>% 
  bind_rows(df_mat_ppa_pnc6 %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Postpartum anemia at PNC6", 
                      count = sum(ppa_pnc6 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_mat_ppa_pnc26 %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Postpartum anemia at PNC26",
                      count = sum(ppa_pnc26 == 1), 
                      hb = first(hb)) %>% 
              ungroup()) %>% 
  bind_rows(df_mat_anc_dpr %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Likelihood of depression", 
                      count = sum(depress == 1), 
                      hb = first(hb)) %>% 
              ungroup())

outcome_order <- c("Postpartum Hemorrhage", "Postpartum anemia at PNC6", "Postpartum anemia at PNC26",
                   "Likelihood of depression")

df_heat_mat$outcome <- factor(df_heat_mat$outcome,levels = rev(outcome_order))

save(df_heat_mat, file = "derived_data/df_heat_mat.rda")