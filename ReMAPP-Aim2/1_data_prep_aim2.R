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

UploadDate = "2026-01-30"

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
derived_rr_dir <- file.path(new_dir, "derived_data_rr")
iso_data_dir <- file.path(new_dir, "iso_results")

# Create Derived_Data directory if it doesn't exist  
if (!dir.exists(derived_data_dir)) {
  dir.create(derived_data_dir)
}

if (!dir.exists(derived_rr_dir)) {
  dir.create(derived_rr_dir)
}

if (!dir.exists(iso_data_dir)) {
  dir.create(iso_data_dir)
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

# Robust date cleaner
clean_date <- function(x) {
  # If already Date, return as-is (keeps NA_Date_)
  if (inherits(x, "Date")) return(x)
  
  # If numeric (Excel serials), convert (as.Date handles NA)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  
  # Coerce to character and clean
  x_chr <- as.character(x)
  x_chr <- stringr::str_trim(x_chr)
  x_chr[x_chr == ""] <- NA_character_
  
  # Parse common date formats (POSIXct)
  parsed <- suppressWarnings(lubridate::parse_date_time(
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
  need_num <- is.na(parsed) &
    !is.na(x_chr) &
    stringr::str_detect(x_chr, "^[0-9]{5}$")
  
  if (any(need_num, na.rm = TRUE)) {
    parsed[need_num] <- as.POSIXct(
      as.Date(as.numeric(x_chr[need_num]), origin = "1899-12-30"),
      tz = "UTC"
    )
  }
  
  as.Date(parsed)
}

#Function to cap HB
cap_hb <- function(x, lo = 5, hi = 18) {
  x <- as.numeric(x)
  if_else(is.na(x), NA_real_,
          if_else(x < lo, lo,
                  if_else(x > hi, hi, x)))
}

#Function to create a clean dataset for RR analysis
outcome_df_rr <- function(df,
                          hb_var,
                          outcome_var,
                          trimester_var,
                          ga_var,
                          id_vars,
                          out_rda_dir,
                          out_xlsx_dir,
                          file_name) {
  
  cleaned <- df %>%
    transmute(
      across(all_of(id_vars)),
      !!ga_var := .data[[ga_var]],
      !!trimester_var := .data[[trimester_var]],
      !!hb_var := cap_hb(.data[[hb_var]]),
      !!outcome_var := .data[[outcome_var]]
    ) %>%
    
    #Fix trimester if missing using GA_WKS
    mutate(
      !!trimester_var := case_when(
        !is.na(.data[[trimester_var]]) &
          .data[[trimester_var]] %in% c(1,2,3) ~ .data[[trimester_var]],
        
        is.na(.data[[trimester_var]]) &
          !is.na(.data[[ga_var]]) &
          .data[[ga_var]] < 14 ~ 1,
        
        is.na(.data[[trimester_var]]) &
          !is.na(.data[[ga_var]]) &
          .data[[ga_var]] >= 14 &
          .data[[ga_var]] < 28 ~ 2,
        
        is.na(.data[[trimester_var]]) &
          !is.na(.data[[ga_var]]) &
          .data[[ga_var]] >= 28 ~ 3,
        
        TRUE ~ NA_real_
      )
    ) %>%
    mutate(
      ANEMIA_4 = case_when(
        .data[[trimester_var]] %in% c(1,3) & .data[[hb_var]] > 0 & .data[[hb_var]] < 7    ~ 3,
        .data[[trimester_var]] %in% c(1,3) & .data[[hb_var]] >= 7 & .data[[hb_var]] < 10  ~ 2,
        .data[[trimester_var]] %in% c(1,3) & .data[[hb_var]] >= 10 & .data[[hb_var]] < 11 ~ 1,
        .data[[trimester_var]] %in% c(1,3) & .data[[hb_var]] >= 11                        ~ 0,
        .data[[trimester_var]] == 2 & .data[[hb_var]] > 0 & .data[[hb_var]] < 7           ~ 3,
        .data[[trimester_var]] == 2 & .data[[hb_var]] >= 7 & .data[[hb_var]] < 9.5        ~ 2,
        .data[[trimester_var]] == 2 & .data[[hb_var]] >= 9.5 & .data[[hb_var]] < 10.5     ~ 1,
        .data[[trimester_var]] == 2 & .data[[hb_var]] >= 10.5                             ~ 0,
        TRUE ~ NA_real_
      ),
      ANEMIA_MODSEV = if_else(ANEMIA_4 %in% c(2,3), 1,
                              if_else(ANEMIA_4 %in% c(0,1), 0, NA_real_)),
      ANEMIA_ANY = if_else(ANEMIA_4 == 0, 0,
                           if_else(ANEMIA_4 %in% c(1,2,3), 1, NA_real_))
    ) %>%
    rename_with(toupper)
  
  # Save as dataset name
  assign(file_name, cleaned)
  save(list = file_name,
       file = file.path(out_rda_dir, paste0(file_name, ".rda")))
  rm(list = file_name)
  
  # Save Excel
  write.xlsx(cleaned,
             file = file.path(out_xlsx_dir, paste0(file_name, ".xlsx")),
             overwrite = TRUE)
  
  return(cleaned)
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
) %>% select(SITE, MOMID, PREGID, PREG_END_DATE, PREG_END, PREG_END_SOURCE, 
             PREG_END_GA, MAT_DEATH_INFAGE, CLOSEOUT_TYPE, MAT_DEATH )

## load MAT_FLOWCHART  ----
MAT_FLOWCHART <- load_or_build(
  "MAT_FLOWCHART",
  "derived_data/MAT_FLOWCHART.rda",   # keep a raw cache, then select below
  function() read.csv(paste0("Z:/Outcome Data/", UploadDate, "/MAT_FLOWCHART.csv")))

## load MNH03 for smoking adjustment of hb  ----
mnh03 <- load_or_build(
  "mnh03",
  "derived_data/mnh03.rda",
  function() read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh03_merged.csv")) %>%
    select(SITE, MOMID, PREGID, M03_SMOKE_OECOCCUR)
)

## load MNH04 for anc visit completion  ----
mnh04 <- load_or_build(
  "mnh04",
  "derived_data/mnh04.rda",
  function() read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh04_merged.csv")) %>%
    select(SITE, MOMID, PREGID, M04_ANC_OBSSTDAT, M04_TYPE_VISIT, 
           M04_MAT_VITAL_MNH04, M04_MAT_VISIT_MNH04, M04_PRG_DSDECOD )
)


## load MNH06 and process ----
# raw_mnh06 for smoking adjustment of hb
## MNH06: compute GA/age + expected type visit 
## then DEDUPLICATE BY DATE within each (SITE, MOMID, PREGID, VISIT DATE)
##Make a wide dataset with a list of the dates, reported type visits and hb reported

raw_mnh06 <- load_or_build(
  "raw_mnh06",
  "derived_data/raw_mnh06.rda",
  function()  read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh06_merged.csv"))
)

# join enroll + endpoints, compute ga_days / age_days
raw_mnh06_2 <- raw_mnh06 %>%
  select(
    SITE, MOMID, PREGID,
    M06_TYPE_VISIT, M06_HB_POC_LBORRES, M06_SPHB_LBORRES, M06_DIAG_VSDAT
  ) %>%
  left_join(
    MAT_ENROLL %>% select(SITE, MOMID, PREGID, BOE_GA_DAYS_ENROLL, PREG_START_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  left_join(
    MAT_FLOWCHART %>% select(SITE, MOMID, PREGID, REMAPP_SCRN, ENROLL_NO_ISSUES),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  left_join(
    MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  mutate(
    M06_DIAG_VSDAT  = clean_date(M06_DIAG_VSDAT),
    PREG_START_DATE = clean_date(PREG_START_DATE),
    PREG_END_DATE   = clean_date(PREG_END_DATE),
    
    M06_HB_POC_LBORRES = as.numeric(M06_HB_POC_LBORRES),
    M06_SPHB_LBORRES   = as.numeric(M06_SPHB_LBORRES),
    M06_HB_POC_LBORRES = if_else(M06_HB_POC_LBORRES < 0, NA_real_, M06_HB_POC_LBORRES),
    M06_SPHB_LBORRES   = if_else(M06_SPHB_LBORRES < 0, NA_real_, M06_SPHB_LBORRES),
    
    date_diff = as.numeric(
      M06_DIAG_VSDAT -
        if_else(M06_TYPE_VISIT %in% c(1:5, 13), PREG_START_DATE, PREG_END_DATE)
    ),
    ga_days  = if_else(M06_TYPE_VISIT %in% c(1:5, 13), date_diff, NA_real_),
    age_days = if_else(M06_TYPE_VISIT %in% c(6:12, 14), date_diff, NA_real_)
  )

prep_mnh06 <- raw_mnh06_2 %>%
  filter(
    REMAPP_SCRN == 1 & ENROLL_NO_ISSUES == 1 & 
      (is.na(M06_DIAG_VSDAT) | M06_DIAG_VSDAT != ymd("1907-07-07")) &
      (
        # Pregnancy visits: within pregnancy window
        (
          M06_TYPE_VISIT %in% c(1:5, 13) &
            !is.na(PREG_START_DATE) &
            M06_DIAG_VSDAT >= PREG_START_DATE &
            (is.na(PREG_END_DATE) | M06_DIAG_VSDAT <= PREG_END_DATE)
        ) |
          # Postnatal visits: after PREG_END_DATE
          (
            M06_TYPE_VISIT %in% c(6:12, 14) &
              !is.na(PREG_END_DATE) &
              M06_DIAG_VSDAT >= PREG_END_DATE
          )
      )
  ) %>%
  mutate(
    EXPECTED_TYPE_VISIT = case_when(
      # Gestational visit logic
      ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL < ga_days ~ 2,
      ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
      ga_days >= 182 & ga_days <= 216 ~ 3,
      ga_days >= 217 & ga_days <= 237 ~ 4,
      ga_days >= 238 & ga_days <= 300 ~ 5,
      
      # Postnatal visit types
      M06_TYPE_VISIT == 6 | age_days < 3 ~ 6,
      age_days >= 3 & age_days <= 5 ~ 7,
      age_days >= 7 & age_days <= 14 ~ 8,
      age_days >= 28 & age_days <= 35 ~ 9,
      age_days >= 42 & age_days <= 104 ~ 10,
      age_days >= 104 & age_days <= 279 ~ 11,
      age_days >= 279 & age_days <= 454 ~ 12,
      
      # .5 visits between protocol windows
      age_days > 5 & age_days < 7 ~ 7.5,
      age_days > 14 & age_days < 28 ~ 8.5,
      age_days > 35 & age_days < 42 ~ 9.5,
      
      TRUE ~ 55
    ),
    EXPECTED_TYPE_VISIT = if_else(
      EXPECTED_TYPE_VISIT %% 1 == 0,
      as.integer(EXPECTED_TYPE_VISIT),
      EXPECTED_TYPE_VISIT
    )
  ) %>%
  filter(ga_days < 300 | age_days < 455)


## store all non-missing hb data
all_poc_hb <- prep_mnh06 %>%
  filter(!(is.na(M06_HB_POC_LBORRES) & is.na(M06_SPHB_LBORRES))) %>%
  select(
    SITE, MOMID, PREGID,
    EXPECTED_TYPE_VISIT,
    M06_DIAG_VSDAT,
    M06_HB_POC_LBORRES, M06_SPHB_LBORRES,
    ga_days, age_days) %>%
  rename(
    M06_TYPE_VISIT = EXPECTED_TYPE_VISIT,
    ga_days_06 = ga_days,
    age_days_06 = age_days)


### function: resolve_mnh06_hb ----
# Purpose: For each MOMID-PREGID-VISITDATE group in the mnh06 dataset,
#          this function selects a single Hb value (M06_HB_POC_LBORRES/SPHB) by:
#          - Calculating the median Hb within the group
#          - Measuring the distance of each Hb from that median
#          - Selecting the value farthest from the median (most "extreme")
#          - Randomly choosing one if there are ties
# Notes:
#   - Groups with all missing Hb values are excluded from processing
#   - This helps resolve duplicates by keeping one meaningful, extreme Hb value

resolve_mnh06_hb <- function(df, hb_var, keep_vars = 
                               c("SITE", "MOMID", "PREGID", "M06_TYPE_VISIT", 
                                 "M06_DIAG_VSDAT","ga_days_06","age_days_06")) {
  hb_var <- rlang::ensym(hb_var)
  
  set.seed(100)  # For reproducibility
  df %>%
    group_by(MOMID, PREGID, M06_DIAG_VSDAT) %>%
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
resolved_poc <- all_poc_hb %>%
  filter(!is.na(M06_HB_POC_LBORRES)) %>%
  resolve_mnh06_hb("M06_HB_POC_LBORRES")

# For SpHb
resolved_sphb <- all_poc_hb %>%
  filter(!is.na(M06_SPHB_LBORRES)) %>%
  resolve_mnh06_hb("M06_SPHB_LBORRES")


mnh06_long <- full_join(resolved_poc, resolved_sphb,
                      by = c("SITE", "MOMID", "PREGID", "M06_TYPE_VISIT", 
                             "M06_DIAG_VSDAT", "ga_days_06","age_days_06")) %>%
  mutate(
    VISIT_PHASE_06 = case_when(
      M06_TYPE_VISIT %in% c(1:5) ~ "ANC",
      M06_TYPE_VISIT %in% c(6:12) ~ "PNC",
      TRUE ~ "OTHER")) %>%
  arrange(SITE, MOMID, PREGID, VISIT_PHASE_06, M06_DIAG_VSDAT)


mnh06_wide <- mnh06_long %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    #POC
    ANC_POC_HB_LIST   = list(M06_HB_POC_LBORRES[VISIT_PHASE_06 == "ANC" & !is.na(M06_HB_POC_LBORRES)]),
    ANC_POC_DATE_LIST = list(M06_DIAG_VSDAT    [VISIT_PHASE_06 == "ANC" & !is.na(M06_HB_POC_LBORRES)]),
    N_ANC_POC_HB      = sum(VISIT_PHASE_06 == "ANC" & !is.na(M06_HB_POC_LBORRES)),
    
    PNC_POC_HB_LIST   = list(M06_HB_POC_LBORRES[VISIT_PHASE_06 == "PNC" & !is.na(M06_HB_POC_LBORRES)]),
    PNC_POC_DATE_LIST = list(M06_DIAG_VSDAT    [VISIT_PHASE_06 == "PNC" & !is.na(M06_HB_POC_LBORRES)]),
    N_PNC_POC_HB      = sum(VISIT_PHASE_06 == "PNC" & !is.na(M06_HB_POC_LBORRES)),
    
    #SPHB
    ANC_SPHB_HB_LIST   = list(M06_SPHB_LBORRES[VISIT_PHASE_06 == "ANC" & !is.na(M06_SPHB_LBORRES)]),
    ANC_SPHB_DATE_LIST = list(M06_DIAG_VSDAT  [VISIT_PHASE_06 == "ANC" & !is.na(M06_SPHB_LBORRES)]),
    N_ANC_SPHB_HB      = sum(VISIT_PHASE_06 == "ANC" & !is.na(M06_SPHB_LBORRES)),
    
    PNC_SPHB_HB_LIST   = list(M06_SPHB_LBORRES[VISIT_PHASE_06 == "PNC" & !is.na(M06_SPHB_LBORRES)]),
    PNC_SPHB_DATE_LIST = list(M06_DIAG_VSDAT  [VISIT_PHASE_06 == "PNC" & !is.na(M06_SPHB_LBORRES)]),
    N_PNC_SPHB_HB      = sum(VISIT_PHASE_06 == "PNC" & !is.na(M06_SPHB_LBORRES)),
    
    .groups = "drop"
  ) 

## load MNH08 and process ----

#mnh08 for hemoglobin 
raw_mnh08 <- load_or_build(
  "raw_mnh08",
  "derived_data/raw_mnh08.rda",
  function()  read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh08_merged.csv")))


raw_mnh08 <- load_or_build(
  "raw_mnh08",
  "derived_data/raw_mnh08.rda",
  function()  read.csv(paste0("Z:/Stacked Data/", UploadDate, "/mnh08_merged.csv"))
)

# join enroll + endpoints, compute ga_days / age_days
raw_mnh08_2 <- raw_mnh08 %>%
  select(
    SITE, MOMID, PREGID, M08_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES
  ) %>%
  left_join(
    MAT_ENROLL %>% select(SITE, MOMID, PREGID, BOE_GA_DAYS_ENROLL, PREG_START_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  left_join(
    MAT_FLOWCHART %>% select(SITE, MOMID, PREGID, REMAPP_SCRN, ENROLL_NO_ISSUES),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  left_join(
    MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  mutate(
    M08_LBSTDAT  = clean_date(M08_LBSTDAT),
    PREG_START_DATE = clean_date(PREG_START_DATE),
    PREG_END_DATE   = clean_date(PREG_END_DATE),
    
    M08_CBC_HB_LBORRES = as.numeric(M08_CBC_HB_LBORRES),
    M08_CBC_HB_LBORRES = if_else(M08_CBC_HB_LBORRES < 0, NA_real_, M08_CBC_HB_LBORRES),
    date_diff = as.numeric(
      M08_LBSTDAT -
        if_else(M08_TYPE_VISIT %in% c(1:5, 13), PREG_START_DATE, PREG_END_DATE)
    ),
    ga_days  = if_else(M08_TYPE_VISIT %in% c(1:5, 13), date_diff, NA_real_),
    age_days = if_else(M08_TYPE_VISIT %in% c(6:12, 14), date_diff, NA_real_)
  )

prep_mnh08 <- raw_mnh08_2 %>%
  filter(
    REMAPP_SCRN == 1 & ENROLL_NO_ISSUES == 1 & 
      (is.na(M08_LBSTDAT) | M08_LBSTDAT != ymd("1907-07-07")) &
      (
        # Pregnancy visits: within pregnancy window
        (
          M08_TYPE_VISIT %in% c(1:5, 13) &
            !is.na(PREG_START_DATE) &
            M08_LBSTDAT >= PREG_START_DATE &
            (is.na(PREG_END_DATE) | M08_LBSTDAT <= PREG_END_DATE)
        ) |
          # Postnatal visits: after PREG_END_DATE
          (
            M08_TYPE_VISIT %in% c(6:12, 14) &
              !is.na(PREG_END_DATE) &
              M08_LBSTDAT >= PREG_END_DATE
          ))) %>%
  mutate(
    EXPECTED_TYPE_VISIT = case_when(
      # Gestational visit logic
      ga_days >= 126 & ga_days <= 181 & BOE_GA_DAYS_ENROLL < ga_days ~ 2,
      ga_days <= 139 | ga_days == BOE_GA_DAYS_ENROLL ~ 1,
      ga_days >= 182 & ga_days <= 216 ~ 3,
      ga_days >= 217 & ga_days <= 237 ~ 4,
      ga_days >= 238 & ga_days <= 300 ~ 5,
      
      # Postnatal visit types
      M08_TYPE_VISIT == 6 | age_days < 3 ~ 6,
      age_days >= 3 & age_days <= 5 ~ 7,
      age_days >= 7 & age_days <= 14 ~ 8,
      age_days >= 28 & age_days <= 35 ~ 9,
      age_days >= 42 & age_days <= 104 ~ 10,
      age_days >= 104 & age_days <= 279 ~ 11,
      age_days >= 279 & age_days <= 454 ~ 12,
      
      # .5 visits between protocol windows
      age_days > 5 & age_days < 7 ~ 7.5,
      age_days > 14 & age_days < 28 ~ 8.5,
      age_days > 35 & age_days < 42 ~ 9.5,
      
      TRUE ~ 55
    ),
    EXPECTED_TYPE_VISIT = if_else(
      EXPECTED_TYPE_VISIT %% 1 == 0,
      as.integer(EXPECTED_TYPE_VISIT),
      EXPECTED_TYPE_VISIT
    )
  ) %>%
  filter(ga_days < 300 | age_days < 455)


## store all non-missing hb data
all_cbc_hb <- prep_mnh08 %>%
  filter(!(is.na(M08_CBC_HB_LBORRES))) %>%
  select(
    SITE, MOMID, PREGID,
    EXPECTED_TYPE_VISIT,
    M08_LBSTDAT,
    M08_CBC_HB_LBORRES,
    ga_days, age_days) %>%
  rename(
    M08_TYPE_VISIT = EXPECTED_TYPE_VISIT,
    ga_days_08 = ga_days,
    age_days_08 = age_days)

### function: resolve_mnh08_hb ----
# Purpose: For each MOMID-PREGID-VISIT Date group in the mnh08 dataset,
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
    group_by(MOMID, PREGID, M08_LBSTDAT) %>%
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

#applying the function to remove the duplicates to have a long dataset for each age
mnh08_long <- resolve_mnh08_hb(all_cbc_hb)

### mnh08 wide dataset creation ----
mnh08_long <- mnh08_long %>%
  mutate(
    VISIT_PHASE_08 = case_when(
      M08_TYPE_VISIT %in% c(1:5) ~ "ANC",
      M08_TYPE_VISIT %in% c(6:12) ~ "PNC",
      TRUE ~ "OTHER")) %>%
  arrange(SITE, MOMID, PREGID, VISIT_PHASE_08, M08_LBSTDAT)


mnh08_wide <- mnh08_long %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    ANC_CBC_HB_LIST   = list(M08_CBC_HB_LBORRES[VISIT_PHASE_08 == "ANC" & !is.na(M08_CBC_HB_LBORRES)]),
    ANC_CBC_DATE_LIST = list(M08_LBSTDAT    [VISIT_PHASE_08 == "ANC" & !is.na(M08_CBC_HB_LBORRES)]),
    N_ANC_CBC_HB      = sum(VISIT_PHASE_08 == "ANC" & !is.na(M08_CBC_HB_LBORRES)),
    
    PNC_CBC_HB_LIST   = list(M08_CBC_HB_LBORRES[VISIT_PHASE_08 == "PNC" & !is.na(M08_CBC_HB_LBORRES)]),
    PNC_CBC_DATE_LIST = list(M08_LBSTDAT    [VISIT_PHASE_08 == "PNC" & !is.na(M08_CBC_HB_LBORRES)]),
    N_PNC_CBC_HB      = sum(VISIT_PHASE_08 == "PNC" & !is.na(M08_CBC_HB_LBORRES)),
    .groups = "drop"
  ) 


## MNH25 processing ----
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
  filter(!is.na(ga_wks_25) & !is.na(epds_score)) %>%
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

## MNH26 processing Prep----
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
  filter(!is.na(ga_wks_26) & !is.na(fatigue_score)) %>%
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
    select(SITE, MOMID, PREGID, HEM_PPH, HEM_PPH_DATE,
           HEM_PPH_AGE_PP_DAYS, HEM_PPH_FORM_SOURCE, HEM_DENOM))


## load MAT_PRETERM (PPROM_OCCUR) ----
# MAT_PRETERM <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_PRETERM.dta")) 
MAT_PRETERM <- load_or_build(
  "MAT_PRETERM",
  "derived_data/MAT_PRETERM.rda",
  function() MAT_PRETERM <- 
    read_dta(paste0("Z:/Outcome Data/", UploadDate,"/MAT_PRETERM.dta")))


## combined df_maternal data set ----
df_maternal <- MAT_FLOWCHART %>% 
  filter (REMAPP_SCRN ==1 & ENROLL_NO_ISSUES == 1) %>%
  select (SITE, MOMID, PREGID, REMAPP_SCRN, ENROLL_NO_ISSUES, PREG_END_ENROLLED) %>%
  left_join(MAT_ENROLL, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(MAT_ENDPOINTS, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh03, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh08_wide, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh06_wide, by = c("SITE", "MOMID", "PREGID"))

## creating a masterlist for all inclusion criteria ----
df_mat_mlist <-  df_maternal   %>% 
  mutate (MISS_HB_ANC = case_when(N_ANC_CBC_HB < 1 | is.na(N_ANC_CBC_HB) ~ 1,
                                  N_ANC_CBC_HB >= 1 ~ 0,
                                  TRUE ~ 55),
          MISS_HB_PNC = case_when(N_PNC_CBC_HB < 1 | is.na(N_PNC_CBC_HB) ~ 1,
                                  N_PNC_CBC_HB >= 1 ~ 0,
                                  TRUE ~ 55),
          MISS_HB_ALL = case_when(MISS_HB_ANC == 1 & MISS_HB_PNC == 1 ~ 1,
                                  MISS_HB_ANC == 0 | MISS_HB_PNC == 0 ~ 0,
                                  TRUE ~ 55),
          HB_COUNT_ANC = N_ANC_CBC_HB,
          HB_COUNT_PNC = N_PNC_CBC_HB) %>% 
  select (MOMID, PREGID, SITE, REMAPP_ENROLL, PREG_END, MISS_HB_ANC, 
          MISS_HB_PNC, MISS_HB_ALL, HB_COUNT_ANC, HB_COUNT_PNC) %>% 
  distinct(SITE, MOMID, PREGID, .keep_all = TRUE)

##maternal analysis masterlist ----
df_mat_analysis <- df_mat_mlist 
  
##infant remapp masterlist ----
df_inf_mlist <- INF_OUTCOMES %>% 
  left_join(df_maternal, by = c("SITE", "MOMID", "PREGID")) %>% 
  filter (REMAPP_SCRN ==1 & ENROLL_NO_ISSUES == 1) %>%
  mutate (MISS_HB_ANC = case_when(N_ANC_CBC_HB < 1 | is.na(N_ANC_CBC_HB) ~ 1,
                                  N_ANC_CBC_HB >= 1 ~ 0,
                                  TRUE ~ 55),
          MISS_HB_PNC = case_when(N_PNC_CBC_HB < 1 | is.na(N_PNC_CBC_HB) ~ 1,
                                  N_PNC_CBC_HB >= 1 ~ 0,
                                  TRUE ~ 55),
          MISS_HB_ALL = case_when(MISS_HB_ANC == 1 & MISS_HB_PNC == 1 ~ 1,
                                  MISS_HB_ANC == 0 | MISS_HB_PNC == 0 ~ 0,
                                  TRUE ~ 55),
          HB_COUNT_ANC = N_ANC_CBC_HB,
          HB_COUNT_PNC = N_PNC_CBC_HB) %>% 
  select (SITE, MOMID, PREGID, INFANTID, REMAPP_ENROLL = REMAPP_SCRN, LIVEBIRTH,
          BIRTH_OUTCOME_REPORTED, STILLBIRTH_DENOM, STILLBIRTH_20WK,
          MISS_HB_ANC, MISS_HB_PNC, MISS_HB_ALL, HB_COUNT_ANC, 
          HB_COUNT_PNC) 


##infant analysis masterlist ----
df_inf_analysis <- df_inf_mlist 


#*****************************************************************************
#1.5 A. Maternal Denominators ----
#*****************************************************************************
#After enrolment, we want to determine the number of women, who had at least one visit
mnh04_sub <- mnh04 %>%
  mutate(
    M04_ANC_OBSSTDAT = clean_date(M04_ANC_OBSSTDAT)
  ) %>%
  filter (M04_TYPE_VISIT %in% c(1:5, 13)) %>%
  # remove not-complete visits 
  filter(M04_MAT_VISIT_MNH04 %in% c(1, 2)) %>%
  # keep mom alive
  filter(M04_MAT_VITAL_MNH04 == 1)  %>%
  #remove duplicate visit dates
  arrange(SITE, MOMID, PREGID, M04_ANC_OBSSTDAT, desc(M04_TYPE_VISIT)) %>%
  distinct(SITE, MOMID, PREGID, M04_ANC_OBSSTDAT, .keep_all = TRUE) 

last_viable_df <- mnh04_sub %>%
  group_by(SITE, MOMID, PREGID) %>%
  mutate(
    enroll_date = (min(M04_ANC_OBSSTDAT[M04_TYPE_VISIT == 1], na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # keep only visits AFTER enrollment date
  filter(!is.na(enroll_date), M04_ANC_OBSSTDAT > enroll_date) %>%
  # viable / continued viable pregnancy
  #filter(M04_PRG_DSDECOD == 1) %>%
  group_by(SITE, MOMID, PREGID, enroll_date) %>%
  arrange(M04_ANC_OBSSTDAT) %>%
  mutate(N_VISITS_AFTER_ENROLL = n()) %>%  # <- count of viable visits after enrollment
  slice_tail(n = 1) %>%   # last viable visit after enrollment
  ungroup() %>%
  transmute(
    SITE, MOMID, PREGID,
    LAST_VIABLE_VISITDATE = M04_ANC_OBSSTDAT,
    LAST_VIABLE_TYPEVISIT = M04_TYPE_VISIT,
    N_VISITS_AFTER_ENROLL,
    ANALYSIS_DENOM_A = 1
  )

#group infant outcome by mom, to see how many moms have 
mat_inf_outcome <- INF_OUTCOMES %>%
  select(
    SITE, MOMID, PREGID, INFANTID,
    LIVEBIRTH,
    BIRTH_OUTCOME_REPORTED,
    STILLBIRTH_20WK,
    STILLBIRTH_28WK,
    STILLBIRTH_DENOM,
    MISSING_MNH09,
    MISSING_MNH11,
    DTH_TIME_MISSING,
    DOB_AFTER_DEATH,
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    N_INFANTS = n_distinct(INFANTID),
    
    N_LIVEBIRTH = sum(LIVEBIRTH == 1, na.rm = TRUE),
    
    N_STILLBIRTH_20WK = sum(
      STILLBIRTH_DENOM == 1 & STILLBIRTH_20WK == 1,
      na.rm = TRUE
    ),
    
    N_STILLBIRTH_28WK = sum(
      STILLBIRTH_DENOM == 1 & STILLBIRTH_28WK == 1,
      na.rm = TRUE
    ), 
    N_BIRTH_OUTCOME_REPORTED = sum(
      BIRTH_OUTCOME_REPORTED == 1,
      na.rm = TRUE
    ),
    # data quality / missingness flags
    N_MISSING_MNH09 = sum(MISSING_MNH09 == 1, na.rm = TRUE),
    N_MISSING_MNH11 = sum(MISSING_MNH11 == 1, na.rm = TRUE),
    N_DTH_TIME_MISSING = sum(DTH_TIME_MISSING == 1, na.rm = TRUE),
    N_DOB_AFTER_DEATH = sum(DOB_AFTER_DEATH == 1, na.rm = TRUE),
    
    .groups = "drop"
  )

mat_remapp_denom <- last_viable_df %>%
  right_join(
    MAT_FLOWCHART %>%
      filter(REMAPP_SCRN == 1, ENROLL_NO_ISSUES == 1) %>%
      select(SITE, MOMID, PREGID, CLOSEOUT_GA, PREG_END_ENROLLED, starts_with("ELIGIBLE_")),
      by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(MAT_ENDPOINTS, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(MAT_ENROLL %>% select(SITE, MOMID, PREGID, PREG_START_DATE),
      by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mat_inf_outcome, by = c("SITE", "MOMID", "PREGID")) %>%
##Denominator A - Enrolled + One Visit ----
  mutate (ONE_VISIT_DENOM = 
            case_when(N_VISITS_AFTER_ENROLL >= 1 ~ 1, #included
                      ELIGIBLE_EXCL_REASON == 1  ~ 2, #mat death
                      ELIGIBLE_EXCL_REASON %in% c(2,4,5) ~ 3, #loss to followup
                      ELIGIBLE_EXCL_REASON == 3 ~ 4, #withdrew from study
                      is.na(N_VISITS_AFTER_ENROLL) ~ 5, #no visit after enrolment/recode to ltfu after dataclose - PJW Todo
                      TRUE ~ 55),
          DENOM_A_INCLUDED   = if_else(ONE_VISIT_DENOM == 1, 1, 0),
          DENOM_A_EXCL_CODE  = if_else(ONE_VISIT_DENOM == 1, NA_real_, as.numeric(ONE_VISIT_DENOM)),
          DENOM_A_EXCL_LABEL = case_when(
            ONE_VISIT_DENOM == 1 ~ "Included: Enrolled + ≥1 visit after enrollment",
            ONE_VISIT_DENOM == 2 ~ "Maternal death",
            ONE_VISIT_DENOM == 3 ~ "Loss to follow-up",
            ONE_VISIT_DENOM == 4 ~ "Withdrew from study",
            ONE_VISIT_DENOM == 5 ~ "Loss to follow-up",
            TRUE                 ~ "Other unspecified reasons"
          ),
          
          GA_AT_LAST_VISIT = as.numeric(difftime(
              ymd(LAST_VIABLE_VISITDATE),
              ymd(PREG_START_DATE),
              units = "days")),
          FOLLOWUP_END_DATE = ymd(UploadDate), #change to when consortium decides is the site closeout
          AGE_INF_TODAY = as.numeric(difftime(
            ymd(FOLLOWUP_END_DATE),
            ymd(PREG_END_DATE),
            units = "days")),
          
 ##Denominator B - Enrolled + One Visit + G A20 Met ----          
          ONE_VISIT_DENOM_20 =  
            case_when(N_VISITS_AFTER_ENROLL >= 1 & 
                      (GA_AT_LAST_VISIT >= 140 |PREG_END_GA >= 140) ~ 1, #included 
                      ONE_VISIT_DENOM %in% c(2:5) ~ 77, #We already know they got kicked out
                      N_VISITS_AFTER_ENROLL >= 1 & PREG_END_GA < 140 ~ 6, #not included if preg end date <20weeks GA
                      N_VISITS_AFTER_ENROLL >= 1 & GA_AT_LAST_VISIT < 140 ~ 7, #not included if last visit seen <20weeks GA
                      TRUE ~ 55),
         DENOM_B_INCLUDED   = if_else(ONE_VISIT_DENOM_20 == 1, 1, 0),
         DENOM_B_EXCL_CODE  = if_else(ONE_VISIT_DENOM_20 == 1, NA_real_, as.numeric(ONE_VISIT_DENOM_20)),
         DENOM_B_EXCL_LABEL = case_when(
           ONE_VISIT_DENOM_20 == 1  ~ "Included: Enrolled + ≥1 visit + reached ≥20 weeks gestation (GA)",
           ONE_VISIT_DENOM_20 == 6  ~ "Pregnancy ended <20 weeks GA",
           ONE_VISIT_DENOM_20 == 7  ~ "Lost to follow-up",
           ONE_VISIT_DENOM_20 == 77 ~ "Not Applicable",
           TRUE                     ~ "Other unspecified reasons"
         ),
 ##Denominator C - Any Pregnancy End (including death) + After GA20 Met ----
          PREG_END_DENOM_20 = 
            case_when(N_VISITS_AFTER_ENROLL >= 1 & 
                      PREG_END == 1 & PREG_END_GA >= 140 ~ 1, #included
                      ONE_VISIT_DENOM_20 %in% c(2:7,77) ~ 77, #We already know they got kicked out
                      ELIGIBLE_EXCL_REASON %in% c(2,4,5) ~ 3, #loss to followup
                      ELIGIBLE_EXCL_REASON == 3 ~ 4, #withdrew from study
                      TRUE ~ 55),
 
         DENOM_C_INCLUDED   = if_else(PREG_END_DENOM_20 == 1, 1, 0),
         DENOM_C_EXCL_CODE  = if_else(PREG_END_DENOM_20 == 1, NA_real_, as.numeric(PREG_END_DENOM_20)),
         DENOM_C_EXCL_LABEL = case_when(
           PREG_END_DENOM_20 == 1  ~ "Included: Pregnancy endpoint recorded AND reached ≥20 weeks GA",
           PREG_END_DENOM_20 == 3  ~ "Loss to follow-up",
           PREG_END_DENOM_20 == 4  ~ "Withdrew from study",
           PREG_END_DENOM_20 == 77 ~ "Not Applicable",
           TRUE                    ~ "Other unspecified reasons"
         ),
 ##Denominator C.2 - Participants with birthoutcome reported + After GA20 Met ----
 #why am i ordering the 1-7 reasons differently, because of the flowchart or how they get taken out
 BIRTHOUT_DENOM_20 = 
   case_when(N_VISITS_AFTER_ENROLL >= 1 & 
             PREG_END == 1 & PREG_END_GA >= 140 & N_BIRTH_OUTCOME_REPORTED >= 1 ~ 1, #included
             ONE_VISIT_DENOM_20 %in% c(2:7,77) ~ 77, #We already know they got kicked out
             N_MISSING_MNH09 >= 1 | N_MISSING_MNH11 >= 1 ~ 5,
             N_DTH_TIME_MISSING >= 1 | N_DOB_AFTER_DEATH >= 1 ~ 6,
             ELIGIBLE_EXCL_REASON == 1  ~ 2, #mat death
             ELIGIBLE_EXCL_REASON %in% c(2,4,5) ~ 3, #loss to followup
             ELIGIBLE_EXCL_REASON == 3 ~ 4, #withdrew from study
             is.na (PREG_END_GA) ~ 7, #preg end GA is unknown (closeout issues or date issues)
             is.na (N_BIRTH_OUTCOME_REPORTED) |  N_BIRTH_OUTCOME_REPORTED == 0 ~ 8,
             TRUE ~ 55),
 ##Denominator D - Participants with birthoutcome reported + have reached 42days PP ----
 #why am i ordering the 1-7 reasons differently, because of the flowchart or how they get taken out
 POSTPARTUM_DENOM_6 = 
   case_when(PREG_END_DENOM_20 == 1 & 
             AGE_INF_TODAY >= 42 &  
               (MAT_DEATH != 1 | is.na(MAT_DEATH_INFAGE) |  # no maternal death
                MAT_DEATH_INFAGE > 42 )  ~ 1, #included if she didn't die before 42days
             
             PREG_END_DENOM_20 %in% c(2:7,55,77) ~ 77, #We already know they got kicked out
             MAT_DEATH == 1 ~ 2, #mat death
             CLOSEOUT_TYPE == 4 ~ 4, #withdrew from study
             CLOSEOUT_TYPE %in% c(1,5, 6, 77) ~ 3, #loss to followup
             CLOSEOUT_TYPE %in% c(2) ~ 5, #followup complete
             AGE_INF_TODAY < 42 ~ 7, #infants at the time of followup had not passed earliest date for 6months visit (recode to ltfu)
             is.na (N_BIRTH_OUTCOME_REPORTED) |  N_BIRTH_OUTCOME_REPORTED == 0 ~ 6, #no infant/baby outcome recorded (possibly induced abortions)
             TRUE ~ 55),
 
             DENOM_D_INCLUDED   = if_else(POSTPARTUM_DENOM_6 == 1, 1, 0),
             DENOM_D_EXCL_CODE  = if_else(POSTPARTUM_DENOM_6 == 1, NA_real_, as.numeric(POSTPARTUM_DENOM_6)),
             DENOM_D_EXCL_LABEL = case_when(
               POSTPARTUM_DENOM_6 == 1  ~ "Included: Reached ≥42 days postpartum and not maternal death before 42 days",
               POSTPARTUM_DENOM_6 == 2  ~ "Maternal death <6weeks postpartum",
               POSTPARTUM_DENOM_6 == 3  ~ "Loss to follow-up",
               POSTPARTUM_DENOM_6 == 4  ~ "Withdrew from study",
               POSTPARTUM_DENOM_6 == 5  ~ "Loss to follow-up",
               POSTPARTUM_DENOM_6 == 6  ~ "No birth outcome recorded",
               POSTPARTUM_DENOM_6 == 7  ~ "<6 weeks by study end",
               POSTPARTUM_DENOM_6 == 77 ~ "Not Applicable",
               TRUE                     ~ "Other unspecified reasons"
             ),
 ##Denominator E - Participants with birthoutcome reported + have reached 182days PP ----
 #why am i ordering the 1-7 reasons differently, because of the flowchart or how they get taken out
 POSTPARTUM_DENOM_26 = 
   case_when(PREG_END_DENOM_20 == 1 & 
              (N_LIVEBIRTH >= 1 | N_STILLBIRTH_20WK >= 1) & 
              AGE_INF_TODAY >= 182 & 
               (MAT_DEATH != 1 | is.na(MAT_DEATH_INFAGE)  | # no maternal death
                MAT_DEATH_INFAGE > 182 )  ~ 1, #included if she didn't die before 182days
              
              POSTPARTUM_DENOM_6 %in% c(2:7,77,55) ~ 77, #We already know they got kicked out
               MAT_DEATH == 1 ~ 2, #mat death
               CLOSEOUT_TYPE == 4 ~ 4, #withdrew from study
               CLOSEOUT_TYPE %in% c(1,5, 6, 77) ~ 3, #loss to followup
               CLOSEOUT_TYPE %in% c(2) |  N_LIVEBIRTH < 1 |  is.na (N_LIVEBIRTH) ~ 5, #followup complete if she never had a livebirth but we expected her in 6weeks 
               AGE_INF_TODAY < 182 ~ 6, #infants at the time of followup had not passed earliest date for 6months visit (recode to ltfu)
             TRUE ~ 55),
         
         DENOM_E_INCLUDED   = if_else(POSTPARTUM_DENOM_26 == 1, 1, 0),
         DENOM_E_EXCL_CODE  = if_else(POSTPARTUM_DENOM_26 == 1, NA_real_, as.numeric(POSTPARTUM_DENOM_26)),
         DENOM_E_EXCL_LABEL = case_when(
           POSTPARTUM_DENOM_26 == 1  ~ "Included: Livebirth or stillbirth ≥20 weeks AND reached ≥182 days postpartum and not maternal death before 182 days",
           POSTPARTUM_DENOM_26 == 2  ~ "Maternal death <6months postpartum",
           POSTPARTUM_DENOM_26 == 3  ~ "Loss to follow-up",
           POSTPARTUM_DENOM_26 == 4  ~ "Withdrew from study",
           POSTPARTUM_DENOM_26 == 5  ~ "Loss to follow-up",
           POSTPARTUM_DENOM_26 == 6  ~ "<6 months by study end",
           POSTPARTUM_DENOM_26 == 77 ~ "Not Applicable",
           TRUE                      ~ "Other unspecified reasons"
         ))


# denom_20_miss <- remapp_denom_a %>% filter (ONE_VISIT_DENOM %in% c(1))  
table (mat_remapp_denom$ONE_VISIT_DENOM)
table (mat_remapp_denom$ONE_VISIT_DENOM_20)
table (mat_remapp_denom$PREG_END_DENOM_20)
table (mat_remapp_denom$BIRTHOUT_DENOM_20)
table (mat_remapp_denom$POSTPARTUM_DENOM_6)
table (mat_remapp_denom$POSTPARTUM_DENOM_26)

mat_remapp_denom_a <- mat_remapp_denom %>% filter (ONE_VISIT_DENOM %in% c(1,2))
mat_remapp_denom_b <-  mat_remapp_denom %>% filter (ONE_VISIT_DENOM_20 %in% c(1))
mat_remapp_denom_c <-  mat_remapp_denom %>% filter(PREG_END_DENOM_20 %in% c(1))
mat_remapp_denom_c2 <-  mat_remapp_denom %>% filter (BIRTHOUT_DENOM_20 %in% c(1))
mat_remapp_denom_d <-  mat_remapp_denom %>% filter (POSTPARTUM_DENOM_6 %in% c(1))
mat_remapp_denom_e <-  mat_remapp_denom %>% filter (POSTPARTUM_DENOM_26 %in% c(1))

base_remapp_exlusions <- mat_remapp_denom %>%
  select(
    # Core IDs
    SITE, MOMID, PREGID,
    # Key pregnancy timing
    PREG_START_DATE, PREG_END_DATE, PREG_END_GA, LAST_VIABLE_VISITDATE,
    GA_AT_LAST_VISIT, FOLLOWUP_END_DATE, AGE_INF_TODAY,
    # Key counts / endpoints
    N_VISITS_AFTER_ENROLL,PREG_END,N_BIRTH_OUTCOME_REPORTED,
    N_LIVEBIRTH,N_STILLBIRTH_20WK, MAT_DEATH, MAT_DEATH_INFAGE, CLOSEOUT_TYPE,
    # Raw denominators (codes)
    ONE_VISIT_DENOM, ONE_VISIT_DENOM_20, PREG_END_DENOM_20,POSTPARTUM_DENOM_6,
    POSTPARTUM_DENOM_26,
    # Flowchart inclusion flags
    DENOM_A_INCLUDED, DENOM_B_INCLUDED, DENOM_C_INCLUDED, DENOM_D_INCLUDED, 
    DENOM_E_INCLUDED,
    # Flowchart exclusion labels
    DENOM_A_EXCL_LABEL, DENOM_B_EXCL_LABEL, DENOM_C_EXCL_LABEL,
    DENOM_D_EXCL_LABEL, DENOM_E_EXCL_LABEL)

save(base_remapp_exlusions, file = "derived_data/base_remapp_exlusions.rda")

#*****************************************************************************
#1.5 B. Infant Denominators ----
#*****************************************************************************
#get all Infant Records
inf_remapp_denom <- INF_OUTCOMES %>%
  left_join(
    mat_remapp_denom_c %>%
      select(SITE, MOMID, PREGID, dplyr::starts_with("N_", ignore.case = TRUE)),
    by = c("SITE", "MOMID", "PREGID")
  ) %>%
  left_join(
    INF_COD %>%
      select(SITE, MOMID, PREGID, INFANTID, AGE_DEATH, DEATH_DATE),
    by = c("SITE", "MOMID", "PREGID", "INFANTID")
  ) %>%
  filter(
    dplyr::coalesce(N_STILLBIRTH_20WK, 0) >= 1 |
    dplyr::coalesce(N_LIVEBIRTH, 0) >= 1
  ) %>%
  ##Denominator F - Stillbirth + Livebirths ----
 mutate (LIVE_STILL_BIRTH_DENOM = 
          case_when(BIRTH_OUTCOME_REPORTED == 1 ~ 1, #included
                    ADJUD_NEEDED == 1  ~ 2, #adjucation needed for stillbirth/livebirth
                    FETAL_LOSS == 1 ~ 3, #fetal loss
                    BIRTH_OUTCOME_REPORTED != 1 ~ 4, #status of infant unknown
                    is.na(INFANTID) ~ 5,  #infantid not allocated 
                    TRUE ~ 55)) %>%
 ##Denominator G - Stillbirth + Livebirths >= 28wks ----
mutate (LIVE_STILL_BIRTH28_DENOM = 
          case_when(LIVE_STILL_BIRTH_DENOM == 1 & GESTAGEBIRTH_ANY_DAYS >= 196 ~ 1, #included
                    LIVE_STILL_BIRTH_DENOM %in% c(2:5) ~ 77, 
                    GESTAGEBIRTH_ANY_DAYS < 196 ~ 2,  #GA did not reach 28days
                    TRUE ~ 55))  %>%
  ##Denominator H - Livebirths ----
mutate (ALL_LIVEBIRTH_DENOM = 
            case_when(LIVE_STILL_BIRTH_DENOM == 1 & LIVEBIRTH ==1 ~ 1, #included
                        LIVE_STILL_BIRTH_DENOM %in% c(2:5) ~ 77, 
                        STILLBIRTH_28WK == 1 | STILLBIRTH_20WK == 1 ~ 2,  #stillbirth at 28wks
                            TRUE ~ 55))

table (inf_remapp_denom$LIVE_STILL_BIRTH_DENOM)
table (inf_remapp_denom$LIVE_STILL_BIRTH28_DENOM)
table (inf_remapp_denom$ALL_LIVEBIRTH_DENOM)

inf_remapp_denom_f <- inf_remapp_denom %>% filter (LIVE_STILL_BIRTH_DENOM %in% c(1))
inf_remapp_denom_g <- inf_remapp_denom %>% filter (LIVE_STILL_BIRTH28_DENOM %in% c(1))
inf_remapp_denom_h <-  inf_remapp_denom %>% filter (ALL_LIVEBIRTH_DENOM %in% c(1))

test <- inf_remapp_denom_f %>% filter (LIVE_STILL_BIRTH28_DENOM != 1) %>% 
  select(INFANTID, SITE, GESTAGEBIRTH_ANY_DAYS, STILLBIRTH_20WK, 
         STILLBIRTH_28WK, STILLBIRTH_TIMING)
save(inf_remapp_denom, file = "derived_data/inf_remapp_exlusions.rda")

#*****************************************************************************
#2. Hemoglobin data processing ----
#*****************************************************************************
#*prepare data
prep_hb1 <- df_maternal %>% 
  filter (REMAPP_SCRN ==1 & ENROLL_NO_ISSUES == 1) %>%
  select(MOMID, PREGID, SITE, PREG_START_DATE, PREG_END_DATE, M03_SMOKE_OECOCCUR, REMAPP_SCRN, ENROLL_NO_ISSUES) %>%
  left_join(mnh08_long %>% rename (VISITDATE = M08_LBSTDAT), by = c("SITE", "MOMID", "PREGID")) 

prep_hb2 <- df_maternal %>% 
  filter (REMAPP_SCRN ==1 & ENROLL_NO_ISSUES == 1) %>%
  select(MOMID, PREGID, SITE, PREG_START_DATE, PREG_END_DATE, M03_SMOKE_OECOCCUR, REMAPP_SCRN, ENROLL_NO_ISSUES) %>%
  left_join(mnh06_long %>% rename (VISITDATE = M06_DIAG_VSDAT), by = c("SITE", "MOMID", "PREGID")) 

all_remapp_hb <- prep_hb1  %>%
  full_join(prep_hb2, by = c("SITE", "MOMID", "PREGID","VISITDATE")) %>%
  mutate(
    PREG_START_DATE   = coalesce(PREG_START_DATE.x,   PREG_START_DATE.y),
    PREG_END_DATE     = coalesce(PREG_END_DATE.x,     PREG_END_DATE.y),
    M03_SMOKE_OECOCCUR = coalesce(M03_SMOKE_OECOCCUR.x, M03_SMOKE_OECOCCUR.y),
    REMAPP_SCRN       = coalesce(REMAPP_SCRN.x,       REMAPP_SCRN.y),
    ENROLL_NO_ISSUES  = coalesce(ENROLL_NO_ISSUES.x,  ENROLL_NO_ISSUES.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>%
  mutate(
    ## HB CBC (MNH08) ----
    visit_type_hb = M08_TYPE_VISIT,
    ga_days_hb    = ga_days_08,
    age_days_hb   = age_days_08,
    
    age_wks_hb = case_when(
      visit_type_hb >= 6 ~ floor(as.numeric(age_days_hb) / 7),
      visit_type_hb < 6  ~ NA_real_,
      TRUE ~ NA_real_
    ),
    
    ga_wks_hb = case_when(
      is.na(ga_days_hb) & visit_type_hb < 6 ~ NA_real_,
      visit_type_hb < 6 ~ ga_days_hb / 7,
      TRUE ~ NA_real_
    ),
    
    trimester_hb = case_when(
      visit_type_hb >= 6 ~ NA_real_,
      ga_wks_hb > 0  & ga_wks_hb < 14 ~ 1,
      ga_wks_hb >= 14 & ga_wks_hb < 28 ~ 2,
      ga_wks_hb >= 28 & ga_wks_hb <= 40 ~ 3,
      TRUE ~ NA_real_
    ),
    
    hb_alti = case_when(
      SITE %in% c("Kenya","Zambia") & !is.na(M08_CBC_HB_LBORRES) ~ M08_CBC_HB_LBORRES - 0.8,
      !is.na(M08_CBC_HB_LBORRES) ~ M08_CBC_HB_LBORRES,
      TRUE ~ NA_real_
    ),
    
    hb_adj = if_else(M03_SMOKE_OECOCCUR == 1, 0.3, 0, missing = 0),
    hb     = hb_alti - hb_adj,
    
    ## HB POC (MNH06) ----
    visit_type_hb_poc = M06_TYPE_VISIT,
    ga_days_hb_poc    = ga_days_06,
    age_days_hb_poc   = age_days_06,
    
    age_wks_hb_poc = case_when(
      visit_type_hb_poc >= 6 ~ floor(as.numeric(age_days_hb_poc) / 7),
      visit_type_hb_poc < 6  ~ NA_real_,
      TRUE ~ NA_real_
    ),
    
    ga_wks_hb_poc = case_when(
      is.na(ga_days_hb_poc) & visit_type_hb_poc < 6 ~ NA_real_,
      visit_type_hb_poc < 6 ~ ga_days_hb_poc / 7,
      TRUE ~ NA_real_
    ),
    
    trimester_hb_poc = case_when(
      visit_type_hb_poc >= 6 ~ NA_real_,
      ga_wks_hb_poc > 0  & ga_wks_hb_poc < 14 ~ 1,
      ga_wks_hb_poc >= 14 & ga_wks_hb_poc < 28 ~ 2,
      ga_wks_hb_poc >= 28 & ga_wks_hb_poc <= 40 ~ 3,
      TRUE ~ NA_real_
    ),
    
    hb_poc_alti = case_when(
      SITE %in% c("Kenya","Zambia") & !is.na(M06_HB_POC_LBORRES) ~ M06_HB_POC_LBORRES - 0.8,
      !is.na(M06_HB_POC_LBORRES) ~ M06_HB_POC_LBORRES,
      TRUE ~ NA_real_
    ),
    
    hb_poc_adj = if_else(M03_SMOKE_OECOCCUR == 1, 0.3, 0, missing = 0),
    hb_poc     = hb_poc_alti - hb_poc_adj,
    
    ##sphb
    sphb_alti = case_when(
      SITE %in% c("Kenya","Zambia") & !is.na(M06_SPHB_LBORRES) ~ M06_SPHB_LBORRES - 0.8,
      !is.na(M06_SPHB_LBORRES) ~ M06_SPHB_LBORRES,
      TRUE ~ NA_real_
    ),
    
    sphb_adj = if_else(M03_SMOKE_OECOCCUR == 1, 0.3, 0, missing = 0),
    sphb     = sphb_alti - sphb_adj,
    
    ## Combined Hb ----
    visit_type_86 = case_when(
      # ANC: force CBC visit type if either indicates ANC
      (visit_type_hb %in% c(1:5, 13) | visit_type_hb_poc %in% c(1:5, 13)) ~ visit_type_hb,
      
      # PNC: prefer CBC visit type, else take POC visit type
      (visit_type_hb %in% c(6:12, 14) | visit_type_hb_poc %in% c(6:12, 14)) ~ coalesce(visit_type_hb, visit_type_hb_poc),
      
      TRUE ~ NA_real_
    ),
    
    ga_wks_86 = ga_wks_hb,
    
    trimester_86 = trimester_hb, #use CBC trimester for ANC
    
    age_wks_86 = case_when(
      visit_type_86 %in% c(1:5, 13)  ~ NA_real_,
      visit_type_86 %in% c(6:12, 14) ~ coalesce(age_wks_hb, age_wks_hb_poc),
      TRUE ~ NA_real_
    ),
    
    hb_cbc_poc = case_when(
      visit_type_86 %in% c(1:5, 13) ~ hb,                    # ANC: CBC only
      visit_type_86 %in% c(6:12,14) ~ coalesce(hb, hb_poc),  # PNC: prefer CBC else POC
      TRUE ~ NA_real_
    ),
    
    hb_source = case_when(
      visit_type_86 %in% c(1:5, 13) & !is.na(hb) ~ "CBC",
      !is.na(hb) ~ "CBC",
      !is.na(hb_poc) ~ "POC",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    hb_level = case_when(
      trimester_86 %in% c(1,3) & hb_cbc_poc > 0 & hb_cbc_poc < 7 ~ 1,
      trimester_86 %in% c(1,3) & hb_cbc_poc >= 7 & hb_cbc_poc < 10 ~ 2,
      trimester_86 %in% c(1,3) & hb_cbc_poc >= 10 & hb_cbc_poc < 11  ~ 3,
      trimester_86 %in% c(1,3) & hb_cbc_poc >= 11 & hb_cbc_poc <= 13 ~ 4,
      trimester_86 %in% c(1,2,3) & hb_cbc_poc > 13 & hb_cbc_poc < 15 ~ 5,
      trimester_86 %in% c(1,2,3) & hb_cbc_poc >= 15 ~ 6,
      
      trimester_86 == 2 & hb_cbc_poc > 0 & hb_cbc_poc < 7 ~ 1,
      trimester_86 == 2 & hb_cbc_poc >= 7 & hb_cbc_poc < 9.5 ~ 2,
      trimester_86 == 2 & hb_cbc_poc >= 9.5 & hb_cbc_poc < 10.5  ~ 3,
      trimester_86 == 2 & hb_cbc_poc >= 10.5 & hb_cbc_poc <= 13 ~ 4,
      
      visit_type_86 <= 9 & hb_cbc_poc > 0 & hb_cbc_poc < 7 ~ 1,
      visit_type_86 <= 9 & hb_cbc_poc >= 7 & hb_cbc_poc < 10 ~ 2,
      visit_type_86 <= 9 & hb_cbc_poc >= 10 & hb_cbc_poc < 11 ~ 3,
      visit_type_86 <= 9 & hb_cbc_poc >= 11 & hb_cbc_poc <= 13 ~ 4,
      visit_type_86 <= 9 & hb_cbc_poc > 13 & hb_cbc_poc < 15 ~ 5,
      visit_type_86 <= 9 & hb_cbc_poc >= 15 ~ 6,
      TRUE ~ NA_real_)) %>% 
  select(-c(hb_alti, hb_poc_alti, sphb_alti, M03_SMOKE_OECOCCUR)) 

#combined Hb dataset (CBC/POC)
df_hb_long <- all_remapp_hb %>% filter (!is.na(hb_cbc_poc)) 

df_hb_long$hb_level <- factor(
  df_hb_long$hb_level, 
  levels = c(1,2,3,4,5,6),
  labels = c("Severe", "Moderate", "Mild", "Normal", "High hb 13-<15g/dl", "High hb >=15g/dl")
)

#Only CBC Hb dataset
df_hb_long2 <- df_hb_long %>% 
  filter(hb > 0) %>%
  transmute (SITE, MOMID, PREGID, hb, hb_level, 
             trimester = trimester_hb, 
             visit_type = visit_type_hb, 
             ga_wks = round(ga_wks_hb, 0), 
             ga_days_hb, 
             age_wks = round(age_wks_hb, 0), 
             age_days_hb)


## median Hb Value for ANC Visits ----
#************calculate median of the hb value for ANC visits per SITE***********************************
med_hb_anc_site <- df_hb_long2 %>% 
  filter(visit_type < 6) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_exm_anc with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_exm_anc <- df_hb_long2 %>% 
  filter(visit_type %in% 1:5) %>%
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
         trimester, ga_wks, ga_days_hb, visit_type) %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_anc, file = "derived_data/df_hb_exm_anc.rda")

table(df_hb_exm_anc$trimester)

 #Test if the randomization is working welll and the right trimeter is being chosen
test_trimester_match <- df_hb_exm_anc %>%
  select(SITE, MOMID, PREGID, hb_exm, hb, hb_max, hb_min, trimester, ga_wks, visit_type) %>%
  left_join(
    df_hb_long2 %>%
      filter(visit_type %in% 1:5) %>%
      select(SITE, MOMID, PREGID, hb_orig = hb, trimester, TYPE_VISIT = visit_type ),
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
         trimester, ga_wks, visit_type) %>% 
  mutate (hb = hb_exm)

save(df_hb_exm_t12, file = "derived_data/df_hb_exm_t12.rda")

table(df_hb_exm_t12$trimester)


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
###anemia at anc ----
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

###anemia at pnc ----
pnc_hb_df_1 <- df_hb_long %>%
  mutate(
    TYPE_VISIT = as.numeric(visit_type_86),
    AGE_WKS    = as.numeric(age_wks_86),
    HB_CBC_PC  = as.numeric(hb_cbc_poc),
    HB_SOURCE  = hb_source
  ) %>%
  filter(TYPE_VISIT %in% c(9, 10, 11, 12)) %>%
  mutate(
    suffix = case_when(
      TYPE_VISIT == 9  ~ "PNC4",
      TYPE_VISIT == 10 ~ "PNC6",
      TYPE_VISIT == 11 ~ "PNC26",
      TYPE_VISIT == 12 ~ "PNC52" )) %>%
  select(SITE, MOMID, PREGID, TYPE_VISIT, suffix, AGE_WKS, HB_CBC_PC, HB_SOURCE, hb_level)



test<- df_hb_long %>% filter (visit_type_hb_poc == 10)

pnc_hb_df_2 <- pnc_hb_df_1 %>%
  group_by(SITE, TYPE_VISIT) %>%
  group_modify(~{
    med_site <- .x %>%
      group_by(MOMID, PREGID) %>%
      summarise(hb_mom_med = median(HB_CBC_PC, na.rm = TRUE), .groups = "drop") %>%
      summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")
    
    .x %>%
      tidyr::crossing(med_site) %>%   # adds hb_med to every row
      mutate(hb_dis = abs(HB_CBC_PC - hb_med)) %>%
      group_by(MOMID, PREGID) %>%     # <-- removed TYPE_VISIT here
      filter(hb_dis == max(hb_dis, na.rm = TRUE)) %>%
      slice_sample(n = 1) %>%
      ungroup()
  }) %>%
  ungroup() %>%
  mutate(
    ANEMIA = case_when(
      HB_CBC_PC > 0  & HB_CBC_PC < 7   ~ 3,
      HB_CBC_PC >= 7 & HB_CBC_PC < 10  ~ 2,
      HB_CBC_PC >= 10 & HB_CBC_PC < 11 ~ 1,
      HB_CBC_PC >= 11                 ~ 0,
      TRUE ~ NA_real_))


pnc_hb_df <- pnc_hb_df_2 %>% 
  select(SITE, MOMID, PREGID, suffix, HB_CBC_PC, HB_SOURCE, ANEMIA) 

pnc_hb_wide <- pnc_hb_df %>%
  pivot_wider(
    id_cols = c(SITE, MOMID, PREGID),
    names_from = suffix,
    values_from = c(HB_CBC_PC, HB_SOURCE, ANEMIA),
    names_glue = "{.value}_{suffix}"
  )

MAT_ANEMIA <- MAT_FLOWCHART %>% 
  filter (REMAPP_SCRN ==1 & ENROLL_NO_ISSUES == 1) %>%
  select(SITE, SCRNID, MOMID, PREGID) %>%
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

mat_dth <- MAT_COD %>%
  filter(DEATH_DATE_MISS == 0, COD_MISS == 0) %>%
  select(SITE, MOMID, PREGID, DEATH_DATE, COD, COD_TEXT) %>%
  mutate(
    PREG_CAUSE = case_when(
      COD %in% c(2, 3, 4,7, 12, 13) ~ 1L,
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

# prep_mat_compo <- mat_remapp_denom_a %>%
#   select(SITE, MOMID, PREGID, ONE_VISIT_DENOM) %>%
#   left_join(df_hb_exm_anc, by = c("MOMID","PREGID","SITE")) %>%
#   left_join(MAT_NEAR_MISS, by = c("SITE", "MOMID","PREGID")) %>% 
#   left_join(mat_dth %>% select(SITE, MOMID, PREGID, PREG_CAUSE),
#             by = c("SITE", "MOMID","PREGID")) %>% 
#   left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE, PREG_END),
#             by = c("SITE", "MOMID","PREGID")) %>%
#   mutate (  UPLOADDATE    = as.Date(UploadDate),
#             DAYS_PP42 = as.numeric(UPLOADDATE - PREG_END_DATE)) %>%
#   filter (PREG_END == 1 & DAYS_PP42 >= 42 & is.na(NEARMISS_MISS_FORMS)) %>%
#   mutate (  mat_out_compo = case_when (
#             #mat_death from preg causes
#             PREG_CAUSE == 1 | 
#             #severe complications: missing sepsis
#             HEM_PPH_SEV==1 | PREECLAMPSIA_SEV == 1 | UTERINE_RUP == 1 |
#             #critical intervention
#             MAT_ICU == 1 | HYSTERECTOMY == 1 | LAPAROTOMY == 1| TRANSFUSION_NEARMISS == 1|
#             #Life threatning conditions: Organ failure/dysfunction
#             ORG_FAIL == 1 | ORG_FAIL_HRT == 1 | ORG_FAIL_RESP == 1 | 
#             ORG_FAIL_RENAL ==1 |ORG_FAIL_OTHR ==1| ORG_FAIL_LIVER == 1 | 
#             ORG_FAIL_NEUR == 1 | ORG_FAIL_UTER == 1 | ORG_FAIL_HEM == 1 ~ 1L,
#             TRUE ~ 0L))  %>%
#   filter(mat_out_compo %in% c(0, 1))


prep_mat_compo <- mat_remapp_denom_a %>%
  select(SITE, MOMID, PREGID, starts_with ("ELIBILITY_"), ONE_VISIT_DENOM) %>%
  left_join(df_hb_exm_anc, by = c("MOMID", "PREGID", "SITE")) %>%
  left_join(MAT_NEAR_MISS, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mat_dth %>% select(SITE, MOMID, PREGID, PREG_CAUSE),
            by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE, PREG_END),
            by = c("SITE", "MOMID", "PREGID")) %>%
  mutate(
    UPLOADDATE = as.Date(UploadDate),
    DAYS_PP42 = as.numeric(UPLOADDATE - PREG_END_DATE)
  ) %>%
  mutate(
    # First, check if any criteria is 1 (same as your original)
    any_criteria_1 = case_when(
      PREG_CAUSE == 1 |
        HEM_PPH_SEV == 1 | PREECLAMPSIA_SEV == 1 | UTERINE_RUP == 1 |
        MAT_ICU == 1 | HYSTERECTOMY == 1 | LAPAROTOMY == 1 | TRANSFUSION_NEARMISS == 1 |
        ORG_FAIL == 1 | ORG_FAIL_HRT == 1 | ORG_FAIL_RESP == 1 |
        ORG_FAIL_RENAL == 1 | ORG_FAIL_OTHR == 1 | ORG_FAIL_LIVER == 1 |
        ORG_FAIL_NEUR == 1 | ORG_FAIL_UTER == 1 | ORG_FAIL_HEM == 1 ~ 1,
      TRUE ~ 0
    ),
    # Check if any criteria is 0
    any_criteria_0 = case_when(
      PREG_CAUSE == 0 |
        HEM_PPH_SEV == 0 | PREECLAMPSIA_SEV == 0 | UTERINE_RUP == 0 |
        MAT_ICU == 0 | HYSTERECTOMY == 0 | LAPAROTOMY == 0 | TRANSFUSION_NEARMISS == 0 |
        ORG_FAIL == 0 | ORG_FAIL_HRT == 0 | ORG_FAIL_RESP == 0 |
        ORG_FAIL_RENAL == 0 | ORG_FAIL_OTHR == 0 | ORG_FAIL_LIVER == 0 |
        ORG_FAIL_NEUR == 0 | ORG_FAIL_UTER == 0 | ORG_FAIL_HEM == 0 ~ 1,
      TRUE ~ 0
    ),
    
    # Check if ALL criteria are NA
    all_criteria_na = case_when(
      is.na(PREG_CAUSE) &
        is.na(HEM_PPH_SEV) & is.na(PREECLAMPSIA_SEV) & is.na(UTERINE_RUP) &
        is.na(MAT_ICU) & is.na(HYSTERECTOMY) & is.na(LAPAROTOMY) & is.na(TRANSFUSION_NEARMISS) &
        is.na(ORG_FAIL) & is.na(ORG_FAIL_HRT) & is.na(ORG_FAIL_RESP) &
        is.na(ORG_FAIL_RENAL) & is.na(ORG_FAIL_OTHR) & is.na(ORG_FAIL_LIVER) &
        is.na(ORG_FAIL_NEUR) & is.na(ORG_FAIL_UTER) & is.na(ORG_FAIL_HEM) ~ 1,
      TRUE ~ 0
    ),
    
    # Now apply the 4-category logic
    mat_out_compo = case_when(
      any_criteria_1 == 1 ~ 1,
      any_criteria_1 == 0 & any_criteria_0 == 1 ~ 0,
      all_criteria_na == 1 ~ 55,
    ))  %>% 
  filter(mat_out_compo %in% c(0, 1))

dim(prep_mat_compo)

df_mat_compo <- prep_mat_compo %>%
  select(SITE, MOMID, PREGID, mat_out_compo, trimester, ga_wks, hb) %>%
  filter(mat_out_compo %in% c(0, 1)) %>%
  filter(!is.na(hb))

dim(df_mat_compo)

remapp_mat_composite <- outcome_df_rr(
  df = df_mat_compo,
  hb_var = "hb",
  outcome_var = "mat_out_compo",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_composite"
)

df_mat_compo_trim1 <- df_hb_exm_trim1 %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, mat_out_compo)) %>% 
  filter(mat_out_compo %in% c(0, 1) & !is.na(hb))

df_mat_compo_trim2 <- df_hb_exm_trim2 %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, mat_out_compo)) %>% 
  filter(mat_out_compo %in% c(0, 1) & !is.na(hb))

df_mat_compo_trim3 <- df_hb_exm_trim3 %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, mat_out_compo)) %>% 
  filter(mat_out_compo %in% c(0, 1) & !is.na(hb))


df_mat_analysis <- df_mat_analysis %>% 
  left_join(prep_mat_compo %>% select(SITE, MOMID, PREGID, MAT_COMPO=mat_out_compo),
            by = c("SITE", "MOMID","PREGID"))

##postpartum hemorrhage----
df_mat_pph <- df_hb_exm_anc %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH),
                                by = c("SITE", "MOMID","PREGID")) %>% 
  filter(HEM_PPH %in% c(0, 1))

pph_prep <- mat_remapp_denom_c %>%
  select(SITE, MOMID, PREGID, PREG_END_DENOM_20) %>%
  left_join(df_hb_exm_anc, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(MAT_HEMORRHAGE, by = c("SITE", "MOMID", "PREGID")) %>%
  mutate(
    miss_pph = case_when(
      HEM_PPH %in% c(0L, 1L) ~ 0L,
      TRUE ~ 1L
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (pph_prep$HEM_PPH, pph_prep$miss_pph )

missing_pph <- pph_prep %>%
  filter(miss_pph == 1) 

table (missing_pph$miss_pph )

df_mat_pph <- pph_prep %>%
  filter(HEM_PPH >= 0) %>% 
  select(SITE, MOMID, PREGID, HEM_PPH, starts_with("hb"), trimester, starts_with("ga_"))
dim (df_mat_pph)

df_mat_pph <- df_mat_pph %>%
  filter(!is.na(hb))
dim (df_mat_pph)

remapp_mat_pph <- outcome_df_rr(
  df = df_mat_pph,
  hb_var = "hb",
  outcome_var = "HEM_PPH",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_pph"
)


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

ppa_pnc6_prep <- mat_remapp_denom_d %>%
  select(SITE, MOMID, PREGID, PREG_END_DENOM_20) %>%
  left_join(df_hb_exm_anc, by = c("MOMID", "PREGID", "SITE")) %>%
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>%
  mutate(
    
    miss_ppa_pnc6 = case_when(
      is.na(ANEMIA_PNC6) ~ 1L,
      TRUE ~ 0
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (ppa_pnc6_prep$ppa_pnc6, ppa_pnc6_prep$miss_ppa_pnc6 )


#Why is anemia at PNC6 missing?
missing_base <- ppa_pnc6_prep %>%
  filter(is.na(ppa_pnc6)) %>%
  left_join(MAT_ENDPOINTS %>% 
              select(SITE, MOMID, PREGID, PREG_END_DATE),
            by = c("SITE", "MOMID","PREGID")) %>%
  mutate(
    DOB = clean_date(PREG_END_DATE),
    pnc6_start = DOB + lubridate::days(42),
    pnc6_end   = DOB + lubridate::days(104)
  ) %>%
  select(SITE, MOMID, PREGID, DOB, pnc6_start, pnc6_end)

mnh06_pnc <- raw_mnh06 %>%
  filter(M06_TYPE_VISIT %in% c(6:12, 14)) %>%
  transmute(
    SITE, MOMID, PREGID,
    visit_date_06 = clean_date(M06_DIAG_VSDAT),
    type_visit_06 = M06_TYPE_VISIT,
    visit_done_06 = M06_MAT_VISIT_MNH06,
    vital_06      = M06_MAT_VITAL_MNH06,
    hb_poc        = suppressWarnings(readr::parse_number(as.character(M06_HB_POC_LBORRES)))
  )

mnh08_pnc <- raw_mnh08 %>%
  filter(M08_TYPE_VISIT %in% c(6:12, 14)) %>%
  transmute(
    SITE, MOMID, PREGID,
    visit_date_08 = clean_date(M08_LBSTDAT),
    type_visit_08 = M08_TYPE_VISIT,
    visit_done_08 = M08_MAT_VISIT_MNH08,
    vital_08      = M08_MAT_VITAL_MNH08,
    hb_cbc        = suppressWarnings(readr::parse_number(as.character(M08_CBC_HB_LBORRES)))
  )

mnh06_pnc6_summary <- missing_base %>%
  left_join(mnh06_pnc, by = c("SITE","MOMID","PREGID")) %>%
  mutate(
    in_pnc6_06   = !is.na(visit_date_06) & visit_date_06 >= pnc6_start & visit_date_06 <= pnc6_end,
    later_06     = !is.na(visit_date_06) & visit_date_06 > pnc6_end,
    done_alive_06 = visit_done_06 %in% c(1,2) & vital_06 == 1,
    hb_valid_06  = !is.na(hb_poc) & hb_poc > 0
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    has_mnh06_in_pnc6        = any(in_pnc6_06, na.rm = TRUE),
    has_done_alive_in_pnc6_06 = any(in_pnc6_06 & done_alive_06, na.rm = TRUE),
    has_hb_in_pnc6_06        = any(in_pnc6_06 & done_alive_06 & hb_valid_06, na.rm = TRUE),
    
    has_later_done_alive_06  = any(later_06 & done_alive_06, na.rm = TRUE),
    
    # optional: capture a representative hb value from window
    hb_poc_best_pnc6 = suppressWarnings(max(hb_poc[in_pnc6_06 & done_alive_06 & hb_valid_06], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    hb_poc_best_pnc6 = ifelse(is.infinite(hb_poc_best_pnc6), NA_real_, hb_poc_best_pnc6)
  )

mnh08_pnc6_summary <- missing_base %>%
  left_join(mnh08_pnc, by = c("SITE","MOMID","PREGID")) %>%
  mutate(
    in_pnc6_08   = !is.na(visit_date_08) & visit_date_08 >= pnc6_start & visit_date_08 <= pnc6_end,
    later_08     = !is.na(visit_date_08) & visit_date_08 > pnc6_end,
    done_alive_08 = visit_done_08 %in% c(1,2) & vital_08 == 1,
    hb_valid_08  = !is.na(hb_cbc) & hb_cbc > 0
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    has_mnh08_in_pnc6         = any(in_pnc6_08, na.rm = TRUE),
    has_done_alive_in_pnc6_08 = any(in_pnc6_08 & done_alive_08, na.rm = TRUE),
    has_hb_in_pnc6_08         = any(in_pnc6_08 & done_alive_08 & hb_valid_08, na.rm = TRUE),
    
    has_later_done_alive_08   = any(later_08 & done_alive_08, na.rm = TRUE),
    
    hb_cbc_best_pnc6 = suppressWarnings(max(hb_cbc[in_pnc6_08 & done_alive_08 & hb_valid_08], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    hb_cbc_best_pnc6 = ifelse(is.infinite(hb_cbc_best_pnc6), NA_real_, hb_cbc_best_pnc6)
  )

missing_ppa6 <- missing_base %>%
  left_join(mnh06_pnc6_summary, by = c("SITE","MOMID","PREGID")) %>%
  left_join(mnh08_pnc6_summary, by = c("SITE","MOMID","PREGID")) %>%
  mutate(
    # any in-window completed alive visit on either form
    has_done_alive_in_pnc6 = coalesce(has_done_alive_in_pnc6_06, FALSE) |
      coalesce(has_done_alive_in_pnc6_08, FALSE),
    
    # any later completed alive visit on either form
    has_later_done_alive = coalesce(has_later_done_alive_06, FALSE) |
      coalesce(has_later_done_alive_08, FALSE),
    
    # Hb exists from either source in window
    has_hb_in_pnc6 = coalesce(has_hb_in_pnc6_06, FALSE) |
      coalesce(has_hb_in_pnc6_08, FALSE),
    
    # harmonized Hb value (CBC preferred)
    hb_pnc6 = dplyr::coalesce(hb_cbc_best_pnc6, hb_poc_best_pnc6),
    
    hb_source = dplyr::case_when(
      !is.na(hb_cbc_best_pnc6) ~ "cbc (mnh08)",
      !is.na(hb_poc_best_pnc6) ~ "poc (mnh06)",
      TRUE ~ NA_character_
    ),
    
    # “form presence in window” (not necessarily completed)
    has_any_form_in_pnc6 = coalesce(has_mnh06_in_pnc6, FALSE) |
      coalesce(has_mnh08_in_pnc6, FALSE)
  ) %>%
  mutate(
    missing_reason_code = case_when(
      # Not missing if Hb present from either source
      has_hb_in_pnc6 ~ 0L,
      
      # 3) Completed alive in-window visit exists, but Hb missing/invalid on both sources
      has_done_alive_in_pnc6 & !has_hb_in_pnc6 ~ 3L,
      
      # 2) No completed in-window visit/Hb, but later completed alive visit exists
      !has_done_alive_in_pnc6 & has_later_done_alive ~ 2L,
      
      # 1) No in-window evidence and no later follow-up
      TRUE ~ 1L
    ),
    
    missing_ppa6 = case_when(
      missing_reason_code == 0L ~ "Not missing: Hb present (MNH06 and/or MNH08)",
      missing_reason_code == 1L ~ "No PNC6 MNH06/MNH08 in window and no later postnatal visits (forms missing / no follow-up)",
      missing_reason_code == 2L ~ "PNC6 missed (no in-window MNH06/MNH08 Hb) but later completed alive visit(s) exist",
      missing_reason_code == 3L ~ "PNC6 visit completed (alive) but Hb missing/<=0 on both MNH06 POC and MNH08 CBC",
      TRUE ~ "Review"
    )
  ) 

table (missing_ppa6$missing_ppa6, missing_ppa6$SITE)
dim (missing_ppa6)

df_mat_ppa_pnc6 <- ppa_pnc6_prep %>%
  filter(ppa_pnc6 >= 0) 
dim (df_mat_ppa_pnc6)

df_mat_ppa_pnc6 <- df_mat_ppa_pnc6 %>%
  filter(!is.na(hb))
dim (df_mat_ppa_pnc6)

remapp_mat_ppa_pnc6 <- outcome_df_rr(
  df = df_mat_ppa_pnc6,
  hb_var = "hb",
  outcome_var = "ppa_pnc6",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_ppa_pnc6"
)

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
ppa_pnc26_prep <- mat_remapp_denom_e %>%
  select(SITE, MOMID, PREGID, PREG_END_DENOM_20, N_STILLBIRTH_20WK, N_LIVEBIRTH) %>%
  left_join(df_hb_exm_anc, by = c("MOMID", "PREGID", "SITE")) %>%
  left_join(df_mat_anemia, by = c("SITE", "MOMID","PREGID")) %>%
  mutate(    
    miss_ppa_pnc26 = case_when(
      is.na(ANEMIA_PNC6) ~ 1L,
      TRUE ~ 0
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (ppa_pnc26_prep$ppa_pnc26, ppa_pnc26_prep$miss_ppa_pnc26 )


#Why is anemia at PNC26 missing?
missing_26_base <- ppa_pnc26_prep %>%
  filter(is.na(ppa_pnc26)) %>%
  left_join(MAT_ENDPOINTS %>% 
              select(SITE, MOMID, PREGID, PREG_END_DATE),
            by = c("SITE", "MOMID", "PREGID")) %>%
  mutate(
    DOB = clean_date(PREG_END_DATE),
    pnc26_start = DOB + lubridate::days(182),
    pnc26_end   = DOB + lubridate::days(279)
  ) %>%
  select(SITE, MOMID, PREGID, DOB, pnc26_start, pnc26_end)

mnh06_pnc26_summary <- missing_26_base %>%
  left_join(mnh06_pnc, by = c("SITE","MOMID","PREGID")) %>%
  mutate(
    in_pnc26_06   = !is.na(visit_date_06) & visit_date_06 >= pnc26_start & visit_date_06 <= pnc26_end,
    later_06     = !is.na(visit_date_06) & visit_date_06 > pnc26_end,
    done_alive_06 = visit_done_06 %in% c(1,2) & vital_06 == 1,
    hb_valid_06  = !is.na(hb_poc) & hb_poc > 0
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    has_mnh06_in_pnc26        = any(in_pnc26_06, na.rm = TRUE),
    has_done_alive_in_pnc26_06 = any(in_pnc26_06 & done_alive_06, na.rm = TRUE),
    has_hb_in_pnc26_06        = any(in_pnc26_06 & done_alive_06 & hb_valid_06, na.rm = TRUE),
    
    has_later_done_alive_06  = any(later_06 & done_alive_06, na.rm = TRUE),
    
    # optional: capture a representative hb value from window
    hb_poc_best_pnc26 = suppressWarnings(max(hb_poc[in_pnc26_06 & done_alive_06 & hb_valid_06], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    hb_poc_best_pnc26 = ifelse(is.infinite(hb_poc_best_pnc26), NA_real_, hb_poc_best_pnc26)
  )

mnh08_pnc26_summary <- missing_26_base %>%
  left_join(mnh08_pnc, by = c("SITE", "MOMID", "PREGID")) %>%
  mutate(
    in_pnc26_08   = !is.na(visit_date_08) & visit_date_08 >= pnc26_start & visit_date_08 <= pnc26_end,
    later_08     = !is.na(visit_date_08) & visit_date_08 > pnc26_end,
    done_alive_08 = visit_done_08 %in% c(1,2) & vital_08 == 1,
    hb_valid_08  = !is.na(hb_cbc) & hb_cbc > 0
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    has_mnh08_in_pnc26         = any(in_pnc26_08, na.rm = TRUE),
    has_done_alive_in_pnc26_08 = any(in_pnc26_08 & done_alive_08, na.rm = TRUE),
    has_hb_in_pnc26_08         = any(in_pnc26_08 & done_alive_08 & hb_valid_08, na.rm = TRUE),
    
    has_later_done_alive_08   = any(later_08 & done_alive_08, na.rm = TRUE),
    
    hb_cbc_best_pnc26 = suppressWarnings(max(hb_cbc[in_pnc26_08 & done_alive_08 & hb_valid_08], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    hb_cbc_best_pnc26 = ifelse(is.infinite(hb_cbc_best_pnc26), NA_real_, hb_cbc_best_pnc26)
  )

missing_ppa26 <- missing_26_base %>%
  left_join(mnh06_pnc26_summary, by = c("SITE","MOMID","PREGID")) %>%
  left_join(mnh08_pnc26_summary, by = c("SITE","MOMID","PREGID")) %>%
  mutate(
    # any in-window completed alive visit on either form
    has_done_alive_in_pnc26 = coalesce(has_done_alive_in_pnc26_06, FALSE) |
      coalesce(has_done_alive_in_pnc26_08, FALSE),
    
    # any later completed alive visit on either form
    has_later_done_alive = coalesce(has_later_done_alive_06, FALSE) |
      coalesce(has_later_done_alive_08, FALSE),
    
    # Hb exists from either source in window
    has_hb_in_pnc26 = coalesce(has_hb_in_pnc26_06, FALSE) |
      coalesce(has_hb_in_pnc26_08, FALSE),
    
    # harmonized Hb value (CBC preferred)
    hb_pnc26 = dplyr::coalesce(hb_cbc_best_pnc26, hb_poc_best_pnc26),
    
    hb_source = dplyr::case_when(
      !is.na(hb_cbc_best_pnc26) ~ "cbc (mnh08)",
      !is.na(hb_poc_best_pnc26) ~ "poc (mnh06)",
      TRUE ~ NA_character_
    ),
    
    # “form presence in window” (not necessarily completed)
    has_any_form_in_pnc26 = coalesce(has_mnh06_in_pnc26, FALSE) |
      coalesce(has_mnh08_in_pnc26, FALSE)
  ) %>%
  mutate(
    missing_reason_code = case_when(
      # Not missing if Hb present from either source
      has_hb_in_pnc26 ~ 0L,
      
      # 3) Completed alive in-window visit exists, but Hb missing/invalid on both sources
      has_done_alive_in_pnc26 & !has_hb_in_pnc26 ~ 3L,
      
      # 2) No completed in-window visit/Hb, but later completed alive visit exists
      !has_done_alive_in_pnc26 & has_later_done_alive ~ 2L,
      
      # 1) No in-window evidence and no later follow-up
      TRUE ~ 1L
    ),
    
    missing_ppa26 = case_when(
      missing_reason_code == 0L ~ "Not missing: Hb present (MNH06 and/or MNH08)",
      missing_reason_code == 1L ~ "No PNC26 MNH06/MNH08 in window and no later postnatal visits (forms missing / no follow-up)",
      missing_reason_code == 2L ~ "PNC26 missed (no in-window MNH06/MNH08 Hb) but later completed alive visit(s) exist",
      missing_reason_code == 3L ~ "PNC26 visit completed (alive) but Hb missing/<=0 on both MNH06 POC and MNH08 CBC",
      TRUE ~ "Review"
    )
  ) 

table (missing_ppa26$missing_ppa26)
dim (missing_ppa26)

df_mat_ppa_pnc26 <- ppa_pnc26_prep %>%
  filter(ppa_pnc26 >= 0) 
dim (df_mat_ppa_pnc26)

df_mat_ppa_pnc26 <- df_mat_ppa_pnc26 %>%
  filter(!is.na(hb))
dim (df_mat_ppa_pnc26)

remapp_mat_ppa_pnc26 <- outcome_df_rr(
  df = df_mat_ppa_pnc26,
  hb_var = "hb",
  outcome_var = "ppa_pnc26",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_ppa_pnc26"
)
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
  left_join(MAT_PRETERM , by = c("SITE", "MOMID","PREGID")) %>%
  mutate(
    pprom = case_when(
      PPROM_OCCUR %in% c(1,0) ~ PPROM_OCCUR,
      TRUE ~ NA_real_
    )) %>%
  filter(pprom >= 0) 

pprom_prep <- mat_remapp_denom_c %>%
  select(SITE, MOMID, PREGID, PREG_END_DENOM_20) %>%
  left_join(df_hb_exm_anc, by = c("MOMID", "PREGID", "SITE")) %>%
  left_join(MAT_PRETERM, by = c("SITE", "MOMID","PREGID")) %>%
  mutate(
    pprom = case_when(
      PPROM_OCCUR %in% c(1,0) ~ PPROM_OCCUR,
      TRUE ~ -5
    ),
    
    miss_pprom = case_when(
     !is.na(PPROM_OCCUR_MISS) ~ PPROM_OCCUR_MISS,
      TRUE ~ 55L
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (pprom_prep$pprom, pprom_prep$miss_pprom )

missing_pprom <- pprom_prep %>%
  filter(is.na(pprom)) 

df_mat_pprom <- pprom_prep %>%
  filter(pprom >= 0) 
dim (df_mat_pprom)

df_mat_pprom <- df_mat_pprom %>%
  filter(!is.na(hb))
dim (df_mat_pprom)

remapp_mat_pprom <- outcome_df_rr(
  df = df_mat_pprom,
  hb_var = "hb",
  outcome_var = "pprom",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_pprom"
)

df_mat_pprom_trim1 <- df_hb_exm_trim1 %>% 
  left_join(MAT_PRETERM %>% select(SITE, MOMID, PREGID, pprom = PPROM_OCCUR), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(pprom %in% c(0, 1))

df_mat_pprom_trim2 <- df_hb_exm_trim2 %>% 
  left_join(MAT_PRETERM %>% select(SITE, MOMID, PREGID, pprom = PPROM_OCCUR), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(pprom %in% c(0, 1))

df_mat_pprom_trim3 <- df_hb_exm_trim3 %>% 
  left_join(MAT_PRETERM %>% select(SITE, MOMID, PREGID, pprom = PPROM_OCCUR), by = c("SITE", "MOMID","PREGID")) %>% 
  filter(pprom %in% c(0, 1))

df_mat_analysis <- df_mat_analysis %>% 
  left_join(df_mat_pprom %>% select(SITE, MOMID, PREGID, PPROM = pprom),
            by = c("SITE", "MOMID","PREGID"))

##preeclampsia----
df_mat_preclamp <- mat_remapp_denom_b %>%
  select(SITE, MOMID, PREGID, ONE_VISIT_DENOM) %>%
  left_join(df_hb_exm_anc, by = c("MOMID", "PREGID", "SITE")) %>%
  left_join(MAT_HDP_RAW, by = c("SITE", "MOMID","PREGID")) %>%
  mutate(
    preclamp = case_when(
      PREECLAMPSIA == 1 & PREECLAMPSIA_GA_WK >= ga_wks ~ 1L,
      
      # outcome before exposure OR explicitly no preeclampsia
      PREECLAMPSIA == 0 |
      (PREECLAMPSIA == 1 & PREECLAMPSIA_GA_WK < ga_wks) ~ 0L,
      
      TRUE ~ -5L
    ),
    
    miss_preclamp = case_when(
      preclamp %in% c(0L, 1L) ~ 0L,
      preclamp == -5L & HDP_GROUP_MISS %in% 1:3 ~ 1L, #HDP Group missing
      preclamp == -5L & PREEC_CASE_REVIEW == 1 ~ 2L, #Preeclampsia case in-review
      preclamp == -5L & is.na(PREECLAMPSIA_GA_WK) ~ 3L, #missing onset GA
      TRUE ~ 55L
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (df_mat_preclamp$preclamp, df_mat_preclamp$miss_preclamp )

missing_preclampsia <- df_mat_preclamp %>%
  filter(preclamp == -5) %>% 
  filter (PREECLAMPSIA == 1) %>%
  select(PREECLAMPSIA_GA_WK , ga_wks)

df_mat_preclamp <- df_mat_preclamp %>%
  filter(preclamp >= 0) 
dim (df_mat_preclamp)

df_mat_preclamp <- df_mat_preclamp %>%
  filter(!is.na(hb))
dim (df_mat_preclamp)

remapp_mat_preclamp <- outcome_df_rr(
  df = df_mat_preclamp,
  hb_var = "hb",
  outcome_var = "preclamp",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_preclamp"
)
#trimester specific analysis
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
  filter(visit_type %in% 1:5, !is.na(hb)) %>%
  select(SITE, MOMID, PREGID, ga_wks, ga_days_hb, hb, trimester) 

#Basically we want to choose the dpr score measurement closest to the hb measurement
df_hb_anc_dpr <- restruc_mnh25 %>%
  transmute(
    SITE, MOMID, PREGID,
    depress, epds_score,
    TYPE_VISIT_25 = EXPECTED_TYPE_VISIT,
    M25_OBSSTDAT,
    ga_days_25 = as.numeric(ga_days)
  ) %>%
  full_join(
    hb_anc_long %>%
      transmute(
        SITE, MOMID, PREGID, trimester,
        ga_days_hb, ga_wks, hb),
    by = c("SITE","MOMID","PREGID")
  ) %>%
  filter(
    is.finite(hb), hb > 0,
    !is.na(ga_days_hb), !is.na(ga_days_25),
    abs(ga_days_hb - ga_days_25) <= 7
  ) %>%
  mutate(ga_day_diff = abs(ga_days_hb - ga_days_25)) %>%
  group_by(SITE, MOMID, PREGID, TYPE_VISIT_25, M25_OBSSTDAT, ga_days_25) %>%
  slice_min(order_by = ga_day_diff, n = 1, with_ties = FALSE) %>%
  ungroup()


#then or the analysis we want one hb per mom
med_hb_anc_dpr <- df_hb_anc_dpr %>%
  filter(is.finite(hb)) %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")


df_hb_anc_dpr2 <- df_hb_anc_dpr %>%
  inner_join(med_hb_anc_dpr, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  ungroup() %>%
  select(SITE, MOMID, PREGID, epds_score, depress, ga_wks ,trimester ,starts_with("hb", ignore.case = TRUE)) %>% 
  filter (PREGID %in% df_maternal$PREGID)

#we want to get the missing IDs and actually have the right denominator
df_mat_anc_dpr <- mat_remapp_denom_b %>%
  select(SITE, MOMID, PREGID, starts_with("ELIGIBILITY_"), ONE_VISIT_DENOM) %>%
  left_join(df_hb_anc_dpr2, by = c("MOMID", "PREGID", "SITE")) %>%
  mutate(
    miss_dpr = case_when(
      depress %in% c(0L, 1L) ~ 0L,
      is.na(depress) | is.na(hb) | is.na(epds_score) ~ 1L,
      TRUE ~ 1L   # anything weird (e.g., 2, -5) treat as missing
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (df_mat_anc_dpr$miss_dpr)

#why is depression missing?
missing_dpr <- df_mat_anc_dpr %>%
  filter(miss_dpr == 1) %>%
  select(SITE, MOMID, PREGID, ONE_VISIT_DENOM) %>%
  # attach MNH25 timing + EPDS fields
  left_join(
    restruc_mnh25 %>%
      filter(EXPECTED_TYPE_VISIT %in% c(1:5, 13)) %>%
      transmute(
        SITE, MOMID, PREGID,
        ga_days_25 = as.numeric(ga_days),
        q_answered,
        epds_score,
        depress
      ),
    by = c("SITE","MOMID","PREGID")
  ) %>%
  
  # attach Hb timing + value
  left_join(
    hb_anc_long %>%
      transmute(
        SITE, MOMID, PREGID,
        ga_days_hb,
        hb
      ),
    by = c("SITE","MOMID","PREGID")
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    ONE_VISIT_DENOM = first(ONE_VISIT_DENOM),
    
    # Presence checks
    has_mnh25 = any(!is.na(ga_days_25) | !is.na(q_answered) | !is.na(epds_score) | !is.na(depress)),
    has_hb    = any(is.finite(hb) | !is.na(ga_days_hb)),
    
    has_hb_value = any(is.finite(hb)),
    has_hb_ga    = any(!is.na(ga_days_hb)),
    has_25_ga    = any(!is.na(ga_days_25)),
    
    # MNH25 completeness / EPDS (overall, across all MNH25 rows)
    any_mnh25_pre20wk  = any(!is.na(ga_days_25) & ga_days_25 < 20 * 7),
    any_q_answered_lt2 = any(!is.na(q_answered) & q_answered < 2),
    all_epds_na        = all(is.na(epds_score)),
    any_epds_nonmiss   = any(!is.na(epds_score)),
    
    # Timeline match (± 1 week == ±7 days)
    has_pair_within_1wk =
      any(!is.na(ga_days_hb) & !is.na(ga_days_25) & abs(ga_days_hb - ga_days_25) <= 7),
    
    # NEW: EPDS present among the MNH25 rows that actually have an Hb within ±7 days
    any_epds_nonmiss_matched =
      any(
        !is.na(epds_score) &
          !is.na(ga_days_hb) & !is.na(ga_days_25) &
          abs(ga_days_hb - ga_days_25) <= 7
      ),
    
    # Depression value validity (if your depress is 0/1 coded)
    has_depress_value = any(depress %in% c(0L, 1L)),
    
    .groups = "drop"
  ) %>%
  mutate(
    missing_reason = case_when(
      
      # 1) MNH25 collected but too early to be eligible
      has_mnh25 & any_mnh25_pre20wk ~
        "MNH25 collected before 20 weeks (ineligible)",
      
      # 2) MNH25 exists but essentially not completed
      has_mnh25 & any_q_answered_lt2 & all_epds_na ~
        "MNH25 present but incomplete (q_answered < 2; EPDS missing)",
      
      # NEW: EPDS exists elsewhere but not on the Hb-matched MNH25 record
      has_pair_within_1wk & any_epds_nonmiss & !any_epds_nonmiss_matched ~
        "EPDS present at other MNH25 visits, but missing for the Hb-matched MNH25 visit (±7 days)",
      
      # 3) MNH25 present but EPDS missing (everywhere)
      has_mnh25 & all_epds_na ~
        "MNH25 present but EPDS score missing",
      
      # 4) No MNH25 at all
      !has_mnh25 ~
        "No MNH25 record",
      
      # 5) No Hb at all
      !has_hb ~
        "No Hb ANC data",
      
      # 6) Both exist but timing mismatch
      has_mnh25 & has_hb & !has_pair_within_1wk ~
        "Hb and MNH25 exist but no pairing within ±7 days",
      
      # 7) Timing OK but depression still missing
      has_pair_within_1wk & (!has_depress_value|any_q_answered_lt2) ~
        "Timing match exists but depression missing/invalid",
      
      TRUE ~
        "Other / review"
    )
  )
table (missing_dpr$missing_reason)


#anc all dpr dataset
df_mat_anc_dpr <- df_mat_anc_dpr %>%
  filter(depress >= 0 & !is.na(hb)) 
dim (df_mat_anc_dpr)

remapp_mat_anc_dpr <- outcome_df_rr(
  df = df_mat_anc_dpr,
  hb_var = "hb",
  outcome_var = "depress",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_anc_dpr"
)

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
  filter(visit_type == 10 , !is.na(hb)) %>%
  select(SITE, MOMID, PREGID, age_wks, age_days_hb, hb) %>%
  mutate(age_wks = (as.numeric(age_wks)))

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
  full_join(
    df_dpr_pnc,
    by = c("SITE", "MOMID", "PREGID")) %>%
  filter(hb > 0, depress >= 0) %>%
  filter(
    !is.na(age_days_hb), !is.na(age_days_25),
    (abs(age_days_hb - age_days_25) <= 7)) 

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
df_hb_anc_fat <- mnh26_df %>%
  transmute(
    SITE, MOMID, PREGID,
    fatigued,
    fatigue_score,
    TYPE_VISIT_26 = EXPECTED_TYPE_VISIT,
    M26_FTGE_OBSTDAT,
    ga_days_26 = as.numeric(ga_days)
  ) %>%
  full_join(
    hb_anc_long %>%
      transmute(
        SITE, MOMID, PREGID, trimester,
        ga_days_hb, ga_wks, hb),
    by = c("SITE","MOMID","PREGID")
  ) %>%
  # keep only plausible values (and keep in-window pairs)
  filter(
    is.finite(hb), hb > 0,
    !is.na(ga_days_hb), !is.na(ga_days_26),
    abs(ga_days_hb - ga_days_26) <= 7,
    fatigued >= 0,
    fatigue_score >= 0
  ) %>%
  mutate(ga_day_diff = abs(ga_days_hb - ga_days_26)) %>%
  # one Hb match per MNH26 record (closest in GA-days)
  group_by(SITE, MOMID, PREGID, TYPE_VISIT_26, M26_FTGE_OBSTDAT, ga_days_26) %>%
  slice_min(order_by = ga_day_diff, n = 1, with_ties = FALSE) %>%
  ungroup()


#then or the analysis we want one hb per mom
med_hb_anc_fat <- df_hb_anc_fat %>%
  filter(is.finite(hb)) %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")


df_hb_anc_fat2 <- df_hb_anc_fat  %>%
  inner_join(med_hb_anc_fat, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  ungroup() %>%
  select(SITE, MOMID, PREGID, fatigue_score, fatigued, trimester, ga_wks,
         starts_with("hb", ignore.case = TRUE))

#we want to get the missing IDs and actually have the right denominator
df_mat_anc_fat <- mat_remapp_denom_b %>%
  select(SITE, MOMID, PREGID, starts_with("ELIGIBILITY_"), ONE_VISIT_DENOM) %>%
  left_join(df_hb_anc_fat2, by = c("MOMID", "PREGID", "SITE")) %>%
  mutate(
    miss_fat = case_when(
      fatigued %in% c(0L, 1L) ~ 0L,
      is.na(fatigued) | is.na(hb) | is.na(fatigue_score) ~ 1L,
      TRUE ~ 1L   # anything weird (e.g., 2, -5) treat as missing
    )
  ) %>%
  arrange(SITE, MOMID, PREGID)

table (df_mat_anc_fat$miss_fat)

df_mat_anc_fat_c <- df_mat_anc_fat %>% filter (!is.na(hb) & !is.na(fatigue_score))
remapp_mat_anc_fat <- outcome_df_rr(
  df = df_mat_anc_fat_c,
  hb_var = "hb",
  outcome_var = "fatigue_score",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_mat_anc_fat"
)

#why is fatigue missing?
missing_fat <- df_mat_anc_fat %>%
  filter(miss_fat == 1) %>%
  select(SITE, MOMID, PREGID, ONE_VISIT_DENOM) %>%
  # attach MNH26 timing + Fatigue fields
  left_join(
    mnh26_df %>%
      filter(EXPECTED_TYPE_VISIT %in% c(1:5, 13)) %>%
      transmute(
        SITE, MOMID, PREGID,
        ga_days_26 = as.numeric(ga_days),
        q_answered = questionsanswered,
        fatigue_score,
        fatigued
      ),
    by = c("SITE","MOMID","PREGID")
  ) %>%
  
  # attach Hb timing + value
  left_join(
    hb_anc_long %>%
      transmute(
        SITE, MOMID, PREGID,
        ga_days_hb,
        hb
      ),
    by = c("SITE","MOMID","PREGID")
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(
    ONE_VISIT_DENOM = first(ONE_VISIT_DENOM),
    
    # Presence checks
    has_mnh26 = any(!is.na(ga_days_26) | !is.na(q_answered) | !is.na(fatigue_score) | !is.na(fatigued)),
    has_hb    = any(is.finite(hb) | !is.na(ga_days_hb)),
    
    has_hb_value = any(is.finite(hb)),
    has_hb_ga    = any(!is.na(ga_days_hb)),
    has_26_ga    = any(!is.na(ga_days_26)),
    
    # MNH26 completeness / Fatigue (overall, across all MNH26 rows)
    any_mnh26_pre20wk  = any(!is.na(ga_days_26) & ga_days_26 < 20 * 7),
    any_q_answered_lt2 = any(!is.na(q_answered) & q_answered < 2),
    all_epds_na        = all(is.na(fatigue_score)),
    any_epds_nonmiss   = any(!is.na(fatigue_score)),
    
    # Timeline match (± 1 week == ±7 days)
    has_pair_within_1wk =
      any(!is.na(ga_days_hb) & !is.na(ga_days_26) & abs(ga_days_hb - ga_days_26) <= 7),
    
    # NEW: Fatigue present among the MNH26 rows that actually have an Hb within ±7 days
    any_epds_nonmiss_matched =
      any(
        !is.na(fatigue_score) &
          !is.na(ga_days_hb) & !is.na(ga_days_26) &
          abs(ga_days_hb - ga_days_26) <= 7
      ),
    
    # fatigue value validity (if your fatigued is 0/1 coded)
    has_fatigued_value = any(fatigued %in% c(0L, 1L)),
    
    .groups = "drop"
  ) %>%
  mutate(
    missing_reason = case_when(
      
      # 1) MNH26 collected but too early to be eligible
      has_mnh26 & any_mnh26_pre20wk ~
        "MNH26 collected before 20 weeks (ineligible)",
      
      # 2) MNH26 exists but essentially not completed
      has_mnh26 & any_q_answered_lt2 & all_epds_na ~
        "MNH26 present but incomplete (q_answered < 2; Fatigue missing)",
      
      # NEW: Fatigue exists elsewhere but not on the Hb-matched MNH26 record
      has_pair_within_1wk & any_epds_nonmiss & !any_epds_nonmiss_matched ~
        "Fatigue present at other MNH26 visits, but missing for the Hb-matched MNH26 visit (±7 days)",
      
      # 3) MNH26 present but Fatigue missing (everywhere)
      has_mnh26 & all_epds_na ~
        "MNH26 present but Fatigue score missing",
      
      # 4) No MNH26 at all
      !has_mnh26 ~
        "No MNH26 record",
      
      # 5) No Hb at all
      !has_hb ~
        "No Hb ANC data",
      
      # 6) Both exist but timing mismatch
      has_mnh26 & has_hb & !has_pair_within_1wk ~
        "Hb and MNH26 exist but no pairing within ±7 days",
      
      # 7) Timing OK but fatigue still missing
      has_pair_within_1wk & (!has_fatigued_value|any_q_answered_lt2) ~
        "Timing match exists but fatigue missing/invalid",
      
      TRUE ~
        "Other / review"
    )
  )
table (missing_fat$missing_reason)


#anc all fat dataset
df_mat_anc_fat <- df_mat_anc_fat %>%
  filter(fatigued >= 0 & !is.na(hb)) 
dim (df_mat_anc_fat)
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
  full_join(
    df_fat_pnc,
    by = c("SITE", "MOMID", "PREGID")) %>%
  filter(hb > 0, fatigued >= 0) %>%
  filter(
    !is.na(age_days_hb), !is.na(age_days_26),
    (abs(age_days_hb - age_days_26) <= 7)) 

#then or the analysis we want one hb per mom
med_hb_pnc_fat <- df_hb_pnc_fat %>%
  group_by(SITE) %>%
  summarise(site_hb_median = median(hb, na.rm = TRUE), .groups = "drop")

df_mat_pnc_fat_pnc6 <- df_hb_pnc_fat %>%
  inner_join(med_hb_pnc_fat, by = "SITE") %>%
  mutate(
    signed_dev = hb - site_hb_median,
    abs_dev    = abs(signed_dev)
  ) %>%
  group_by(SITE, MOMID, PREGID) %>%
  slice_max(order_by = abs_dev, n = 1, with_ties = FALSE) %>%  # pick the single farthest record
  ungroup() %>%
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
filter (PREGID %in% mat_remapp_denom_c$PREGID) %>% 
left_join(                                      # be explicit: we only keep rows with an hb match
    df_hb_exm_anc, by = c("SITE", "MOMID", "PREGID")
  ) 

##compo data (preterm37, lbw2500, sga10)----
prep_compo <- inf_remapp_denom_h %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant, by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
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
      preterm37 == 0 | lbw2500 == 0 | sga10 == 0 ~ 0, 
      TRUE ~ NA_real_
    ),
    miss_inf_comp = case_when(
      compo_pre_lbw_sga %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1 
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table (prep_compo$miss_inf_comp, prep_compo$compo_pre_lbw_sga )

miss_inf_comp <- prep_compo %>%
  filter(miss_inf_comp == 1) 

table (miss_inf_comp$miss_inf_comp )
dim (miss_inf_comp)

df_inf_compo <- prep_compo %>% 
  filter(compo_pre_lbw_sga >= 0 & hb > 0)
dim (df_inf_compo)

remapp_df_inf_composite <- outcome_df_rr(
  df = df_inf_compo,
  hb_var = "hb",
  outcome_var = "compo_pre_lbw_sga",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_composite"
)

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
df_inf_preterm37_prep <- prep_compo %>% 
  mutate(
  miss_preterm37 = case_when(
    preterm37 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
    TRUE ~ 1)) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table (df_inf_preterm37_prep$preterm37, df_inf_preterm37_prep$miss_preterm37 )

miss_preterm37 <- df_inf_preterm37_prep %>%
  filter(miss_preterm37 == 1) 

table (miss_preterm37$miss_preterm37 )
dim (miss_preterm37)

df_inf_preterm37 <- df_inf_preterm37_prep %>% 
  filter(preterm37 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, preterm37, starts_with("hb"), trimester, ga_wks)
dim (df_inf_preterm37)

remapp_df_inf_preterm37 <- outcome_df_rr(
  df = df_inf_preterm37,
  hb_var = "hb",
  outcome_var = "preterm37",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_preterm37"
)

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
df_inf_lbw2500_prep <- prep_compo %>% 
  mutate(
    miss_lbw2500 = case_when(
      lbw2500 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_lbw2500_prep$lbw2500, df_inf_lbw2500_prep$miss_lbw2500)

miss_lbw2500 <- df_inf_lbw2500_prep %>%
  filter(miss_lbw2500 == 1) %>%
  mutate(missing_lbw2500_rsn = 
           case_when(MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 1, #missing form
                     MISSING_BOTH == 1 ~ 2, #missing birthweight
                     MISSING_TIME == 1 ~ 3, #missing birth time
                     TRUE ~ 55)) %>%
  mutate(missing_lbw2500_label = case_when(
    missing_lbw2500_rsn == 1 ~ "Missing birth outcome form",
    missing_lbw2500_rsn == 2 ~ "Missing birthweight",
    missing_lbw2500_rsn == 3 ~ "Missing time of birth",
    TRUE ~ "Other / unclassified"
  ))

table(miss_lbw2500$missing_lbw2500_rsn)
dim(miss_lbw2500)

df_inf_lbw2500 <- df_inf_lbw2500_prep %>% 
  filter(lbw2500 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, lbw2500, starts_with("hb"),trimester, ga_wks)
dim (df_inf_lbw2500)

remapp_df_inf_lbw2500 <- outcome_df_rr(
  df = df_inf_lbw2500,
  hb_var = "hb",
  outcome_var = "lbw2500",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_lbw2500"
)


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

df_inf_sga10_prep <- prep_compo %>% 
  mutate(
    miss_sga10 = case_when(
      sga10 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_sga10_prep$sga10, df_inf_sga10_prep$miss_sga10)

miss_sga10 <- df_inf_sga10_prep %>%
  filter(miss_sga10 == 1) %>%
  mutate(missing_sga10_rsn = 
           case_when(MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 1, #missing form
                     MISSING_BOTH == 1 ~ 2, #missing birthweight
                     MISSING_TIME == 1 ~ 3, #missing birth time
                     !(SEX %in% (1:2)) | is.na (SEX) ~ 4, #missing sex
                     is.na (SGA_CENTILE) ~ 5, #missing centile
                     TRUE ~ 55)) %>%
  mutate(missing_sga10_label = case_when(
    missing_sga10_rsn == 1 ~ "Missing birth outcome form",
    missing_sga10_rsn == 2 ~ "Missing birthweight",
    missing_sga10_rsn == 3 ~ "Missing time of birth",
    missing_sga10_rsn == 4 ~ "Missing/invalid infant sex",
    missing_sga10_rsn == 5 ~ "Missing SGA centile",
    TRUE ~ "Other / unclassified"
  ))

table(miss_sga10$missing_sga10_rsn)
dim(miss_sga10)

df_inf_sga10 <- df_inf_sga10_prep %>% 
  filter(sga10 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, sga10, starts_with("hb"),trimester,ga_wks)
dim (df_inf_sga10)

remapp_df_inf_sga10 <- outcome_df_rr(
  df = df_inf_sga10,
  hb_var = "hb",
  outcome_var = "sga10",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_sga10"
)

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
df_inf_preterm34_prep <- prep_compo %>% 
  mutate(
    preterm34 = case_when(
      PRETERMBIRTH_CAT %in% c(13, 14, 15) ~ 1,
      PRETERMBIRTH_CAT %in% c(10, 11, 12) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    miss_preterm34 = case_when(
      preterm34 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_preterm34_prep$preterm34, df_inf_preterm34_prep$miss_preterm34)

miss_preterm34 <- df_inf_preterm34_prep %>%
  filter(miss_preterm34 == 1)

table(miss_preterm34$miss_preterm34)
dim(miss_preterm34)

df_inf_preterm34 <- df_inf_preterm34_prep %>% 
  filter(preterm34 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, preterm34, starts_with("hb"), trimester, ga_wks)
dim (df_inf_preterm34)

remapp_df_inf_preterm34 <- outcome_df_rr(
  df = df_inf_preterm34,
  hb_var = "hb",
  outcome_var = "preterm34",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_preterm34"
)


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

df_inf_lbw1500_prep <- prep_compo %>% 
  mutate(lbw1500 = case_when(
    LBW_CAT_ANY %in% c(11) ~ 1,
    LBW_CAT_ANY %in% c(12,13,14) ~ 0, 
    TRUE ~ NA_real_
  )) %>% 
  mutate(
    miss_lbw1500 = case_when(
      lbw1500 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_lbw1500_prep$lbw1500, df_inf_lbw1500_prep$miss_lbw1500)

miss_lbw1500 <- df_inf_lbw1500_prep %>%
  filter(miss_lbw1500 == 1) %>%
  mutate(missing_lbw1500_rsn = 
           case_when(MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 1, #missing form
                     MISSING_BOTH == 1 ~ 2, #missing birthweight
                     MISSING_TIME == 1 ~ 3, #missing birth time
                     TRUE ~ 55)) %>%
  mutate(missing_lbw1500_label = case_when(
    missing_lbw1500_rsn == 1 ~ "Missing birth outcome form",
    missing_lbw1500_rsn == 2 ~ "Missing birthweight",
    missing_lbw1500_rsn == 3 ~ "Missing time of birth",
    TRUE ~ "Other / unclassified"
  ))

table(miss_lbw1500$missing_lbw1500_rsn)
dim(miss_lbw1500)

df_inf_lbw1500 <- df_inf_lbw1500_prep %>% 
  filter(lbw1500 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, lbw1500, starts_with("hb"), trimester, ga_wks)
dim (df_inf_lbw1500)

remapp_df_inf_lbw1500 <- outcome_df_rr(
  df = df_inf_lbw1500,
  hb_var = "hb",
  outcome_var = "lbw1500",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_lbw1500"
)


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

df_inf_sga3_prep <- prep_compo %>% 
  mutate(
    sga3 = case_when(
      SGA_CAT %in% c(11) ~ 1,
      SGA_CAT %in% c(12,13,14) ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>% 
  mutate(
    miss_sga3 = case_when(
      sga3 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_sga3_prep$sga3, df_inf_sga3_prep$miss_sga3)

miss_sga3 <- df_inf_sga3_prep %>%
  filter(miss_sga3 == 1) %>%
  mutate(missing_sga3_rsn = 
           case_when(MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 1, #missing form
                     MISSING_BOTH == 1 ~ 2, #missing birthweight
                     MISSING_TIME == 1 ~ 3, #missing birth time
                     !(SEX %in% (1:2)) | is.na (SEX) ~ 4, #missing sex
                     is.na (SGA_CENTILE) ~ 5, #missing centile
                     TRUE ~ 55)) %>%
  mutate(missing_sga3_label = case_when(
    missing_sga3_rsn == 1 ~ "Missing birth outcome form",
    missing_sga3_rsn == 2 ~ "Missing birthweight",
    missing_sga3_rsn == 3 ~ "Missing time of birth",
    missing_sga3_rsn == 4 ~ "Missing/invalid infant sex",
    missing_sga3_rsn == 5 ~ "Missing SGA centile",
    TRUE ~ "Other / unclassified"
  )) 

table(miss_sga3$missing_sga3_rsn)
dim(miss_sga3)

df_inf_sga3 <- df_inf_sga3_prep %>% 
  filter(sga3 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, sga3, starts_with("hb"), trimester, ga_wks)
dim (df_inf_sga3)

remapp_df_inf_sga3 <- outcome_df_rr(
  df = df_inf_sga3,
  hb_var = "hb",
  outcome_var = "sga3",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID","INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_sga3"
)


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

df_inf_psbi_prep <- inf_remapp_denom_h %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant, by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
  mutate(
    age_days_today = as.integer(ymd(UploadDate) - ymd(DOB)),
    
    inf_psbi = case_when(
      if_any(
        c(INF_PSBI_IPC, INF_PSBI_PNC0, INF_PSBI_PNC1, INF_PSBI_HOSPITAL,
          INF_PSBI_PNC4, INF_PSBI_PNC6, INF_PSBI_UNSCHED),
        ~ .x == 1
      ) ~ 1,
      
      age_days_today >= 59 & if_any(
        c(INF_PSBI_IPC, INF_PSBI_PNC0, INF_PSBI_PNC1, INF_PSBI_HOSPITAL,
          INF_PSBI_PNC4, INF_PSBI_PNC6, INF_PSBI_UNSCHED),
        ~ .x == 0
      ) ~ 0,
      
      TRUE ~ NA_real_)) %>% 
  mutate(
    miss_psbi = case_when(
      inf_psbi %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_psbi_prep$inf_psbi, df_inf_psbi_prep$miss_psbi)

miss_psbi <- df_inf_psbi_prep %>%
  filter(miss_psbi == 1) %>%
  mutate(missing_psbi_rsn = 
           case_when(age_days_today < 59 ~ 1, #participant hasn't passed risk period <59days
                     TRUE ~ 55)) 

table(miss_psbi$missing_psbi_rsn)
dim(miss_psbi)

df_inf_psbi <- df_inf_psbi_prep %>% 
  filter(inf_psbi >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, inf_psbi, starts_with("hb"), trimester, ga_wks)
dim (df_inf_psbi)

remapp_df_inf_psbi <- outcome_df_rr(
  df = df_inf_psbi,
  hb_var = "hb",
  outcome_var = "inf_psbi",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID", "INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_psbi"
)


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
df_inf_asph_prep <- inf_remapp_denom_f %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant, by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
  mutate(inf_asph = ifelse(INF_ASPH %in% c(0,1), INF_ASPH, NA_real_)) %>% 
  mutate(
    miss_inf_asph = case_when(
      inf_asph %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_asph_prep$inf_asph, df_inf_asph_prep$miss_inf_asph)

miss_inf_asph <- df_inf_asph_prep %>%
  filter(miss_inf_asph == 1) %>% 
  mutate(
    miss_inf_asph_rsn = case_when(
      STILLBIRTH_TIMING %in% c(99,55) ~ 1, #stillbirth timing information is missing
      STILLBIRTH_SIGNS_LIFE %in% c(55) ~ 2, #signs of life information is missing
      MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 3, #missing form
      is.na(INFANTID) ~ 4,  #infantid not allocated 
      TRUE ~ 55
    )
  ) 

table(miss_inf_asph$miss_inf_asph_rsn)
dim(miss_inf_asph)

df_inf_asph <- df_inf_asph_prep %>% 
  filter(inf_asph >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, inf_asph, starts_with("hb"), trimester, ga_wks)
dim (df_inf_asph)

remapp_df_inf_asph <- outcome_df_rr(
  df = df_inf_asph,
  hb_var = "hb",
  outcome_var = "inf_asph",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID", "INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_asph"
)


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
  filter(ga_days_hb > 0 & ga_days_hb <= 140) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth20 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth20 <- df_hb_long2 %>% 
  filter(ga_days_hb > 0 & ga_days_hb <= 140) %>%
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
         trimester, ga_wks) %>% 
  mutate (hb = hb_exm) 

save(df_hb_stbirth20, file = "derived_data/df_hb_stbirth20.rda")



###median Hb Value for Stillbirth 20weeks Trimester 1 Visits ----
#*****calculate median hb value for Stillbirth Trimester 1 visits per SITE******
med_hb_stbirth20_trim1 <- df_hb_long2 %>% 
  filter(ga_days_hb > 0 & ga_days_hb < 98) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth20_trim1 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth20_trim1 <- df_hb_long2 %>% 
  filter(ga_days_hb > 0 & ga_days_hb < 98) %>%
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
         trimester, ga_wks) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth20_trim1, file = "derived_data/df_hb_stbirth20_trim1.rda")

table(df_hb_stbirth20_trim1$trimester)

###median Hb Value for Stillbirth 20weeks Trimester 2 Visits ----
#*****calculate median hb value for Stillbirth Trimester 2 visits per SITE******
med_hb_stbirth20_trim2 <- df_hb_long2 %>% 
  filter(ga_days_hb >= 98 & ga_days_hb < 140) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth20_trim1 with one hb_exm value per mom ************
set.seed(100)
df_hb_stbirth20_trim2 <- df_hb_long2 %>% 
  filter(ga_days_hb >= 98 & ga_days_hb < 140) %>%
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
         trimester, ga_wks) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth20_trim2, file = "derived_data/df_hb_stbirth20_trim2.rda")

table(df_hb_stbirth20_trim2$trimester)



#Stillbirth 20weeks outcome dataset generation
df_inf_stillbirth20_prep <- inf_remapp_denom_f %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant %>% 
              select (-(starts_with("hb")), -trimester, -ga_wks), 
            by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
  left_join(       # be explicit: we only keep rows with an hb match
    df_hb_stbirth20, by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  mutate(inf_stillbirth20 = case_when(
    STILLBIRTH_20WK == 1 ~ 1,
    STILLBIRTH_20WK == 0 ~ 0, 
    TRUE ~ NA_real_
  )) %>%
  mutate(
    miss_inf_stillbirth20 = case_when(
      inf_stillbirth20 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1)) 

miss_inf_stillbirth20 <- df_inf_stillbirth20_prep %>%
  filter(miss_inf_stillbirth20 == 1) %>%
  mutate(
    miss_inf_stillbirth20_rsn = case_when(
      is.na(hb) ~ 5,  # Missing <20weeks hb data
      STILLBIRTH_TIMING %in% c(99, 55) ~ 1,  # stillbirth timing info missing
      STILLBIRTH_SIGNS_LIFE %in% c(55) ~ 2,  # signs of life info missing
      MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 3,  # missing form
      is.na(INFANTID) ~ 4,  # infantid not allocated
      TRUE ~ 55
    )
  )

table(miss_inf_stillbirth20$miss_inf_stillbirth20_rsn)
dim(miss_inf_stillbirth20)

df_inf_stillbirth20 <- df_inf_stillbirth20_prep %>% 
  filter(inf_stillbirth20 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, inf_stillbirth20, starts_with("hb"), trimester, ga_wks)
dim (df_inf_stillbirth20)

remapp_df_inf_stillbirth20 <- outcome_df_rr(
  df = df_inf_stillbirth20,
  hb_var = "hb",
  outcome_var = "inf_stillbirth20",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID", "INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_stillbirth20"
)

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
  filter(ga_days_hb > 0 & ga_days_hb <= 196) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth28 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth28 <- df_hb_long2 %>% 
  filter(ga_days_hb > 0 & ga_days_hb <= 196) %>%
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
         trimester, ga_wks) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth28, file = "derived_data/df_hb_stbirth28.rda")

table(df_hb_stbirth28$trimester)

###median Hb Value for Stillbirth 28weeks Trimester 1 Visits ----
#*calculate median hb value for Stillbirth Trimester 1 visits per SITE
med_hb_stbirth28_trim1 <- df_hb_long2 %>% 
  filter(ga_days_hb > 0 & ga_days_hb <= 98) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth28_trim1 with one hb_exm value per mom ---- ************
set.seed(100)
df_hb_stbirth28_trim1 <- df_hb_long2 %>% 
  filter(ga_days_hb > 0 & ga_days_hb <= 98) %>%
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
         trimester, ga_wks) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth28_trim1, file = "derived_data/df_hb_stbirth28_trim1.rda")

table(df_hb_stbirth28_trim1$trimester)

###median Hb Value for Stillbirth 28weeks Trimester 2 Visits ----
#*****calculate median hb value for Stillbirth Trimester 2 visits per SITE******
med_hb_stbirth28_trim2 <- df_hb_long2 %>% 
  filter(ga_days_hb >= 98 & ga_days_hb < 196) %>%
  group_by(SITE, MOMID, PREGID) %>%
  summarise(hb_mom_med = median(hb, na.rm = TRUE), .groups = "drop") %>%
  group_by(SITE) %>%
  summarise(hb_med = median(hb_mom_med, na.rm = TRUE), .groups = "drop")

#*********df_hb_stbirth28_trim1 with one hb_exm value per mom ************
set.seed(100)
df_hb_stbirth28_trim2 <- df_hb_long2 %>% 
  filter(ga_days_hb >= 98 & ga_days_hb < 196) %>%
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
         trimester, ga_wks) %>% 
  mutate (hb = hb_exm)

save(df_hb_stbirth28_trim2, file = "derived_data/df_hb_stbirth28_trim2.rda")

table(df_hb_stbirth28_trim2$trimester)


#Stillbirth 28weeks outcome dataset generation
df_inf_stillbirth28_prep <- inf_remapp_denom_g %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant %>% 
              select (-(starts_with("hb")), -trimester, -ga_wks), 
            by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
  left_join(       # be explicit: we only keep rows with an hb match
    df_hb_stbirth28, by = c("SITE", "MOMID", "PREGID")
  ) %>% 
  mutate(inf_stillbirth28 = case_when(
    STILLBIRTH_28WK == 1 ~ 1,
    STILLBIRTH_28WK == 0 ~ 0, 
    TRUE ~ NA_real_
  )) %>%
  mutate(
    miss_inf_stillbirth28 = case_when(
      inf_stillbirth28 %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) 

miss_inf_stillbirth28 <- df_inf_stillbirth28_prep %>%
  filter(miss_inf_stillbirth28 == 1) %>%
  mutate(
    miss_inf_stillbirth28_rsn = case_when(
      is.na(hb) ~ 5,  # Missing <28weeks hb data
      STILLBIRTH_TIMING %in% c(99, 55) ~ 1,  # stillbirth timing info missing
      STILLBIRTH_SIGNS_LIFE %in% c(55) ~ 2,  # signs of life info missing
      MISSING_MNH09 == 1 | MISSING_MNH11 == 1 ~ 3,  # missing form
      is.na(INFANTID) ~ 4,  # infantid not allocated
      TRUE ~ 55
    )
  )

table(miss_inf_stillbirth28$miss_inf_stillbirth28_rsn)
dim(miss_inf_stillbirth28)

df_inf_stillbirth28 <- df_inf_stillbirth28_prep %>% 
  filter(inf_stillbirth28 >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, inf_stillbirth28, starts_with("hb"),
         trimester, ga_wks)
dim (df_inf_stillbirth28)
table(df_inf_stillbirth28$inf_stillbirth28)

remapp_df_inf_stillbirth28 <- outcome_df_rr(
  df = df_inf_stillbirth28,
  hb_var = "hb",
  outcome_var = "inf_stillbirth28",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID", "INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_stillbirth28"
)

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
df_inf_hyperbili_prep <- inf_remapp_denom_h %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant, by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
  mutate(
    age_days_today = as.integer(ymd(UploadDate) - ymd(DOB)),
    hyperbili = case_when(
    INF_HYPERBILI_AAP_24HR == 1 | INF_HYPERBILI_AAP_5DAY == 1 | INF_HYPERBILI_AAP_14DAY == 1 ~ 1,
    (age_days_today >= 14) & (INF_HYPERBILI_AAP_24HR == 0 | INF_HYPERBILI_AAP_5DAY == 0 | INF_HYPERBILI_AAP_14DAY == 0) ~ 0,
    TRUE ~ NA_real_
  )) %>% 
  mutate(
    miss_hyperbili = case_when(
      hyperbili %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_hyperbili_prep$hyperbili, df_inf_hyperbili_prep$miss_hyperbili)

miss_hyperbili <- df_inf_hyperbili_prep %>%
  filter(miss_hyperbili == 1) %>%
  mutate(missing_hyperbili_rsn = 
           case_when(age_days_today < 14 ~ 1, #participant hasn't passed risk period <14days
                     DENOM_HYPERBILI_ANY == 0 ~ 2, #no tcb measurement
                     TRUE ~ 55)) 

table(miss_hyperbili$missing_hyperbili_rsn)
dim(miss_hyperbili)

df_inf_hyperbili <- df_inf_hyperbili_prep %>% 
  filter(hyperbili >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, hyperbili, starts_with("hb"), trimester, ga_wks)
dim (df_inf_hyperbili)

remapp_df_inf_hyperbili <- outcome_df_rr(
  df = df_inf_hyperbili,
  hb_var = "hb",
  outcome_var = "hyperbili",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID", "INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_hyperbili"
)

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

df_inf_mortality_prep <- inf_remapp_denom_h %>%
  select(SITE, MOMID, PREGID, INFANTID) %>%
  left_join(df_infant%>%
              select(SITE, MOMID, PREGID, INFANTID, DOB, hb, hb_exm,
                     hb_max, hb_min, trimester, ga_wks), 
            by = c("SITE", "MOMID", "PREGID", "INFANTID"))  %>% 
  left_join(
    dplyr::select(inf_mortality, SITE, MOMID, PREGID, INFANTID, NEONATAL_DTH),
    by = c("SITE", "MOMID", "PREGID", "INFANTID")
  ) %>%
  mutate(
    age_days_today = as.integer(ymd(UploadDate) - ymd(DOB)),
    neo_mortality = case_when(
      NEONATAL_DTH == 1 ~ 1,
      age_days_today < 28 ~ NA_real_,
      TRUE ~ 0
    )) %>% 
  mutate(
    miss_neo_mortality = case_when(
      neo_mortality %in% c(0, 1) & !is.na(hb) & hb > 0 ~ 0,
      TRUE ~ 1
    )
  ) %>%
  arrange(SITE, MOMID, PREGID, INFANTID)

table(df_inf_mortality_prep$neo_mortality, df_inf_mortality_prep$miss_neo_mortality)

miss_neo_mortality <- df_inf_mortality_prep %>%
  filter(miss_neo_mortality == 1) %>%
  mutate(missing_neo_mortality_rsn = 
           case_when(age_days_today < 28 ~ 1, #participant hasn't passed risk period <28days
                     TRUE ~ 55)) 

table(miss_neo_mortality$missing_neo_mortality_rsn)
dim(miss_neo_mortality)

df_inf_mortality <- df_inf_mortality_prep %>% 
  filter(neo_mortality >= 0 & hb > 0) %>% 
  select(SITE, MOMID, PREGID, INFANTID, neo_mortality, starts_with("hb"), trimester, ga_wks)
dim (df_inf_mortality)

remapp_df_inf_mortality <- outcome_df_rr(
  df = df_inf_mortality,
  hb_var = "hb",
  outcome_var = "neo_mortality",
  trimester_var = "trimester",
  ga_var = "ga_wks",
  id_vars = c("SITE", "MOMID", "PREGID", "INFANTID"),
  out_rda_dir = derived_rr_dir,
  out_xlsx_dir = paste0("Z:/Precious_working_files/ReMAPP datasets/", UploadDate),
  file_name = "remapp_df_inf_mortality"
)

df_inf_mortality_trim1 <- df_inf_mortality %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_mortality_trim2 <- df_inf_mortality %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2, by = c("SITE", "MOMID","PREGID")) %>% 
  filter(hb > 0)

df_inf_mortality_trim3 <- df_inf_mortality %>% 
  select(-starts_with("hb")) %>% 
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
    POC_HB = hb_poc,
    SPHB = sphb
  ) %>%
  pivot_longer(c("CBC_HB", "POC_HB", "SPHB"),
               names_to = c("hb_type"),
               values_to = "hb_value") %>%
  #redefine ga_wks by using different ga_wks for different hb type
  mutate(ga_wks = case_when(
    hb_type == "CBC_HB" ~ ga_wks_hb,
    hb_type %in% c("POC_HB", "SPHB") ~ round(ga_wks_hb_poc,1),
    TRUE ~ NA_real_
  )) %>% 
  #redefine visit_type by using different visit_type for different hb type
  mutate(visit_type = case_when(
    hb_type == "CBC_HB" ~ visit_type_hb,
    hb_type %in% c("POC_HB", "SPHB") ~ visit_type_hb_poc,
    TRUE ~ NA_real_
  )) %>% 
  select(MOMID, PREGID, SITE, hb_type, hb_value, visit_type, ga_wks)


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
save(INF_OUTCOMES, file = "derived_data/INF_OUTCOMES.rda")

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