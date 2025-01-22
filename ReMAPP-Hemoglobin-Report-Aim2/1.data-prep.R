#*****************************************************************************
#Hemoglobin report - ReMAPP Aim2
#Author: Xiaoyan
#Email: xyh@gwu.edu
#*****************************************************************************
library(tidyverse)
library(lubridate)
library(naniar)
library(haven)

UploadDate = "2024-06-28"

#*****************************************************************************
#*1. Load and merge data
#*****************************************************************************
#load MAT_ENROLL 
MAT_ENROLL <- read.csv(paste0("Z:/Outcome Data/",UploadDate,"/MAT_ENROLL.csv"))
save(MAT_ENROLL, file = "derived_data/MAT_ENROLL.rda")

#load INF_OUTCOMES 
INF_OUTCOMES <- read.csv(paste0("Z:/Outcome Data/",UploadDate,"/INF_OUTCOMES.csv")) 
save(INF_OUTCOMES, file = "derived_data/INF_OUTCOMES.rda")

#load MAT_anmeia
MAT_ANEMIA <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_ANEMIA.dta")) %>% 
  select(SITE, MOMID, PREGID, ANEMIA_T1, ANEMIA_T2, ANEMIA_T3,   
         ANEMIA_ANC, ANEMIA_PNC6, ANEMIA_PNC26)
save(MAT_ANEMIA, file = "derived_data/MAT_ANEMIA.rda")

#load MAT_ENDPOINT --> PREG_END
MAT_ENDPOINTS <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_ENDPOINTS.dta")) %>% 
  select(SITE, MOMID, PREGID, PREG_END)
save(MAT_ENDPOINTS, file = "derived_data/MAT_ENDPOINTS.rda")

#mnh03 for smoking adjustment for hb
mnh03 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh03_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, M03_SMOKE_OECOCCUR)

#mnh06
mnh06 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh06_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, 
         M06_TYPE_VISIT, M06_HB_POC_LBORRES, M06_SPHB_LBORRES, M06_DIAG_VSDAT) %>% 
  filter(M06_TYPE_VISIT < 13) %>% 
  group_by(SITE, MOMID, PREGID, M06_TYPE_VISIT) %>% 
  mutate(n=n()) %>% 
  filter(n == 1) %>% 
  ungroup() %>% 
  select(-n) %>% 
  pivot_wider(names_from = M06_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID))

#mnh08 data
mnh08 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh08_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, 
         M08_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES) %>% 
  filter(M08_TYPE_VISIT < 13) %>% 
  group_by(SITE, MOMID, PREGID, M08_TYPE_VISIT) %>% 
  mutate(n=n()) %>% 
  filter(n == 1) %>% 
  ungroup() %>% 
  select(-n) %>% 
  pivot_wider(names_from = M08_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID)) 

#mnh09
mnh09 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh09_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, M09_MAT_LD_OHOSTDAT,
         matches("_INF1"), matches("_INF2"),
         matches("_INF3"), matches("_INF4")) 

df_maternal <- MAT_ENROLL %>% 
  mutate(remapp = case_when(
    (SITE == "Ghana" & M02_SCRN_OBSSTDAT >= "2022-12-28") |
      (SITE == "Kenya" & M02_SCRN_OBSSTDAT >= "2023-04-14") |
      (SITE == "Zambia" & M02_SCRN_OBSSTDAT >= "2022-12-15") |
      (SITE == "Pakistan" & M02_SCRN_OBSSTDAT >= "2022-09-22" & M02_SCRN_OBSSTDAT <= "2024-04-05") |
      (SITE == "India-CMC" & M02_SCRN_OBSSTDAT >= "2023-06-20") |
      (SITE == "India-SAS" & M02_SCRN_OBSSTDAT >= "2023-08-15") ~ 1, 
    TRUE ~ NA_real_)) %>% 
  filter(remapp == 1) %>% 
  left_join(MAT_ANEMIA) %>% 
  left_join(mnh03, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh06, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh08, by = c("SITE", "MOMID", "PREGID")) %>% 
  left_join(mnh09 %>% select(SITE, MOMID, PREGID, M09_MAT_LD_OHOSTDAT), by = c("SITE", "MOMID", "PREGID"))

#save data
save(df_maternal, file = "derived_data/df_maternal.rda")

#*****************************************************************************
#*2. hb data - long, wide and other datasets
#*****************************************************************************
#*prepare data
prep_hb <- df_maternal %>% 
  dplyr:: select("SCRNID", "MOMID", "PREGID", "SITE",
                 EST_CONCEP_DATE,
                 M03_SMOKE_OECOCCUR,
                 num_range("M06_HB_POC_LBORRES_", 1:12),
                 num_range("M06_SPHB_LBORRES_", 1:12),
                 num_range("M06_DIAG_VSDAT_",1:12),
                 num_range("M08_CBC_HB_LBORRES_",1:12),
                 num_range("M08_LBSTDAT_",1:12),
                 M09_MAT_LD_OHOSTDAT, 
                 ANEMIA_T1, ANEMIA_T2, ANEMIA_T3, 
                 ANEMIA_ANC, ANEMIA_PNC6, ANEMIA_PNC26
  ) %>%
  #replace 7s and 5s with NA
  replace_with_na_all(condition = ~.< 0) %>%
  replace_with_na_all(condition = ~.== "1907-07-07") %>% 
  replace_with_na_all(condition = ~.== "1905-05-05") %>% 
  distinct()

save(prep_hb, file = "derived_data/prep_hb.rda")

#long hb data - basic hb data (no filter)
df_hb_long <- prep_hb %>% 
  #to long format
  pivot_longer(
    -c("SCRNID","MOMID","PREGID","SITE", EST_CONCEP_DATE, M03_SMOKE_OECOCCUR, 
       M09_MAT_LD_OHOSTDAT, ANEMIA_T1, ANEMIA_T2, ANEMIA_T3,  
       ANEMIA_ANC, ANEMIA_PNC6, ANEMIA_PNC26
       ),
    names_to = c(".value", "visit_type"), 
    names_pattern = "^M\\d{2}_(.+)_(\\d+)"
  ) %>% 
  mutate(
    visit_type = as.numeric(visit_type),
    #ga_wks is calculated from MNH08 for general use of hemoglobin (cbc hb)
    ga_wks = case_when(
      visit_type >= 6 ~ NA_real_,
      visit_type < 6 ~ as.numeric(ymd(LBSTDAT) - ymd(EST_CONCEP_DATE))/7
    ),
    #ga_wks_06 is calculated from MNH06 for pochb and sphb
    ga_wks_06 = case_when(
      visit_type >= 6 ~ NA_real_,
      visit_type < 6 ~ as.numeric(ymd(DIAG_VSDAT) - ymd(EST_CONCEP_DATE))/7
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
    hb_alti = case_when(
      SITE == "Kenya" ~ CBC_HB_LBORRES - 0.8,
      SITE == "Zambia" ~ CBC_HB_LBORRES - 0.8,
      !is.na(CBC_HB_LBORRES) ~ CBC_HB_LBORRES,
      TRUE ~ NA_real_
    ), 
    #adjust for both smoke and altitude
    hb = case_when(
      M03_SMOKE_OECOCCUR == 1 ~ hb_alti - 0.3,
      M03_SMOKE_OECOCCUR == 0 ~ hb_alti,
      TRUE ~ NA_real_
    ),
    poc_hb_alti = case_when(
      SITE == "Kenya" ~ HB_POC_LBORRES - 0.8,
      SITE == "Zambia" ~ HB_POC_LBORRES - 0.8,
      !is.na(HB_POC_LBORRES) ~ HB_POC_LBORRES,
      TRUE ~ NA_real_
    ), 
    #adjust for both smoke and altitude
    poc_hb = case_when(
      M03_SMOKE_OECOCCUR == 1 ~ poc_hb_alti - 0.3,
      M03_SMOKE_OECOCCUR == 0 ~ poc_hb_alti,
      TRUE ~ NA_real_
    ),
    sphb_alti = case_when(
      SITE == "Kenya" ~ SPHB_LBORRES - 0.8,
      SITE == "Zambia" ~ SPHB_LBORRES - 0.8,
      !is.na(SPHB_LBORRES) ~ SPHB_LBORRES,
      TRUE ~ NA_real_
    ), 
    #adjust for both smoke and altitude
    sphb = case_when(
      M03_SMOKE_OECOCCUR == 1 ~ sphb_alti - 0.3,
      M03_SMOKE_OECOCCUR == 0 ~ sphb_alti,
      TRUE ~ NA_real_
    ),
    hb_cbc_poc = case_when(
      !is.na(hb) ~ hb, 
      !is.na(poc_hb) ~ poc_hb,
      TRUE ~ NA_real_
    ), 
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
  select(-c(hb_alti, poc_hb_alti, sphb_alti, M03_SMOKE_OECOCCUR)) 

#factor hb_level
df_hb_long$hb_level <- factor(
  df_hb_long$hb_level, 
  levels = c(1,2,3,4,5,6),
  labels = c("Severe", "Moderate", "Mild", "Normal", "High hb 13-<15g/dl", "High hb >=15g/dl")
  )


#wide hb data 
df_hb_wide <- df_hb_long %>% 
  pivot_wider(
    names_from = visit_type,
    values_from = -c("SCRNID","MOMID","PREGID","SITE", EST_CONCEP_DATE, visit_type, M09_MAT_LD_OHOSTDAT)
  ) 

save(df_hb_long, file = "derived_data/df_hb_long.rda")
save(df_hb_wide, file = "derived_data/df_hb_wide.rda")


#*********************long data including hb_type ******************************
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
  select(MOMID, PREGID, SITE, time, hb_type, hb_value, visit_type, ga_wks)

save(df_hb_3type, file = "derived_data/df_hb_3type.rda")

#*******anemia status in each trimeseter
df_anemia_trim <- MAT_ANEMIA %>%
  filter(MOMID %in% df_maternal$MOMID) %>% 
  select(-c(ANEMIA_ANC, ANEMIA_PNC6, ANEMIA_PNC26)) %>% 
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

save(df_anemia_trim, file = "derived_data/df_anemia_trim.rda")

#*******Anemia for moms who complete pregnancy
df_anemia_complet_preg <- MAT_ANEMIA %>%
  filter(MOMID %in% df_maternal$MOMID) %>% 
  left_join(MAT_ENDPOINTS) %>% 
  filter(PREG_END == 1) 

save(df_anemia_complet_preg, file = "derived_data/df_anemia_complet_preg.rda")





