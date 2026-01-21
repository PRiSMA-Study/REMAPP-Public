#*****************************************************************************
#ReMAPP Aim2 Analysis
#Author: Xiaoyan Hu
#Email; xyh@gwu.edu
#*****************************************************************************

library(tidyverse)
library(lubridate)
library(growthstandards) ## INTERGROWTH PACKAGE
library(haven)

UploadDate = "2024-09-20"

#*****************************************************************************
#*1. Load data
#*****************************************************************************
#load MAT_ENROLL 
MAT_ENROLL <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_ENROLL.dta"))
save(MAT_ENROLL, file = "derived_data/MAT_ENROLL")

#load INF_OUTCOMES(composite, preterm37, lbw2500, sga10)
INF_OUTCOMES <- read.csv(paste0("Z:/Outcome Data/",UploadDate,"/INF_OUTCOMES.csv")) 
save(INF_OUTCOMES, file = "derived_data/INF_OUTCOMES")

# load MAT_ENDPOINS
MAT_ENDPOINTS <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_ENDPOINTS.dta")) 
save(MAT_ENDPOINTS, file = "derived_data/MAT_ENDPOINTS")
MAT_ENDPOINTS <- MAT_ENDPOINTS %>% 
  select(SITE, MOMID, PREGID, PREG_END_DATE)

#mnh03 for smoking adjustment of hb
mnh03 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh03_merged.csv")) %>% 
  select(SITE, MOMID, PREGID, M03_SMOKE_OECOCCUR)

#mnh08 for hemoglobin 
mnh08 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh08_merged.csv")) %>% 
  group_by(SITE, MOMID, M08_TYPE_VISIT) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  ungroup() %>% 
  select(SITE, MOMID, PREGID, 
         M08_TYPE_VISIT, M08_LBSTDAT, M08_CBC_HB_LBORRES) %>% 
  filter(M08_TYPE_VISIT < 13) %>% 
  pivot_wider(names_from = M08_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID)) 

#mnh25 for depression
prep_mnh25 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh25_merged.csv")) %>% 
  filter(M25_TYPE_VISIT < 13) %>% 
  left_join(MAT_ENROLL %>% select(SITE, MOMID, PREGID, EST_CONCEP_DATE)) %>% 
  left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE)) %>% 
  group_by(SITE, MOMID, PREGID, M25_TYPE_VISIT) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  mutate_all(~ if_else(. < 0, NA, .)) %>% 
  mutate_all(~ if_else(. %in% c(55,77), NA, .)) %>% 
  ungroup() 
mnh25 <- prep_mnh25 %>% 
  mutate(
    ga_wks_25 = case_when(
      M25_TYPE_VISIT >= 6 ~ NA_real_,
      M25_TYPE_VISIT < 6 ~ as.numeric(ymd(M25_OBSSTDAT) - ymd(EST_CONCEP_DATE))/7
    ), 
    pst_wks_25 =case_when(
      M25_TYPE_VISIT <= 6 ~ NA_real_,
      M25_TYPE_VISIT > 6 ~ as.numeric(ymd(M25_OBSSTDAT) - PREG_END_DATE)/7
    ),
    score1 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0101 == 1 ~ case_when(
          M25_EPDS0101_Y == 2 ~ 0,
          M25_EPDS0101_Y == 1 ~ 1),
        M25_EPDS0101 %in% c(0,2) ~ case_when(
          M25_EPDS0101_N == 1 ~ 2,
          M25_EPDS0101_N == 2 ~ 3),
      ),
      SITE != "Ghana" ~ M25_EPDS0101 - 1), 
    score2 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0102 == 1 ~ case_when(
          M25_EPDS0102_Y == 2 ~ 0,
          M25_EPDS0102_Y == 1 ~ 1),
        M25_EPDS0102 %in% c(0,2) ~ case_when(
          M25_EPDS0102_N == 1 ~ 2,
          M25_EPDS0102_N == 2 ~ 3),
      ), 
      SITE != "Ghana" ~ M25_EPDS0102 - 1),
    score3 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0103 %in% c(0,2) ~ case_when(
          M25_EPDS0103_Y == 2 ~ 0,
          M25_EPDS0103_Y == 1 ~ 1),
        M25_EPDS0103 == 1 ~ case_when(
          M25_EPDS0103_N == 1 ~ 2,
          M25_EPDS0103_N == 2 ~ 3),
      ),
      SITE %in% c("Kenya") ~ M25_EPDS0103 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0103
    ),
    score4 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0104 %in% c(0,2) ~ case_when(
          M25_EPDS0104_Y == 2 ~ 0,
          M25_EPDS0104_Y == 1 ~ 1),
        M25_EPDS0104 == 1 ~ case_when(
          M25_EPDS0104_N == 1 ~ 2,
          M25_EPDS0104_N == 2 ~ 3),
      ),
      SITE != "Ghana" ~ M25_EPDS0104 - 1),
    score5 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0105 %in% c(0,2) ~ case_when(
          M25_EPDS0105_Y == 2 ~ 0,
          M25_EPDS0105_Y == 1 ~ 1),
        M25_EPDS0105 == 1 ~ case_when(
          M25_EPDS0105_N == 1 ~ 2,
          M25_EPDS0105_N == 2 ~ 3),
      ), 
      SITE == "Kenya" ~ M25_EPDS0105 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0105),
    score6 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0106 %in% c(0,2) ~ case_when(
          M25_EPDS0106_Y == 2 ~ 0,
          M25_EPDS0106_Y == 1 ~ 1),
        M25_EPDS0106 == 1 ~ case_when(
          M25_EPDS0106_N == 1 ~ 2,
          M25_EPDS0106_N == 2 ~ 3),
      ),
      SITE == "Kenya" ~ M25_EPDS0106 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0106),
    score7 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0107 %in% c(0,2) ~ case_when(
          M25_EPDS0107_Y == 2 ~ 0,
          M25_EPDS0107_Y == 1 ~ 1),
        M25_EPDS0107 == 1 ~ case_when(
          M25_EPDS0107_N == 1 ~ 2,
          M25_EPDS0107_N == 2 ~ 3),
      ),
      SITE == "Kenya" ~ M25_EPDS0107 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0107),
    score8 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0108 %in% c(0,2) ~ case_when(
          M25_EPDS0108_Y == 2 ~ 0,
          M25_EPDS0108_Y == 1 ~ 1),
        M25_EPDS0108 == 1 ~ case_when(
          M25_EPDS0108_N == 1 ~ 2,
          M25_EPDS0108_N == 2 ~ 3),
      ),
      SITE == "Kenya" ~ M25_EPDS0108 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0108),
    score9 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0109 %in% c(0,2) ~ case_when(
          M25_EPDS0109_Y == 2 ~ 0,
          M25_EPDS0109_Y == 1 ~ 1),
        M25_EPDS0109 == 1 ~ case_when(
          M25_EPDS0109_N == 1 ~ 2,
          M25_EPDS0109_N == 2 ~ 3),
      ),
      SITE == "Kenya" ~ M25_EPDS0109 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0109),
    score10 = case_when(
      SITE == "Ghana" ~ case_when(
        M25_EPDS0110 %in% c(0,2) ~ case_when(
          M25_EPDS0110_Y == 2 ~ 0,
          M25_EPDS0110_Y == 1 ~ 1),
        M25_EPDS0110 == 1 ~ case_when(
          M25_EPDS0110_N == 1 ~ 2,
          M25_EPDS0110_N == 2 ~ 3),
      ),
      SITE == "Kenya" ~ M25_EPDS0110 - 1,
      !SITE %in% c("Ghana", "Kenya") ~ 4 - M25_EPDS0110),
  ) %>% 
  rowwise() %>% 
  mutate(
    n_answered = sum(across(starts_with("score"), ~.x %in% c(0:3)), na.rm = TRUE),
    epds_score = case_when(
      n_answered != 0 ~ round(sum(c_across(starts_with("score")), na.rm = TRUE)*10/n_answered), 
      n_answered == 0 ~ NA_real_),
    depress = case_when(
      SITE == "Ghana" & epds_score >= 11 ~ 1, 
      SITE == "India-CMC" & epds_score >= 8 ~ 1,
      SITE == "India-SAS" & epds_score >= 10 ~ 1,
      SITE == "Kenya" & epds_score >= 13 ~ 1,
      SITE == "Pakistan" & epds_score >= 14 ~ 1,
      SITE == "Zambia" & epds_score >= 10 ~ 1,
      epds_score >= 0 ~ 0,
      TRUE ~ NA_real_)
  ) %>% 
  ungroup() %>% 
  select(SITE, MOMID, PREGID, M25_TYPE_VISIT, ga_wks_25, pst_wks_25, epds_score, depress,
         M25_EPDS01_SCORRES, M25_EPDS01_CAT_SCORRES, M25_EPDS0110_SCORRES) %>% 
  pivot_wider(names_from = M25_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID))

#mnh26 for fatigue
prep_mnh26 <- read.csv(paste0("Z:/Stacked Data/",UploadDate,"/mnh26_merged.csv")) %>% 
  filter(M26_TYPE_VISIT < 13) %>%
  left_join(MAT_ENROLL %>% select(SITE, MOMID, PREGID, EST_CONCEP_DATE)) %>% 
  left_join(MAT_ENDPOINTS %>% select(SITE, MOMID, PREGID, PREG_END_DATE)) %>% 
  group_by(SITE, MOMID, PREGID, M26_TYPE_VISIT) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  ungroup() %>% 
  mutate_all(~ if_else(. < 0, NA, .)) %>% 
  mutate_all(~ if_else(. %in% c(55,66,77), NA, .)) 
mnh26 <- prep_mnh26 %>% 
  mutate(
    ga_wks_26 = case_when(
      M26_TYPE_VISIT >= 6 ~ NA_real_,
      M26_TYPE_VISIT < 6 ~ as.numeric(ymd(M26_FTGE_OBSTDAT) - ymd(EST_CONCEP_DATE))/7
    ), 
    pst_wks_26 =case_when(
      M26_TYPE_VISIT <= 6 ~ NA_real_,
      M26_TYPE_VISIT > 6 ~ as.numeric(ymd(M26_FTGE_OBSTDAT) - PREG_END_DATE)/7
    )) %>% 
  select(-M26_FTGE_ASSIST, -M26_FTGE_OBSTDAT) %>% 
  rowwise() %>% 
  mutate(
    #reverse two variables with different logic than others
    M26_FTGE_AN5 = 4 - M26_FTGE_AN5, 
    M26_FTGE_AN7 = 4 - M26_FTGE_AN7, 
    n_answered = sum(across(starts_with("M26_FTGE_"), ~.x %in% c(0:4)), na.rm = TRUE),
    fatigue_score = case_when(
      n_answered != 0 ~ round(sum(4-c_across(starts_with("M26_FTGE_")), na.rm = TRUE)*13/n_answered), 
      n_answered == 0 ~ NA_real_)
  ) %>% 
  ungroup() %>% 
  select(SITE, MOMID, PREGID, M26_TYPE_VISIT, ga_wks_26, pst_wks_26, fatigue_score) %>% 
  pivot_wider(names_from = M26_TYPE_VISIT,
              values_from = -c(SITE, MOMID, PREGID))

#load MAT_HEMORRHAGE (HEM_PPH)
MAT_HEMORRHAGE <- read.csv(paste0("Z:/Outcome Data/",UploadDate,"/MAT_HEMORRHAGE.csv")) %>% 
  select(SITE, MOMID, PREGID, HEM_PPH)

#loda MAT_anmeia(ANEMIA_PNC26, ANEMIA_PNC6)
MAT_ANEMIA <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_ANEMIA.dta")) %>% 
  select(SITE, MOMID, PREGID, ANEMIA_PNC6, ANEMIA_PNC26)

# load MAT_HDP
MAT_HDP <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_HDP.dta")) 

# load MAT_PRETERM(PPROM_OCCUR)
MAT_PRETERM <- read_dta(paste0("Z:/Outcome Data/",UploadDate,"/MAT_PRETERM.dta")) 

#maternal data
df_maternal <- MAT_ENROLL %>%
  mutate(remapp = case_when(
    (SITE == "Ghana" & M02_SCRN_OBSSTDAT >= "2022-12-28") |
      (SITE == "Kenya" & M02_SCRN_OBSSTDAT >= "2023-04-03") |
      (SITE == "Zambia" & M02_SCRN_OBSSTDAT >= "2022-12-15") |
      (SITE == "Pakistan" & M02_SCRN_OBSSTDAT >= "2022-09-22" & M02_SCRN_OBSSTDAT <= "2024-04-05") |
      (SITE == "India-CMC" & M02_SCRN_OBSSTDAT >= "2023-06-20") |
      (SITE == "India-SAS" & M02_SCRN_OBSSTDAT >= "2023-08-15") ~ 1, 
    TRUE ~ NA_real_)) %>% 
  filter(remapp == 1) %>% 
  left_join(MAT_ENDPOINTS, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh03, by = c("SITE", "MOMID", "PREGID")) %>%
  left_join(mnh08, by = c("SITE", "MOMID", "PREGID")) 

#*****************************************************************************
#*2. hb data 
#*****************************************************************************
#*prepare data
prep_hb2 <- df_maternal %>% 
  dplyr:: select("SCRNID", "MOMID", "PREGID", "SITE",
                 EST_CONCEP_DATE, PREG_END_DATE,
                 M03_SMOKE_OECOCCUR,
                 num_range("M08_TYPE_VISIT_",1:12), #update to include all visits
                 num_range("M08_CBC_HB_LBORRES_",1:12), #update to include all visits
                 num_range("M08_LBSTDAT_",1:12) #update to include all visits
  ) %>% 
  mutate(across(where(is.integer), ~ as.integer(.))) %>%
  #replace 7s and 5s with NA hb can't be 0
  mutate(across(where(is.numeric), ~ ifelse(. < 0, NA, .))) %>%
  mutate(across(where(is.character), ~ case_when(
    . %in% c("1907-07-07", "1905-05-05") ~ NA_character_,
    TRUE ~ .
  ))) 


#long hb data - basic hb data (no filter)
df_hb_long2 <- prep_hb2 %>% 
  #to long format
  pivot_longer(
    -c("SCRNID","MOMID","PREGID","SITE", EST_CONCEP_DATE, PREG_END_DATE, M03_SMOKE_OECOCCUR),
    names_to = c(".value", "visit_type"), 
    names_pattern = "^M\\d{2}_(.+)_(\\d+)"
  ) %>% 
  mutate(
    ga_wks = case_when(
      TYPE_VISIT >= 6 ~ NA_real_,
      TYPE_VISIT < 6 ~ as.numeric(ymd(LBSTDAT) - ymd(EST_CONCEP_DATE))/7
    ),
    pst_wks =case_when(
      TYPE_VISIT <= 6 ~ NA_real_,
      TYPE_VISIT > 6 ~ as.numeric(ymd(LBSTDAT) - PREG_END_DATE)/7
    ),
    trimester = case_when(
      TYPE_VISIT >= 6 ~ NA_real_,
      ga_wks > 0 & ga_wks < 14 ~ 1,
      ga_wks >= 14 & ga_wks < 28 ~ 2,
      ga_wks >= 28 & ga_wks <= 43 ~ 3,
      TRUE ~ NA_real_
    ), 
    hb_alti = case_when(
      #new adjustment
      SITE %in% c("Kenya", "Zambia") ~ CBC_HB_LBORRES - 0.8,
      !is.na(CBC_HB_LBORRES) ~ CBC_HB_LBORRES,
      TRUE ~ NA_real_
    ), 
    #adjust for both smoke and altitude
    hb = case_when(
      M03_SMOKE_OECOCCUR == 1 ~ hb_alti - 0.3,
      M03_SMOKE_OECOCCUR == 0 ~ hb_alti,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-c(hb_alti, M03_SMOKE_OECOCCUR, CBC_HB_LBORRES, visit_type)) %>% 
  filter(hb > 0)

#wide hb data 
df_hb_wide2 <- df_hb_long2 %>%
  pivot_wider(
    names_from = TYPE_VISIT,
    values_from = -c("SCRNID","MOMID","PREGID","SITE", EST_CONCEP_DATE, PREG_END_DATE)
  )

#************calculate median of the hb value for ANC visits***********************************
med_hb_anc <- df_hb_long2 %>% 
  filter(TYPE_VISIT < 6) %>% 
  group_by(MOMID, PREGID, SITE) %>%
  mutate(
    hb_mom_med = median(hb, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  select(MOMID, PREGID, SITE, hb_mom_med) %>% 
  distinct() %>% 
  mutate(hb_med = median(hb_mom_med, na.rm = TRUE)) %>% 
  select(hb_med) %>% 
  distinct() %>% 
  as.numeric()

#********************df_hb_exm_anc with one hb_exm value per mom********************
set.seed(100)
df_hb_exm_anc <- df_hb_long2 %>% 
  #use ANC hb value only for current task
  filter(TYPE_VISIT < 6) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  mutate(hb_dis = abs(hb - med_hb_anc)) %>% 
  group_by(MOMID, PREGID, hb_dis) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(MOMID, PREGID) %>% 
  #calculate one extreme hb value for each mom
  mutate(
    hb_exm = case_when(
      #if the maximum distance is unique
      n == 1 & hb_dis == max(hb_dis) ~ hb,
      #if the maximum distance is not unique
      n > 1 & hb_dis == max(hb_dis) ~ sample(hb[hb_dis == max(hb_dis)],1), 
      TRUE ~ NA_real_
    ),
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  #remove the rows if the hb value is not extreme value
  filter(!is.na(hb_exm)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb_max, hb_min) %>% 
  distinct() %>% 
  mutate(hb = round(hb_exm, 1))

#************calculate median of the hb value for ANC and PNC0-4 visits***********************************
med_hb_anc_pnc0_4 <- df_hb_long2 %>% 
  filter(TYPE_VISIT < 10) %>% 
  group_by(MOMID, PREGID, SITE) %>%
  mutate(
    hb_mom_med = median(hb, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  select(MOMID, PREGID, SITE, hb_mom_med) %>% 
  distinct() %>% 
  mutate(hb_med = median(hb_mom_med, na.rm = TRUE)) %>% 
  select(hb_med) %>% 
  distinct() %>% 
  as.numeric()

#********************df_hb_exm_pnc6 with one hb_exm value per mom***************
df_hb_exm_anc_pnc0_4 <- df_hb_long2 %>% 
  #use ANC and pnc0-pnc4 hb value only 
  filter(TYPE_VISIT < 10) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  mutate(hb_dis = abs(hb - med_hb_anc_pnc0_4)) %>% 
  group_by(MOMID, PREGID, hb_dis) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(MOMID, PREGID) %>% 
  #calculate one extreme hb value for each mom
  mutate(
    hb_exm = case_when(
      #if the maximum distance is unique
      n == 1 & hb_dis == max(hb_dis) ~ hb,
      #if the maximum distance is not unique
      n > 1 & hb_dis == max(hb_dis) ~ sample(hb[hb_dis == max(hb_dis)],1), 
      TRUE ~ NA_real_
    ),
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  #remove the rows if the hb value is not extreme value
  filter(!is.na(hb_exm)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb_max, hb_min) %>% 
  distinct() %>% 
  mutate(hb = round(hb_exm, 1))

#************calculate median of the hb value for ANC and PNC0-6 visits***********************************
med_hb_anc_pnc0_6 <- df_hb_long2 %>% 
  filter(TYPE_VISIT < 11) %>% 
  group_by(MOMID, PREGID, SITE) %>%
  mutate(
    hb_mom_med = median(hb, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  select(MOMID, PREGID, SITE, hb_mom_med) %>% 
  distinct() %>% 
  mutate(hb_med = median(hb_mom_med, na.rm = TRUE)) %>% 
  select(hb_med) %>% 
  distinct() %>% 
  as.numeric()

#********************df_hb_exm_pnc6 with one hb_exm value per mom***************
df_hb_exm_anc_pnc0_6 <- df_hb_long2 %>% 
  #use ANC and pnc0-pnc4 hb value only 
  filter(TYPE_VISIT < 11) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  mutate(hb_dis = abs(hb - med_hb_anc_pnc0_6)) %>% 
  group_by(MOMID, PREGID, hb_dis) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(MOMID, PREGID) %>% 
  #calculate one extreme hb value for each mom
  mutate(
    hb_exm = case_when(
      #if the maximum distance is unique
      n == 1 & hb_dis == max(hb_dis) ~ hb,
      #if the maximum distance is not unique
      n > 1 & hb_dis == max(hb_dis) ~ sample(hb[hb_dis == max(hb_dis)],1), 
      TRUE ~ NA_real_
    ),
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  #remove the rows if the hb value is not extreme value
  filter(!is.na(hb_exm)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb_max, hb_min) %>% 
  distinct() %>% 
  mutate(hb = round(hb_exm, 1))

#************calculate median of the hb value for trimester 1***********************************
med_hb_trim1 <- df_hb_long2 %>% 
  filter(trimester == 1) %>% 
  group_by(MOMID, PREGID, SITE) %>%
  mutate(
    hb_mom_med = median(hb, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  select(MOMID, PREGID, SITE, hb_mom_med) %>% 
  distinct() %>% 
  mutate(hb_med = median(hb_mom_med, na.rm = TRUE)) %>% 
  select(hb_med) %>% 
  distinct() %>% 
  as.numeric()

#********************df_hb_exm_trim1 with one hb_exm value per mom********************
df_hb_exm_trim1 <- df_hb_long2 %>%
  filter(trimester == 1) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  mutate(hb_dis = abs(hb - med_hb_trim1)) %>% 
  group_by(MOMID, PREGID, hb_dis) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(MOMID, PREGID) %>% 
  #calculate one extreme hb value for each mom
  mutate(
    hb_exm = case_when(
      #if the maximum distance is unique
      n == 1 & hb_dis == max(hb_dis) ~ hb,
      #if the maximum distance is not unique
      n > 1 & hb_dis == max(hb_dis) ~ sample(hb[hb_dis == max(hb_dis)],1), 
      TRUE ~ NA_real_
    ),
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  #remove the rows if the hb value is not extreme value
  filter(!is.na(hb_exm)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb_max, hb_min) %>% 
  distinct() %>% 
  mutate(hb = round(hb_exm, 1))

#************calculate median of the hb value for trimester2***********************************
med_hb_trim2 <- df_hb_long2 %>% 
  filter(trimester == 2) %>% 
  # filter(hb < 30 & hb > 3) %>% 
  group_by(MOMID, PREGID, SITE) %>%
  mutate(
    hb_mom_med = median(hb, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  select(MOMID, PREGID, SITE, hb_mom_med) %>% 
  distinct() %>% 
  mutate(hb_med = median(hb_mom_med, na.rm = TRUE)) %>% 
  select(hb_med) %>% 
  distinct() %>% 
  as.numeric()

#********************df_hb_exm_trim2 with one hb_exm value per mom********************
df_hb_exm_trim2 <- df_hb_long2 %>% 
  #use ANC hb value only for current task
  filter(trimester == 2) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  mutate(hb_dis = abs(hb - med_hb_anc)) %>% 
  group_by(MOMID, PREGID, hb_dis) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(MOMID, PREGID) %>% 
  #calculate one extreme hb value for each mom
  mutate(
    hb_exm = case_when(
      #if the maximum distance is unique
      n == 1 & hb_dis == max(hb_dis) ~ hb,
      #if the maximum distance is not unique
      n > 1 & hb_dis == max(hb_dis) ~ sample(hb[hb_dis == max(hb_dis)],1), 
      TRUE ~ NA_real_
    ),
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  #remove the rows if the hb value is not extreme value
  filter(!is.na(hb_exm)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb_max, hb_min) %>% 
  distinct() %>% 
  mutate(hb = round(hb_exm, 1))

#************calculate median of the hb value for trimester 3***********************************
med_hb_trim3 <- df_hb_long2 %>% 
  filter(trimester == 3) %>% 
  group_by(MOMID, PREGID, SITE) %>%
  mutate(
    hb_mom_med = median(hb, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  select(MOMID, PREGID, SITE, hb_mom_med) %>% 
  distinct() %>% 
  mutate(hb_med = median(hb_mom_med, na.rm = TRUE)) %>% 
  select(hb_med) %>% 
  distinct() %>% 
  as.numeric()

#********************df_hb_exm_trim3 with one hb_exm value per mom********************
df_hb_exm_trim3 <- df_hb_long2 %>% 
  #use ANC hb value only for current task
  filter(trimester == 3) %>% 
  select(SITE, MOMID, PREGID, hb) %>% 
  mutate(hb_dis = abs(hb - med_hb_anc)) %>% 
  group_by(MOMID, PREGID, hb_dis) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(MOMID, PREGID) %>% 
  #calculate one extreme hb value for each mom
  mutate(
    hb_exm = case_when(
      #if the maximum distance is unique
      n == 1 & hb_dis == max(hb_dis) ~ hb,
      #if the maximum distance is not unique
      n > 1 & hb_dis == max(hb_dis) ~ sample(hb[hb_dis == max(hb_dis)],1), 
      TRUE ~ NA_real_
    ),
    hb_max = max(hb, na.rm = TRUE),
    hb_min = min(hb, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  #remove the rows if the hb value is not extreme value
  filter(!is.na(hb_exm)) %>%
  select(SITE, MOMID, PREGID, hb_exm, hb_max, hb_min) %>% 
  distinct() %>% 
  mutate(hb = round(hb_exm, 1))

#*****************************************************************************
#*3. Maternal outcome data
#*****************************************************************************
prep_dpr_ftg <- df_hb_wide2 %>% 
  select(MOMID, PREGID, SITE, EST_CONCEP_DATE, starts_with("hb_"), starts_with("ga_wks_")) %>% 
  left_join(mnh25) %>% 
  left_join(mnh26) %>% 
  select(MOMID, PREGID, SITE, EST_CONCEP_DATE,
         starts_with("hb_"), 
         starts_with("epds_score_"),
         starts_with("depress_"),
         starts_with("fatigue_score_"),
         starts_with("ga_wks"),
         starts_with("pst_wks")
  ) %>%
  mutate(across(where(is.character), ~ na_if(., "n/a"))) %>%
  mutate(across(where(is.numeric), 
                ~ ifelse(. < 0, NA, .))) %>%
  mutate(across(where(is.character), 
                ~ ifelse(. %in% c("1907-07-07", "1905-05-05"), NA, .)))

#******Depress and fatigue data*************************************************
#maternal outcome data using one hb value at ANC20 (using part of enrollment data)
dpr_ftg_anc20 <- prep_dpr_ftg %>% 
  mutate(
    dpr_score = case_when(
      !is.na(epds_score_2) & !is.na(hb_2) ~ epds_score_2,
      !is.na(epds_score_1) & !is.na(hb_1) ~ epds_score_1,
      abs(ga_wks_25_3 - 20) <= 2 & abs(ga_wks_3 - 20) <= 2 & abs(ga_wks_25_3 - ga_wks_3) <= 2 &
        !is.na(epds_score_3) & !is.na(hb_3) ~ epds_score_3,
      TRUE ~ NA_real_),
    dpr = case_when(
      !is.na(epds_score_2) & !is.na(hb_2) ~ depress_2,
      !is.na(epds_score_1) & !is.na(hb_1) ~ depress_1,
      abs(ga_wks_25_3 - 20) <= 2 & abs(ga_wks_3 - 20) <= 2 & abs(ga_wks_25_3 - ga_wks_3) <= 2 &
        !is.na(epds_score_3) & !is.na(hb_3) ~ depress_3,
      TRUE ~ NA_real_),
    ftg_score = case_when(
      !is.na(fatigue_score_2) & !is.na(hb_2) ~ fatigue_score_2,
      !is.na(fatigue_score_1) & !is.na(hb_1) ~ fatigue_score_1,
      abs(ga_wks_26_3 - 20) <= 2 & abs(ga_wks_3 - 20) <= 2 & abs(ga_wks_26_3 - ga_wks_3) <= 2 &
        !is.na(fatigue_score_3) & !is.na(hb_3) ~ fatigue_score_3,
      TRUE ~ NA_real_), 
    hb_dpr = case_when(
      !is.na(epds_score_2) & !is.na(hb_2) ~ hb_2,
      !is.na(epds_score_1) & !is.na(hb_1) ~ hb_1,
      abs(ga_wks_25_3 - 20) <= 2 & abs(ga_wks_3 - 20) <= 2 & abs(ga_wks_25_3 - ga_wks_3) <= 2 &
        !is.na(epds_score_3) & !is.na(hb_3) ~ hb_3,
      TRUE ~ NA_real_),
    hb_ftg = case_when(
      !is.na(fatigue_score_2) & !is.na(hb_2) ~ hb_2,
      !is.na(fatigue_score_1) & !is.na(hb_1) ~ hb_1,
      abs(ga_wks_26_3 - 20) <= 2 & abs(ga_wks_3 - 20) <= 2 & abs(ga_wks_26_3 - ga_wks_3) <= 2 &
        !is.na(fatigue_score_3) & !is.na(hb_3) ~ hb_3,
      TRUE ~ NA_real_),
  ) %>% 
  select(SITE, MOMID, PREGID, dpr_score, dpr, hb_dpr, ftg_score, hb_ftg)

#maternal outcome data using one hb value at ANC32 
dpr_ftg_anc32 <- prep_dpr_ftg %>% 
  mutate(
    dpr_score = case_when(
      !is.na(epds_score_4) & !is.na(hb_4) ~ epds_score_4,
      abs(ga_wks_25_3 - 32) <= 2 & abs(ga_wks_3 - 32) <= 2 & abs(ga_wks_25_3 - ga_wks_3) <= 2 &
        !is.na(epds_score_3) & !is.na(hb_3) ~ epds_score_3,
      abs(ga_wks_25_5 - 32) <= 2 & abs(ga_wks_5 - 32) <= 2 & abs(ga_wks_25_5 - ga_wks_5) <= 2 &
        !is.na(epds_score_5) & !is.na(hb_5) ~ epds_score_5,
      TRUE ~ NA_real_),
    dpr = case_when(
      !is.na(epds_score_4) & !is.na(hb_4) ~ depress_4,
      abs(ga_wks_25_3 - 32) <= 2 & abs(ga_wks_3 - 32) <= 2 & abs(ga_wks_25_3 - ga_wks_3) <= 2 &
        !is.na(epds_score_3) & !is.na(hb_3) ~ depress_3,
      abs(ga_wks_25_5 - 32) <= 2 & abs(ga_wks_5 - 32) <= 2 & abs(ga_wks_25_5 - ga_wks_5) <= 2 &
        !is.na(epds_score_5) & !is.na(hb_5) ~ depress_5,
      TRUE ~ NA_real_),
    ftg_score = case_when(
      !is.na(fatigue_score_4) & !is.na(hb_4) ~ fatigue_score_4,
      abs(ga_wks_26_3 - 32) <= 2 & abs(ga_wks_3 - 32) <= 2 & abs(ga_wks_26_3 - ga_wks_3) <= 2 &
        !is.na(fatigue_score_3) & !is.na(hb_3) ~ fatigue_score_3,
      abs(ga_wks_26_5 - 32) <= 2 & abs(ga_wks_5 - 32) <= 2 & abs(ga_wks_26_5 - ga_wks_5) <= 2 &
        !is.na(fatigue_score_5) & !is.na(hb_5) ~ fatigue_score_5,
      TRUE ~ NA_real_), 
    hb_dpr = case_when(
      !is.na(epds_score_4) & !is.na(hb_4) ~ hb_4,
      abs(ga_wks_25_3 - 32) <= 2 & abs(ga_wks_3 - 32) <= 2 & abs(ga_wks_25_3 - ga_wks_3) <= 2 &
        !is.na(epds_score_3) & !is.na(hb_3) ~ hb_3,
      abs(ga_wks_25_5 - 32) <= 2 & abs(ga_wks_5 - 32) <= 2 & abs(ga_wks_25_5 - ga_wks_5) <= 2 &
        !is.na(epds_score_5) & !is.na(hb_5) ~ hb_5,
      TRUE ~ NA_real_),
    hb_ftg = case_when(
      !is.na(fatigue_score_4) & !is.na(hb_4) ~ hb_4,
      abs(ga_wks_26_3 - 32) <= 2 & abs(ga_wks_3 - 32) <= 2 & abs(ga_wks_26_3 - ga_wks_3) <= 2 &
        !is.na(fatigue_score_3) & !is.na(hb_3) ~ hb_3,
      abs(ga_wks_26_5 - 32) <= 2 & abs(ga_wks_5 - 32) <= 2 & abs(ga_wks_26_5 - ga_wks_5) <= 2 &
        !is.na(fatigue_score_5) & !is.na(hb_5) ~ hb_5,
      TRUE ~ NA_real_),
  ) %>% 
  select(SITE, MOMID, PREGID, dpr_score, dpr, hb_dpr, ftg_score, hb_ftg)

#maternal outcome data using one hb value at PNC6
dpr_ftg_pnc6 <- prep_dpr_ftg %>% 
  mutate(
    dpr_score = case_when(
      !is.na(epds_score_10) & !is.na(hb_10) ~ epds_score_10,
      TRUE ~ NA_real_),
    dpr = case_when(
      !is.na(epds_score_10) & !is.na(hb_10) ~ depress_10,
      TRUE ~ NA_real_),
    ftg_score = case_when(
      !is.na(fatigue_score_10) & !is.na(hb_10) ~ fatigue_score_10,
      TRUE ~ NA_real_), 
    hb_dpr = case_when(
      !is.na(epds_score_10) & !is.na(hb_10) ~ hb_10,
      TRUE ~ NA_real_),
    hb_ftg = case_when(
      !is.na(fatigue_score_10) & !is.na(hb_10) ~ hb_10,
      TRUE ~ NA_real_),
  ) %>% 
  select(SITE, MOMID, PREGID, dpr_score, dpr, hb_dpr, ftg_score, hb_ftg)

df_mat_dpr_ftg <- dpr_ftg_anc20 %>% mutate(visit = "ANC20") %>% 
  bind_rows(dpr_ftg_anc32 %>% mutate(visit = "ANC32")) %>% 
  bind_rows(dpr_ftg_pnc6 %>% mutate(visit = "PNC6")) 

df_mat_dpr <- df_mat_dpr_ftg %>% 
  filter(dpr_score >= 0 & hb_dpr >= 0) %>% 
  mutate(hb = round(hb_dpr,1)) %>% 
  select(-c(ftg_score, hb_ftg))

df_mat_ftg <- df_mat_dpr_ftg %>%
  filter(ftg_score >= 0 & hb_ftg >= 0) %>% 
  mutate(hb = round(hb_ftg,1)) %>% 
  select(-c(dpr_score, hb_dpr))

#******postpartum hemorrhage****************************************************
df_mat_pph <- df_hb_exm_anc %>% 
  left_join(MAT_HEMORRHAGE %>% select(SITE, MOMID, PREGID, HEM_PPH)) %>% 
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

#******Maternal postpartum anemia at PNC6***************************************
df_mat_anemia <- MAT_ANEMIA %>% 
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

df_mat_ppa_pnc6 <- df_hb_exm_anc_pnc0_4 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

df_mat_ppa_pnc6_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

df_mat_ppa_pnc6_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

df_mat_ppa_pnc6_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc6 %in% c(0, 1))

#******Maternal postpartum anemia at PNC26***************************************
df_mat_ppa_pnc26 <- df_hb_exm_anc_pnc0_6 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

df_mat_ppa_pnc26_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_anemia) %>% 
  filter(ppa_pnc26 %in% c(0, 1))

#******Preterm premature rupture of membranes***************************************
df_mat_pprom <- df_hb_exm_anc %>%
  left_join(MAT_PRETERM %>% select(SITE, MOMID, PREGID, PPROM_OCCUR)) %>%
  mutate(
    pprom = case_when(
      PPROM_OCCUR %in% c(1,0) ~ PPROM_OCCUR,
      TRUE ~ NA_real_
    )) %>%
  filter(pprom >= 0) 

df_mat_pprom_trim1 <- df_hb_exm_trim1 %>% 
  left_join(df_mat_pprom) %>% 
  filter(pprom %in% c(0, 1))

df_mat_pprom_trim2 <- df_hb_exm_trim2 %>% 
  left_join(df_mat_pprom) %>% 
  filter(pprom %in% c(0, 1))

df_mat_pprom_trim3 <- df_hb_exm_trim3 %>% 
  left_join(df_mat_pprom) %>% 
  filter(pprom %in% c(0, 1))

#*****************************************************************************
#*4. Infant outcome data 
#*****************************************************************************
df_infant <- INF_OUTCOMES %>% 
  filter(!is.na(INFANTID)) %>% 
  select(SITE, MOMID, PREGID, INFANTID,
         LBW2500_PRISMA, LBW1500_PRISMA, 
         PRETERMBIRTH_LT37, PRETERMBIRTH_LT34, 
         SGA_CENTILE, 
         INF_ASPH,
         INF_PSBI_IPC, INF_PSBI_PNC0, INF_PSBI_PNC1, INF_PSBI_PNC4, INF_PSBI_DENOM, 
         STILLBIRTH_20WK,  STILLBIRTH_22WK, STILLBIRTH_24WK, STILLBIRTH_28WK,
         INF_HYPERBILI_TCB15_24HR, INF_HYPERBILI_TCB15_5DAY, INF_HYPERBILI_TCB15_14DAY,
         INF_HYPERBILI_AAP_24HR, INF_HYPERBILI_AAP_5DAY, INF_HYPERBILI_AAP_14DAY
         ) %>% 
  group_by(SITE, MOMID, PREGID) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  ungroup() %>% 
  filter(MOMID %in% df_maternal$MOMID) %>% 
  left_join(df_hb_exm_anc) 


#******compo data (preterm37, lbw2500, sga10)************************************
prep_compo <- df_infant %>% 
  mutate(
    preterm37 = ifelse(PRETERMBIRTH_LT37 %in% c(0, 1), PRETERMBIRTH_LT37, NA_real_), 
    lbw2500 = ifelse(LBW2500_PRISMA %in% c(0, 1), LBW2500_PRISMA, NA_real_), 
    sga10 = case_when(
      SGA_CENTILE >= 0 & SGA_CENTILE < 10 ~ 1,
      SGA_CENTILE >= 10 & SGA_CENTILE <= 100 ~ 0, 
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
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_compo_trim2 <- df_inf_compo %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_compo_trim3 <- df_inf_compo %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

df_inf_preterm37 <- prep_compo %>% 
  filter(preterm37 >= 0)

df_inf_preterm37_trim1 <- df_inf_preterm37 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_preterm37_trim2 <- df_inf_preterm37 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_preterm37_trim3 <- df_inf_preterm37 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

df_inf_lbw2500 <- prep_compo %>% 
  filter(lbw2500 >= 0)

df_inf_lbw2500_trim1 <- df_inf_lbw2500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_lbw2500_trim2 <- df_inf_lbw2500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_lbw2500_trim3 <- df_inf_lbw2500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

df_inf_sga10 <- prep_compo %>% 
  filter(sga10 >= 0)

df_inf_sga10_trim1 <- df_inf_sga10 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_sga10_trim2 <- df_inf_sga10 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_sga10_trim3 <- df_inf_sga10 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)


#******preterm34********************
df_inf_preterm34 <- df_infant %>% 
  mutate(preterm34 = ifelse(PRETERMBIRTH_LT34 %in% c(1,0), PRETERMBIRTH_LT34, NA)) %>% 
  filter(preterm34 >= 0 & hb > 0)

df_inf_preterm34_trim1 <- df_inf_preterm34 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_preterm34_trim2 <- df_inf_preterm34 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_preterm34_trim3 <- df_inf_preterm34 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

#******lbw1500********************
df_inf_lbw1500 <- df_infant %>% 
  mutate(lbw1500 = ifelse(LBW1500_PRISMA %in% c(1,0), LBW1500_PRISMA, NA)) %>% 
  filter(lbw1500 >= 0 & hb > 0)

df_inf_lbw1500_trim1 <- df_inf_lbw1500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_lbw1500_trim2 <- df_inf_lbw1500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_lbw1500_trim3 <- df_inf_lbw1500 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

#******sga3********************
df_inf_sga3 <- df_infant %>% 
  mutate(
    sga3 = case_when(
      SGA_CENTILE >= 0 & SGA_CENTILE < 3 ~ 1,
      SGA_CENTILE >= 3 & SGA_CENTILE <= 100 ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>% 
  filter(sga3 >= 0 & hb > 0)

df_inf_sga3_trim1 <- df_inf_sga3 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_sga3_trim2 <- df_inf_sga3 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_sga3_trim3 <- df_inf_sga3 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

#******neonatal psbi (defined at IPC and PNC0 visits)********************
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
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_psbi_trim2 <- df_inf_psbi %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_psbi_trim3 <- df_inf_psbi %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

#******neonatal asphyxia (defined at IPC and PNC0 visits)********************
df_inf_asph <- df_infant %>% 
  mutate(inf_asph = ifelse(INF_ASPH %in% c(0,1), INF_ASPH, NA_real_)) %>% 
  filter(hb > 0 & inf_asph >= 0)

df_inf_asph_trim1 <- df_inf_asph %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_asph_trim2 <- df_inf_asph %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_asph_trim3 <- df_inf_asph %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

#******stillbirth: death before delivery at GA >=20********************
df_inf_stillbirth20 <- df_infant %>% 
  mutate(inf_stillbirth20 = ifelse(STILLBIRTH_20WK %in% c(0,1), STILLBIRTH_20WK, NA_real_)) %>% 
  filter(hb > 0 & inf_stillbirth20 >= 0)

df_inf_stillbirth20_trim1 <- df_inf_stillbirth20 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_stillbirth20_trim2 <- df_inf_stillbirth20 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_stillbirth20_trim3 <- df_inf_stillbirth20 %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)

#******Neonatal hyperbilirubinemia********************
df_inf_hyperbili <- df_infant %>% 
  mutate(hyperbili = case_when(
    INF_HYPERBILI_AAP_24HR == 1 | INF_HYPERBILI_AAP_5DAY == 1 | INF_HYPERBILI_AAP_14DAY == 1 ~ 1,
    INF_HYPERBILI_AAP_24HR == 0 | INF_HYPERBILI_AAP_5DAY == 0 | INF_HYPERBILI_AAP_14DAY == 0 ~ 0,
    TRUE ~ NA_real_
  )) %>% 
  filter(hb > 0 & hyperbili >= 0)

df_inf_hyperbili_trim1 <- df_inf_hyperbili %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim1) %>% 
  filter(hb > 0)

df_inf_hyperbili_trim2 <- df_inf_hyperbili %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim2) %>% 
  filter(hb > 0)

df_inf_hyperbili_trim3 <- df_inf_hyperbili %>% 
  select(-starts_with("hb")) %>% 
  left_join(df_hb_exm_trim3) %>% 
  filter(hb > 0)
#*****************************************************************************
#*5. save data
#*****************************************************************************
#maternal data(ReMAPP)
save(df_maternal, file = "derived_data/df_maternal.rda")

#hb data
save(df_hb_long2, file = "derived_data/df_hb_long2.rda")
save(df_hb_wide2, file = "derived_data/df_hb_wide2.rda")
save(df_hb_exm_anc, file = "derived_data/df_hb_exm_anc.rda")

#matoutcome data
save(df_mat_dpr, file = "derived_data/df_mat_dpr.rda")
save(df_mat_ftg, file = "derived_data/df_mat_ftg.rda")

save(df_mat_pph, file = "derived_data/df_mat_pph.rda")

save(df_mat_ppa_pnc6, file = "derived_data/df_mat_ppa_pnc6.rda")
save(df_mat_ppa_pnc26, file = "derived_data/df_mat_ppa_pnc26.rda")

save(df_mat_pprom, file = "derived_data/df_mat_pprom.rda")

#trim1
save(df_mat_pph_trim1, file = "derived_data/df_mat_pph_trim1.rda")

save(df_mat_ppa_pnc6_trim1, file = "derived_data/df_mat_ppa_pnc6_trim1.rda")
save(df_mat_ppa_pnc26_trim1, file = "derived_data/df_mat_ppa_pnc26_trim1.rda")

save(df_mat_pprom_trim1, file = "derived_data/df_mat_pprom_trim1.rda")

#trim2
save(df_mat_pph_trim2, file = "derived_data/df_mat_pph_trim2.rda")

save(df_mat_ppa_pnc6_trim2, file = "derived_data/df_mat_ppa_pnc6_trim2.rda")
save(df_mat_ppa_pnc26_trim2, file = "derived_data/df_mat_ppa_pnc26_trim2.rda")

save(df_mat_pprom_trim2, file = "derived_data/df_mat_pprom_trim2.rda")

#trim3
save(df_mat_pph_trim3, file = "derived_data/df_mat_pph_trim3.rda")

save(df_mat_ppa_pnc6_trim3, file = "derived_data/df_mat_ppa_pnc6_trim3.rda")
save(df_mat_ppa_pnc26_trim3, file = "derived_data/df_mat_ppa_pnc26_trim3.rda")

save(df_mat_pprom_trim3, file = "derived_data/df_mat_pprom_trim3.rda")

#infant data
save(df_infant, file = "derived_data/df_infant.rda")

#infoutcome data
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

save(df_inf_hyperbili, file = "derived_data/df_inf_hyperbili.rda")

#trim1
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

save(df_inf_hyperbili_trim1, file = "derived_data/df_inf_hyperbili_trim1.rda")

#trim2
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

save(df_inf_hyperbili_trim2, file = "derived_data/df_inf_hyperbili_trim2.rda")

#trim3
save(df_inf_compo_trim3, file = "derived_data/df_inf_compo_trim3.rda")
save(df_inf_preterm37_trim3, file = "derived_data/df_inf_preterm37_trim3.rda")
save(df_inf_lbw2500_trim3, file = "derived_data/df_inf_lbw2500_trim3.rda")
save(df_inf_sga10_trim3, file = "derived_data/df_inf_sga10_trim3.rda")

save(df_inf_preterm34_trim3, file = "derived_data/df_inf_preterm34_trim3.rda")
save(df_inf_lbw1500_trim3, file = "derived_data/df_inf_lbw1500_trim3.rda")
save(df_inf_sga3_trim3, file = "derived_data/df_inf_sga3_trim3.rda")

save(df_inf_psbi_trim3, file = "derived_data/df_inf_psbi_trim3.rda")

save(df_inf_asph_trim3, file = "derived_data/df_inf_asph_trim3.rda")

save(df_inf_stillbirth20_trim3, file = "derived_data/df_inf_stillbirth20_trim3.rda")

save(df_inf_hyperbili_trim3, file = "derived_data/df_inf_hyperbili_trim3.rda")

#*****************************************************************************
#*6. heatmap data - added 2024-07-29
#*****************************************************************************
#Heat map data for infant outcome
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

#Heat map data for maternal outcome
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
  bind_rows(df_mat_dpr %>% 
              mutate(hb = round(hb)) %>% 
              group_by(hb) %>% 
              reframe(outcome = "Likelihood of depression", 
                      count = sum(dpr == 1), 
                      hb = first(hb)) %>% 
              ungroup())

outcome_order <- c("Postpartum Hemorrhage", "Postpartum anemia at PNC6", "Postpartum anemia at PNC26",
                   "Likelihood of depression")

df_heat_mat$outcome <- factor(df_heat_mat$outcome,levels = rev(outcome_order))

save(df_heat_mat, file = "derived_data/df_heat_mat.rda")
