#*****************************************************************************
#*derived variables for statistical analysis plan
#*
#*input: matData_Wide.RData, 
#*       matData_Anc_Visits.RData, 
#*       infData_Wide. Rdata
     
#*output: df_maternal.rda (keep PRiSMA enrolled only)
#*        df_infant.rda
#*        df_baseline.rda
#*        df_hb_long.rda
#*        df_hb_wide.rda
#*        df_sensitvie.rda
    
#*****************************************************************************
rm(list = ls())

library(tidyverse)
library(lubridate)
library(naniar)

#*****************************************************************************
#*set directory and prepare data 
#*****************************************************************************
#******read data 
uploadDate = "2023-06-09"
# today's date
todayDate <- Sys.Date()

setwd(paste0("~/github/ReMAPP-SAP/", uploadDate))

load(paste0("Z:/Merged Data/", uploadDate, "/MatData_Wide_", uploadDate, ".RData"))
load(paste0("Z:/Merged Data/", uploadDate, "/MatData_Anc_Visits_", uploadDate, ".RData"))
load(paste0("Z:/Merged Data/", uploadDate, "/InfData_Wide_", uploadDate, ".RData"))

#******prepare maternal and infant data 
#df_maternal
df_maternal <- MatData_Wide %>% 
  distinct() %>% 
  #keep PRiSMA enrolled only
  filter(M02_CONSENT_IEORRES_1 == 1) %>% 
  #merge EDD variables in
  left_join(MatData_Anc_Visits %>% dplyr::select(MOMID, PREGID, EDD), by = c("MOMID", "PREGID")) %>% 
  distinct() %>% 
  #temp solution of Ghana data (Some EDD are earlier than US visit date). Remove this step if we have final data
  mutate(
    bestedd = case_when(
      SITE != "Ghana" ~ EDD,
      SITE == "Ghana" & M01_US_OHOSTDAT_1 < EDD ~ EDD, 
      SITE == "Ghana" & M01_US_OHOSTDAT_1 >= EDD ~ dmy("1907-07-07")
    )
)

df_infant <- InfData_Wide %>% 
  #add all M09 birthoutcome variables
  left_join(df_maternal %>% dplyr::select("MOMID", "PREGID", "SITE", matches("_INF1_")), 
            by = c("MOMID", "PREGID", "INFANTID"="M09_INFANTID_INF1_6", "SITE")) %>% 
  left_join(df_maternal %>% dplyr::select("MOMID", "PREGID", "SITE", matches("_INF2_")), 
            by = c("MOMID", "PREGID", "INFANTID"="M09_INFANTID_INF2_6", "SITE")) %>% 
  left_join(df_maternal %>% dplyr::select("MOMID", "PREGID", "SITE", matches("_INF3_")), 
            by = c("MOMID", "PREGID", "INFANTID"="M09_INFANTID_INF3_6", "SITE")) %>% 
  left_join(df_maternal %>% dplyr::select("MOMID", "PREGID", "SITE", matches("_INF4_")), 
            by = c("MOMID", "PREGID", "INFANTID"="M09_INFANTID_INF4_6", "SITE")) %>% 
  #add M01 ultrasound GA, date (using GA_US_DAYS_1 should be better than using individual ga days?)
  left_join(df_maternal %>% dplyr::select("MOMID", "PREGID", "SITE", "M01_US_OHOSTDAT_1", "GA_US_DAYS_1"), 
            by = c("MOMID", "PREGID", "SITE")) 

save(df_maternal, file = "derived_data/df_maternal.rda")
save(df_infant, file = "derived_data/df_infant.rda")

#*****************************************************************************
#*Baseline characters data
#*age, bmi_enroll, ga_wks_enroll, school_yrs, married, nulliparous
#*****************************************************************************

#******prepare data
prep_baseline <- df_maternal %>% 
  dplyr::select("SCRNID", "MOMID", "PREGID", "SITE",
         M00_BRTHDAT_1, M00_ESTIMATED_AGE_1,
         M00_SCHOOL_YRS_SCORRES_1, M00_SCHOOL_SCORRES_1,
         num_range("M01_US_EDD_BRTHDAT_FTS",1:4,"_1"),bestedd,GA_US_DAYS_1,
         M02_SCRN_OBSSTDAT_1, 
         M03_MARITAL_SCORRES_1,
         M04_PH_PREV_RPORRES_1,
         M05_WEIGHT_PERES_1, M05_HEIGHT_PERES_1
         ) %>% 
  #replace default value 7s with NA
  replace_with_na_all(condition = ~.== -7) %>%
  replace_with_na_all(condition = ~.== "1907-07-07") %>% 
  distinct()  

#other baseline variable
df_baseline <- prep_baseline %>%
  mutate(
    #Age
    age = case_when(
      !is.na(M02_SCRN_OBSSTDAT_1) & !is.na(M00_BRTHDAT_1) ~ as.numeric(M02_SCRN_OBSSTDAT_1 - M00_BRTHDAT_1) %/% 365,
      !is.na(M00_ESTIMATED_AGE_1) ~ as.numeric(M00_ESTIMATED_AGE_1),
      TRUE ~ NA_real_
    ),
    #!!!replace age == 0 with NA temp solution
    age = ifelse(age == 0, NA, age),
    #Body mass index 
    bmi_enroll = case_when(
      M05_WEIGHT_PERES_1 > 0 & M05_HEIGHT_PERES_1 > 0 ~ 
      M05_WEIGHT_PERES_1 / M05_HEIGHT_PERES_1 / M05_HEIGHT_PERES_1 * 10000, 
      TRUE ~ NA_real_
    ), 
    #GA
    ga_wks_enroll = GA_US_DAYS_1 / 7,
    #Years of formal education 
    school_yrs = case_when(
      M00_SCHOOL_YRS_SCORRES_1 >= 0 ~ as.numeric(M00_SCHOOL_YRS_SCORRES_1), 
      M00_SCHOOL_SCORRES_1 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    #Married or cohabiting
    married = case_when(
      M03_MARITAL_SCORRES_1 %in% c(1,2) ~ 1, #married/cohabiting
      M03_MARITAL_SCORRES_1 %in% c(3:5) ~ 0, #divorced/permanently separated/widowed/single,never married
      TRUE ~ NA_real_
    ),
    #Nulliparous
    nulliparous = case_when(
      M04_PH_PREV_RPORRES_1 == 0 ~ 1, #never pregnant before
      M04_PH_PREV_RPORRES_1 == 1 ~ 0, 
      TRUE ~ NA_real_
    )) %>% 
  dplyr::select(-c(matches("M\\d{2}_"),GA_US_DAYS_1, bestedd))

save(df_baseline, file = "derived_data/df_baseline.rda")

#*****************************************************************************
#*HB data
#*****************************************************************************
#*prepare data
prep_hb <- df_maternal %>% 
  dplyr:: select("SCRNID", "MOMID", "PREGID", "SITE",
                 bestedd, GA_US_DAYS_1,
                 num_range("M08_CBC_HB_LBORRES_",1:5),
                 #NOTE: we should use MAT_SPEC_COLLECT_DAT to calculate GA once data is cleaned. 
                 # num_range("M07_MAT_SPEC_COLLECT_DAT_",1:5),
                 #Use LBSTDAT as temp solution for issues in MAT_SPEC_COLLECT_DAT
                 num_range("M08_LBSTDAT_",1:5)
  ) %>%
  #replace 7s with NA
  replace_with_na_all(condition = ~.== -7) %>%
  replace_with_na_all(condition = ~.== "1907-07-07") %>% 
  distinct() 


#long hb data
df_hb_long_base <- prep_hb %>% 
  #to long format
  pivot_longer(
    -c("SCRNID","MOMID","PREGID","SITE","bestedd", "GA_US_DAYS_1"),#variables to be longer
    names_to = c(".value", "visit_type"), 
    names_pattern = "^M\\d{2}_(.+)_(\\d+)"
  ) %>% 
  #rename var names (replace CBC_HB_LBORRES with MAT_SPEC_COLLECT_DAT once data is cleaned)
  rename(hb = CBC_HB_LBORRES, collect_date = LBSTDAT) %>% 
  #remove NAs
  filter(!is.na(hb) & !is.na(bestedd) & !is.na(collect_date)) %>% 
  #derive variables
  mutate(
    ga_wks = (280 - as.numeric(bestedd - collect_date)) / 7, 
    trimester = case_when(
      ga_wks > 3 & ga_wks < 14 ~ 1,
      ga_wks >= 14 & ga_wks < 27 ~ 2,
      ga_wks >= 27 & ga_wks < 43 ~ 3,
      TRUE ~ -5
    ),
    anemia_level = case_when(
      hb >= 10 & hb < 11 ~ 1, #mild
      hb >= 7 & hb < 10 ~ 2, #moderate
      hb < 7 ~ 3, #severe
      TRUE ~ -5
    ),
    anemia = case_when(
      hb < 11 ~ 1, 
      hb >= 11 ~ 0,
      TRUE ~ -5
    )
  ) %>% 
  #Remove rows that have gestational weeks less than 10 or greater than 50 (42-week pregnancy + 8-week postpartum).
  #Remove rows that have hemoglobin less than 5 or greater than 18.
  #Remove rows that have missing values on gestational weeks.
  filter(ga_wks >= 10 & ga_wks <= 50 ) %>% 
  filter(hb >= 5 & hb <= 18) %>% 
  filter(!is.na(ga_wks) & !is.na(hb))


#wide hb data 
df_hb_wide <- df_hb_long_base %>% 
  group_by(SCRNID, MOMID, PREGID, SITE) %>% 
  #find the last visit data for the mom
  mutate(lastVisit = max(as.numeric(visit_type))) %>% 
  ungroup() %>% 
  #to wide dat
  pivot_wider(
    names_from = visit_type,
    values_from = c("collect_date", "hb", "ga_wks", "trimester", "anemia", "anemia_level")
  ) %>% 
  rowwise() %>% 
  #latest hb status for the mom
  mutate(
    trimester = eval(parse(text = paste0("trimester_", lastVisit))),
    hb = eval(parse(text = paste0("hb_", lastVisit))),
    anemia = eval(parse(text = paste0("anemia_", lastVisit))),
    anemia_level = eval(parse(text = paste0("anemia_level_", lastVisit))),
  ) %>% 
  ungroup()

#save data
save(df_hb_wide, file = "derived_data/df_hb_wide.rda")

#*****************************************************************************
#*Sensitivity Analysis Data
#*****************************************************************************
#*prepare data
prep_sensitive_long <- df_hb_long_base %>% 
  left_join(df_maternal %>% 
              dplyr::select(M01_US_OHOSTDAT_1, #use us ga for each fetus?
                     num_range("M09_INFANTID_INF",1:4,"_6"), 
                     num_range("M09_DELIV_DSSTDAT_INF",1:4,"_6"), 
                     num_range("M09_BIRTH_DSTERM_INF",1:4,"_6")), 
            by = c("SCRNID", "MOMID", "PREGID", "SITE")) %>% 
  replace_with_na_all(condition = ~.== "1907-07-07") %>%
  replace_with_na_all(condition = ~.x %in% c("n/a")) %>%
  replace_with_na_at(.vars = c("M09_BIRTH_DSTERM_INF1_6", "M09_BIRTH_DSTERM_INF2_6",
                               "M09_BIRTH_DSTERM_INF3_6", "M09_BIRTH_DSTERM_INF4_6"),
                     condition = ~.x == 77) %>%
  distinct() 

#long data
df_sensitive_long <- prep_sensitive_long %>% 
  mutate(
    preterm = case_when(
      ((M09_DELIV_DSSTDAT_INF1_6 - M01_US_OHOSTDAT_1) + GA_US_DAYS_1)/7 < 37 &
        M09_BIRTH_DSTERM_INF1_6 == 1 ~ 1, #preterm (as long as there is one preterm, woman has preterm birth)
      ((M09_DELIV_DSSTDAT_INF2_6 - M01_US_OHOSTDAT_1) + GA_US_DAYS_1)/7 < 37 &
        M09_BIRTH_DSTERM_INF2_6 == 1 ~ 1, 
      ((M09_DELIV_DSSTDAT_INF3_6 - M01_US_OHOSTDAT_1) + GA_US_DAYS_1)/7 < 37 &
        M09_BIRTH_DSTERM_INF3_6 == 1 ~ 1, 
      ((M09_DELIV_DSSTDAT_INF4_6 - M01_US_OHOSTDAT_1) + GA_US_DAYS_1)/7 < 37 &
        M09_BIRTH_DSTERM_INF4_6 == 1 ~ 1, 
      M09_BIRTH_DSTERM_INF1_6 == 1 ~ 2, #term birth
      M09_BIRTH_DSTERM_INF2_6 == 1 ~ 2,
      M09_BIRTH_DSTERM_INF3_6 == 1 ~ 2,
      M09_BIRTH_DSTERM_INF4_6 == 1 ~ 2,
      M09_BIRTH_DSTERM_INF1_6 == 2 ~ 3, #fetal death
      M09_BIRTH_DSTERM_INF2_6 == 2 ~ 3,
      M09_BIRTH_DSTERM_INF3_6 == 2 ~ 3,
      M09_BIRTH_DSTERM_INF4_6 == 2 ~ 3,
      TRUE ~ 55)) 

#*****************************************************************************
#*save hb data with all derived variables 
#*****************************************************************************
#add preterm variable
df_hb_long <- df_hb_long_base %>% 
  left_join(df_sensitive_long %>% 
              dplyr::select(SCRNID, MOMID, PREGID, SITE, visit_type, preterm), 
            by = c("SCRNID", "MOMID", "PREGID", "SITE", "visit_type")) %>% 
  distinct()
#save data
save(df_hb_long, file = "derived_data/df_hb_long.rda")

#*****************************************************************************
#*Infant outcome data
#*****************************************************************************
#*prepare data
prep_infoutcome_long <- df_infant %>% 
  dplyr::select(MOMID, PREGID, INFANTID, SITE, 
         M01_US_OHOSTDAT_1, GA_US_DAYS_1,
         num_range("M09_DELIV_DSSTDAT_INF",1:4,"_6"), 
         num_range("M09_BIRTH_DSTERM_INF",1:4,"_6"),
         M11_BW_FAORRES, M11_BW_EST_FAORRES
         ) %>% 
  replace_with_na_all(condition = ~.< 0) %>%
  replace_with_na_all(condition = ~.== "1907-07-07") %>%
  replace_with_na_all(condition = ~.x %in% c("n/a")) %>%
  replace_with_na_at(.vars = c("M09_BIRTH_DSTERM_INF1_6", "M09_BIRTH_DSTERM_INF2_6",
                               "M09_BIRTH_DSTERM_INF3_6", "M09_BIRTH_DSTERM_INF4_6"),
                     condition = ~.x == 77) %>% distinct() 

#derive variables
df_infoutcome_long <- prep_infoutcome_long %>% 
mutate(
  delivery_date = case_when(
    !is.na(M09_DELIV_DSSTDAT_INF1_6) ~ M09_DELIV_DSSTDAT_INF1_6, 
    !is.na(M09_DELIV_DSSTDAT_INF2_6) ~ M09_DELIV_DSSTDAT_INF2_6, 
    !is.na(M09_DELIV_DSSTDAT_INF3_6) ~ M09_DELIV_DSSTDAT_INF3_6, 
    !is.na(M09_DELIV_DSSTDAT_INF4_6) ~ M09_DELIV_DSSTDAT_INF4_6
  ), 
  birth_outcome = case_when(
    !is.na(M09_BIRTH_DSTERM_INF1_6) ~ M09_BIRTH_DSTERM_INF1_6, 
    !is.na(M09_BIRTH_DSTERM_INF2_6) ~ M09_BIRTH_DSTERM_INF2_6, 
    !is.na(M09_BIRTH_DSTERM_INF3_6) ~ M09_BIRTH_DSTERM_INF3_6, 
    !is.na(M09_BIRTH_DSTERM_INF4_6) ~ M09_BIRTH_DSTERM_INF4_6
  ),
  preterm = case_when(
    ((delivery_date - M01_US_OHOSTDAT_1) + GA_US_DAYS_1)/7 < 37 & birth_outcome == 1 ~ 1, #preterm
    ((delivery_date - M01_US_OHOSTDAT_1) + GA_US_DAYS_1)/7 >= 37 & birth_outcome == 1 ~ 0, #term birth
    birth_outcome == 2 ~ 3, #fetal death
    TRUE ~ 55
  ),
  lbw = case_when(
    #temp code for Pakistan unit problem
    SITE == "Pakistan" & M11_BW_FAORRES < 2.5 & M11_BW_EST_FAORRES < 72 ~ 1, #lbw
    SITE == "Pakistan" & M11_BW_FAORRES >= 2.5 & M11_BW_EST_FAORRES < 72 ~ 1, #lbw
    M11_BW_FAORRES < 2500 & M11_BW_EST_FAORRES < 72 ~ 1, #lbw
    M11_BW_FAORRES >= 2500 & M11_BW_EST_FAORRES < 72 ~ 0, #not lbw
    TRUE ~ 55 
    )) %>% 
  left_join(df_hb_wide %>% dplyr::select(MOMID, PREGID, SITE, starts_with("hb_")),
            by = c("MOMID", "PREGID", "SITE")) %>% 
  pivot_longer(starts_with("hb_"),
    names_to = c(".value", "visit_type"), 
    names_pattern = "(.+)_(\\d+)"
  ) %>% 
  filter(hb >= 5 & hb <= 18 & !is.na(hb)) 

#df_inf_preterm data
df_inf_preterm_prep <- df_infoutcome_long %>%  
  dplyr::select(MOMID, PREGID, INFANTID, SITE, preterm, hb) %>%
  filter(preterm == 1 | preterm == 0)  #keep preterm and term delivery

df_inf_preterm <- df_inf_preterm_prep %>% 
  arrange(hb) %>% 
  mutate(group1 = round(hb,1), #round to one decimal place
         group2 = ceiling(hb/0.5)*0.5, #use the larger boundary every 0.5 interval
         group_ref = ntile(x = hb, n = round(nrow(df_inf_preterm_prep)/25))) %>%  #each group have 25 and above people
  group_by(group_ref) %>% 
  mutate(group3 = median(hb)) %>%
  ungroup()

save(df_inf_preterm, file = "derived_data/df_inf_preterm.rda")


#******************************primary cohort maternal outcomes****************
#*Please use this part as reference ONLY, it was developed at early stage, the definition and 
#*also the variables might changed.
#******Composite maternal morbidity and mortality through 42 days postpartum 
#******maternal death (maternal mortality)
#*ANC
# var_mat_death_anc <- matData %>% 
#   select(all_of(matID), 
#          matches("M23_CLOSE_DSDECOD|M23_ACC_DDORRES"),
#          matches("M09_MAT_LD_OHOSTDAT"), 
#          matches("MAT_VITAL_MNH0\\d{1}_V[1-5]"),
#          matches("M16_MAT_VITAL_MNH16"),
#          matches("M19_MAT_ARRIVAL_DSDECOD_V[1-5]"),
#          matches("M19_VISIT_FAORRES_V[1-5]"),
#          matches("M21_AETERM")
#          ) %>% 
#   mutate(
#     #maternal mortality at ANC
#     mat_death_anc = case_when(
#       # M23_CLOSE_DSDECOD == 3 & M23_ACC_DDORRES == 0 &
#       #   is.na(M09_MAT_LD_OHOSTDAT) #close out due to mother's death, but not caused by incident or accident
#       # ~ 1,
#       # M23_CLOSE_DSDECOD == 3 & M23_ACC_DDORRES == 1 &
#       #   (is.na(M09_MAT_LD_OHOSTDAT) #close out due to mother's death which caused by incident or accident
#       # ~ 0, 
#       if_any(matches("M04_MAT_VITAL_MNH0\\d{1}_V\\d{1}"), ~. == 2) & 
#         if_any(matches("M04_ACC_DDORRES\\d{1}_V\\d{1}"), ~. == 2) #ANC
#       ~ 1, 
#       if_any(matches("M04_MAT_VITAL_MNH0\\d{1}_V\\d{1}"), ~. == 2) & 
#         if_any(matches("M04_ACC_DDORRES\\d{1}_V\\d{1}"), ~. == 1) #ANC death due to accident
#       ~ 0, 
#       #code below to be considered whether accident or not
#       if_any(matches("MAT_VITAL_MNH0\\d{1}_V[1-5]"), ~. == 2) #| #any death during ANC 
#         # M16_MAT_VITAL_MNH16 == 2 | #ANC exit
#         # if_any(matches("^M19_MAT_ARRIVAL_DSDECOD_V[1-5]"), ~.x == 2)  #hospital death before care at ANC
#         # if_any(matches("^M19_VISIT_FAORRES_V[1-5]"), ~.x == 5) |#hospital death after care at ANC
#         # M21_AETERM == 1 | is.na(M09_MAT_LD_OHOSTDAT) #adverse event - maternal death
#       ~ 1, 
#       if_all(matches("MAT_VITAL_MNH0\\d{1}_V[1-5]"), ~. %in% c(1,NA)) #| #any death during ANC 
#       # M16_MAT_VITAL_MNH16 == 2 | #ANC exit
#       # if_any(matches("^M19_MAT_ARRIVAL_DSDECOD_V[1-5]"), ~.x != 2)  
#       # if_any(matches("^M19_VISIT_FAORRES_V[1-5]"), ~.x != 5) 
#       # M21_AETERM != 1 | is.na(M09_MAT_LD_OHOSTDAT) #adverse event - not maternal death
#       ~ 0, 
#       TRUE ~ -5
#     )
#   ) %>% 
#   select(-c(matches("M\\d{2}_")))
# 
# #IPC
# var_mat_death_ipc <- matData %>% 
#   select(all_of(matID), 
#          matches("M23_CLOSE_DSDECOD|M23_ACC_DDORRES"),
#          matches("M09_MAT_LD_OHOSTDAT|M09_MAT_VITAL_MNH09"), 
#          matches("M10_MAT_VITAL_MNH10|M10_MAT_DSTERM"), 
#          matches("M17_MAT_VITAL_MNH17"),
#          matches("M19_MAT_ARRIVAL_DSDECOD_V6|M19_VISIT_FAORRES_V6"),
#          matches("M21_AETERM")
#   ) %>% 
#   mutate(
#     mat_death_ipc = case_when(
#       #maternal mortality at IPC
#       # (M23_CLOSE_DSDECOD == 3  & M23_ACC_DDORRES == 0 & 
#       #    (M09_MAT_VITAL_MNH09 == 2 | M10_MAT_VITAL_MNH10 == 2 | M10_MAT_DSTERM == 3)) #IPC maternal death, not accident
#       # ~ 1,
#       # (M23_CLOSE_DSDECOD == 3  & M23_ACC_DDORRES == 1 & 
#       #    (M09_MAT_VITAL_MNH09 == 2 | M10_MAT_VITAL_MNH10 == 2 | M10_MAT_DSTERM == 3)) #IPC maternal death, accident
#       #  ~ 0,
#       #code below can be considered
#       #   M17_MAT_VITAL_MNH17 == 2 #IPC exit
#       # ~ 1, 
#       # # if_any(matches("^M19_MAT_ARRIVAL_DSDECOD_V6"), ~.x ==2)  
#       # # if_any(matches("^M19_VISIT_FAORRES_V6"), ~.x == 5) ~ 1, 
#       # # if_any(matches("^M19_MAT_ARRIVAL_DSDECOD_V6"), ~.x !=2)  
#       # # if_any(matches("^M19_VISIT_FAORRES_V6"), ~.x != 5) ~ 0,         
#       # M21_AETERM == 1 & is.na(M09_MAT_LD_OHOSTDAT) #adverse event #revise
#       # ~ 0, 
#       TRUE ~ -5
#     )
#   ) %>% 
#   select(-c(matches("M\\d{2}_")))
# 
# #PNC
# var_mat_death_pnc42 <- matData %>% 
#   select(all_of(matID), 
#          matches("M23_CLOSE_DSDECOD|M23_ACC_DDORRES"),
#          matches("M09_MAT_LD_OHOSTDAT|M09_MAT_VITAL_MNH09"), 
#          matches("MAT_VITAL_MNH\\d{2}_V[7-9,10]"),
#          matches("M12_MATERNAL_DSDECOD_V[7-9,10]"),
#          matches("M18_MAT_VITAL_MNH18"),
#          matches("M19_MAT_ARRIVAL_DSDECOD_V[7-9,10]"),
#          matches("19_VISIT_FAORRES_V[7-9,10]"),
#          matches("M21_AETERM")
#   ) %>% 
#   mutate(
#     mat_death_pnc42 = case_when(
#       #maternal mortality at PNC < 42 days type = 10
#       # M23_CLOSE_DSDECOD == 3 & M23_ACC_DDORRES == 0 & 
#       #   (is.na(M09_MAT_LD_OHOSTDAT) & (dmy(M23_CLOSE_DSSTDAT) - dmy(M09_MAT_LD_OHOSTDAT)) <= 42) #close out due to mother's death, not accident, also before 42 days PNC
#       # ~ 1, 
#       # M23_CLOSE_DSDECOD == 3 & M23_ACC_DDORRES == 1 &
#       #   (is.na(M09_MAT_LD_OHOSTDAT) & (dmy(M23_CLOSE_DSSTDAT) - dmy(M09_MAT_LD_OHOSTDAT)) <= 42) #close out due to mother's death which caused by incident or accident, also before 42 days PNC
#       # ~ 0, 
#       # if_any(matches("MAT_VITAL_MNH\\d{2}_V[7-10]"), ~. == 2) #| #any death during ANC or IPC
#       #   # if_any(matches("M12_MATERNAL_DSDECOD_V[7-10]", ~. == 3) 
#       # ~ 1, 
#       # # (M18_MAT_VITAL_MNH18 == 2 & dmy(M18_VISDAT) - dmy(M09_MAT_LD_OHOSTDAT) <= 42) |#PNC exit
#       # # if_any(matches("^M19_MAT_ARRIVAL_DSDECOD_V\\d{1}"), ~.x == 2)  | M19_MAT_ARRIVAL_DSDECOD_V10 == 2 | 
#       # # if_any(matches("^M19_VISIT_FAORRES_V\\d{1}"), ~.x == 5)  | M19_VISIT_FAORRES_V10 == 5 |#hospital until 42 days PNC
#       # # M21_AETERM == 1  #adverse event
#       # # ~ 1, 
#       # if_any(matches("MAT_VITAL_MNH\\d{2}_V[7-10]"), ~. %in% c(1,NA)) #| #any death during pnc
#       #   # if_any(matches("M12_MATERNAL_DSDECOD_V[7-10]", ~. != 3) 
#       #   ~ 0, 
#       # # (M18_MAT_VITAL_MNH18 == 2 & dmy(M18_VISDAT) - dmy(M09_MAT_LD_OHOSTDAT) <= 42) |#PNC exit
#       # # if_any(matches("^M19_MAT_ARRIVAL_DSDECOD_V\\d{1}"), ~.x == 2)  | M19_MAT_ARRIVAL_DSDECOD_V10 == 2 | 
#       # # if_any(matches("^M19_VISIT_FAORRES_V\\d{1}"), ~.x == 5)  | M19_VISIT_FAORRES_V10 == 5 |#hospital until 42 days PNC
#       # # M21_AETERM == 1  #adverse event #revise
#       # # ~ 0, 
#       TRUE ~ -5
#     )
#   ) %>% 
#   select(-c(matches("M\\d{2}_")))
# 
# vars_mat_death <- var_mat_death_anc %>% 
#   left_join(var_mat_death_ipc, by=c(all_of(matID))) %>% 
#   left_join(var_mat_death_pnc42, by=c(all_of(matID))) %>% 
#   mutate(
#   mat_death = case_when(
#     mat_death_anc == 1 | mat_death_ipc == 1 | mat_death_pnc42 == 1 ~ 1, 
#     mat_death_pnc42 == 0 ~ 0,
#     mat_death_ipc == 0 & mat_death_pnc42 == -5 ~ 0,
#     mat_death_anc == 0 & mat_death_ipc == -5 & mat_death_pnc42 == -5 ~ 0,
#     TRUE ~ -5
#    )
#   )
# 
# #******WHO maternal near-miss
#   vars_near_miss <- matData %>% 
#     select(all_of(matID), 
#            num_range("M19_ORGAN_FAIL_MHTERM_", c(1:7,88)),
#            num_range("M09_ORG_FAIL_MHOCCUR_", c(1:7,88)), 
#            num_range("M12_ORG_FAIL_MHOCCUR_", c(1:4,88))
#     ) %>% 
#     left_join(vars_mat_death %>% select(all_of(matID), mat_death), by=c(all_of(matID))) %>% 
#     mutate(
#       near_miss = case_when(
#     if_any(matches("M19_ORGAN_FAIL_MHTERM_"), ~. == 1) & mat_death == 0 ~ 1, #hospital
#     if_any(matches("M09_ORG_FAIL_MHOCCUR_"), ~. == 1) & mat_death == 0 ~ 1, #IPC
#     if_any(matches("M12_ORG_FAIL_MHOCCUR_"), ~. == 1) & mat_death == 0 ~ 1, #PNC
#     is.na(mat_death) ~ 0,
#     TRUE ~ -5
#     )
#   ) %>% 
#     select(-c(matches("M\\d{2}_")))
# 
#   #******evere complication(potential life_threatening complications)
#   #******postpartum hemorrhage at Labor
#   var_hemorr_ld <- matData %>% 
#     select(all_of(matID), M09_PPH_CEOCCUR, starts_with("M09_DELIV_PRROUTE_INF"), M09_PPH_ESTIMATE_FASTAT, M09_PPH_ESTIMATE_FAORRES) %>% 
#     mutate(
#       hemorr_ld = case_when(
#         # hemorrhage
#         M09_PPH_CEOCCUR == 1 ~ 1, #Did mother experience postpartum hemorrhage (PPH)?
#         M09_PPH_ESTIMATE_FASTAT == 1 & M09_PPH_ESTIMATE_FAORRES > 1000 ~ 1, # estimated blood loss
#         if_any(starts_with("M09_DELIV_PRROUTE_INF"), ~. == 2) & M09_PPH_ESTIMATE_FASTAT == 1 & M09_PPH_ESTIMATE_FAORRES > 1000 ~ 1, #Cesarean
#         if_any(starts_with("M09_DELIV_PRROUTE_INF"), ~. == 1) & M09_PPH_ESTIMATE_FASTAT == 1 & M09_PPH_ESTIMATE_FAORRES > 500 ~ 1, #Vaginal
#         # no hemorrhage
#         M09_PPH_CEOCCUR == 0 ~ 0,
#         M09_PPH_ESTIMATE_FASTAT == 1 & M09_PPH_ESTIMATE_FAORRES <= 500 ~ 0,
#         if_any(starts_with("M09_DELIV_PRROUTE_INF"), ~. == 2) & M09_PPH_ESTIMATE_FASTAT == 1 & M09_PPH_ESTIMATE_FAORRES <= 1000 ~ 0,
#         if_any(starts_with("M09_DELIV_PRROUTE_INF"), ~. == 1) & M09_PPH_ESTIMATE_FASTAT == 1 & M09_PPH_ESTIMATE_FAORRES <= 500 ~ 0,
#         # don't know
#         M09_PPH_CEOCCUR == 99 ~ 99,
#         # missing
#         TRUE ~ -5
#       )
#     ) %>% 
#     select(-c(matches("M\\d{2}_")))
#   
#   #postpartum hemorrhage at PNC
#   var_hemorr_pnc <- matData %>%
#     # # keep subjects with M12 visits
#     # filter(!is.na(M12_VISIT_TOT)) %>%
#     select(all_of(matID), M12_VISIT_TOT, matches("^M12_(BIRTH_COMPL_MHTERM_1)_V\\d+")) %>%
#     pivot_longer(!c(all_of(matID), M12_VISIT_TOT),
#                  names_to = c(".value", "pnc_order"),
#                  names_pattern = "^M12_(.+)_V(\\d+)"
#     ) %>%
#     mutate(
#       hemorr_pnc_t1 = case_when(
#         # hemorrhage at PNC visit
#         BIRTH_COMPL_MHTERM_1 == 1 ~ 1,
#         # no hemorrhage at PNC visit
#         BIRTH_COMPL_MHTERM_1 == 0 ~ 0,
#         TRUE ~ -5
#       )
#     ) %>%
#     mutate(pnc_order = as.numeric(pnc_order)) %>%
#     group_by(SCRNID, MOMID, PREGID, SITE) %>%
#     # arrange(pnc_order) %>%
#     # slice(c(1:M12_VISIT_TOT[1])) %>%
#     summarise(
#       hemorr_pnc = case_when(
#         # if any PNC visit has hemorrhage diagnosed
#         1 %in% hemorr_pnc_t1 ~ 1,
#         # if none PNC visit has hemorrhage
#         0 %in% hemorr_pnc_t1 ~ 0, ###*** 0 if all reported 0
#         TRUE ~ -5
#       )
#     ) %>%
#     ungroup() %>% 
#     select(-c(matches("M\\d{2}_")))
#   
#   #postpartum hemorrhage
#   vars_pos_hemorr <- vars_hemorr_ld %>%
#     full_join(vars_hemorr_pnc, by = c(all_of(matID))) %>%
#     distinct() %>%
#     mutate(
#       pos_hemorr = case_when(
#         # hemorrhage occurs
#         hemorr_ld == 1 | hemorr_pnc == 1 ~ 1,
#         # hemorrhage not occurs
#         # if hemorrhage not occurs at LD or no visit at LD
#         # And if hemorrhage not occurs at PNC or no visit at PNC
#         (hemorr_ld == 0 | is.na(hemorr_ld)) &
#           (hemorr_pnc == 0 | is.na(hemorr_pnc)) ~ 0,
#         # missing
#         TRUE ~ -5
#       )) %>%
#     select(-c(matches("M\\d{2}_")))
#  
#   #******preclampsia/eclampsia
#   # preeclampsia/eclampsia at ANC
#   var_clampsia_anc <- matData %>%
#     # # keep subjects with M04 visits
#     # filter(!is.na(M04_VISIT_TOT)) %>%
#     select(all_of(matID), M04_VISIT_TOT, matches("^M04_(HPD_MHOCCUR|HPD_MHTERM)_V\\d+")) %>%
#     pivot_longer(!c(all_of(matID), M04_VISIT_TOT),
#                  names_to = c(".value", "anc_order"),
#                  names_pattern = "^M04_(.+)_V(\\d+)"
#     ) %>%
#     mutate(
#       # for each ANC visit
#       mat_clampsia_anc_t1 = case_when(
#         # eclampsia at ANC
#         HPD_MHOCCUR == 1 & HPD_MHTERM == 2 ~ 1,
#         # no eclampsia
#         HPD_MHOCCUR == 0 | HPD_MHOCCUR == 1 & HPD_MHTERM %in% c(1,3) ~ 0,
#         # don't know
#         HPD_MHOCCUR == 99 | HPD_MHOCCUR == 1 & HPD_MHTERM == 99 ~ 99,
#         # missing
#         TRUE ~ -5
#       )) %>%
#     mutate(anc_order = as.numeric(anc_order)) %>%
#     group_by(SCRNID, MOMID, PREGID, SITE) %>% 
#     # arrange(ANC_ORDER) %>%
#     # slice(c(1:M04_VISIT_TOT[1])) %>%
#     summarise(
#       mat_clampsia_anc = case_when(
#         # clampsia if any clampsia at ANC visits
#         1 %in% mat_clampsia_anc_t1 ~ 1,
#         # no clampsia is none at ANC visits
#         all(mat_clampsia_anc_t1 == 0) ~ 0, ###*** 0 if all reported 0
#         # don't know if any is don't know
#         any(mat_clampsia_anc_t1 == 99) ~ 99,
#         # missing
#         TRUE ~ -5
#       )) %>%
#     ungroup() %>%
#     select(-c(matches("M\\d{2}_")))
#   
#   # clampsia at Labor
#   var_clampsia_ld <- matData %>%
#     # filter(!is.na(M09_FORMCOMPLDAT_MNH09)) %>%
#     select(all_of(matID), M09_HDP_HTN_MHOCCUR_3,
#            matches("M09_PREECLAMPSIA_CEOCCUR_")) %>%
#     mutate(
#       mat_clampsia_ld = case_when(
#         # clampsia
#         M09_HDP_HTN_MHOCCUR_3 == 1 ~ 1,
#         # no clampsia
#         M09_HDP_HTN_MHOCCUR_3 == 0 ~ 0,
#         TRUE ~ -5
#       ),
#       mat_clampsia_ld_severe = case_when(
#         #severe clampsia
#         mat_clampsia_ld == 1 & M09_PREECLAMPSIA_CEOCCUR_77 == 0 ~ 1, 
#         mat_clampsia_ld == 1 & M09_PREECLAMPSIA_CEOCCUR_77 == 1 ~ 0, 
#         mat_clampsia_ld != 1 ~ mat_clampsia_ld
#       )
#     ) %>%
#     select(-c(matches("M\\d{2}_")))
#   
#   # clampsia at PNC
#   var_clampsia_pnc <- matData %>%
#     # filter(!is.na(M12_VISIT_TOT)) %>%
#     select(all_of(matID), M12_VISIT_TOT, matches("^M12_BIRTH_COMPL_MHTERM_2_V\\d+")) %>%
#     pivot_longer(!c(all_of(matID), M12_VISIT_TOT),
#                  names_to = c(".value", "pnc_order"),
#                  names_pattern = "^M12_(.+)_V(\\d+)"
#     ) %>%
#     mutate(
#       # at each PNC visit
#       mat_clampsia_pnc_t1 = case_when(
#         BIRTH_COMPL_MHTERM_2 == 1 ~ 1,
#         BIRTH_COMPL_MHTERM_2 == 0 ~ 0,
#         TRUE ~ -5
#       )) %>%
#     mutate(pnc_order = as.numeric(pnc_order)) %>%
#     group_by(SCRNID, MOMID, PREGID, SITE) %>%
#     # arrange(PNC_ORDER) %>%
#     # # !!! we only keep the first total number of visits (must to distinguish NO visit and MISSING visit)
#     # slice(c(1:M12_VISIT_TOT[1])) %>%
#     summarise(
#       mat_clampsia_pnc = case_when(
#         1 %in% mat_clampsia_pnc_t1 ~ 1,
#         all(mat_clampsia_pnc_t1 == 0) ~ 0, ###*** 0 if all reported 0
#         TRUE ~ -5
#       )) %>%
#     ungroup() %>%
#     select(-c(matches("M\\d{2}_")))
#   
#   # combine
#   vars_clampsia <- var_clampsia_anc %>%
#     left_join(var_clampsia_ld, by = c(all_of(matID))) %>%
#     full_join(var_clampsia_pnc, by = c(all_of(matID))) %>%
#     distinct() %>%
#     mutate(
#       mat_clampsia = case_when(
#         # clampisa if occurs at ANC, LD or PNC
#         mat_clampsia_anc == 1 | mat_clampsia_ld == 1 | mat_clampsia_pnc == 1 ~ 1,
#         # no clampsia if
#         # if not occurs at ANC or no ANC visit yet AND
#         # if not occurs at LD or no LD visit yet AND
#         # if not occurs at PNC or no PNC visit yet
#         (mat_clampsia_anc == 0 | mat_clampsia_anc == -5) &
#           (mat_clampsia_ld == 0 | mat_clampsia_ld == -5) &
#           (mat_clampsia_pnc == 0 | mat_clampsia_pnc == -5) ~ 0,
#         # missing
#         TRUE ~ -5
#       )
#     ) %>%
#     select(-c(matches("M\\d{2}_")))
  
# 
#   spesis_severe_sysinfect = case_when(
#     
#   ), 
#   reptured_uterus = case_when(
#     
#   ), 
#   severe_complication = case_when(
#     pos_hemorr == 1 | preeclampsia == 1 | eclampsia == 1 |
#       spesis_sysinfect == 1 | reptured_uterus == 1 ~ 1, 
#     pos_hemorr == 0 & preeclampsia == 0 & eclampsia == 0 &
#       spesis_sysinfect == 0 & reptured_uterus == 0 ~ 0,
#     TRUE ~ -5
#   )
#   ) %>% 
#   select(all_of(matID), starts_with("mat_death"),) # near_miss 




#*severe complication

#******Perinatal depression

#******Postpartum anemia (6 weeks)
#*******Postpartum anemia (6 months)
#*
#******************************primary cohort fetal/neonatal outcomes****************
#******preterm birth
#******small-for-gestational-age
#******stillbirth

#******************************secondary cohort maternal outcomes****************
#******Late maternal mortality

#******Preeclampsia/eclampsia
#******Postpartum hemorrhage
#******Preterm premature rupture of membranes
#******Maternal sepsis
#******Maternal fatigue

#******************************secondary cohort fetal/neonatal outcomes****************
#******Low birth weight
#******Hyperbilirubinemia
#******Possible severe bacterial infection/sepsis
#******Birth asphyxia
#******Infant neurodevelopment

#******************************Hemoglobin in ANC & PNC visits***********************
# POC at ANC
# vars_hemoglobin_pocANC <- matData %>% 
#   # only keep participants with at least one M06 visit
#   filter(!is.na(M06_VISIT_TOT)) %>% 
#   select(all_of(matID), M06_VISIT_TOT, 
#          matches("^M06_(DIAG_VSDAT|HB_POC_LBORRES)_V\\d+")) %>% 
#   left_join(matOutcome %>% select(MOMID, PREGID, bestedd), 
#             by = c("MOMID", "PREGID")) %>% 
#   # wide to long
#   pivot_longer(!c(MOMID, PREGID, M06_VISIT_TOT, bestedd),
#                names_to = c(".value", "VISIT_ORDER"),
#                names_pattern = "^M06_(.+)_V(\\d+)"
#   ) %>% 
#   mutate(VISIT_ORDER = as.numeric(VISIT_ORDER)) %>% 
#   group_by(MOMID, PREGID) %>%
#   arrange(VISIT_ORDER) %>% 
#   slice(c(1:M06_VISIT_TOT[1])) %>% 
#   # keep ANC visits
#   filter(VISIT_ORDER >= 1 & VISIT_ORDER <= 5) %>% 
#   mutate(
#     # gestational age in weeks at the visit time
#     GA_WKS = as.numeric(dmy(DIAG_VSDAT) - (bestedd - 280))/7,
#     # set outliers when GA <=0 or GA >= 42
#     GA_OUTLIER = case_when(GA_WKS <= 0 | GA_WKS >= 42 ~ 1,
#                            !is.na(GA_WKS) ~ 0,
#                            TRUE ~ 99)) %>% 
#   ungroup() %>% 
#   # select HB measurements and GA
#   select(all_of(matID), HB_POC_LBORRES, GA_WKS, GA_OUTLIER) %>% 
#   left_join(vars_miscarrSB, by = c("MOMID", "PREGID")) %>% 
#   distinct()
# 
# # CBC at ANC
# vars_hemoglobin_cbcANC <- matData %>% 
#   # only keep participants with at least one M06 visit
#   filter(!is.na(M08_VISIT_TOT)) %>% 
#   select(all_of(matID), M08_VISIT_TOT, 
#          matches("^M08_(CBC_LBTSTDAT|CBC_HB_LBORRES)_V\\d+")) %>% 
#   left_join(matOutcome %>% select(MOMID, PREGID, bestedd), 
#             by = c("MOMID", "PREGID")) %>% 
#   # wide to long
#   pivot_longer(!c(MOMID, PREGID, M08_VISIT_TOT, bestedd),
#                names_to = c(".value", "VISIT_ORDER"),
#                names_pattern = "^M08_(.+)_V(\\d+)"
#   ) %>% 
#   mutate(VISIT_ORDER = as.numeric(VISIT_ORDER)) %>% 
#   group_by(MOMID, PREGID) %>%
#   arrange(VISIT_ORDER) %>% 
#   slice(c(1:M08_VISIT_TOT[1])) %>% 
#   # keep ANC visits
#   filter(VISIT_ORDER >= 1 & VISIT_ORDER <= 5) %>% 
#   mutate(
#     # gestational age in weeks at the visit time
#     GA_WKS = as.numeric(CBC_LBTSTDAT - (bestedd - 280))/7,
#     # set outliers when GA <=0 or GA >= 42
#     GA_OUTLIER = case_when(GA_WKS <= 0 | GA_WKS >= 42 ~ 1,
#                            !is.na(GA_WKS) ~ 0,
#                            TRUE ~ 99)) %>% 
#   ungroup() %>% 
#   # select HB measurements and GA
#   select(all_of(matID), CBC_HB_LBORRES, GA_WKS, GA_OUTLIER) %>% 
#   left_join(vars_miscarrSB, by = c("MOMID", "PREGID")) %>% 
#   distinct()
# 
# #******************************Hemoglobin in ANC & PNC visits***********************
# # POC hemoglobin at ANC & PNC
# vars_hemoglobin_poc <- matData %>% 
#   # only keep participants with at least one M06 visit
#   filter(!is.na(M06_VISIT_TOT)) %>% 
#   select(all_of(matID), M06_VISIT_TOT, 
#          matches("^M06_(DIAG_VSDAT|HB_POC_LBORRES)_V\\d+")) %>% 
#   left_join(matOutcome %>% select(all_of(matID), bestedd), 
#             by = c(all_of(matID))) %>% 
#   # wide to long
#   pivot_longer(!c(all_of(matID), M06_VISIT_TOT, bestedd),
#                names_to = c(".value", "visit_order"),
#                names_pattern = "^M06_(.+)_V(\\d+)"
#   ) %>% 
#   mutate(visit_oder = as.numeric(visit_order)) %>% 
#   group_by(SCRNID, MOMID, PREGID, SITE) %>% 
#   arrange(visit_order) %>% 
#   slice(c(1:M06_VISIT_TOT[1])) %>% 
#   mutate(
#     # gestational age in weeks at the visit time
#     ga_wks = as.numeric(dmy(DIAG_VSDAT) - (bestedd - 280))/7
#   ) %>% 
#   ungroup() %>% 
#   # select HB measurements and GA
#   select(all_of(matID), HB_POC_LBORRES, ga_wks) %>% #HB_POC_LBPERF
#   rename(hb = HB_POC_LBORRES) %>% 
#   mutate(hbtype = "POC")
# 
# # CBC hemoglobin at ANC & PNC
# vars_hemoglobin_cbc <- matData %>% 
#   # only keep participants with at least one M06 visit
#   filter(!is.na(M08_VISIT_TOT)) %>% 
#   select(all_of(matID), M08_VISIT_TOT, 
#          matches("^M08_(CBC_LBTSTDAT|CBC_HB_LBORRES)_V\\d+")) %>% 
#   left_join(matOutcome %>% select(all_of(matID), bestedd), 
#             by = c(all_of(matID))) %>% 
#   # wide to long
#   pivot_longer(!c(all_of(matID), M08_VISIT_TOT, bestedd),
#                names_to = c(".value", "visit_order"),
#                names_pattern = "^M08_(.+)_V(\\d+)"
#   ) %>% 
#   mutate(visit_order = as.numeric(visit_order)) %>% 
#   group_by(SCRNID, MOMID, PREGID, SITE) %>%
#   arrange(visit_order) %>% 
#   slice(c(1:M08_VISIT_TOT[1])) %>% 
#   mutate(
#     # gestational age in weeks at the visit time
#     ga_wks = as.numeric(dmy(CBC_LBTSTDAT) - (bestedd - 280))/7
#   ) %>% 
#   ungroup() %>% 
#   # select HB measurements and GA
#   select(all_of(matID), CBC_HB_LBORRES, ga_wks) %>%  #CBC_LBPERF_1, 
#   rename(hb = CBC_HB_LBORRES) %>% 
#   mutate(hbtype = "CBC")
# 
# # merge
# hemoglobinData <- vars_hemoglobin_poc %>% 
#   rbind(vars_hemoglobin_cbc) %>% 
#   filter(!is.na(hb)) %>% 
#   distinct()


