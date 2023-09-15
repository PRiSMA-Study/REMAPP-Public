#*****************************************************************************
#*Sensitivity Analysis
#*Criteria update dates: 
#*2023-05-05
#*2023-05-12
#*2023-05-26
#*****************************************************************************

library(tidyverse)
library(labelled)  
library(dplyr)
library(naniar)

#load data created for healthy outcome 
#use healthyoutcomeFull if it's available
load("Z:/Merged Data/2023-05-26/vars_criteria_FullData_2023-05-26.RData") 
load("Z:/Merged Data/2023-05-26/MatData_Anc_Visits_2023-05-26.RData")

#*****************************************************************************
#*Create data, add label, derive variables
#*****************************************************************************
df_sensitive <- vars_criteria %>% 
  filter(SITE != "Zambia") %>% #temp remove zambia because all analysis has 0 obs
  left_join(MatData_Anc_Visits %>% 
              dplyr::select(MOMID, PREGID, starts_with("M04_PPH_COMP_RPORRES_")), 
            by = c("MOMID", "PREGID")) %>% #temp solution. Once healthy data udpated, this step can be removed
  distinct() %>% 
  mutate_at(c("CRIT_AGE", "CRIT_GA","CRIT_BMI", "CRIT_MUAC", "CRIT_HEIGHT",
              "CRIT_SINGLEPREG", "CRIT_BP", "CRIT_PREV_MISCARR", "CRIT_PREV_PRETERM_LBW", "CRIT_FETALDEATH",
              "CRIT_COMPLICATION", "CRIT_SMOKE", "CRIT_DRINK", "CRIT_HIV", "CRIT_CHRONIC", 
              "CRIT_MALARIA", "CRIT_HEPATITISB", "CRIT_HEPATITISC"), 
            ~replace_na(.,0)) %>% 
  mutate(
    #additional derived variables
    #preterm
    CRIT_PRETERM = case_when(
      M04_PRETERM_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PRETERM_RPORRES_1 == 0 ~ 1, 
      M04_PRETERM_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #low birth weight
    CRIT_LBW = case_when(
      M04_LOWBIRTHWT_RPORRES_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 | M04_LOWBIRTHWT_RPORRES_1 == 0 ~ 1,
      M04_LOWBIRTHWT_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #Unplanned Cesarean delivery or Gestational diabetes
    CRIT_UNPL_CESARIAN_DIAB = case_when(
      M04_UNPL_CESARIAN_PROCCUR_1 == 1 | M04_GEST_DIAB_RPORRES_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 | (M04_UNPL_CESARIAN_PROCCUR_1 == 0 & M04_GEST_DIAB_RPORRES_1 == 0) ~ 1,
      M04_UNPL_CESARIAN_PROCCUR_1 == 99 | M04_GEST_DIAB_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #Unplanned Cesarean 
    CRIT_UNPL_CESARIAN = case_when(
      M04_UNPL_CESARIAN_PROCCUR_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 | M04_UNPL_CESARIAN_PROCCUR_1 == 0 ~ 1,
      M04_UNPL_CESARIAN_PROCCUR_1 == 99 ~ 0,
      TRUE ~ 55 
    ),
    #Gestational diabetes
    CRIT_TEMP_DIAB = case_when(
      M04_GEST_DIAB_RPORRES_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 | M04_GEST_DIAB_RPORRES_1 == 0 ~ 1,
      M04_GEST_DIAB_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #Preeclampsia/eclampsia
    CRIT_TEMP_ECLAMPSIA = case_when(
      M04_PREECLAMPSIA_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PREECLAMPSIA_RPORRES_1 == 0 ~ 1,
      M04_PREECLAMPSIA_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #premature rupture of membranes
    CRIT_TEMP_PROM = case_when(
      M04_PREMATURE_RUPTURE_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PREMATURE_RUPTURE_RPORRES_1 == 0 ~ 1,
      M04_PREMATURE_RUPTURE_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #macrosomia (>4000g)
    CRIT_TEMP_MACROSOMIA = case_when(
      M04_MACROSOMIA_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_MACROSOMIA_RPORRES_1 == 0 ~ 1,
      M04_MACROSOMIA_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #oligohydramnios
    CRIT_TEMP_OLIGOHYDRAMNIOS = case_when(
      M04_OLIGOHYDRAMNIOS_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_OLIGOHYDRAMNIOS_RPORRES_1 == 0 ~ 1,
      M04_OLIGOHYDRAMNIOS_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #antepartum hemorrhage
    CRIT_TEMP_APH = case_when(
      M04_APH_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_APH_RPORRES_1 == 0 ~ 1,
      M04_APH_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #postpartum hemorrhage
    CRIT_TEMP_PPH = case_when(
      M04_PPH_RPORRES_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PPH_RPORRES_1 == 0 ~ 1,
      M04_PPH_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55
    ),
    #uterine atony contribute to postpartum hemorrhage
    CRIT_TEMP_UTERINE_ATONY = case_when(
      M04_PPH_COMP_RPORRES_1_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PPH_RPORRES_1 == 0 | M04_PPH_COMP_RPORRES_1_1 == 0 ~ 1,
      TRUE ~ 55
    ),
    #Retained placenta contribute to postpartum hemorrhage
    CRIT_TEMP_RETAINED_PLACENTA = case_when(
      M04_PPH_COMP_RPORRES_2_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PPH_RPORRES_1 == 0 | M04_PPH_COMP_RPORRES_2_1 == 0 ~ 1,
      TRUE ~ 55
    ),
    #Uterine rupture contribute to postpartum hemorrhage
    CRIT_TEMP_UTERINE_RUPTURE = case_when(
      M04_PPH_COMP_RPORRES_3_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PPH_RPORRES_1 == 0 | M04_PPH_COMP_RPORRES_3_1 == 0 ~ 1,
      TRUE ~ 55
    ),
    #Cervical/vaginal laceration contribute to postpartum hemorrhage, other 
    CRIT_TEMP_LACERATION = case_when(
      M04_PPH_COMP_RPORRES_4_1 == 1 | M04_PPH_COMP_RPORRES_88_1 == 1 ~ 0,
      M04_PH_PREV_RPORRES_1 == 0 | M04_PPH_RPORRES_1 == 0 | 
        M04_PPH_COMP_RPORRES_4_1 == 0 | M04_PPH_COMP_RPORRES_88_1 == 0 ~ 1,
      TRUE ~ 55
    ),
    #still birth
    CRIT_STILLBIRTH = case_when(
      M04_STILLBIRTH_RPORRES_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 | M04_STILLBIRTH_RPORRES_1 == 0 ~ 1,
      M04_STILLBIRTH_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55),
    #new height criteria
    CRIT_HEIGHT_NEW = case_when(
      M05_HEIGHT_PERES_1 < 150 ~ 0, 
      M05_HEIGHT_PERES_1 >= 150 ~ 1,
      TRUE ~ 55),
    #new complication criteria
    CRIT_COMPLICATION_NEW = case_when(
      M04_PREECLAMPSIA_RPORRES_1== 1 |  
        M04_GEST_DIAB_RPORRES_1 == 1 | M04_PREMATURE_RUPTURE_RPORRES_1 == 1 |  
        M04_MACROSOMIA_RPORRES_1 == 1 | M04_OLIGOHYDRAMNIOS_RPORRES_1 == 1 |  
        M04_APH_RPORRES_1 == 1 |  M04_PPH_RPORRES_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 |
        (M04_PREECLAMPSIA_RPORRES_1 == 0 & 
           M04_GEST_DIAB_RPORRES_1 == 0 & M04_PREMATURE_RUPTURE_RPORRES_1 == 0 & 
           M04_MACROSOMIA_RPORRES_1 == 0 & M04_OLIGOHYDRAMNIOS_RPORRES_1 == 0 & 
           M04_APH_RPORRES_1 == 0 & M04_PPH_RPORRES_1 == 0) ~ 1,
      M04_PREECLAMPSIA_RPORRES_1 == 99 |  
        M04_GEST_DIAB_RPORRES_1 == 99 | M04_PREMATURE_RUPTURE_RPORRES_1== 99 |    
        M04_MACROSOMIA_RPORRES_1== 99 | M04_OLIGOHYDRAMNIOS_RPORRES_1 == 99 |    
        M04_APH_RPORRES_1 == 99 | M04_PPH_RPORRES_1 == 99 ~ 0,
      TRUE ~ 55),
    # unplanned cesarean delivery
    CRIT_UNPL_CESARIAN = case_when(
      M04_UNPL_CESARIAN_PROCCUR_1 == 1 ~ 0, 
      M04_PH_PREV_RPORRES_1 == 0 | M04_UNPL_CESARIAN_PROCCUR_1 == 0 ~ 1,
      M04_UNPL_CESARIAN_PROCCUR_1 == 99 ~ 0,
      TRUE ~ 55),
#all 
    eligibleall = case_when(CRIT_AGE == 1 &
                       CRIT_GA == 1 &
                       CRIT_BMI == 1 & 
                       CRIT_MUAC == 1 &
                       CRIT_HEIGHT == 1 &
                       CRIT_SINGLEPREG == 1 &
                       CRIT_BP == 1 &
                       CRIT_PREV_MISCARR == 1 &
                       CRIT_PREV_PRETERM_LBW == 1 &
                       CRIT_STILLBIRTH == 1 &
                       CRIT_COMPLICATION_NEW == 1 &
                       CRIT_UNPL_CESARIAN == 1 &
                       CRIT_SMOKE == 1 &
                       ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                       CRIT_HIV == 1 &
                       CRIT_CHRONIC == 1 &
                       CRIT_MALARIA == 1 &
                       CRIT_HEPATITISB == 1 &
                       CRIT_HEPATITISC == 1 &
                       CRIT_HEMOGLOBINOPATHIES == 1 &
                       CRIT_IRON == 1 &
                       CRIT_INFLAM == 1 ~ 1, 
                       TRUE ~ 0),
    #1.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM 
    eligible1 = case_when(CRIT_AGE == 1 &
                           CRIT_GA == 1 &
                           CRIT_BMI == 1 & 
                           CRIT_MUAC == 1 &
                           CRIT_HEIGHT == 1 &
                           CRIT_SINGLEPREG == 1 &
                           CRIT_BP == 1 &
                           CRIT_PREV_MISCARR == 1 &
                           CRIT_PREV_PRETERM_LBW == 1 & #confirmed results are same
                           # CRIT_PRETERM == 1 &
                           # CRIT_LBW == 1 &
                           CRIT_STILLBIRTH == 1 &
                           CRIT_COMPLICATION_NEW == 1 &
                           CRIT_UNPL_CESARIAN == 1 &
                           CRIT_SMOKE == 1 &
                           ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                           CRIT_HIV == 1 &
                           CRIT_CHRONIC == 1 &
                           CRIT_MALARIA == 1 &
                           CRIT_HEPATITISB == 1 &
                           CRIT_HEPATITISC == 1 ~ 1, 
                       TRUE ~ 0),
#2.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM and CRIT_PREV_MISCARR.
    eligible2 = ifelse(CRIT_AGE == 1 &
                         CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         # CRIT_PREV_MISCARR == 1 &
                         CRIT_PREV_PRETERM_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
    #3.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM and CRIT_PREV_PRETERM_LBW
    eligible3 = ifelse(CRIT_AGE == 1 &
                         CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         # CRIT_PREV_PRETERM_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
    #4.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM and CRIT_STILLBIRTH
    eligible4 = ifelse(CRIT_AGE == 1 &
                         CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         CRIT_PREV_PRETERM_LBW == 1 &
                         # CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
    #5.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM and CRIT_COMPLICATION_NEW.
    eligible5 = ifelse(CRIT_AGE == 1 &
                         CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         CRIT_PREV_PRETERM_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         # CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
    #6.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM and CRIT_PREV_MISCARR, CRIT_PREV_PRETERM_LBW, CRIT_STILLBIRTH, CRIT_COMPLICATION_NEW.
    eligible6 = ifelse(CRIT_AGE == 1 &
                         CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         # CRIT_PREV_MISCARR == 1 &
                         # CRIT_PREV_PRETERM_LBW == 1 &
                         # CRIT_STILLBIRTH == 1 &
                         # CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
    #7.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM and CRIT_GA.
    eligible7 = ifelse(CRIT_AGE == 1 &
                         # CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         CRIT_PREV_PRETERM_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
    #8.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM CRIT_PREV_MISCARR, 
    #CRIT_PREV_PRETERM_LBW, CRIT_STILLBIRTH, CRIT_COMPLICATION_NEW and CRIT_GA.
    eligible8 = ifelse(CRIT_AGE == 1 &
                         # CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         # CRIT_PREV_MISCARR == 1 &
                         # CRIT_PREV_PRETERM_LBW == 1 &
                         # CRIT_STILLBIRTH == 1 &
                         # CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1,
                           1, 0 ),
#9. exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and revise height criteria to >=150cm
    eligible9 = ifelse(CRIT_AGE == 1 &
                         # CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         # CRIT_HEIGHT == 1 &
                         M05_HEIGHT_PERES_1 >=150 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         CRIT_PREV_PRETERM_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1,1, 0 ),
#10. exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON,CRIT_INFLAM, CRIT_PREV_MISCARR, 
#CRIT_PREV_PRETERM_LBW, CRIT_STILLBIRTH,CRIT_COMPLICATION_NEW, CRIT_GA and revise height criteria to >=150cm
    eligible10 = ifelse(CRIT_AGE == 1 &
                         # CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         # CRIT_HEIGHT == 1 &
                         M05_HEIGHT_PERES_1 >=150 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         # CRIT_PREV_MISCARR == 1 &
                         # CRIT_PREV_PRETERM_LBW == 1 &
                         # CRIT_STILLBIRTH == 1 &
                         # CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1,1,0),
#11.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_PREV_MISCARR, CRIT_GA
    eligible11 = ifelse(CRIT_AGE == 1 &
                          # CRIT_GA == 1 &
                          CRIT_BMI == 1 & 
                          CRIT_MUAC == 1 &
                          CRIT_HEIGHT == 1 &
                          # M05_HEIGHT_PERES_1 >=150 &
                          CRIT_SINGLEPREG == 1 &
                          CRIT_BP == 1 &
                          # CRIT_PREV_MISCARR == 1 &
                          CRIT_PREV_PRETERM_LBW == 1 &
                          CRIT_STILLBIRTH == 1 &
                          CRIT_COMPLICATION_NEW == 1 &
                          CRIT_UNPL_CESARIAN == 1 &
                          CRIT_SMOKE == 1 &
                          ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                          CRIT_HIV == 1 &
                          CRIT_CHRONIC == 1 &
                          CRIT_MALARIA == 1 &
                          CRIT_HEPATITISB == 1 &
                          CRIT_HEPATITISC == 1,1, 0 ),
#12.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_PREV_PRETERM_LBW, CRIT_GA
    eligible12 = ifelse(CRIT_AGE == 1 &
                          # CRIT_GA == 1 &
                          CRIT_BMI == 1 & 
                          CRIT_MUAC == 1 &
                          CRIT_HEIGHT == 1 &
                          # M05_HEIGHT_PERES_1 >=150 &
                          CRIT_SINGLEPREG == 1 &
                          CRIT_BP == 1 &
                          CRIT_PREV_MISCARR == 1 &
                          # CRIT_PREV_PRETERM_LBW == 1 &
                          CRIT_STILLBIRTH == 1 &
                          CRIT_COMPLICATION_NEW == 1 &
                          CRIT_UNPL_CESARIAN == 1 &
                          CRIT_SMOKE == 1 &
                          ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                          CRIT_HIV == 1 &
                          CRIT_CHRONIC == 1 &
                          CRIT_MALARIA == 1 &
                          CRIT_HEPATITISB == 1 &
                          CRIT_HEPATITISC == 1,1, 0 ),
#13.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_STILLBIRTH, CRIT_GA
    eligible13 = ifelse(CRIT_AGE == 1 &
                          # CRIT_GA == 1 &
                          CRIT_BMI == 1 & 
                          CRIT_MUAC == 1 &
                          CRIT_HEIGHT == 1 &
                          # M05_HEIGHT_PERES_1 >=150 &
                          CRIT_SINGLEPREG == 1 &
                          CRIT_BP == 1 &
                          CRIT_PREV_MISCARR == 1 &
                          CRIT_PREV_PRETERM_LBW == 1 &
                          # CRIT_STILLBIRTH == 1 &
                          CRIT_COMPLICATION_NEW == 1 &
                          CRIT_UNPL_CESARIAN == 1 &
                          CRIT_SMOKE == 1 &
                          ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                          CRIT_HIV == 1 &
                          CRIT_CHRONIC == 1 &
                          CRIT_MALARIA == 1 &
                          CRIT_HEPATITISB == 1 &
                          CRIT_HEPATITISC == 1,1, 0 ),
#14.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_COMPLICATION_NEW, CRIT_GA
    eligible14 = ifelse(CRIT_AGE == 1 &
                          # CRIT_GA == 1 &
                          CRIT_BMI == 1 & 
                          CRIT_MUAC == 1 &
                          CRIT_HEIGHT == 1 &
                          # M05_HEIGHT_PERES_1 >=150 &
                          CRIT_SINGLEPREG == 1 &
                          CRIT_BP == 1 &
                          CRIT_PREV_MISCARR == 1 &
                          CRIT_PREV_PRETERM_LBW == 1 &
                          CRIT_STILLBIRTH == 1 &
                          # CRIT_COMPLICATION_NEW == 1 &
                          CRIT_UNPL_CESARIAN == 1 &
                          CRIT_SMOKE == 1 &
                          ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                          CRIT_HIV == 1 &
                          CRIT_CHRONIC == 1 &
                          CRIT_MALARIA == 1 &
                          CRIT_HEPATITISB == 1 &
                          CRIT_HEPATITISC == 1,1, 0 ),   
#15.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and
#revise height criteria to >=150cm
    eligible15 = ifelse(CRIT_AGE == 1 &
                          # CRIT_GA == 1 &
                          CRIT_BMI == 1 & 
                          CRIT_MUAC == 1 &
                          # CRIT_HEIGHT == 1 &
                          M05_HEIGHT_PERES_1 >=150 &
                          CRIT_SINGLEPREG == 1 &
                          CRIT_BP == 1 &
                          CRIT_PREV_MISCARR == 1 &
                          CRIT_PREV_PRETERM_LBW == 1 &
                          CRIT_STILLBIRTH == 1 &
                          CRIT_COMPLICATION_NEW == 1 &
                          CRIT_UNPL_CESARIAN == 1 &
                          CRIT_SMOKE == 1 &
                          ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                          CRIT_HIV == 1 &
                          CRIT_CHRONIC == 1 &
                          CRIT_MALARIA == 1 &
                          CRIT_HEPATITISB == 1 &
                          CRIT_HEPATITISC == 1,1, 0 ),
    #16.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and CRIT_PRETERM
    eligible16 = ifelse(CRIT_AGE == 1 &
                         # CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         # CRIT_PREV_PRETERM_LBW == 1 &
                         # CRIT_PRETERM == 1 &
                         CRIT_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
#17.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and CRIT_LBW
    eligible17 = ifelse(CRIT_AGE == 1 &
                         # CRIT_GA == 1 &
                         CRIT_BMI == 1 & 
                         CRIT_MUAC == 1 &
                         CRIT_HEIGHT == 1 &
                         CRIT_SINGLEPREG == 1 &
                         CRIT_BP == 1 &
                         CRIT_PREV_MISCARR == 1 &
                         # CRIT_PREV_PRETERM_LBW == 1 &
                         CRIT_PRETERM == 1 &
                         # CRIT_LBW == 1 &
                         CRIT_STILLBIRTH == 1 &
                         CRIT_COMPLICATION_NEW == 1 &
                         CRIT_UNPL_CESARIAN == 1 &
                         CRIT_SMOKE == 1 &
                         ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                         CRIT_HIV == 1 &
                         CRIT_CHRONIC == 1 &
                         CRIT_MALARIA == 1 &
                         CRIT_HEPATITISB == 1 &
                         CRIT_HEPATITISC == 1, 1, 0 ),
#18.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and CRIT_STILLBIRTH
eligible18 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      # CRIT_STILLBIRTH == 1 &
                      CRIT_COMPLICATION_NEW == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#19.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and CRIT_PREV_MISCARR
eligible19 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      CRIT_STILLBIRTH == 1 &
                      CRIT_COMPLICATION_NEW == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#20.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and preterm, low birth weight, fetal death, miscarriage
eligible20 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      # CRIT_PRETERM == 1 &
                      # CRIT_LBW == 1 &
                      # CRIT_STILLBIRTH == 1 &
                      CRIT_COMPLICATION_NEW == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#21.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and preeclampsia/eclampsia
eligible21 = ifelse(CRIT_AGE == 1 &
                     # CRIT_GA == 1 &
                     CRIT_BMI == 1 & 
                     CRIT_MUAC == 1 &
                     CRIT_HEIGHT == 1 &
                     CRIT_SINGLEPREG == 1 &
                     CRIT_BP == 1 &
                     CRIT_PREV_MISCARR == 1 &
                     # CRIT_PREV_PRETERM_LBW == 1 &
                     CRIT_PRETERM == 1 &
                     CRIT_LBW == 1 &
                     CRIT_STILLBIRTH == 1 &
                     # CRIT_COMPLICATION_NEW == 1 &
                     CRIT_UNPL_CESARIAN_DIAB == 1 &
                     # CRIT_TEMP_ECLAMPSIA == 1 & 
                     CRIT_TEMP_PROM == 1 & 
                     CRIT_TEMP_MACROSOMIA == 1 & 
                     CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                     CRIT_TEMP_APH == 1 & 
                     CRIT_TEMP_PPH == 1 &
                     CRIT_SMOKE == 1 &
                     ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                     CRIT_HIV == 1 &
                     CRIT_CHRONIC == 1 &
                     CRIT_MALARIA == 1 &
                     CRIT_HEPATITISB == 1 &
                     CRIT_HEPATITISC == 1, 1, 0 ),
#22.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and premature rupture of membranes (before labor began)
eligible22 = ifelse(CRIT_AGE == 1 &
                     # CRIT_GA == 1 &
                     CRIT_BMI == 1 & 
                     CRIT_MUAC == 1 &
                     CRIT_HEIGHT == 1 &
                     CRIT_SINGLEPREG == 1 &
                     CRIT_BP == 1 &
                     CRIT_PREV_MISCARR == 1 &
                     # CRIT_PREV_PRETERM_LBW == 1 &
                     CRIT_PRETERM == 1 &
                     CRIT_LBW == 1 &
                     CRIT_STILLBIRTH == 1 &
                     # CRIT_COMPLICATION_NEW == 1 &
                     CRIT_UNPL_CESARIAN_DIAB == 1 &
                     CRIT_TEMP_ECLAMPSIA == 1 &
                     # CRIT_TEMP_PROM == 1 & 
                     CRIT_TEMP_MACROSOMIA == 1 & 
                     CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                     CRIT_TEMP_APH == 1 & 
                     CRIT_TEMP_PPH == 1 &
                     CRIT_SMOKE == 1 &
                     ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                     CRIT_HIV == 1 &
                     CRIT_CHRONIC == 1 &
                     CRIT_MALARIA == 1 &
                     CRIT_HEPATITISB == 1 &
                     CRIT_HEPATITISC == 1, 1, 0 ),
#23.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and macrosomia (>4000g)
eligible23 = ifelse(CRIT_AGE == 1 &
                     # CRIT_GA == 1 &
                     CRIT_BMI == 1 & 
                     CRIT_MUAC == 1 &
                     CRIT_HEIGHT == 1 &
                     CRIT_SINGLEPREG == 1 &
                     CRIT_BP == 1 &
                     CRIT_PREV_MISCARR == 1 &
                     # CRIT_PREV_PRETERM_LBW == 1 &
                     CRIT_PRETERM == 1 &
                     CRIT_LBW == 1 &
                     CRIT_STILLBIRTH == 1 &
                     # CRIT_COMPLICATION_NEW == 1 &
                     CRIT_UNPL_CESARIAN_DIAB == 1 &
                     CRIT_TEMP_ECLAMPSIA == 1 &
                     CRIT_TEMP_PROM == 1 & 
                     # CRIT_TEMP_MACROSOMIA == 1 & 
                     CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                     CRIT_TEMP_APH == 1 & 
                     CRIT_TEMP_PPH == 1 &
                     CRIT_SMOKE == 1 &
                     ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                     CRIT_HIV == 1 &
                     CRIT_CHRONIC == 1 &
                     CRIT_MALARIA == 1 &
                     CRIT_HEPATITISB == 1 &
                     CRIT_HEPATITISC == 1, 1, 0 ),
#24.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and oligohydramnios
eligible24 = ifelse(CRIT_AGE == 1 &
                     # CRIT_GA == 1 &
                     CRIT_BMI == 1 & 
                     CRIT_MUAC == 1 &
                     CRIT_HEIGHT == 1 &
                     CRIT_SINGLEPREG == 1 &
                     CRIT_BP == 1 &
                     CRIT_PREV_MISCARR == 1 &
                     # CRIT_PREV_PRETERM_LBW == 1 &
                     CRIT_PRETERM == 1 &
                     CRIT_LBW == 1 &
                     CRIT_STILLBIRTH == 1 &
                     # CRIT_COMPLICATION_NEW == 1 &
                     CRIT_UNPL_CESARIAN_DIAB == 1 &
                     CRIT_TEMP_ECLAMPSIA == 1 &
                     CRIT_TEMP_PROM == 1 &
                     CRIT_TEMP_MACROSOMIA == 1 & 
                     # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                     CRIT_TEMP_APH == 1 & 
                     CRIT_TEMP_PPH == 1 &
                     CRIT_SMOKE == 1 &
                     ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                     CRIT_HIV == 1 &
                     CRIT_CHRONIC == 1 &
                     CRIT_MALARIA == 1 &
                     CRIT_HEPATITISB == 1 &
                     CRIT_HEPATITISC == 1, 1, 0 ),
#25.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and antepartum hemorrhage
eligible25 = ifelse(CRIT_AGE == 1 &
                     # CRIT_GA == 1 &
                     CRIT_BMI == 1 & 
                     CRIT_MUAC == 1 &
                     CRIT_HEIGHT == 1 &
                     CRIT_SINGLEPREG == 1 &
                     CRIT_BP == 1 &
                     CRIT_PREV_MISCARR == 1 &
                     # CRIT_PREV_PRETERM_LBW == 1 &
                     CRIT_PRETERM == 1 &
                     CRIT_LBW == 1 &
                     CRIT_STILLBIRTH == 1 &
                     # CRIT_COMPLICATION_NEW == 1 &
                     CRIT_UNPL_CESARIAN_DIAB == 1 &
                     CRIT_TEMP_ECLAMPSIA == 1 &
                     CRIT_TEMP_PROM == 1 & 
                     CRIT_TEMP_MACROSOMIA == 1 & 
                     CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                     # CRIT_TEMP_APH == 1 & 
                     CRIT_TEMP_PPH == 1 &
                     CRIT_SMOKE == 1 &
                     ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                     CRIT_HIV == 1 &
                     CRIT_CHRONIC == 1 &
                     CRIT_MALARIA == 1 &
                     CRIT_HEPATITISB == 1 &
                     CRIT_HEPATITISC == 1, 1, 0 ),
#26.exclude CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA and postpartum hemorrhage
eligible26 = ifelse(CRIT_AGE == 1 &
                     # CRIT_GA == 1 &
                     CRIT_BMI == 1 & 
                     CRIT_MUAC == 1 &
                     CRIT_HEIGHT == 1 &
                     CRIT_SINGLEPREG == 1 &
                     CRIT_BP == 1 &
                     CRIT_PREV_MISCARR == 1 &
                     # CRIT_PREV_PRETERM_LBW == 1 &
                     CRIT_PRETERM == 1 &
                     CRIT_LBW == 1 &
                     CRIT_STILLBIRTH == 1 &
                     # CRIT_COMPLICATION_NEW == 1 &
                     CRIT_UNPL_CESARIAN_DIAB == 1 &
                     CRIT_TEMP_ECLAMPSIA == 1 &
                     CRIT_TEMP_PROM == 1 & 
                     CRIT_TEMP_MACROSOMIA == 1 & 
                     CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                     CRIT_TEMP_APH == 1 & 
                     # CRIT_TEMP_PPH == 1 &
                     CRIT_SMOKE == 1 &
                     ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                     CRIT_HIV == 1 &
                     CRIT_CHRONIC == 1 &
                     CRIT_MALARIA == 1 &
                     CRIT_HEPATITISB == 1 &
                     CRIT_HEPATITISC == 1, 1, 0 ),
#27.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
#CRIT_PRETERM, CRIT_LBW, CRIT_STILLBIRTH, CRIT_PREV_MISCARR and gestational diabetes, 
#preeclampsia/eclampsia, premature rupture of membranes (before labor began), 
#macrosomia (>4000g), oligohydramnios, antepartum hemorrhage, postpartum hemorrhage
eligible27 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      # CRIT_PRETERM == 1 &
                      # CRIT_LBW == 1 &
                      # CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      # CRIT_UNPL_CESARIAN_DIAB == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      # CRIT_TEMP_DIAB == 1 &
                      # CRIT_TEMP_ECLAMPSIA == 1 &
                      # CRIT_TEMP_PROM == 1 & 
                      # CRIT_TEMP_MACROSOMIA == 1 & 
                      # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                      # CRIT_TEMP_APH == 1 & 
                      # CRIT_TEMP_APH == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#28.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
#CRIT_PRETERM, CRIT_STILLBIRTH, CRIT_PREV_MISCARR and gestational diabetes, 
# preeclampsia/eclampsia, premature rupture of membranes (before labor began), 
# macrosomia (>4000g), oligohydramnios, antepartum hemorrhage, postpartum hemorrhage
eligible28 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      # CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      # CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      # CRIT_UNPL_CESARIAN_DIAB == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      # CRIT_TEMP_DIAB == 1 &
                      # CRIT_TEMP_ECLAMPSIA == 1 &
                      # CRIT_TEMP_PROM == 1 & 
                      # CRIT_TEMP_MACROSOMIA == 1 & 
                      # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                      # CRIT_TEMP_APH == 1 & 
                      # CRIT_TEMP_APH == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#29.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
# CRIT_STILLBIRTH, CRIT_PREV_MISCARR and gestational diabetes, preeclampsia/eclampsia, 
# premature rupture of membranes (before labor began), macrosomia (>4000g), oligohydramnios, 
# antepartum hemorrhage, postpartum hemorrhage
eligible29 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      # CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      # CRIT_UNPL_CESARIAN_DIAB == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      # CRIT_TEMP_DIAB == 1 &
                      # CRIT_TEMP_ECLAMPSIA == 1 &
                      # CRIT_TEMP_PROM == 1 & 
                      # CRIT_TEMP_MACROSOMIA == 1 & 
                      # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                      # CRIT_TEMP_APH == 1 & 
                      # CRIT_TEMP_APH == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#30.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
# CRIT_PRETERM, CRIT_PREV_MISCARR and gestational diabetes, preeclampsia/eclampsia, 
# premature rupture of membranes (before labor began), macrosomia (>4000g), 
# oligohydramnios, antepartum hemorrhage, postpartum hemorrhage
eligible30 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      # CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      # CRIT_UNPL_CESARIAN_DIAB == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      # CRIT_TEMP_DIAB == 1 &
                      # CRIT_TEMP_ECLAMPSIA == 1 &
                      # CRIT_TEMP_PROM == 1 & 
                      # CRIT_TEMP_MACROSOMIA == 1 & 
                      # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                      # CRIT_TEMP_APH == 1 & 
                      # CRIT_TEMP_APH == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#31.eexclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
# CRIT_PREV_MISCARR and gestational diabetes, preeclampsia/eclampsia, 
# premature rupture of membranes (before labor began), macrosomia (>4000g), 
# oligohydramnios, antepartum hemorrhage, postpartum hemorrhage
eligible31 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      CRIT_HEIGHT == 1 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      # CRIT_UNPL_CESARIAN_DIAB == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      # CRIT_TEMP_DIAB == 1 &
                      # CRIT_TEMP_ECLAMPSIA == 1 &
                      # CRIT_TEMP_PROM == 1 & 
                      # CRIT_TEMP_MACROSOMIA == 1 & 
                      # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                      # CRIT_TEMP_APH == 1 & 
                      # CRIT_TEMP_APH == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#32.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
# CRIT_PREV_MISCARR and gestational diabetes, preeclampsia/eclampsia, 
# premature rupture of membranes (before labor began), macrosomia (>4000g), 
# oligohydramnios, antepartum hemorrhage, postpartum hemorrhage and revise height criteria to >=150cm
eligible32 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      # CRIT_HEIGHT == 1 &
                      M05_HEIGHT_PERES_1 >=150 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PREV_PRETERM_LBW == 1 &
                      CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      # CRIT_UNPL_CESARIAN_DIAB == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      # CRIT_TEMP_DIAB == 1 &
                      # CRIT_TEMP_ECLAMPSIA == 1 &
                      # CRIT_TEMP_PROM == 1 & 
                      # CRIT_TEMP_MACROSOMIA == 1 & 
                      # CRIT_TEMP_OLIGOHYDRAMNIOS == 1 &
                      # CRIT_TEMP_APH == 1 & 
                      # CRIT_TEMP_APH == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
#33.exclude criteria CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, CRIT_INFLAM, CRIT_GA, 
# CRIT_PREV_MISCARR, CRIT_PRETERM, CRIT_COMPLICATION_NEW and revise height criteria to >=150cm
eligible33 = ifelse(CRIT_AGE == 1 &
                      # CRIT_GA == 1 &
                      CRIT_BMI == 1 & 
                      CRIT_MUAC == 1 &
                      # CRIT_HEIGHT == 1 &
                      M05_HEIGHT_PERES_1 >=150 &
                      CRIT_SINGLEPREG == 1 &
                      CRIT_BP == 1 &
                      # CRIT_PREV_MISCARR == 1 &
                      # CRIT_PRETERM == 1 &
                      CRIT_LBW == 1 &
                      CRIT_STILLBIRTH == 1 &
                      # CRIT_COMPLICATION_NEW == 1 &
                      CRIT_UNPL_CESARIAN == 1 &
                      CRIT_SMOKE == 1 &
                      ((SITE != "Pakistan" & CRIT_DRINK == 1) | (SITE == "Pakistan")) &
                      CRIT_HIV == 1 &
                      CRIT_CHRONIC == 1 &
                      CRIT_MALARIA == 1 &
                      CRIT_HEPATITISB == 1 &
                      CRIT_HEPATITISC == 1, 1, 0 ),
)

#*****************************************************************************
#*hb data
#*****************************************************************************
prep_hball <- df_sensitive %>% distinct() %>% 
  left_join(MatData_Anc_Visits %>% 
              dplyr::select(all_of(matID), EDD, matches("M08_(LBSTDAT|CBC_HB_LBORRES)_\\d+")),
            by = c("SCRNID", "MOMID", "PREGID", "SITE", "EDD")) %>% 
  distinct() %>% 
  dplyr::select(SCRNID, MOMID, PREGID, SITE, EDD, matches("M08_(LBSTDAT|CBC_HB_LBORRES)_\\d+"), matches("eligible\\d+")) %>% 
  filter(SITE == "Pakistan") %>% #because it's the only site with CBC hemoglobin data across gestational ages
  mutate(bestEdd = EDD) %>% 
  replace_with_na_all(condition = ~.== -7) %>%
  replace_with_na_all(condition = ~.== "1907-07-07") %>% 
  distinct()

#long data for eligible 1
df_hb_long_sensitive <- prep_hball %>% 
  pivot_longer(
    starts_with("M08_"), 
    names_to = c(".value", "visit_type"), 
    names_pattern = "^M\\d{2}_(.+)_(\\d+)",
    values_drop_na = TRUE
  ) %>% 
  rename(hb = CBC_HB_LBORRES) %>% 
  rename(collect_date = LBSTDAT) %>% 
  filter(!is.na(hb) & !is.na(bestEdd)) %>% 
  mutate(
    ga_wks = (280 - as.numeric(bestEdd - collect_date)) / 7
  ) %>% 
  filter(ga_wks >=10 & ga_wks <=50) %>% 
  filter(hb >= 5 & hb <= 18) 


#save data
# setwd("D:/Users/xyh/Documents/github/Sensitivity-Analysis/2023-05-26")
# setwd("D:/Users/xyh/Documents/github/Sensitivity-Analysis/Transfer")
save(df_sensitive, file= "derived_data/df_sensitive.rda")
save(df_hb_long_sensitive, file = "derived_data/df_hb_long_sensitive.rda")

