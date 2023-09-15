#*******************************************************************************
#*Frequency table for healthy cohort criteria
#*input: healthyOutcome.rda
#*Note: Add export to excel code if prefer. 
#*crosstab and histogram codes start at line 200. 
#*#*******************************************************************************
library(crosstable)
library(tidyverse)
library(ggplot2)
library(labelled)  
library(naniar)

# making a test change.

#use healthyOutcome_FullCohort data is it's available
load("Z:/Merged Data/2023-05-26/vars_criteria_FullData_2023-05-26.RData")
load("Z:/Merged Data/2023-05-26/MatData_Anc_Visits_2023-05-26.RData")
# load("Z:/Merged Data/2023-05-26/healthyOutcome_FullCohort_2023-05-26.RData")

#*#*******************************************************************************
#create data with labels for healthy cohort                       
healthy_labels <- c(
  SCRNID = "Screening ID", 
  MOMID = "Mother ID", 
  PREGID = "Pregnancy ID", 
  SITE = " Site", 
  ENROLL = "PRiSMA participants",
  M00_KNOWN_DOBYN_SCORRES_1 = "Mom's birthday known", 
  M00_BRTHDAT_1 = "Provide date of birth",
  M00_ESTIMATED_AGE_1 = "Estimate age in years", 
  M02_SCRN_OBSSTDAT_1 = "Screening date",
  M02_CONSENT_IEORRES_1 = "PRiSMA enrolled",
  M03_SMOKE_OECOCCUR_1 = "Smoked cigarettes/cigars/pipe",
  M03_CHEW_BNUT_OECOCCUR_1 = "Chewed betel nut in last month",
  M03_CHEW_OECOCCUR_1 = "Chewed tobacoo in last month",    
  M03_DRINK_OECOCCUR_1 = "Had any drink containing alcohol in last month", 
  M05_ANT_PEDAT_1 = "Date of anthropometric assessment", 
  M05_WEIGHT_PERES_1 = "Record weight: kg",
  M05_HEIGHT_PERES_1 = "Record height: cm",               
  M05_MUAC_PERES_1 = "Record Mid-Upper Arm Circumference (MUAC): cm",
  GA_US_DAYS_1 = "Gestational age",
  M04_PRETERM_RPORRES_1 = "Previous pre-term delivery",
  M04_PH_PREV_RPORRES_1 = "Has ever been pregnant",
  M04_PH_PREVN_RPORRES_1 = "Total number of previous pregnancies", 
  M04_PH_LIVE_RPORRES_1 = "Total number of previous pregnancies resulted in a live birth",
  M04_MISCARRIAGE_RPORRES_1 = "Experienced spontaneous miscarriage",
  M04_MISCARRIAGE_CT_RPORRES_1 = "Total number of spontaneous miscarriage", 
  M04_PH_OTH_RPORRES_1 = "Total number of pregnancy ended in a loss (stillbirth, miscarriage, or abortion)",
  M04_STILLBIRTH_RPORRES_1 = "Experienced stillbirth previously", 
  M04_LOWBIRTHWT_RPORRES_1 = "Experienced low birth weight previously",
  M04_MALARIA_EVER_MHOCCUR_1 = "Had malaria with current pregnancy",
  M04_CANCER_EVER_MHOCCUR_1 = "Ever been diagnosed with cancer",
  M04_KIDNEY_EVER_MHOCCUR_1 = "Ever been diagnosed with kidney disease",
  M04_CARDIAC_EVER_MHOCCUR_1 = "Ever been diagnosed with cardiac disease",
  M04_HIV_MHOCCUR_1 = "Had HIV with current pregnancy",
  M04_HIV_EVER_MHOCCUR_1 = "Ever been diagnosed with HIV",
  M04_UNPL_CESARIAN_PROCCUR_1 = "Previous unplanned Cesarean delivery", 
  M04_PREECLAMPSIA_RPORRES_1 = "Previous preeclampsia/eclampsia", 
  M04_GEST_DIAB_RPORRES_1 = "Previous gestational diabetes", 
  M04_PREMATURE_RUPTURE_RPORRES_1 = "Previous premature membranes rupture", 
  M04_MACROSOMIA_RPORRES_1 = "Previous macrosomia", 
  M04_OLIGOHYDRAMNIOS_RPORRES_1 = "Previous oligohydramnios", 
  M04_APH_RPORRES_1 = "Previous antepartum hemorrhage", 
  M04_PPH_RPORRES_1 = "Previous postpartum hemorrhage",
  M06_SINGLETON_PERES_1 = "Singleton pregnancy", 
  M06_BP_SYS_VSORRES_1_1 = "Record 1st BP measurement: (systolic) mmHg", 
  M06_BP_SYS_VSORRES_2_1 = "Record 2nd BP measurement: (systolic) mmHg", 
  M06_BP_SYS_VSORRES_3_1 = "Record 3rd BP measurement: (systolic) mmHg", 
  M06_BP_DIA_VSORRES_1_1 = "Record 1st BP measurement: (diastolic) mmHg", 
  M06_BP_DIA_VSORRES_2_1 = "Record 2nd BP measurement: (diastolic) mmHg", 
  M06_BP_DIA_VSORRES_3_1 = "Record 3rd BP measurement: (diastolic) mmHg", 
  M06_HBV_POC_LBORRES_1 = "HBV results",
  M06_HBV_POC_LBPERF_1 = "Point-of-care Hepatitis B (HBV) test performed", 
  M06_HCV_POC_LBORRES_1 = "HCV results",
  M06_HCV_POC_LBPERF_1 = "Point-of-care Hepatitis C (HCV) test performed", 
  M06_HIV_POC_LBORRES_1 = "HIV results",
  M06_HIV_POC_LBPERF_1 = "Point-of-care Hepatitis HIV test performed", 
  M08_MN_LBPERF_8_1 = "Ferritin test performed",
  M08_FERRITIN_LBORRES_1 = "Ferritin test results",
  M08_RBC_LBPERF_2_1 = "Hemoglobinopathies & Thalassemias test performed", 
  M08_RBC_THALA_LBORRES_1 = "Hemoglobinopathies & Thalassemias results", 
  M08_RBC_LBPERF_3_1 = "Glucose-6-Phosphate Dehydrogenase test performed", 
  M08_RBC_GLUC6_LBORRES_1 = "G6PD results", 
  M08_MN_LBPERF_12_1 = "C-reactive protein (CRP) test performed",
  M08_CRP_LBORRES_1 = "C-reactive protein (CRP) test results",
  M08_MN_LBPERF_13_1 = "Alpha 1-acid glycoprotein (AGP) test performed",
  M08_AGP_LBORRES_1 = "Alpha 1-acid glycoprotein (AGP) test results",
  EDD = "Estimated due date (EDD)",
  BASELINEDATE_1 = "Screening Date",
  # REMAPP_LAUNCH = "Remapp launch",
  AGE_ENROLL = "Age at enrollment",
  CRIT_AGE = "Aged 18 to 34 years", 
  BASELINE_GA_WKS = "GA at screening: weeks", 
  CRIT_GA = "Gestational age <14 weeks at enrollment", 
  BMI = paste0("Body mass index (BMI) kg/", expression(m^2)), 
  CRIT_BMI = paste0("BMI >18.5 and <30kg/", expression(m^2)), 
  CRIT_MUAC = "Mid-upper arm circumference (MUAC) > 23cm", 
  CRIT_HEIGHT = "Height >153 cm", 
  CRIT_SINGLEPREG = "Singleton pregnancy", 
  M06_BP_SYS_1 = "mean BP measurement:  (systolic)",
  M06_BP_DIA_1 = "mean BP measurement: (diastolic)",
  CRIT_BP = "Systolic blood pressure <140 mmHg and diastolic blood pressure <90 mmHg", 
  CRIT_PREV_MISCARR = "< 1 miscarriage in two consecutive pregnancies", 
  # CRIT_PREV_PRETERM_LBW = "No previous preterm or low birth weight delivery",
  CRIT_FETALDEATH = "No previous fetal death", 
  CRIT_COMPLICATION = "No history of pregnancy complications", 
  CRIT_SMOKE = "Non-cigarette smoking, tobacco chewing, or betel nut use ", 
  CRIT_DRINK = "No alcohol consumption during pregnancy", 
  CRIT_HIV = "No known history or current HIV",
  CRIT_CHRONIC = "No known history or current chronic disease including cancer, kidney disease, and cardiac conditions", 
  CRIT_MALARIA = "No current malaria infection (per RDT) at baseline", 
  CRIT_HEPATITISB = "No Hepatitis B virus infection", 
  CRIT_HEPATITISC = "No Hepatitis C virus infection", 
  CRIT_HEMOGLOBINOPATHIES = "No hemoglobinopathies", 
  CRIT_IRON = "Not iron deficient (serum ferritin >15 mcg/L-adjusted for inflammation)", 
  CRIT_INFLAM = "No subclinical inflammation (CRP > 5 mg/L and/or AGP > 1 g/L)", 
  #for new variables merged in, not in healthyoutcome data
  M04_APH_COMP_RPTEST_1_1 = "Have Placental abruption that contributed to antepartum hemorrhage",
  M04_PPH_COMP_RPORRES_1_1 = "Have Uterine atony that contributed to postpartum hemorrhage",
  M04_PPH_COMP_RPORRES_2_1 = "Have Retained placenta that contributed to postpartum hemorrhage",
  M04_PPH_COMP_RPORRES_3_1 = "Have Uterine rupture that contributed to postpartum hemorrhage",
  M04_PPH_COMP_RPORRES_4_1 = "Have Cervical/vaginal laceration that contributed to postpartum hemorrhage",
  CRIT_PRETERM = "No previous preterm",
  CRIT_LBW = "No previous low birth weight delivery",
  CRIT_STILLBIRTH = "No previous stillbirth",
  CRIT_HEIGHT_NEW = "Height >= 150 cm",
  CRIT_COMPLICATION_NEW = "No history of pregnancy complications",
  CRIT_UNPL_CESARIAN = "No previous unplanned cesarean delivery",
  HEALTHY_CHECK = "Num of criteria questions answered",
  HEALTHY_ELIGIBLE = "Eligible status"
) 

healthyOutcomeLabel <- vars_criteria %>% 
  #merge new variables in
  left_join(MatData_Anc_Visits %>% select(MOMID, PREGID, M04_APH_COMP_RPTEST_1_1, 
                                          num_range("M04_PPH_COMP_RPORRES_",1:4,"_1")), 
            by = c("MOMID", "PREGID")) %>% 
  distinct() %>% 
  select(-"CRIT_PREV_PRETERM_LBW") %>% 
  mutate(#additional derived variables
         #preterm
         CRIT_PRETERM = case_when(
           M04_PRETERM_RPORRES_1 == 1 ~ 0,
           M04_PH_PREV_RPORRES_1 == 0 | M04_PRETERM_RPORRES_1 == 0 ~ 1, 
           M04_PRETERM_RPORRES_1 == 99 ~ 0,
           TRUE ~ 55),
         #low birth weight
         CRIT_LBW = case_when(
           M04_LOWBIRTHWT_RPORRES_1 == 1 ~ 0, 
           M04_PH_PREV_RPORRES_1 == 0 | M04_LOWBIRTHWT_RPORRES_1 == 0 ~ 1,
           M04_LOWBIRTHWT_RPORRES_1 == 99 ~ 0,
           TRUE ~ 55),
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
        ) %>% 
  rowwise() %>%
  mutate(HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE)) %>%
  mutate(
    HEALTHY_ELIGIBLE = case_when(
      if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
      if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
      HEALTHY_CHECK < 22 ~ 3 #pending
    ) )%>%
  ungroup() %>%
  mutate(HEALTHY_CHECK = as.factor(HEALTHY_CHECK),
         CRIT_DRINK = as.character(CRIT_DRINK),
         CRIT_MALARIA = as.character(CRIT_MALARIA))

healthyOutcomeLabel <- set_variable_labels(healthyOutcomeLabel, .labels = healthy_labels) 

#crosstabs on newly added criteria
generalTB1 = crosstable(as.data.frame(healthyOutcomeLabel),
                        c(CRIT_HEIGHT_NEW, CRIT_LBW, CRIT_STILLBIRTH, CRIT_UNPL_CESARIAN, CRIT_COMPLICATION_NEW),
                        by = c(SITE),
                        percent_pattern="{n} ({p_col})", #{p_col}
                        percent_digits=0,
                        showNA = "ifany",
                        # funs = c(n, p_col), 
                        total = "row",
                        label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
generalTB1

#crosstab overall
#old criteria
healthyOutcomeLabel_old <- healthyOutcomeLabel %>% 
  select(-c("CRIT_HEIGHT_NEW", "CRIT_STILLBIRTH", "CRIT_UNPL_CESARIAN", 
            "CRIT_COMPLICATION_NEW", "CRIT_HEMOGLOBINOPATHIES", "CRIT_IRON") ) %>% 
  rowwise() %>%
  mutate(HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE)) %>% 
  mutate(
    HEALTHY_ELIGIBLE = case_when(
      if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
      if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
      HEALTHY_CHECK < 22 ~ 3 #pending 
    ) )%>%
  ungroup() %>% 
  mutate(HEALTHY_CHECK = as.factor(HEALTHY_CHECK)) 

#all old criteria 
tb_criteria_old = crosstable(as.data.frame(healthyOutcomeLabel_old),
                       c(CRIT_AGE, CRIT_GA, CRIT_BMI, 
                         CRIT_MUAC, CRIT_HEIGHT, CRIT_SINGLEPREG, CRIT_BP, CRIT_PREV_MISCARR, CRIT_PRETERM, CRIT_LBW,
                         CRIT_FETALDEATH, CRIT_COMPLICATION, CRIT_SMOKE, CRIT_DRINK, CRIT_HIV, CRIT_CHRONIC, CRIT_MALARIA, CRIT_HEPATITISB, CRIT_HEPATITISC,
                         CRIT_INFLAM, HEALTHY_CHECK, HEALTHY_ELIGIBLE #CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, 
                         ),
                       by = c(SITE),
                       percent_pattern="{n} ({p_col})",
                       percent_digits=0,
                       showNA = "ifany",
                       total = "row",
                       label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) 
tb_criteria_old
#old criteria eligible
tb_criteria_old_1 = crosstable(as.data.frame(healthyOutcomeLabel_old),
                             c(HEALTHY_ELIGIBLE),
                             by = c(SITE),
                             percent_pattern="{n} ({p_col})", 
                             percent_digits=0,
                             showNA = "ifany",
                             total = "row",
                             label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) 
tb_criteria_old_1


#new criterias
healthyOutcomeLabel_new <- healthyOutcomeLabel %>% 
  select(-c("CRIT_HEIGHT", "CRIT_FETALDEATH", "CRIT_COMPLICATION", 
            "CRIT_COMPLICATION_NEW", "CRIT_HEMOGLOBINOPATHIES", "CRIT_IRON") ) %>% 
  rowwise() %>%
  mutate(HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE)) %>% 
  mutate(
    HEALTHY_ELIGIBLE = case_when(
      if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
      if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
      HEALTHY_CHECK < 22 ~ 3 #pending 
    ) )%>%
  ungroup() %>% 
  mutate(HEALTHY_CHECK = as.factor(HEALTHY_CHECK))

#crosstab overall on new criteria
tb_criteria_new = crosstable(as.data.frame(healthyOutcomeLabel_new),
                       c(CRIT_AGE, CRIT_GA, CRIT_BMI, 
                         CRIT_MUAC, CRIT_HEIGHT_NEW, CRIT_SINGLEPREG, CRIT_BP, CRIT_PREV_MISCARR, CRIT_PRETERM, CRIT_LBW,
                         CRIT_STILLBIRTH, CRIT_UNPL_CESARIAN, CRIT_SMOKE, CRIT_DRINK, CRIT_HIV, CRIT_CHRONIC, CRIT_MALARIA, CRIT_HEPATITISB, CRIT_HEPATITISC,
                         CRIT_INFLAM, HEALTHY_CHECK, HEALTHY_ELIGIBLE #CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, 
                       ),
                       by = c(SITE),
                       percent_pattern="{n} ({p_col})", 
                       percent_digits=0,
                       showNA = "ifany",
                       # funs = c(n, p_col), 
                       total = "row",
                       label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) 
tb_criteria_new

#new criteria eligible
tb_criteria_new_1 = crosstable(as.data.frame(healthyOutcomeLabel_new),
                             c(HEALTHY_ELIGIBLE),
                             by = c(SITE),
                             percent_pattern="{n} ({p_col})", 
                             percent_digits=0,
                             showNA = "ifany",
                             total = "row",
                             label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) 
tb_criteria_new_1

#crosstab overall
#old criteria but remove CRIT_GA
healthyOutcomeLabel_old_noga <- healthyOutcomeLabel %>% 
  select(-c("CRIT_HEIGHT_NEW", "CRIT_STILLBIRTH", "CRIT_UNPL_CESARIAN", 
            "CRIT_COMPLICATION_NEW", "CRIT_GA", "HEALTHY_CHECK", "HEALTHY_ELIGIBLE", 
            "CRIT_HEMOGLOBINOPATHIES", "CRIT_IRON") ) %>% 
  rowwise() %>%
  mutate(HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE)) %>% 
  mutate(
    HEALTHY_ELIGIBLE = case_when(
      if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
      if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
      HEALTHY_CHECK < 21 ~ 3 #pending 
    ) )%>%
  ungroup() %>% 
  mutate(HEALTHY_CHECK = as.factor(HEALTHY_CHECK)) 
#old criteria remove CRIT_GA
tb_criteria_old_noga = crosstable(as.data.frame(healthyOutcomeLabel_old_noga),
                             c(CRIT_AGE, CRIT_BMI, 
                               CRIT_MUAC, CRIT_HEIGHT, CRIT_SINGLEPREG, CRIT_BP, CRIT_PREV_MISCARR, CRIT_PRETERM, CRIT_LBW,
                               CRIT_FETALDEATH, CRIT_COMPLICATION, CRIT_SMOKE, CRIT_DRINK, CRIT_HIV, CRIT_CHRONIC, CRIT_MALARIA, CRIT_HEPATITISB, CRIT_HEPATITISC,
                               CRIT_INFLAM, HEALTHY_CHECK, HEALTHY_ELIGIBLE #CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, 
                             ),
                             by = c(SITE),
                             percent_pattern="{n} ({p_col})", 
                             percent_digits=0,
                             showNA = "ifany",
                             total = "row",
                             label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) 
tb_criteria_old_noga
#old criteria remove CRIT_GA eligible
tb_criteria_old_noga_1 = crosstable(as.data.frame(healthyOutcomeLabel_old_noga),
                                  c(HEALTHY_ELIGIBLE
                                  ),
                                  by = c(SITE),
                                  percent_pattern="{n} ({p_col})", 
                                  percent_digits=0,
                                  showNA = "ifany",
                                  total = "row",
                                  label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
tb_criteria_old_noga_1

#new criterias
healthyOutcomeLabel_new_noga <- healthyOutcomeLabel %>% 
  select(-c("CRIT_HEIGHT", "CRIT_FETALDEATH", "CRIT_COMPLICATION", 
            "CRIT_COMPLICATION_NEW", "CRIT_GA", "HEALTHY_CHECK", 
            "HEALTHY_ELIGIBLE", "CRIT_HEMOGLOBINOPATHIES", "CRIT_IRON") ) %>% 
  # rename(old_height = CRIT_HEIGHT, 
  #        fetaldeath = CRIT_FETALDEATH, 
  #        old_complication = CRIT_COMPLICATION,
  #        new_complication = CRIT_COMPLICATION_NEW) %>% 
  rowwise() %>%
  mutate(HEALTHY_CHECK = sum(across(starts_with("CRIT_"), ~ .x %in% c(1, 0, 666)), na.rm = TRUE)) %>% 
  mutate(
    HEALTHY_ELIGIBLE = case_when(
      if_all(starts_with("CRIT_"), ~.x %in% c(1, 666)) ~ 1, #eligible
      if_any(starts_with("CRIT_"), ~.x == 0) ~ 0, #Not eligible
      HEALTHY_CHECK < 21 ~ 3 #pending 
    ) )%>%
  ungroup() %>% 
  mutate(HEALTHY_CHECK = as.factor(HEALTHY_CHECK),
         CRIT_DRINK = as.character(CRIT_DRINK), 
         CRIT_MALARIA = as.character(CRIT_MALARIA))

#crosstab overall on new criteria
tb_criteria_new_noga = crosstable(as.data.frame(healthyOutcomeLabel_new_noga),
                             c(CRIT_AGE, CRIT_BMI, 
                               CRIT_MUAC, CRIT_HEIGHT_NEW, CRIT_SINGLEPREG, CRIT_BP, CRIT_PREV_MISCARR, CRIT_PRETERM, CRIT_LBW,
                               CRIT_STILLBIRTH, CRIT_UNPL_CESARIAN, CRIT_SMOKE, CRIT_DRINK, CRIT_HIV, CRIT_CHRONIC, CRIT_MALARIA, CRIT_HEPATITISB, CRIT_HEPATITISC,
                               CRIT_INFLAM, HEALTHY_CHECK, HEALTHY_ELIGIBLE #CRIT_HEMOGLOBINOPATHIES, CRIT_IRON, 
                             ),
                             by = c(SITE),
                             percent_pattern="{n} ({p_col})", #{p_col}
                             percent_digits=0,
                             showNA = "ifany",
                             # funs = c(n, p_col), 
                             total = "row",
                             label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
tb_criteria_new_noga
#crosstab overall on new criteria eligible
tb_criteria_new_noga_1 = crosstable(as.data.frame(healthyOutcomeLabel_new_noga),
                                  c(HEALTHY_ELIGIBLE
                                  ),
                                  by = c(SITE),
                                  percent_pattern="{n} ({p_col})", #{p_col}
                                  percent_digits=0,
                                  showNA = "ifany",
                                  # funs = c(n, p_col), 
                                  total = "row",
                                  label = TRUE) %>%
  as_flextable(keep_id = FALSE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
tb_criteria_new_noga_1

#******************************************************************************
#function for floor decimal
myFloor <- function(x, Decimals=1) {
  x2<-x*10^Decimals
  floor(x2)/10^Decimals
}

#table on original vars
dataOriginal <- healthyOutcomeLabel %>% 
  mutate(AGE_ENROLL = floor(AGE_ENROLL),
         BASELINE_GA_WKS = floor(BASELINE_GA_WKS), 
         BMI = myFloor(BMI,1), 
         M05_MUAC_PERES_1 = floor(M05_MUAC_PERES_1),
         M05_HEIGHT_PERES_1 = floor(M05_HEIGHT_PERES_1), 
         M06_BP_SYS_1 = floor(M06_BP_SYS_1), 
         M06_BP_DIA_1 = floor(M06_BP_DIA_1), 
         )
dataOriginal[] <- lapply(dataOriginal, function(x) if(is.numeric(x)) 
  as.character(x) else x)
dataOriginal <-set_variable_labels(dataOriginal, .labels = healthy_labels)
criteriaTB = crosstable(dataOriginal,
                       c(AGE_ENROLL, #1.age
                         BASELINE_GA_WKS, #2.ga
                         BMI, #3.BMI
                         M05_MUAC_PERES_1, #4.MUAC
                         M05_HEIGHT_PERES_1, #5.height
                         M06_SINGLETON_PERES_1, #6.singleton
                         M06_BP_SYS_1, M06_BP_DIA_1, #7.blood pressure
                         M04_PH_PREV_RPORRES_1, M04_PH_OTH_RPORRES_1, M04_MISCARRIAGE_RPORRES_1, M04_MISCARRIAGE_CT_RPORRES_1, #8.miscarriage 
                         M04_PH_PREV_RPORRES_1, M04_PRETERM_RPORRES_1, M04_LOWBIRTHWT_RPORRES_1, #9. preterm & low birthweight
                         M04_PH_PREV_RPORRES_1, M04_PH_PREVN_RPORRES_1, M04_PH_LIVE_RPORRES_1, M04_PH_OTH_RPORRES_1, #10.neonetal/fetal death
                         M04_UNPL_CESARIAN_PROCCUR_1, M04_PREECLAMPSIA_RPORRES_1, M04_GEST_DIAB_RPORRES_1, 
                         M04_PREMATURE_RUPTURE_RPORRES_1, M04_MACROSOMIA_RPORRES_1, M04_OLIGOHYDRAMNIOS_RPORRES_1, 
                         M04_APH_RPORRES_1, M04_PPH_RPORRES_1, #11.complication
                         M03_SMOKE_OECOCCUR_1, M03_CHEW_BNUT_OECOCCUR_1, M03_CHEW_OECOCCUR_1, #12.smoke
                         M03_DRINK_OECOCCUR_1, #13.drink
                         M06_HIV_POC_LBORRES_1, M04_HIV_EVER_MHOCCUR_1, M04_HIV_MHOCCUR_1, #14.HIV
                         M04_CANCER_EVER_MHOCCUR_1, M04_KIDNEY_EVER_MHOCCUR_1, M04_CARDIAC_EVER_MHOCCUR_1, #15.cancer, kidney, cardiac disease
                         M04_MALARIA_EVER_MHOCCUR_1, #16.malaria
                         M06_HBV_POC_LBORRES_1, #17.HBV
                         M06_HCV_POC_LBORRES_1, #18.HCV
                         M08_RBC_LBPERF_2_1, M08_RBC_THALA_LBORRES_1, M08_RBC_LBPERF_3_1, M08_RBC_GLUC6_LBORRES_1, #19.hemoglobianpathies
                         M08_MN_LBPERF_8_1, M08_FERRITIN_LBORRES_1, #20.iron 
                         M08_MN_LBPERF_12_1, M08_CRP_LBORRES_1, M08_MN_LBPERF_13_1, M08_AGP_LBORRES_1 #21. inflamation
                         ),
                       by = c(SITE),
                       percent_pattern="{n} ({p_col})", #{p_col}
                       percent_digits=0,
                       showNA = "ifany",
                       # funs = (n),
                       total = "row",
                       label = TRUE) %>%
  as_flextable(keep_id = TRUE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
criteriaTB

#crosstab for preterm, low birth weight, stillbirth, miscarrige
table_fetal = crosstable(as.data.frame(dataOriginal),
                         c(M04_PH_PREV_RPORRES_1, M04_PRETERM_RPORRES_1, M04_LOWBIRTHWT_RPORRES_1,
                           M04_STILLBIRTH_RPORRES_1, 
                           M04_MISCARRIAGE_RPORRES_1, M04_MISCARRIAGE_CT_RPORRES_1
                           ),
                         by = c(SITE),
                         percent_pattern="{n} ({p_col})", #{p_col}
                         percent_digits=0,
                         showNA = "ifany",
                         # funs = c(n, p_col), 
                         total = "row",
                         label = TRUE) %>%
  as_flextable(keep_id = TRUE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
table_fetal


#original var for new criteria
table_newvar = crosstable(as.data.frame(dataOriginal),
                         c(M04_LOWBIRTHWT_RPORRES_1,
                           M04_STILLBIRTH_RPORRES_1, 
                           M04_UNPL_CESARIAN_PROCCUR_1, 
                           M05_HEIGHT_PERES_1
                         ),
                         by = c(SITE),
                         percent_pattern="{n} ({p_col})", #{p_col}
                         percent_digits=0,
                         showNA = "ifany",
                         # funs = c(n, p_col), 
                         total = "row",
                         label = TRUE) %>%
  as_flextable(keep_id = TRUE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
table_newvar


#crosstab for history of pregnancy complications
data_complication <- dataOriginal %>% filter(M04_PH_PREV_RPORRES_1 == 1)
table_complication = crosstable(as.data.frame(data_complication),
                         c(#M04_UNPL_CESARIAN_PROCCUR_1, 
                           M04_PREECLAMPSIA_RPORRES_1, M04_GEST_DIAB_RPORRES_1,
                           M04_PREMATURE_RUPTURE_RPORRES_1, M04_MACROSOMIA_RPORRES_1,
                           M04_OLIGOHYDRAMNIOS_RPORRES_1, 
                           M04_APH_RPORRES_1, M04_APH_COMP_RPTEST_1_1,
                           M04_PPH_RPORRES_1, M04_PPH_COMP_RPORRES_1_1,
                           M04_PPH_COMP_RPORRES_2_1, M04_PPH_COMP_RPORRES_3_1,
                           M04_PPH_COMP_RPORRES_4_1
                         ),
                         by = c(SITE),
                         percent_pattern="{n} ({p_col})", #{p_col}
                         percent_digits=0,
                         showNA = "ifany",
                         # funs = c(n, p_col), 
                         total = "row",
                         label = TRUE) %>%
  as_flextable(keep_id = TRUE, header_show_n = 1:2,
               header_show_n_pattern=list(cell="{.col} (N={.n})", total="Total\n(N={.n})")) #compact = TRUE, header_show_n = 1:2
table_complication
#******************************************************************************
#histogram
#1.age
# table(vars_criteria$AGE_ENROLL)
healthyOutcomeLabel %>% 
  filter(!is.na(AGE_ENROLL)) %>% 
  ggplot(aes(x = floor(AGE_ENROLL), fill = SITE)) +
  geom_histogram(color = "black",
                          alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=seq(14,50,1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Age (year)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Age by Site") +
  theme_bw() +
  geom_vline(aes(xintercept = 18), colour="green") +
  geom_vline(aes(xintercept = 34), colour="purple") 

#2.baseline ga
# table(vars_criteria$BASELINE_GA_WKS)
healthyOutcomeLabel %>% 
  filter(!is.na(BASELINE_GA_WKS)) %>% 
  ggplot(aes(x = floor(BASELINE_GA_WKS), fill = SITE)) +
  geom_histogram(color = "black",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=seq(0,20,1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Gestational Age (week)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of GA at Baseline by Site") +
  theme_bw() +
  geom_vline(aes(xintercept = 14), colour="purple") 

#3.BMI
# table(vars_criteria$BMI)
plot_BMI <- healthyOutcomeLabel %>% 
  filter(!is.na(BMI)) 
plot_BMI %>% 
  ggplot(aes(x = myFloor(BMI,1), fill = SITE)) +
  geom_histogram(color = "black",
                 alpha = 0.5, position = "stack", binwidth =0.1) +
  scale_x_continuous(breaks=seq(floor(min(plot_BMI$BMI)),floor(max(plot_BMI$BMI)+1),0.5))+
  scale_y_continuous(breaks=seq(0,500,1))+
  xlab("Body Mass Index (kg/m^2)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of BMI at Baseline by Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = 18.5), colour="coral") +
  geom_vline(aes(xintercept = 30), colour="purple") 

#4.MUAC
table(healthyOutcomeLabel$M05_MUAC_PERES_1)
plot_MUAC <- healthyOutcomeLabel %>% 
  filter(!is.na(M05_MUAC_PERES_1)) 
plot_MUAC %>% 
  ggplot(aes(x = floor(M05_MUAC_PERES_1), fill = SITE)) +
  geom_histogram(color = "black",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=seq(floor(min(plot_MUAC$M05_MUAC_PERES_1)),floor(max(plot_MUAC$M05_MUAC_PERES_1)+1),1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Mid-Upper Arm Circumference (cm)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of MUAC at Baseline by Site") +
  theme_bw() +
  geom_vline(aes(xintercept = 23), colour="purple") 

#5.height
table(healthyOutcomeLabel$M05_HEIGHT_PERES_1)
plot_height <- healthyOutcomeLabel %>% 
  filter(!is.na(M05_HEIGHT_PERES_1)) 
plot_height %>% 
  ggplot(aes(x = floor(M05_HEIGHT_PERES_1), fill = SITE)) +
  geom_histogram(color = "black",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=
    seq(floor(min(plot_height$M05_HEIGHT_PERES_1)),floor(max(plot_height$M05_HEIGHT_PERES_1)+1),1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Height (cm)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Height at Baseline by Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 153), colour="purple") 

#6.BP measurement (systolic)
table(healthyOutcomeLabel$M06_BP_SYS_1)
plot_bpSYS <- healthyOutcomeLabel %>% 
  filter(!is.na(M06_BP_SYS_1)) 
plot_bpSYS %>% 
  ggplot(aes(x = floor(M06_BP_SYS_1), fill = SITE)) +
  geom_histogram(color = "black",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=
                       seq(floor(min(plot_bpSYS$M06_BP_SYS_1)),floor(max(plot_bpSYS$M06_BP_SYS_1)+1),2))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Systolic blood pressure (mmHg)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Systolic blood pressure at Baseline by Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = 140), colour="purple") 

#7. BP measurement: (diastolic)
table(healthyOutcomeLabel$M06_BP_DIA_1, healthyOutcomeLabel$SITE)
plot_bpDIA <- healthyOutcomeLabel %>% 
  filter(!is.na(M06_BP_DIA_1)) 
plot_bpDIA %>% 
  ggplot(aes(x = floor(M06_BP_DIA_1), fill = SITE)) +
  geom_histogram(color = "black",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=
                       seq(floor(min(plot_bpDIA$M06_BP_DIA_1)),floor(max(plot_bpDIA$M06_BP_DIA_1)+1),1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Diastolic Blood Pressure (mmHg)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Diastolic Blood Pressure by Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 90), colour="purple") 

#*****
#1. age
healthyOutcomeLabel %>% 
filter(!is.na(AGE_ENROLL)) %>% 
  ggplot(aes(x = floor(AGE_ENROLL))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=seq(14,50,1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Age (year)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Age by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  geom_vline(aes(xintercept = 18), colour="green4") +
  geom_vline(aes(xintercept = 34), colour="purple")

#2.baseline ga
# table(vars_criteria$BASELINE_GA_WKS)
healthyOutcomeLabel %>% 
  filter(!is.na(BASELINE_GA_WKS)) %>% 
  ggplot(aes(x = floor(BASELINE_GA_WKS))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=seq(0,20,1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Gestational Age (week)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of GA at Baseline by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  geom_vline(aes(xintercept = 14), colour="purple")



#4.MUAC
table(healthyOutcomeLabel$M05_MUAC_PERES_1)
plot_MUAC <- healthyOutcomeLabel %>% 
  filter(!is.na(M05_MUAC_PERES_1)) %>% 
  #!!!remove negative values
  filter(M05_MUAC_PERES_1 > 0 & M05_MUAC_PERES_1 <100)  
plot_MUAC %>% 
  ggplot(aes(x = floor(M05_MUAC_PERES_1))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=seq(floor(min(plot_MUAC$M05_MUAC_PERES_1)),floor(max(plot_MUAC$M05_MUAC_PERES_1)+1),1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Mid-Upper Arm Circumference (cm)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of MUAC at Baseline by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  geom_vline(aes(xintercept = 23), colour="purple")

#5.height
table(healthyOutcomeLabel$M05_HEIGHT_PERES_1)
plot_height <- healthyOutcomeLabel %>% 
  filter(!is.na(M05_HEIGHT_PERES_1)) %>% 
  #!!!remove negative values
  filter(M05_HEIGHT_PERES_1 > 120 ) 
plot_height %>% 
  ggplot(aes(x = floor(M05_HEIGHT_PERES_1))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=
                       seq(floor(min(plot_height$M05_HEIGHT_PERES_1)),floor(max(plot_height$M05_HEIGHT_PERES_1)+1),1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Height (cm)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Height at Baseline by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = 150), colour="green4") +
  geom_vline(aes(xintercept = 153), colour="purple")

#6.BP measurement (systolic)
table(healthyOutcomeLabel$M06_BP_SYS_1)
plot_bpSYS <- healthyOutcomeLabel %>% 
  filter(!is.na(M06_BP_SYS_1) & M06_BP_SYS_1>0) 
plot_bpSYS %>% 
  ggplot(aes(x = floor(M06_BP_SYS_1))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=
                       seq(floor(min(plot_bpSYS$M06_BP_SYS_1)),floor(max(plot_bpSYS$M06_BP_SYS_1)+1),2))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Systolic blood pressure (mmHg)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Systolic blood pressure at Baseline by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 140), colour="purple") 

#7. BP measurement: (diastolic)
table(healthyOutcomeLabel$M06_BP_DIA_1, healthyOutcomeLabel$SITE)
plot_bpDIA <- healthyOutcomeLabel %>% 
  filter(!is.na(M06_BP_DIA_1) & M06_BP_DIA_1>0) 
plot_bpDIA %>% 
  ggplot(aes(x = floor(M06_BP_DIA_1))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =1) +
  scale_x_continuous(breaks=
                       seq(floor(min(plot_bpDIA$M06_BP_DIA_1)),floor(max(plot_bpDIA$M06_BP_DIA_1)+1),1))+
  scale_y_continuous(breaks=seq(0,500,10))+
  xlab("Diastolic Blood Pressure (mmHg)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of Diastolic Blood Pressure by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 90), colour="purple") 

#3.BMI
# table(vars_criteria$BMI)
plot_BMI <- healthyOutcomeLabel %>% 
  filter(!is.na(BMI)) %>% 
  filter(BMI > 0 & BMI < 50) 
plot_BMI %>% 
  ggplot(aes(x = myFloor(BMI,1))) +
  geom_histogram(color = "black", fill = "#56B4E9",
                 alpha = 0.5, position = "stack", binwidth =0.1) +
  scale_x_continuous(breaks=seq(floor(min(plot_BMI$BMI)),floor(max(plot_BMI$BMI)+1),0.5))+
  scale_y_continuous(breaks=seq(0,500,1))+
  xlab("Body Mass Index (kg/m^2)") +
  ylab("Frequency Count") + 
  ggtitle("Histogram of BMI at Baseline by Site") +
  facet_grid(rows = vars(SITE), scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1),
        legend.position="bottom",
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 18.5), colour="green4") +
  geom_vline(aes(xintercept = 30), colour="purple") 
