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
library(zoo)
library(openxlsx)

transform_to_NA<-function(variable){
  variable[variable==-7]<-NA
  variable[variable==-5]<-NA
  variable[variable==-6]<-NA
  variable[variable==-9]<-NA
  variable[variable==77]<-NA
  variable[variable==55]<-NA
  variable[variable==66]<-NA
  variable[variable==88]<-NA
  variable[variable==99]<-NA
  return(variable)
}

data.aim3.risk.part_B <- read_dta("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-10-18/ReMAPP-Aim3.dta")
data.anemia <- read_dta("D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_ANEMIA.dta")

data.aim3.partB <- data.aim3.risk.part_B %>% left_join(data.anemia,by = c("SITE", "MOMID", "PREGID"))

# Function to collect RR and CI for GLM models, skipping missing predictors and errors
collect_glm_rr_results <- function(data, outcome_var, predictors) {
  # Initialize an empty data frame to store results
  results <- data.frame(
    Predictor = character(),
    RR = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each predictor variable, fit the GLM, and extract results
  for (predictor in predictors) {
    # Check if predictor exists in the data
    if (!(predictor %in% names(data))) {
      message(paste("Skipping predictor", predictor, "- variable not found in dataset"))
      next
    }
    
    # Define the formula
    formula <- as.formula(paste(outcome_var, "~", predictor))
    
    # Fit the GLM model and handle errors
    model <- tryCatch({
      glm(formula, data = data, family = binomial(link = "log"))
    }, error = function(e) {
      message(paste("Skipping predictor", predictor, "due to error:", e$message))
      return(NULL)
    })
    
    # If model fitting was successful, extract the RR and CI
    if (!is.null(model) && nrow(summary(model)$coefficients) >= 2) {
      coef_summary <- summary(model)$coefficients
      estimate <- coef_summary[2, "Estimate"]
      se <- coef_summary[2, "Std. Error"]
      
      # Calculate RR and CI
      rr <- round(exp(estimate),2)
      ci_lower <- round(exp(estimate - 1.96 * se),2)
      ci_upper <- round(exp(estimate + 1.96 * se),2)
      
      # Append to results data frame
      results <- rbind(results, data.frame(Predictor = predictor, RR = rr, CI_Lower = ci_lower, CI_Upper = ci_upper))
    }
  }
  
  return(results)
}

data.aim3.partB.T1 <- data.aim3.partB %>% 
  filter(REMAPP_T1 == 1) %>% 
  select(SITE, MOMID, PREGID, matches("T1$")) %>% 
  mutate(
    ANEMIA_OUTCOME1 = case_when(
      ANEMIA_T1 %in% c(1, 2, 3) ~ 1,
      ANEMIA_T1 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    HEP_T1 = case_when(
      HEP_T1 == 1 ~ "low",
      HEP_T1 == 2 ~ "normal",
      HEP_T1 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    TIBC_T1 = case_when(
      TIBC_T1 == 1 ~ "low",
      TIBC_T1 == 2 ~ "normal",
      TIBC_T1 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    VITA_T1 = case_when(
      VITA_T1 == 1 ~ "severe",
      VITA_T1 == 2 ~ "moderate",
      VITA_T1 == 3 ~ "mild",
      TRUE ~ NA_character_
    ),
    ZINC_T1 = case_when(
      ZINC_T1 == 1 ~ "low",
      ZINC_T1 == 2 ~ "normal",
      ZINC_T1 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    LEAD5_T1 = case_when(
      LEAD5_T1 == 1 ~ "normal",
      LEAD5_T1 == 2 ~ "high",
      TRUE ~ NA_character_
    ),
    LEAD10_T1 = case_when(
      LEAD10_T1 == 1 ~ "normal",
      LEAD10_T1 == 2 ~ "high",
      TRUE ~ NA_character_
    ),
    HELM_T1 = case_when(
      HELM_T1 == 0 ~ "normal",
      HELM_T1 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    HELM2_T1 = case_when(
      HELM2_T1 == 0 ~ "normal",
      HELM2_T1 == 1 ~ "moderate/severe",
      TRUE ~ NA_character_
    ),
    SCHISTO_STOOL_T1 = case_when(
      SCHISTO_STOOL_T1 == 0 ~ "normal",
      SCHISTO_STOOL_T1 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    SCHISTO_URINE_T1 = case_when(
      SCHISTO_URINE_T1 == 0 ~ "normal",
      SCHISTO_URINE_T1 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    MALBL_T1 = case_when(
      MALBL_T1 == 0 ~ "normal",
      MALBL_T1 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    RBC_MORPH_T1 = case_when(
      RBC_MORPH_T1 == 0 ~ "normal",
      RBC_MORPH_T1 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    WBC_MORPH_T1 = case_when(
      WBC_MORPH_T1 == 0 ~ "normal",
      WBC_MORPH_T1 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    PL_MORPH_T1 = case_when(
      PL_MORPH_T1 == 0 ~ "normal",
      PL_MORPH_T1 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    PARA_MORPH_T1 = case_when(
      PARA_MORPH_T1 == 0 ~ "normal",
      PARA_MORPH_T1 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
  ) %>% dummy_cols(select_columns = c(
    "HEP_T1", "TIBC_T1", "VITA_T1", "ZINC_T1", "LEAD5_T1",
    "LEAD10_T1", "HELM_T1", "HELM2_T1", "SCHISTO_STOOL_T1",
    "SCHISTO_URINE_T1", "MALBL_T1", "RBC_MORPH_T1",
    "WBC_MORPH_T1", "PL_MORPH_T1", "PARA_MORPH_T1"
  ))

data.aim3.partB.T1.Ghana<-data.aim3.partB.T1%>%filter(SITE=="Ghana")
data.aim3.partB.T1.Kenya<-data.aim3.partB.T1%>%filter(SITE=="Kenya")
data.aim3.partB.T1.Pakistan<-data.aim3.partB.T1%>%filter(SITE=="Pakistan")
data.aim3.partB.T1.Zambia<-data.aim3.partB.T1%>%filter(SITE=="Zambia")
data.aim3.partB.T1.India.CMC<-data.aim3.partB.T1%>%filter(SITE=="India-CMC")
data.aim3.partB.T1.India.SAS<-data.aim3.partB.T1%>%filter(SITE=="India-SAS")

# Define the list of predictors
predictors.T1 <- c(
  "HEP_T1_normal","HEP_T1_high", "HEP_T1_low","TIBC_T1_normal","TIBC_T1_high", "TIBC_T1_low",
  "VITA_T1_mild","VITA_T1_severe", "VITA_T1_moderate","ZINC_T1_normal", "ZINC_T1_high", "ZINC_T1_low",
  "LEAD5_T1_normal","LEAD5_T1_high","LEAD10_T1_normal", "LEAD10_T1_high", "HELM_T1_normal","HELM_T1_positive", "HELM2_T1_normal","HELM2_T1_moderate/severe",
  "SCHISTO_STOOL_T1_normal","SCHISTO_STOOL_T1_positive", "SCHISTO_URINE_T1_normal","SCHISTO_URINE_T1_positive", "MALBL_T1_normal", "MALBL_T1_positive",
  "RBC_MORPH_T1_normal","RBC_MORPH_T1_abnormal", "WBC_MORPH_T1_normal","WBC_MORPH_T1_abnormal","PL_MORPH_T1_normal","PL_MORPH_T1_abnormal",
  "PARA_MORPH_T1_normal","PARA_MORPH_T1_abnormal"
)

# Collect results for each site
results.T1.Pakistan <- collect_glm_rr_results(data = data.aim3.partB.T1.Pakistan, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T1)
results.T1.Ghana <- collect_glm_rr_results(data = data.aim3.partB.T1.Ghana, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T1)
results.T1.Kenya <- collect_glm_rr_results(data = data.aim3.partB.T1.Kenya, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T1)
results.T1.Zambia <- collect_glm_rr_results(data = data.aim3.partB.T1.Zambia, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T1)
results.T1.India.CMC <- collect_glm_rr_results(data = data.aim3.partB.T1.India.CMC, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T1)
results.T1.India.SAS <- collect_glm_rr_results(data = data.aim3.partB.T1.India.SAS, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T1)

data.aim3.partB.T2 <- data.aim3.partB %>% 
  filter(REMAPP_T2 == 1) %>% 
  select(SITE, MOMID, PREGID, matches("T2$")) %>% 
  mutate(
    ANEMIA_OUTCOME1 = case_when(
      ANEMIA_T2 %in% c(1, 2, 3) ~ 1,
      ANEMIA_T2 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    HEP_T2 = case_when(
      HEP_T2 == 1 ~ "low",
      HEP_T2 == 2 ~ "normal",
      HEP_T2 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    TIBC_T2 = case_when(
      TIBC_T2 == 1 ~ "low",
      TIBC_T2 == 2 ~ "normal",
      TIBC_T2 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    VITA_T2 = case_when(
      VITA_T2 == 1 ~ "severe",
      VITA_T2 == 2 ~ "moderate",
      VITA_T2 == 3 ~ "mild",
      TRUE ~ NA_character_
    ),
    ZINC_T2 = case_when(
      ZINC_T2 == 1 ~ "low",
      ZINC_T2 == 2 ~ "normal",
      ZINC_T2 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    LEAD5_T2 = case_when(
      LEAD5_T2 == 1 ~ "normal",
      LEAD5_T2 == 2 ~ "high",
      TRUE ~ NA_character_
    ),
    LEAD10_T2 = case_when(
      LEAD10_T2 == 1 ~ "normal",
      LEAD10_T2 == 2 ~ "high",
      TRUE ~ NA_character_
    ),
    HELM_T2 = case_when(
      HELM_T2 == 0 ~ "normal",
      HELM_T2 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    HELM2_T2 = case_when(
      HELM2_T2 == 0 ~ "normal",
      HELM2_T2 == 1 ~ "moderate/severe",
      TRUE ~ NA_character_
    ),
    SCHISTO_STOOL_T2 = case_when(
      SCHISTO_STOOL_T2 == 0 ~ "normal",
      SCHISTO_STOOL_T2 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    SCHISTO_URINE_T2 = case_when(
      SCHISTO_URINE_T2 == 0 ~ "normal",
      SCHISTO_URINE_T2 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    MALBL_T2 = case_when(
      MALBL_T2 == 0 ~ "normal",
      MALBL_T2 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    RBC_MORPH_T2 = case_when(
      RBC_MORPH_T2 == 0 ~ "normal",
      RBC_MORPH_T2 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    WBC_MORPH_T2 = case_when(
      WBC_MORPH_T2 == 0 ~ "normal",
      WBC_MORPH_T2 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    PL_MORPH_T2 = case_when(
      PL_MORPH_T2 == 0 ~ "normal",
      PL_MORPH_T2 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    PARA_MORPH_T2 = case_when(
      PARA_MORPH_T2 == 0 ~ "normal",
      PARA_MORPH_T2 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
  )%>% dummy_cols(select_columns = c(
    "HEP_T2", "TIBC_T2", "VITA_T2", "ZINC_T2", "LEAD5_T2",
    "LEAD10_T2", "HELM_T2", "HELM2_T2", "SCHISTO_STOOL_T2",
    "SCHISTO_URINE_T2", "MALBL_T2", "RBC_MORPH_T2",
    "WBC_MORPH_T2", "PL_MORPH_T2", "PARA_MORPH_T2"
  ))

data.aim3.partB.T2.Ghana<-data.aim3.partB.T2%>%filter(SITE=="Ghana")
data.aim3.partB.T2.Kenya<-data.aim3.partB.T2%>%filter(SITE=="Kenya")
data.aim3.partB.T2.Pakistan<-data.aim3.partB.T2%>%filter(SITE=="Pakistan")
data.aim3.partB.T2.Zambia<-data.aim3.partB.T2%>%filter(SITE=="Zambia")
data.aim3.partB.T2.India.CMC<-data.aim3.partB.T2%>%filter(SITE=="India-CMC")
data.aim3.partB.T2.India.SAS<-data.aim3.partB.T2%>%filter(SITE=="India-SAS")

# Define the list of predictors
predictors.T2 <- c(
  "HEP_T2_normal","HEP_T2_high", "HEP_T2_low","TIBC_T2_normal","TIBC_T2_high", "TIBC_T2_low",
  "VITA_T2_mild","VITA_T2_severe", "VITA_T2_moderate","ZINC_T2_normal", "ZINC_T2_high", "ZINC_T2_low",
  "LEAD5_T2_normal","LEAD5_T2_high","LEAD10_T2_normal", "LEAD10_T2_high", "HELM_T2_normal","HELM_T2_positive", "HELM2_T2_normal","HELM2_T2_moderate/severe",
  "SCHISTO_STOOL_T2_normal","SCHISTO_STOOL_T2_positive", "SCHISTO_URINE_T2_normal","SCHISTO_URINE_T2_positive", "MALBL_T2_normal", "MALBL_T2_positive",
  "RBC_MORPH_T2_normal","RBC_MORPH_T2_abnormal", "WBC_MORPH_T2_normal","WBC_MORPH_T2_abnormal","PL_MORPH_T2_normal","PL_MORPH_T2_abnormal",
  "PARA_MORPH_T2_normal","PARA_MORPH_T2_abnormal"
)
# Collect results for each site
results.T2.Pakistan <- collect_glm_rr_results(data = data.aim3.partB.T2.Pakistan, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T2)
results.T2.Ghana <- collect_glm_rr_results(data = data.aim3.partB.T2.Ghana, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T2)
results.T2.Kenya <- collect_glm_rr_results(data = data.aim3.partB.T2.Kenya, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T2)
results.T2.Zambia <- collect_glm_rr_results(data = data.aim3.partB.T2.Zambia, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T2)
results.T2.India.CMC <- collect_glm_rr_results(data = data.aim3.partB.T2.India.CMC, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T2)
results.T2.India.SAS <- collect_glm_rr_results(data = data.aim3.partB.T2.India.SAS, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T2)


data.aim3.partB.T3 <- data.aim3.partB %>% 
  filter(REMAPP_T3 == 1) %>% 
  select(SITE, MOMID, PREGID, matches("T3$")) %>% 
  mutate(
    ANEMIA_OUTCOME1 = case_when(
      ANEMIA_T3 %in% c(1, 2, 3) ~ 1,
      ANEMIA_T3 == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    HEP_T3 = case_when(
      HEP_T3 == 1 ~ "low",
      HEP_T3 == 2 ~ "normal",
      HEP_T3 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    TIBC_T3 = case_when(
      TIBC_T3 == 1 ~ "low",
      TIBC_T3 == 2 ~ "normal",
      TIBC_T3 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    VITA_T3 = case_when(
      VITA_T3 == 1 ~ "severe",
      VITA_T3 == 2 ~ "moderate",
      VITA_T3 == 3 ~ "mild",
      TRUE ~ NA_character_
    ),
    ZINC_T3 = case_when(
      ZINC_T3 == 1 ~ "low",
      ZINC_T3 == 2 ~ "normal",
      ZINC_T3 == 3 ~ "high",
      TRUE ~ NA_character_
    ),
    LEAD5_T3 = case_when(
      LEAD5_T3 == 1 ~ "normal",
      LEAD5_T3 == 2 ~ "high",
      TRUE ~ NA_character_
    ),
    LEAD10_T3 = case_when(
      LEAD10_T3 == 1 ~ "normal",
      LEAD10_T3 == 2 ~ "high",
      TRUE ~ NA_character_
    ),
    HELM_T3 = case_when(
      HELM_T3 == 0 ~ "normal",
      HELM_T3 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    HELM2_T3 = case_when(
      HELM2_T3 == 0 ~ "normal",
      HELM2_T3 == 1 ~ "moderate/severe",
      TRUE ~ NA_character_
    ),
    SCHISTO_STOOL_T3 = case_when(
      SCHISTO_STOOL_T3 == 0 ~ "normal",
      SCHISTO_STOOL_T3 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    SCHISTO_URINE_T3 = case_when(
      SCHISTO_URINE_T3 == 0 ~ "normal",
      SCHISTO_URINE_T3 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    MALBL_T3 = case_when(
      MALBL_T3 == 0 ~ "normal",
      MALBL_T3 == 1 ~ "positive",
      TRUE ~ NA_character_
    ),
    RBC_MORPH_T3 = case_when(
      RBC_MORPH_T3 == 0 ~ "normal",
      RBC_MORPH_T3 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    WBC_MORPH_T3 = case_when(
      WBC_MORPH_T3 == 0 ~ "normal",
      WBC_MORPH_T3 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    PL_MORPH_T3 = case_when(
      PL_MORPH_T3 == 0 ~ "normal",
      PL_MORPH_T3 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
    PARA_MORPH_T3 = case_when(
      PARA_MORPH_T3 == 0 ~ "normal",
      PARA_MORPH_T3 == 1 ~ "abnormal",
      TRUE ~ NA_character_
    ),
  )%>% dummy_cols(select_columns = c(
    "HEP_T3", "TIBC_T3", "VITA_T3", "ZINC_T3", "LEAD5_T3",
    "LEAD10_T3", "HELM_T3", "HELM2_T3", "SCHISTO_STOOL_T3",
    "SCHISTO_URINE_T3", "MALBL_T3", "RBC_MORPH_T3",
    "WBC_MORPH_T3", "PL_MORPH_T3", "PARA_MORPH_T3"
  ))

data.aim3.partB.T3.Ghana<-data.aim3.partB.T3%>%filter(SITE=="Ghana")
data.aim3.partB.T3.Kenya<-data.aim3.partB.T3%>%filter(SITE=="Kenya")
data.aim3.partB.T3.Pakistan<-data.aim3.partB.T3%>%filter(SITE=="Pakistan")
data.aim3.partB.T3.Zambia<-data.aim3.partB.T3%>%filter(SITE=="Zambia")
data.aim3.partB.T3.India.CMC<-data.aim3.partB.T3%>%filter(SITE=="India-CMC")
data.aim3.partB.T3.India.SAS<-data.aim3.partB.T3%>%filter(SITE=="India-SAS")

# Define the list of predictors
predictors.T3 <- c(
  "HEP_T3_normal","HEP_T3_high", "HEP_T3_low","TIBC_T3_normal","TIBC_T3_high", "TIBC_T3_low",
  "VITA_T3_mild","VITA_T3_severe", "VITA_T3_moderate","ZINC_T3_normal", "ZINC_T3_high", "ZINC_T3_low",
  "LEAD5_T3_normal","LEAD5_T3_high","LEAD10_T3_normal", "LEAD10_T3_high", "HELM_T3_normal","HELM_T3_positive", "HELM2_T3_normal","HELM2_T3_moderate/severe",
  "SCHISTO_STOOL_T3_normal","SCHISTO_STOOL_T3_positive", "SCHISTO_URINE_T3_normal","SCHISTO_URINE_T3_positive", "MALBL_T3_normal", "MALBL_T3_positive",
  "RBC_MORPH_T3_normal","RBC_MORPH_T3_abnormal", "WBC_MORPH_T3_normal","WBC_MORPH_T3_abnormal","PL_MORPH_T3_normal","PL_MORPH_T3_abnormal",
  "PARA_MORPH_T3_normal","PARA_MORPH_T3_abnormal"
)

# Collect results for each site
results.T3.Pakistan <- collect_glm_rr_results(data = data.aim3.partB.T3.Pakistan, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T3)
results.T3.Ghana <- collect_glm_rr_results(data = data.aim3.partB.T3.Ghana, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T3)
results.T3.Kenya <- collect_glm_rr_results(data = data.aim3.partB.T3.Kenya, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T3)
results.T3.Zambia <- collect_glm_rr_results(data = data.aim3.partB.T3.Zambia, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T3)
results.T3.India.CMC <- collect_glm_rr_results(data = data.aim3.partB.T3.India.CMC, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T3)
results.T3.India.SAS <- collect_glm_rr_results(data = data.aim3.partB.T3.India.SAS, outcome_var = "ANEMIA_OUTCOME1", predictors = predictors.T3)

collect_glmer_rr_results <- function(data, outcome_var, predictors, random_effect) {
  # Initialize an empty data frame to store results
  results <- data.frame(
    Predictor = character(),
    Numerator = numeric(),
    Total = numeric(),
    Percentage = character(),
    RR = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    Model_Type = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each predictor variable
  for (predictor in predictors) {
    # Check if predictor exists in the data
    if (!(predictor %in% names(data))) {
      message(paste("Skipping predictor", predictor, "- variable not found in dataset"))
      next
    }
    
    # Calculate summary statistics for the predictor
    Numerator <- sum(data[[predictor]] == 1 & data[[outcome_var]] == 1, na.rm = TRUE)
    Total <- sum(data[[predictor]] == 1, na.rm = TRUE)
    Percentage <- if (Total > 0) {
      paste0(round(Numerator / Total * 100, 2), "%")
    } else {
      "0%"
    }
    
    # If predictor ends with "normal," skip model fitting and set RR = 1
    if (grepl("_normal$", predictor)) {
      results <- rbind(results, data.frame(
        Predictor = predictor,
        Numerator = Numerator,
        Total = Total,
        Percentage = Percentage,
        RR = 1,
        CI_Lower = NA,
        CI_Upper = NA,
        Model_Type = "Normal (Default)",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Skip predictors with no occurrences
    if (Total == 0) {
      message(paste("Skipping predictor", predictor, "- no data for this variable"))
      results <- rbind(results, data.frame(
        Predictor = predictor,
        Numerator = Numerator,
        Total = Total,
        Percentage = Percentage,
        RR = NA,
        CI_Lower = NA,
        CI_Upper = NA,
        Model_Type = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Define the formula
    formula <- as.formula(paste(outcome_var, "~", predictor, "+ (1 |", random_effect, ")"))
    
    # Fit the glmer model and handle errors
    model <- tryCatch({
      glmer(formula, data = data, family = binomial(link = "log"))
    }, error = function(e) {
      tryCatch({
        glmer(formula, data = data, family = poisson(link = "log"))
      }, error = function(e) {
        message(paste("Skipping predictor", predictor, "due to model fitting error:", e$message))
        return(NULL)
      })
    })
    
    # If model fitting was unsuccessful, append NA results
    if (is.null(model)) {
      results <- rbind(results, data.frame(
        Predictor = predictor,
        Numerator = Numerator,
        Total = Total,
        Percentage = Percentage,
        RR = NA,
        CI_Lower = NA,
        CI_Upper = NA,
        Model_Type = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Extract coefficients and calculate RR and CI
    coef_summary <- summary(model)$coefficients
    estimate <- coef_summary[2, "Estimate"]
    se <- coef_summary[2, "Std. Error"]
    
    threshold <- 1 * 10^5
    # Calculate RR and CI
    rr <- exp(estimate)
    ci_lower <- exp(estimate - 1.96 * se)
    ci_upper <- exp(estimate + 1.96 * se)
    
    # Set values above the threshold to Inf
    rr <- ifelse(rr > threshold, Inf, round(rr, 2))
    ci_lower <- ifelse(ci_lower > threshold, Inf, round(ci_lower, 2))
    ci_upper <- ifelse(ci_upper > threshold, Inf, round(ci_upper, 2))
    
    # Determine model type
    model_type <- ifelse(family(model)$family == "binomial", "Binomial", "Poisson")
    
    # Append to results data frame
    results <- rbind(results, data.frame(
      Predictor = predictor,
      Numerator = Numerator,
      Total = Total,
      Percentage = Percentage,
      RR = rr,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      Model_Type = model_type,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

outcome_var <- "ANEMIA_OUTCOME1"
random_effect <- "SITE"

data.aim3.partB.T1.noZambia<-data.aim3.partB.T1%>%filter(SITE!="Zambia")
data.aim3.partB.T2.noZambia<-data.aim3.partB.T2%>%filter(SITE!="Zambia")
data.aim3.partB.T3.noZambia<-data.aim3.partB.T3%>%filter(SITE!="Zambia")

results.T1.glmer <- collect_glmer_rr_results(
  data = data.aim3.partB.T1.noZambia,
  outcome_var = outcome_var,
  predictors = predictors.T1,
  random_effect = random_effect
)

results.T2.glmer <- collect_glmer_rr_results(
  data = data.aim3.partB.T2.noZambia,
  outcome_var = outcome_var,
  predictors = predictors.T2,
  random_effect = random_effect
)

results.T3.glmer <- collect_glmer_rr_results(
  data = data.aim3.partB.T3.noZambia,
  outcome_var = outcome_var,
  predictors = predictors.T3,
  random_effect = random_effect
)

predictor_mapping_T1 <-c(
  "HEP_T1_normal"= "Hepcidin (normal)",
  "HEP_T1_high"= "Hepcidin (high)",
  "HEP_T1_low"= "Hepcidin (low)",
  "TIBC_T1_normal"= "Total iron binding capacity (TIBC, normal)",
  "TIBC_T1_high"= "Total iron binding capacity (TIBC, high)",
  "TIBC_T1_low"= "Total iron binding capacity (TIBC, low)",
  "VITA_T1_mild"= "Vitamin A (mild deficiency)",
  "VITA_T1_severe"= "Vitamin A (severe deficiency)",
  "VITA_T1_moderate"= "Vitamin A (moderate deficiency)",
  "ZINC_T1_normal"= "Zinc (normal)",
  "ZINC_T1_high"= "Zinc (high)",
  "ZINC_T1_low"= "Zinc (low)",
  "LEAD5_T1_normal"= "Lead (5 ug/dL, normal)",
  "LEAD5_T1_high"= "Lead (5 ug/dL, high)",
  "LEAD10_T1_normal"= "Lead (10 ug/dL, normal)",
  "LEAD10_T1_high"= "Lead (10 ug/dL, high)",
  "HELM_T1_normal"= "Helminth (normal)",
  "HELM_T1_positive"= "Helminth (positive)",
  "HELM2_T1_normal"= "Helminth (normal)",
  "HELM2_T1_moderate/severe"= "Helminth (moderate/severe)",
  "SCHISTO_STOOL_T1_normal"= "Schistosomiasis (stool, normal)",
  "SCHISTO_STOOL_T1_positive"= "Schistosomiasis (stool, positive)",
  "SCHISTO_URINE_T1_normal"= "Schistosomiasis (urine, normal)",
  "SCHISTO_URINE_T1_positive"= "Schistosomiasis (urine, positive)",
  "MALBL_T1_normal"= "Malaria (normal)",
  "MALBL_T1_positive"= "Malaria (positive)",
  "RBC_MORPH_T1_normal"= "Red cell morphology (normal)",
  "RBC_MORPH_T1_abnormal"= "Red cell morphology (abnormal)",
  "WBC_MORPH_T1_normal"= "White cell morphology (normal)",
  "WBC_MORPH_T1_abnormal"= "White cell morphology (abnormal)",
  "PL_MORPH_T1_normal"= "Platelet morphology (normal)",
  "PL_MORPH_T1_abnormal"= "Platelet morphology (abnormal)",
  "PARA_MORPH_T1_normal"= "Parasite morphology (normal)",
  "PARA_MORPH_T1_abnormal"= "Parasite morphology (abnormal)"
)

predictor_mapping_T2 <-c(
  "HEP_T2_normal"= "Hepcidin (normal)",
  "HEP_T2_high"= "Hepcidin (high)",
  "HEP_T2_low"= "Hepcidin (low)",
  "TIBC_T2_normal"= "Total iron binding capacity (TIBC, normal)",
  "TIBC_T2_high"= "Total iron binding capacity (TIBC, high)",
  "TIBC_T2_low"= "Total iron binding capacity (TIBC, low)",
  "VITA_T2_mild"= "Vitamin A (mild deficiency)",
  "VITA_T2_severe"= "Vitamin A (severe deficiency)",
  "VITA_T2_moderate"= "Vitamin A (moderate deficiency)",
  "ZINC_T2_normal"= "Zinc (normal)",
  "ZINC_T2_high"= "Zinc (high)",
  "ZINC_T2_low"= "Zinc (low)",
  "LEAD5_T2_normal"= "Lead (5 ug/dL, normal)",
  "LEAD5_T2_high"= "Lead (5 ug/dL, high)",
  "LEAD10_T2_normal"= "Lead (10 ug/dL, normal)",
  "LEAD10_T2_high"= "Lead (10 ug/dL, high)",
  "HELM_T2_normal"= "Helminth (normal)",
  "HELM_T2_positive"= "Helminth (positive)",
  "HELM2_T2_normal"= "Helminth (normal)",
  "HELM2_T2_moderate/severe"= "Helminth (moderate/severe)",
  "SCHISTO_STOOL_T2_normal"= "Schistosomiasis (stool, normal)",
  "SCHISTO_STOOL_T2_positive"= "Schistosomiasis (stool, positive)",
  "SCHISTO_URINE_T2_normal"= "Schistosomiasis (urine, normal)",
  "SCHISTO_URINE_T2_positive"= "Schistosomiasis (urine, positive)",
  "MALBL_T2_normal"= "Malaria (normal)",
  "MALBL_T2_positive"= "Malaria (positive)",
  "RBC_MORPH_T2_normal"= "Red cell morphology (normal)",
  "RBC_MORPH_T2_abnormal"= "Red cell morphology (abnormal)",
  "WBC_MORPH_T2_normal"= "White cell morphology (normal)",
  "WBC_MORPH_T2_abnormal"= "White cell morphology (abnormal)",
  "PL_MORPH_T2_normal"= "Platelet morphology (normal)",
  "PL_MORPH_T2_abnormal"= "Platelet morphology (abnormal)",
  "PARA_MORPH_T2_normal"= "Parasite morphology (normal)",
  "PARA_MORPH_T2_abnormal"= "Parasite morphology (abnormal)"
)

predictor_mapping_T3 <-c(
  "HEP_T3_normal"= "Hepcidin (normal)",
  "HEP_T3_high"= "Hepcidin (high)",
  "HEP_T3_low"= "Hepcidin (low)",
  "TIBC_T3_normal"= "Total iron binding capacity (TIBC, normal)",
  "TIBC_T3_high"= "Total iron binding capacity (TIBC, high)",
  "TIBC_T3_low"= "Total iron binding capacity (TIBC, low)",
  "VITA_T3_mild"= "Vitamin A (mild deficiency)",
  "VITA_T3_severe"= "Vitamin A (severe deficiency)",
  "VITA_T3_moderate"= "Vitamin A (moderate deficiency)",
  "ZINC_T3_normal"= "Zinc (normal)",
  "ZINC_T3_high"= "Zinc (high)",
  "ZINC_T3_low"= "Zinc (low)",
  "LEAD5_T3_normal"= "Lead (5 ug/dL, normal)",
  "LEAD5_T3_high"= "Lead (5 ug/dL, high)",
  "LEAD10_T3_normal"= "Lead (10 ug/dL, normal)",
  "LEAD10_T3_high"= "Lead (10 ug/dL, high)",
  "HELM_T3_normal"= "Helminth (normal)",
  "HELM_T3_positive"= "Helminth (positive)",
  "HELM2_T3_normal"= "Helminth (normal)",
  "HELM2_T3_moderate/severe"= "Helminth (moderate/severe)",
  "SCHISTO_STOOL_T3_normal"= "Schistosomiasis (stool, normal)",
  "SCHISTO_STOOL_T3_positive"= "Schistosomiasis (stool, positive)",
  "SCHISTO_URINE_T3_normal"= "Schistosomiasis (urine, normal)",
  "SCHISTO_URINE_T3_positive"= "Schistosomiasis (urine, positive)",
  "MALBL_T3_normal"= "Malaria (normal)",
  "MALBL_T3_positive"= "Malaria (positive)",
  "RBC_MORPH_T3_normal"= "Red cell morphology (normal)",
  "RBC_MORPH_T3_abnormal"= "Red cell morphology (abnormal)",
  "WBC_MORPH_T3_normal"= "White cell morphology (normal)",
  "WBC_MORPH_T3_abnormal"= "White cell morphology (abnormal)",
  "PL_MORPH_T3_normal"= "Platelet morphology (normal)",
  "PL_MORPH_T3_abnormal"= "Platelet morphology (abnormal)",
  "PARA_MORPH_T3_normal"= "Parasite morphology (normal)",
  "PARA_MORPH_T3_abnormal"= "Parasite morphology (abnormal)"
)

# Function to rename predictors based on mapping
rename_predictors <- function(results, mapping) {
  results$Predictor <- sapply(results$Predictor, function(x) {
    if (x %in% names(mapping)) {
      mapping[x]
    } else {
      x
    }
  })
  return(results)
}

# Apply renaming function for T1, T2, and T3
results.T1.glmer <- rename_predictors(results.T1.glmer, predictor_mapping_T1)
results.T2.glmer <- rename_predictors(results.T2.glmer, predictor_mapping_T2)
results.T3.glmer <- rename_predictors(results.T3.glmer, predictor_mapping_T3)

# View updated results
print(results.T1.glmer)
print(results.T2.glmer)
print(results.T3.glmer)

# Create a workbook
wb <- createWorkbook()

# Add worksheets and data
addWorksheet(wb, "T1 Results")
writeData(wb, "T1 Results", results.T1.glmer)

addWorksheet(wb, "T2 Results")
writeData(wb, "T2 Results", results.T2.glmer)

addWorksheet(wb, "T3 Results")
writeData(wb, "T3 Results", results.T3.glmer)

# Save workbook
saveWorkbook(wb, "D:/Users/yipeng_wei/Documents/Output/ReMAPP aim3 output/aim3_partB.xlsx", overwrite = TRUE)
