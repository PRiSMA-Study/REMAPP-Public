library(gee)
library(dplyr)
library(openxlsx)

# Function to transform specified values to NA
transform_to_NA <- function(variable){
  variable[variable %in% c(-7, -5, -6, -9, 77, 55, 66, 88, 99)] <- NA
  return(variable)
}

# Assuming `long_data.aim3` and `data.aim3` are already loaded and transformed appropriately
analysis_data.aim3 <- long_data.aim3 %>% select("MOMID", "SITE", "ANEMIA_TRIMESTER", "ANEMIA_OUTCOME1", "age", "SCHOOL_MORE10", "water_improved", "toilet_improved", "WEALTH_QUINT_5", "WEALTH_QUINT_4", "WEALTH_QUINT_3", "WEALTH_QUINT_2", "WEALTH_QUINT_1") %>%
  mutate(
    MOMID = as.factor(MOMID),
    SITE = as.factor(SITE),
    ANEMIA_TRIMESTER = as.factor(ANEMIA_TRIMESTER),
    age = as.numeric(age),
    SCHOOL_MORE10 = as.factor(SCHOOL_MORE10),
    water_improved = as.factor(water_improved),
    toilet_improved = as.factor(toilet_improved),
    WEALTH_QUINT_5 = as.factor(WEALTH_QUINT_5),
    WEALTH_QUINT_4 = as.factor(WEALTH_QUINT_4),
    WEALTH_QUINT_3 = as.factor(WEALTH_QUINT_3),
    WEALTH_QUINT_2 = as.factor(WEALTH_QUINT_2),
    WEALTH_QUINT_1 = as.factor(WEALTH_QUINT_1),
    SITE = case_when(
      SITE == "Ghana" ~ "2_Ghana",
      SITE == "Kenya" ~ "1_Kenya",
      SITE == "Pakistan" ~ "0_Pakistan",
      SITE == "Zambia" ~ "3_Zambia",
      SITE == "India-CMC" ~ "4_India.CMC",
      SITE == "India-SAS" ~ "5_India.SAS",
      TRUE ~ NA_character_
    )
  )

# GEE function
gee_function1 <- function(variable) {
  analysis_data <- cbind(analysis_data.aim3, long_data.aim3[, variable])
  analysis_data <- analysis_data %>% na.omit()
  
  formula_str <- paste("ANEMIA_OUTCOME1 ~", variable, "+ANEMIA_TRIMESTER+ SITE + age + SCHOOL_MORE10 + water_improved + toilet_improved + WEALTH_QUINT_4 + WEALTH_QUINT_3 + WEALTH_QUINT_2 + WEALTH_QUINT_1")
  
  # Calculate Numerator, Total, and Percentage
  Numerator <- sum(data.aim3[, variable] == 1 & data.aim3[, "ANEMIA_IND1"] == 1, na.rm = TRUE)
  Total <- sum(data.aim3[, variable] == 1, na.rm = TRUE)
  Percentage <- paste0(round(Numerator / Total * 100, 2), "%")
  
  # Define the result in case of error
  default_result <- c(Numerator, Total, Percentage, "incalculable", NA, NA)
  
  # Try running the GEE model and handle any errors
  result <- tryCatch({
    capture.output(gee_model <- gee(
      formula = as.formula(formula_str),
      id = MOMID,
      data = analysis_data,
      family = poisson(link = "log"),
      corstr = "exchangeable"
    ))
    model_summary <- summary(gee_model)
    
    variable_coeff <- round(model_summary$coefficients[variable, "Estimate"], 3)
    variable_se <- round(model_summary$coefficients[variable, "Robust S.E."], 3)
    
    RR <- round(exp(variable_coeff), 3)
    CI_lower <- round(exp(variable_coeff - 1.96 * variable_se), 3)
    CI_upper <- round(exp(variable_coeff + 1.96 * variable_se), 3)
    
    c(Numerator, Total, Percentage, RR, CI_lower, CI_upper)
  }, error = function(e) {
    # Return the default result in case of error
    default_result
  })
  
  return(result)
}

gee_function_distal1 <- function(variable) {
  analysis_data <- analysis_data.aim3 %>% mutate(
    SCHOOL_MORE10 = as.numeric(SCHOOL_MORE10),
    water_improved = as.numeric(water_improved),
    toilet_improved = as.numeric(toilet_improved),
    WEALTH_QUINT_5 = as.numeric(WEALTH_QUINT_5),
    WEALTH_QUINT_4 = as.numeric(WEALTH_QUINT_4),
    WEALTH_QUINT_3 = as.numeric(WEALTH_QUINT_3),
    WEALTH_QUINT_2 = as.numeric(WEALTH_QUINT_2),
    WEALTH_QUINT_1 = as.numeric(WEALTH_QUINT_1),
  ) %>% na.omit()
  
  formula_str <- paste("ANEMIA_OUTCOME1 ~ ANEMIA_TRIMESTER+ SITE + age + SCHOOL_MORE10 + water_improved + toilet_improved + WEALTH_QUINT_4 + WEALTH_QUINT_3 + WEALTH_QUINT_2 + WEALTH_QUINT_1")
  
  # Calculate Numerator, Total, and Percentage
  Numerator <- sum(data.aim3[, variable] == 1 & data.aim3[, "ANEMIA_IND1"] == 1, na.rm = TRUE)
  Total <- sum(data.aim3[, variable] == 1, na.rm = TRUE)
  Percentage <- paste0(round(Numerator / Total * 100, 2), "%")
  
  # Define the result in case of error
  default_result <- c(Numerator, Total, Percentage, "incalculable", NA, NA)
  
  # Try running the GEE model and handle any errors
  result <- tryCatch({
    capture.output(gee_model <- gee(
      formula = as.formula(formula_str),
      id = MOMID,
      data = analysis_data,
      family = poisson(link = "log"),
      corstr = "exchangeable"
    ))
    model_summary <- summary(gee_model)
    
    variable_coeff <- round(model_summary$coefficients[variable, "Estimate"], 3)
    variable_se <- round(model_summary$coefficients[variable, "Robust S.E."], 3)
    
    RR <- round(exp(variable_coeff), 3)
    CI_lower <- round(exp(variable_coeff - 1.96 * variable_se), 3)
    CI_upper <- round(exp(variable_coeff + 1.96 * variable_se), 3)
    
    c(Numerator, Total, Percentage, RR, CI_lower, CI_upper)
  }, error = function(e) {
    # Return the default result in case of error
    default_result
  })
  
  return(result)
}

# Summary function
summary_function0 <- function(variable){
  Numerator <- sum(data.aim3[, variable] == 0 & data.aim3[, "ANEMIA_IND1"] == 1, na.rm = TRUE)
  Total <- sum(data.aim3[, variable] == 0, na.rm = TRUE)
  Percentage <- paste0(round(Numerator / Total * 100, 2), "%")
  result <- c(Numerator, Total, Percentage, 1, NA, NA)
  return(result)
}

summary_function1<-function(variable){
  Numerator<-sum(data.aim3[,variable] == 1 & data.aim3[,"ANEMIA_IND1"] == 1, na.rm = TRUE)
  Total<-sum(data.aim3[,variable] == 1, na.rm = TRUE)
  Percentage<-paste0(round(sum(data.aim3[,variable] == 1 & data.aim3[,"ANEMIA_IND1"] == 1, na.rm = TRUE)/sum(data.aim3[,variable] == 1, na.rm = TRUE)*100,2),"%")
  result<-rbind(c(Numerator,Total,Percentage,1,NA,NA))
  return(result)
}

# List of variables and their corresponding row names
variables <- list(
  "### Nutrition" = "Nutrition",
  "## Ferritin at ANC20" = "Ferritin at ANC20",
  "FERRITIN70_ANC20_neg" = "Ferritin > 70ng/mL at ANC20",
  "FERRITIN70_ANC20_pos" = "Ferritin < 70ng/mL at ANC20",
  "FERRITIN15_ANC20_neg" = "Ferritin > 15ng/mL at ANC20",
  "FERRITIN15_ANC20_pos" = "Ferritin < 15ng/mL at ANC20",
  "## Serum transferrin receptor at ANC20 (high= likely iron deficiency)" = "Serum transferrin receptor at ANC20 (high= likely iron deficiency)",
  "STFR_ANC20_Normal_neg1" = "Normal Serum transferrin receptor at ANC20",
  "STFR_ANC20_Low_pos" = "Low Serum transferrin receptor at ANC20",
  "STFR_ANC20_High_pos" = "High Serum transferrin receptor at ANC20",
  "## High iodine thyroglobulin at ANC20" = "High iodine thyroglobulin at ANC20",
  "HIGH_TG_44_ANC20_neg" = "High thyroglobulin < 43.5 ug/L Tg at ANC20",
  "HIGH_TG_44_ANC20_pos" = "High thyroglobulin > 43.5 ug/L Tg at ANC20",
  "## RBP4 level at ANC20 (low indicates vitamin A deficiency)" = "RBP4 level at ANC20 (low indicates vitamin A deficiency)",
  "RBP4_ANC20_No_deficiency_neg1" = "No deficiency of RBP4 level at ANC20",
  "RBP4_ANC20_Mild_deficiency_pos" = "Mild deficiency of RBP4 level at ANC20",
  "RBP4_ANC20_Moderate_deficiency_pos" = "Moderate deficiency of RBP4 level at ANC20",
  "RBP4_ANC20_Severe_deficiency_pos" = "Severe deficiency of RBP4 level at ANC20",
  "## Serum B12 at ANC20" = "Serum B12 at ANC20",
  "VITB12_COB_ANC20_Sufficient_neg1" = "Sufficient serum B12 at ANC20",
  "VITB12_COB_ANC20_Deficient_pos" = "Deficient serum B12 at ANC20",
  "VITB12_COB_ANC20_Insufficient_pos" = "Insufficient serum B12 at ANC20",
  "## HoloTC at ANC20 (active form of vitamin B12)" = "HoloTC at ANC20 (active form of vitamin B12)",
  "VITB12_HOL_ANC20_neg" = "Normal holoTC at ANC20",
  "VITB12_HOL_ANC20_pos" = "Low holoTC at ANC20",
  "### RBC" = "RBC",
  "## Sickle cell" = "Sickle cell",
  "RBC_SICKLE_Normal_neg1" = "Sickle cell normal",
  "RBC_SICKLE_Disease_pos" = "Sickle cell disease",
  "RBC_SICKLE_Trait_pos" = "Sickle cell trait",
  "## Beta / alpha thalassemias" = "Beta / alpha thalassemias",
  "RBC_Thalassemia_Normal_neg1" = "Beta / alpha thalassemias Normal",
  "RBC_Thalassemia_Disease_pos" = "Beta / alpha thalassemias disease",
  "## G6PD" = "G6PD",
  "M08_RBC_G6PD_LBORRES_ind_neg" = "Normal G6PD (G6PD > 6.1 U/g Hb )",
  "M08_RBC_G6PD_LBORRES_ind_pos" = "Abnormal G6PD (G6PD < 6.1 U/g Hb) ",
  "## Mean cell volume (MCV)" = "Mean cell volume (MCV)",
  "MCV_ANC20_Normal_neg1" = "Normal MCV",
  "MCV_ANC20_Microcytic_pos" = "Microcytic MCV",
  "MCV_ANC20_Macrocytic_pos" = "Macrocytic MCV",
  "### Communicable disease" = "Communicable disease",
  "## HIV" = "HIV",
  "HIV_POSITIVE_ENROLL_neg" = "HIV Negative",
  "HIV_POSITIVE_ENROLL_pos" = "HIV Positive",
  "## Syphilis" = "Syphilis",
  "SYPH_POSITIVE_ENROLL_neg" = "Syphilis Negative",
  "SYPH_POSITIVE_ENROLL_pos" = "Syphilis Positive",
  "## Malaria" = "Malaria",
  "MAL_POSITIVE_ENROLL_neg" = "Malaria Negative",
  "MAL_POSITIVE_ENROLL_pos" = "Malaria Positive",
  "## Hep B" = "Hep B",
  "HBV_POSITIVE_ENROLL_neg" = "Hep B Negative",
  "HBV_POSITIVE_ENROLL_pos" = "Hep B Positive",
  "## Hep C" = "Hep C",
  "HCV_POSITIVE_ENROLL_neg" = "Hep C Negative",
  "HCV_POSITIVE_ENROLL_pos" = "Hep C Positive",
  "## TB" = "TB",
  "TB_SYMP_POSITIVE_ENROLL_neg" = "TB Negative",
  "TB_SYMP_POSITIVE_ENROLL_pos" = "TB Positive",
  "### Supplement" = "Supplement",
  "## Iron supplement" = "Iron supplement",
  "IRON_Supplement_neg" = "No Iron supplement",
  "IRON_Supplement_pos" = "With Iron supplement",
  "M04_IRON_ORAL_CMOCCUR_neg" = "No oral iron supplement",
  "M04_IRON_ORAL_CMOCCUR_pos" = "With oral iron supplement",
  "## Folic acid supplement" = "Folic acid supplement",
  "M04_IFA_CMOCCUR_neg" = "No folic acid supplement",
  "M04_IFA_CMOCCUR_pos" = "With folic acid supplement",
  "## Calcium supplement" = "Calcium supplement",
  "M04_CALCIUM_CMOCCUR_neg" = "No calcium supplement",
  "M04_CALCIUM_CMOCCUR_pos" = "With calcium supplement",
  "## Vitamin A supplement" = "Vitamin A supplement",
  "M04_VITAMIN_A_CMOCCUR_neg" = "No Vitamin A supplement",
  "M04_VITAMIN_A_CMOCCUR_pos" = "With Vitamin A supplement",
  "## Multiple micronutrient supplement" = "Multiple micronutrient supplement",
  "M04_MICRONUTRIENT_CMOCCUR_neg" = "No multiple micronutrient supplement",
  "M04_MICRONUTRIENT_CMOCCUR_pos" = "With multiple micronutrient supplement",
  "## Anthelmintic treatment" = "Anthelmintic treatment",
  "M04_ANTHELMINTHIC_CMOCCUR_neg" = "No anthelmintic treatment",
  "M04_ANTHELMINTHIC_CMOCCUR_pos" = "With anthelmintic treatment",
  "### Ambient air pollution" = "Ambient air pollution",
  "## Cooking fuel" = "Cooking fuel",
  "M03_STOVE_FCORRESR_ind_neg" = "Use of polluted cooking fuel",
  "M03_STOVE_FCORRESR_ind_pos" = "Use of clean cooking fuel",
  "## Smoking" = "Smoking",
  "hh_smoke_neg" = "No smoking inside the house",
  "hh_smoke_pos" = "Smoking inside the house",
  "### Inflammation" = "Inflammation",
  "## CRP level (high=inflammation)" = "CRP level (high=inflammation)",
  "CRP_ANC20_neg" = "Normal CRP level",
  "CRP_ANC20_pos" = "High CRP level (>5 mg/L)",
  "## AGP level at ANC20" = "AGP level at ANC20",
  "AGP_ANC20_neg" = "Normal AGP level",
  "AGP_ANC20_pos" = "High AGP level (>1 g/L)",
  "### Demographics" = "Demographics",
  "## Education" = "Education",
  "SCHOOL_MORE10_neg" = "Attended school less than 10 years",
  "SCHOOL_MORE10_pos1" = "Attended school more than 10 years",
  "## Water" = "Water",
  "water_improved_neg" = "Use of unimproved water",
  "water_improved_pos1" = "Use of improved water",
  "## Toilet" = "Toilet",
  "toilet_improved_neg" = "Use of unimproved toilet",
  "toilet_improved_pos1" = "Use of improved toilet",
  "## Wealth" = "Wealth",
  "WEALTH_QUINT_5_neg1" = "Wealth: Quintile 5",
  "WEALTH_QUINT_4_pos1" = "Wealth: Quintile 4",
  "WEALTH_QUINT_3_pos1" = "Wealth: Quintile 3",
  "WEALTH_QUINT_2_pos1" = "Wealth: Quintile 2",
  "WEALTH_QUINT_1_pos1" = "Wealth: Quintile 1"
)

# Initialize result table with different levels of headings
result_table <- data.frame(
  Heading = character(),
  Numerator = integer(),
  Total = integer(),
  Percentage = character(),
  RR = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

# Collect results
current_heading1 <- ""
current_heading2 <- ""
for (variable in names(variables)) {
  if (startsWith(variable, "###")) {
    current_heading1 <- variables[[variable]]
    print(paste("Adding Heading 1:", current_heading1))
    result_table <- rbind(result_table, data.frame(Heading = current_heading1, Numerator = NA, Total = NA, Percentage = NA, RR = NA, CI_Lower = NA, CI_Upper = NA, stringsAsFactors = FALSE))
  } else if (startsWith(variable, "##")) {
    current_heading2 <- variables[[variable]]
    print(paste("Adding Heading 2:", current_heading2))
    result_table <- rbind(result_table, data.frame(Heading = current_heading2, Numerator = NA, Total = NA, Percentage = NA, RR = NA, CI_Lower = NA, CI_Upper = NA, stringsAsFactors = FALSE))
  } else if (endsWith(variable, "_neg")) {
    result <- summary_function0(sub("_neg$", "", variable))
    print(paste("Adding Summary for:", variables[[variable]]))
    result_table <- rbind(result_table, data.frame(Heading = variables[[variable]], Numerator = result[1], Total = result[2], Percentage = result[3], RR = result[4], CI_Lower = result[5], CI_Upper = result[6], stringsAsFactors = FALSE))
  } else if (endsWith(variable, "_neg1")) {
    result <- summary_function1(sub("_neg1$", "", variable))
    print(paste("Adding Summary for:", variables[[variable]]))
    result_table <- rbind(result_table, data.frame(Heading = variables[[variable]], Numerator = result[1], Total = result[2], Percentage = result[3], RR = result[4], CI_Lower = result[5], CI_Upper = result[6], stringsAsFactors = FALSE))
  } else if (endsWith(variable, "_pos1")) {
    result <- gee_function_distal1(sub("_pos1$", "", variable))
    print(paste("Adding Summary for:", variables[[variable]]))
    result_table <- rbind(result_table, data.frame(Heading = variables[[variable]], Numerator = result[1], Total = result[2], Percentage = result[3], RR = result[4], CI_Lower = result[5], CI_Upper = result[6], stringsAsFactors = FALSE))
  } 
  else {
    result <- gee_function1(sub("_pos$", "", variable))
    print(paste("Adding GEE for:", variables[[variable]]))
    result_table <- rbind(result_table, data.frame(Heading = variables[[variable]], Numerator = result[1], Total = result[2], Percentage = result[3], RR = result[4], CI_Lower = result[5], CI_Upper = result[6], stringsAsFactors = FALSE))
  }
}

# Convert any NA values in the heading columns to empty strings
result_table$Heading[is.na(result_table$Heading)] <- ""

# Write the result table to an Excel file with appropriate formatting
wb <- createWorkbook()
addWorksheet(wb, "Results")
writeData(wb, "Results", result_table, startRow = 1, startCol = 1, colNames = TRUE)

# Apply formatting
headerStyle1 <- createStyle(fontSize = 14, textDecoration = "bold")
headerStyle2 <- createStyle(fontSize = 12, textDecoration = "bold")
headerStyle3 <- createStyle(fontSize = 10, textDecoration = "bold")

for (i in 1:nrow(result_table)) {
  if (!is.na(result_table$Numerator[i])) {
    next
  }
  if (nchar(result_table$Heading[i]) > 0 && substr(result_table$Heading[i], 1, 1) == " ") {
    if (startsWith(result_table$Heading[i], "##")) {
      addStyle(wb, sheet = "Results", headerStyle2, rows = i + 1, cols = 1, gridExpand = TRUE)
    } else if (startsWith(result_table$Heading[i], "#")) {
      addStyle(wb, sheet = "Results", headerStyle1, rows = i + 1, cols = 1, gridExpand = TRUE)
    } else {
      addStyle(wb, sheet = "Results", headerStyle3, rows = i + 1, cols = 1, gridExpand = TRUE)
    }
  } else {
    addStyle(wb, sheet = "Results", headerStyle3, rows = i + 1, cols = 1, gridExpand = TRUE)
  }
}

saveWorkbook(wb, "D:/Users/yipeng_wei/Documents/Output/ReMAPP aim3 output/aim3_output1.xlsx", overwrite = TRUE)

# Print the result table for verification
print(result_table)