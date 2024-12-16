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

N_perc_non_missing_fun<-function(data,variable){
  paste0(nrow(data),"(",round(sum(!is.na(variable))/nrow(data),4)*100,"%",")")
}


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

# Define the function to calculate frequency and percentage of non-missing values
calculate_frequencies <- function(data, column_name, site_column, site_names) {
  # Frequency table
  Fre_N_disease <- table(data[[column_name]], data[[site_column]])
  
  # Calculate non-missing percentages for each site
  Non_Missing_Perc <- sapply(site_names, function(site) {
    site_data <- data[data[[site_column]] == site,]
    N_perc_non_missing_fun(site_data, site_data[[column_name]])
  })
  
  # Combine percentages and frequency table
  result <- rbind(Non_Missing_Perc, Fre_N_disease)
  return(result)
}

#Site names
site_names <- c("Ghana", "India-CMC", "India-SAS", "Kenya", "Pakistan", "Zambia")

###Import data
UploadDate = "2024-06-28"
path_to_data <- paste0('D:/Users/yipeng_wei/Documents/Stacked data/',UploadDate)# - for AWS data

uploadDate <- as.Date(UploadDate, format = "%Y-%m-%d")
setwd(path_to_data)

###Load data
########################################################################################
recreate_mnh_list <- FALSE #if I want to recreate the big data list then set this to TRUE
########################################################################################
tic <- Sys.time()
if(recreate_mnh_list == TRUE){
  mnh_names <- c()
  list_mnh <- dir(path = path_to_data, pattern = "*.csv", full.names = TRUE) #creates a list of all the csv files in the directory
  for (data_file in list_mnh[]) { #can test by just bringing in a small number (add 1:2 inside the bracket to do so)
    form_name <- substr(basename(data_file), 1,5) #substr pulls out the 1:5 spaces in a char (will pull out "mnh00" etc);
    #basename() pulls out just the name of the file from the entire directory/path.
    print(paste("Reading", form_name))
    assign(form_name, read.csv(data_file))
    mnh_names <- c(mnh_names, form_name) # it's buidling up a vector of obj names.
  }
  save(list=mnh_names, file = paste0("MNH_dfs_", basename(path_to_data), ".RDS"))
} else{
  load(file = paste0('MNH_dfs_',  basename(path_to_data), ".RDS"))
}

toc <- Sys.time()
print(toc - tic)

### Hard code TYPE_VISIT for M00=1, M02=1, M03=1, M09=6, M10=6, M11=6, M17=6, M18=12 ***
if (exists("mnh00")) {
  mnh00 <- mnh00 %>% 
    mutate(M00_TYPE_VISIT = 1)
} else{
  "MNH00 dataframe doesn't exist"
}

if (exists("mnh02")) {
  mnh02 <- mnh02 %>% 
    mutate(M02_TYPE_VISIT = 1)
} else{
  "MNH02 dataframe doesn't exist"
}

if (exists("mnh03")) {
  mnh03 <- mnh03 %>% 
    mutate(M03_TYPE_VISIT = 1)
} else{
  "MNH03 dataframe doesn't exist"
}

if (exists("mnh09")) {
  mnh09 <- mnh09 %>% 
    mutate(M09_TYPE_VISIT = 6)
} else{
  "MNH09 dataframe doesn't exist"
}

if (exists("mnh10")) {
  mnh10 <- mnh10 %>% 
    mutate(M10_TYPE_VISIT = 6)
} else{
  "MNH10 dataframe doesn't exist"
}

if (exists("mnh11")) {
  mnh11 <- mnh11 %>% 
    mutate(M11_TYPE_VISIT = 6)
} else{
  "MNH11 dataframe doesn't exist"
}

if (exists("mnh16")) {
  mnh16 <- mnh16 %>% 
    mutate(M16_TYPE_VISIT = 5)
} else{
  "MNH16 dataframe doesn't exist"
}

if (exists("mnh17")) {
  mnh17 <- mnh17 %>% 
    mutate(M17_TYPE_VISIT = 6)
} else{
  "MNH17 dataframe doesn't exist"
}

if (exists("mnh18")) {
  mnh18 <- mnh18 %>% 
    mutate(M18_TYPE_VISIT = 12)
} else{
  "MNH18 dataframe doesn't exist"
}

if (exists("mnh19")) {
  mnh19 <- mnh19 %>% 
    mutate(M19_TYPE_VISIT = 13)
} else{
  "MNH19 dataframe doesn't exist"
}

###Create a vector of women who are enrolled. 
mnh02 <- mnh02 %>% 
  mutate(ENROLLED_FF = ifelse(M02_AGE_IEORRES==1 & M02_PC_IEORRES==1 & 
                                M02_CATCHMENT_IEORRES==1 & M02_CATCH_REMAIN_IEORRES==1 & 
                                M02_CONSENT_IEORRES==1, 1,0),
         remapp = case_when((SITE == "Ghana" & M02_SCRN_OBSSTDAT >= "2022-12-28") |
    (SITE == "Kenya" & M02_SCRN_OBSSTDAT >= "2023-04-14") |
             (SITE == "Zambia" & M02_SCRN_OBSSTDAT >= "2022-12-15") |
    (SITE == "Pakistan" & M02_SCRN_OBSSTDAT >= "2022-09-22" & M02_SCRN_OBSSTDAT <= "2024-04-05") |
             (SITE == "India-CMC" & M02_SCRN_OBSSTDAT >= "2023-06-20") |
             (SITE == "India-SAS" & M02_SCRN_OBSSTDAT >= "2023-08-15") ~ 1,
           TRUE ~ NA_real_)) %>% filter(ENROLLED_FF ==1 & remapp==1) # THIS DOES DROP SOME WOMEN!

unique <- mnh02 %>% group_by(SITE, PREGID) %>% unique()

enrolled_ids_vec <- as.vector(mnh02$PREGID)


###Rename TYPE_VISIT variables:
mnh02 <- mnh02 %>% rename(TYPE_VISIT = M02_TYPE_VISIT)
mnh03 <- mnh03 %>% rename(TYPE_VISIT = M03_TYPE_VISIT)
mnh08 <- mnh08 %>% rename(TYPE_VISIT = M08_TYPE_VISIT)
mnh04 <- mnh04 %>% rename(TYPE_VISIT = M04_TYPE_VISIT)

###Clean up the data frames: 
mnh02 <- mnh02 %>%
  filter(MOMID!="") %>% #Moved filter to the end so we can create all the variables with TYPE_VISIT at the end of the var name before rows get dropped out
  drop_na(MOMID, PREGID) %>%
  distinct (MOMID, PREGID, SITE, TYPE_VISIT, .keep_all = TRUE) %>% select(-TYPE_VISIT)

mnh03 <- mnh03 %>%
  filter(MOMID!="") %>% #Moved filter to the end so we can create all the variables with TYPE_VISIT at the end of the var name before rows get dropped out
  drop_na(MOMID, PREGID) %>%
  distinct (MOMID, PREGID, SITE, TYPE_VISIT, .keep_all = TRUE)

mnh08 <- mnh08 %>%
  filter(MOMID!="") %>% #Moved filter to the end so we can create all the variables with TYPE_VISIT at the end of the var name before rows get dropped out
  drop_na(MOMID, PREGID) %>%
  distinct (MOMID, PREGID, SITE, TYPE_VISIT, .keep_all = TRUE)

mnh04 <- mnh04 %>%
  filter(MOMID!="") %>% #Moved filter to the end so we can create all the variables with TYPE_VISIT at the end of the var name before rows get dropped out
  drop_na(MOMID, PREGID) %>%
  distinct (MOMID, PREGID, SITE, TYPE_VISIT, .keep_all = TRUE)

###Merge data
merged_df08 <- left_join(mnh02, mnh08, by = c("SITE", "MOMID", "PREGID"))
merged_df03 <- left_join(mnh02, mnh03, by = c("SITE", "MOMID", "PREGID"))
merged_df04 <- left_join(mnh02, mnh04, by = c("SITE", "MOMID", "PREGID"))

merged_df08 <- merged_df08 %>% 
  mutate(
    #RBC_SICKLE
    M08_RBC_SICKLE_LBORRES = transform_to_NA(M08_RBC_SICKLE_LBORRES),
    M08_RBC_THALA_1 = transform_to_NA(M08_RBC_THALA_1),
    M08_RBC_THALA_3 = transform_to_NA(M08_RBC_THALA_3),
    M08_RBC_THALA_5 = transform_to_NA(M08_RBC_THALA_5),
    M08_RBC_THALA_15 = transform_to_NA(M08_RBC_THALA_15),
    M08_RBC_THALA_16 = transform_to_NA(M08_RBC_THALA_16),
    #RBC_THALA
    M08_RBC_THALA_6 = transform_to_NA(M08_RBC_THALA_6),
    M08_RBC_THALA_7 = transform_to_NA(M08_RBC_THALA_7),
    M08_RBC_THALA_8 = transform_to_NA(M08_RBC_THALA_8),
    M08_RBC_THALA_11 = transform_to_NA(M08_RBC_THALA_11),
    M08_RBC_THALA_12 = transform_to_NA(M08_RBC_THALA_12),
    M08_RBC_THALA_13 = transform_to_NA(M08_RBC_THALA_13),
    M08_RBC_THALA_14 = transform_to_NA(M08_RBC_THALA_14),
    #RBC_G6PD
    M08_RBC_G6PD_LBORRES = transform_to_NA(M08_RBC_G6PD_LBORRES),
    M08_RBC_G6PD_LBORRES_ind = case_when(
      M08_RBC_G6PD_LBORRES < 6.1 ~ 1,
      M08_RBC_G6PD_LBORRES >= 6.1 ~ 0,
      TRUE ~ NA_real_))

merged_df04 <- merged_df04 %>% 
  mutate(
    #Iron supplement
    M04_IRON_CMOCCUR = transform_to_NA(M04_IRON_CMOCCUR),
    M04_IRON_ORAL_CMOCCUR = transform_to_NA(M04_IRON_ORAL_CMOCCUR),
    M04_IRON_IV_CMOCCUR = transform_to_NA(M04_IRON_IV_CMOCCUR),
    #Folic acid supplement
    M04_IFA_CMOCCUR = transform_to_NA(M04_IFA_CMOCCUR),
    #Calcium supplement
    M04_CALCIUM_CMOCCUR = transform_to_NA(M04_CALCIUM_CMOCCUR),
    #Vitamin A supplement
    M04_VITAMIN_A_CMOCCUR = transform_to_NA(M04_VITAMIN_A_CMOCCUR),
    #Multiple micronutrient supplement
    M04_MICRONUTRIENT_CMOCCUR = transform_to_NA(M04_MICRONUTRIENT_CMOCCUR),
    #Anthelmintic treatment
    M04_ANTHELMINTHIC_CMOCCUR = transform_to_NA(M04_ANTHELMINTHIC_CMOCCUR),
  )

merged_df03 <- merged_df03 %>% 
  mutate(
    #clean fuel
    M03_STOVE_FCORRESR = transform_to_NA(M03_STOVE_FCORRES),
    M03_STOVE_FCORRESR_ind = case_when(
      M03_STOVE_FCORRESR %in% c(1,2,3,4,5) ~ 1,
      M03_STOVE_FCORRESR %in% c(7,8,9,11) ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_1==1 ~ 1,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_2==1 ~ 1,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_3==1 ~ 1,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_4==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_5==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_6==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_7==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_8==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_9==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_10==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_11==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_12==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_13==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_14==1 ~ 0,
      M03_STOVE_FCORRESR==6 & M03_STOVE_FUEL_FCORRES_15==1 ~ 1,
      TRUE ~ NA_real_),
  )

###Build up dataset
#supplement
df_supplement<-merged_df04 %>% select("SCRNID", "MOMID", "PREGID", "SITE",TYPE_VISIT,
  #Iron supplement
  M04_IRON_CMOCCUR,
  M04_IRON_ORAL_CMOCCUR,
  M04_IRON_IV_CMOCCUR,
  #Folic acid supplement
  M04_IFA_CMOCCUR,
  #Calcium supplement
  M04_CALCIUM_CMOCCUR,
  #Vitamin A supplement
  M04_VITAMIN_A_CMOCCUR,
  #Multiple micronutrient supplement
  M04_MICRONUTRIENT_CMOCCUR,
  #Anthelmintic treatment
  M04_ANTHELMINTHIC_CMOCCUR
  )%>%filter(TYPE_VISIT==1)

df_RBC<-merged_df08 %>% select(
  "SCRNID", "MOMID", "PREGID", "SITE","TYPE_VISIT",
  #RBC_SICKLE
  M08_RBC_SICKLE_LBORRES,
  M08_RBC_THALA_1,
  M08_RBC_THALA_3,
  M08_RBC_THALA_5,
  M08_RBC_THALA_15,
  M08_RBC_THALA_16,
  #RBC_THALA
  M08_RBC_THALA_6,
  M08_RBC_THALA_7,
  M08_RBC_THALA_8,
  M08_RBC_THALA_11,
  M08_RBC_THALA_12,
  M08_RBC_THALA_13,
  M08_RBC_THALA_14,
  #RBC_G6PD
  M08_RBC_G6PD_LBORRES,
  M08_RBC_G6PD_LBORRES_ind,
)%>%filter(TYPE_VISIT==1)

df_fuel<-merged_df03 %>% select(
  "SCRNID", "MOMID", "PREGID", "SITE",M03_STOVE_FCORRESR,M03_STOVE_FCORRESR_ind
)

write.csv(df_supplement,"D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_SUPPLEMENT.csv")
write.csv(df_RBC,"D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_RBC.csv")
write.csv(df_fuel,"D:/Users/yipeng_wei/Documents/ReMAPP aim 3 data/2024-06-28/MAT_FUEL.csv")
