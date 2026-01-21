#****************************************************************************
#Maternal Depression Likelihood
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_dpr.rda")

# Create Antenatal (ANC) and Postnatal (PNC) datasets
df_anc <- df_mat_dpr %>%
  filter(VISIT %in% c("ANC20", "ANC32")) %>%
  arrange(ID, VISIT)

df_pnc <- df_mat_dpr %>%
  filter(VISIT == "PNC6")

# ****************************************************************************
# Antenatal Models (ANC)

# Spline model
spline_dpr_anc <- knot_fun(df_anc, "hb", "dpr")
saveRDS(spline_dpr_anc, "iso_results/spline_dpr_anc.rds")

# Isotonic model
iso_dpr_anc <- flexstepreg_glmer(df_anc$dpr, df_anc$hb, df_anc$SITE,
                                 covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_dpr_anc, "iso_results/iso_dpr_anc.rds")



# ****************************************************************************
# Postnatal Models (PNC)

# Spline model
spline_dpr_pnc <- knot_fun(df_pnc, "hb", "dpr")
saveRDS(spline_dpr_pnc, "iso_results/spline_dpr_pnc.rds")

# Isotonic model
iso_dpr_pnc <- flexstepreg_glmer(df_pnc$dpr, df_pnc$hb, df_pnc$SITE,
                                 covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_dpr_pnc, "iso_results/iso_dpr_pnc.rds")



# Load dataset
load("derived_data/df_mat_dpr_trim1.rda")

# Split antenatal (ANC20 + ANC32) and postpartum (PNC6)
df_anc_trim1 <- df_mat_dpr_trim1 %>% filter(VISIT %in% c("ANC20", "ANC32"))
df_pnc_trim1 <- df_mat_dpr_trim1 %>% filter(VISIT == "PNC6")

# ****************************************************************************
# Antenatal (ANC20 + ANC32) models -- NO special suffix

spline_dpr_trim1 <- knot_fun(df_anc_trim1, "hb", "dpr")
saveRDS(spline_dpr_trim1, "iso_results/spline_dpr_trim1.rds")

iso_dpr_trim1 <- flexstepreg_glmer(df_anc_trim1$dpr, df_anc_trim1$hb, df_anc_trim1$SITE,
                                   covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_dpr_trim1, "iso_results/iso_dpr_trim1.rds")

# ****************************************************************************
# Postpartum (PNC6) models -- PNC suffix

spline_dpr_trim1_pnc <- knot_fun(df_pnc_trim1, "hb", "dpr")
saveRDS(spline_dpr_trim1_pnc, "iso_results/spline_dpr_trim1_pnc.rds")

iso_dpr_trim1_pnc <- flexstepreg_glmer(df_pnc_trim1$dpr, df_pnc_trim1$hb, df_pnc_trim1$SITE,
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_dpr_trim1_pnc, "iso_results/iso_dpr_trim1_pnc.rds")


# Load dataset
load("derived_data/df_mat_dpr_trim2.rda")

# Split antenatal (ANC20 + ANC32) and postpartum (PNC6)
df_anc_trim2 <- df_mat_dpr_trim2 %>% filter(VISIT %in% c("ANC32"))
df_pnc_trim2 <- df_mat_dpr_trim2 %>% filter(VISIT == "PNC6")

# ****************************************************************************
# Antenatal (ANC20 + ANC32) models -- NO special suffix

spline_dpr_trim2 <- knot_fun(df_anc_trim2, "hb", "dpr")
saveRDS(spline_dpr_trim2, "iso_results/spline_dpr_trim2.rds")

iso_dpr_trim2 <- flexstepreg_glmer(df_anc_trim2$dpr, df_anc_trim2$hb, df_anc_trim2$SITE,
                                   covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_dpr_trim2, "iso_results/iso_dpr_trim2.rds")

# ****************************************************************************
# Postpartum (PNC6) models -- PNC suffix

spline_dpr_trim2_pnc <- knot_fun(df_pnc_trim2, "hb", "dpr")
saveRDS(spline_dpr_trim2_pnc, "iso_results/spline_dpr_trim2_pnc.rds")

iso_dpr_trim2_pnc <- flexstepreg_glmer(df_pnc_trim2$dpr, df_pnc_trim2$hb, df_pnc_trim2$SITE,
                                       covar2 = NULL, random_effect = NULL, alpha = 0.01)
saveRDS(iso_dpr_trim2_pnc, "iso_results/iso_dpr_trim2_pnc.rds")






# #run spline model
# spline_dpr <- knot_fun(df_mat_dpr, "hb", "dpr")
# saveRDS(spline_dpr, "iso_results/spline_dpr.rds")
# 
# #run isotonic model
# iso_dpr <- flexstepreg_glmer(df_mat_dpr$dpr, df_mat_dpr$hb, df_mat_dpr$SITE, 
#                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
# saveRDS(iso_dpr, "iso_results/iso_dpr.rds")
# 
# #****************************************************************************
# #*******trim1
# #prepare data
# load("derived_data/df_mat_dpr_trim1.rda")
# 
# #run spline model
# spline_dpr_trim1 <- knot_fun(df_mat_dpr_trim1, "hb", "dpr")
# saveRDS(spline_dpr_trim1, "iso_results/spline_dpr_trim1.rds")
# 
# #run isotonic model
# iso_dpr_trim1 <- flexstepreg_glmer(df_mat_dpr_trim1$dpr, df_mat_dpr_trim1$hb, df_mat_dpr_trim1$SITE,
#                              covar2=NULL, random_effect=NULL, alpha = 0.01)
# saveRDS(iso_dpr_trim1, "iso_results/iso_dpr_trim1.rds")
# 
# #****************************************************************************
# #*******trim2
# #prepare data
# load("derived_data/df_mat_dpr_trim2.rda")
# 
# #run spline model
# spline_dpr_trim2 <- knot_fun(df_mat_dpr_trim2, "hb", "dpr")
# saveRDS(spline_dpr_trim2, "iso_results/spline_dpr_trim2.rds")
# 
# #run isotonic model
# iso_dpr_trim2 <- flexstepreg_glmer(df_mat_dpr_trim2$dpr, df_mat_dpr_trim2$hb, df_mat_dpr_trim2$SITE,
#                                    covar2=NULL, random_effect=NULL, alpha = 0.01)
# saveRDS(iso_dpr_trim2, "iso_results/iso_dpr_trim2.rds")
# 
# #****************************************************************************
# #*******trim3
# #prepare data
# load("derived_data/df_mat_dpr_trim3.rda")
# 
# #run spline model
# spline_dpr_trim3 <- knot_fun(df_mat_dpr_trim3, "hb", "dpr")
# saveRDS(spline_dpr_trim3, "iso_results/spline_dpr_trim3.rds")
# 
# #run isotonic model
# iso_dpr_trim3 <- flexstepreg_glmer(df_mat_dpr_trim3$dpr, df_mat_dpr_trim3$hb, df_mat_dpr_trim3$SITE,
#                                    covar2=NULL, random_effect=NULL, alpha = 0.01)
# saveRDS(iso_dpr_trim3, "iso_results/iso_dpr_trim3.rds")
