#****************************************************************************
#Possible severe bacterial infection
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_psbi.rda")

#run spline model
spline_psbi <- knot_fun(df_inf_psbi, "hb", "inf_psbi")
saveRDS(spline_psbi, "iso_results/spline_psbi.rds")

#run isotonic model
iso_psbi <- flexstepreg_glmer(df_inf_psbi$inf_psbi, df_inf_psbi$hb, df_inf_psbi$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_psbi, "iso_results/iso_psbi.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_psbi_trim1.rda")

#run spline model
spline_psbi_trim1 <- knot_fun(df_inf_psbi_trim1, "hb", "inf_psbi")
saveRDS(spline_psbi_trim1, "iso_results/spline_psbi_trim1.rds")

#run isotonic model
iso_psbi_trim1 <- flexstepreg_glmer(df_inf_psbi_trim1$inf_psbi, df_inf_psbi_trim1$hb, df_inf_psbi_trim1$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_psbi_trim1, "iso_results/iso_psbi_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_psbi_trim2.rda")

#run spline model
spline_psbi_trim2 <- knot_fun(df_inf_psbi_trim2, "hb", "inf_psbi")
saveRDS(spline_psbi_trim2, "iso_results/spline_psbi_trim2.rds")

#run isotonic model
iso_psbi_trim2 <- flexstepreg_glmer(df_inf_psbi_trim2$inf_psbi, df_inf_psbi_trim2$hb, df_inf_psbi_trim2$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_psbi_trim2, "iso_results/iso_psbi_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_psbi_trim3.rda")

#run spline model
spline_psbi_trim3 <- knot_fun(df_inf_psbi_trim3, "hb", "inf_psbi")
saveRDS(spline_psbi_trim3, "iso_results/spline_psbi_trim3.rds")

#run isotonic model
iso_psbi_trim3 <- flexstepreg_glmer(df_inf_psbi_trim3$inf_psbi, df_inf_psbi_trim3$hb, df_inf_psbi_trim3$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_psbi_trim3, "iso_results/iso_psbi_trim3.rds")

