#****************************************************************************
#Preterm <37 weeks
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_preterm37.rda")

#run spline model
spline_preterm37 <- knot_fun(df_inf_preterm37, "hb", "preterm37")
saveRDS(spline_preterm37, "iso_results/spline_preterm37.rds")

#run isotonic model
iso_preterm37 <- flexstepreg_glmer(df_inf_preterm37$preterm37, df_inf_preterm37$hb, df_inf_preterm37$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm37, "iso_results/iso_preterm37.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_preterm37_trim1.rda")

#run spline model
spline_preterm37_trim1 <- knot_fun(df_inf_preterm37_trim1, "hb", "preterm37")
saveRDS(spline_preterm37_trim1, "iso_results/spline_preterm37_trim1.rds")

#run isotonic model
iso_preterm37_trim1 <- flexstepreg_glmer(df_inf_preterm37_trim1$preterm37, df_inf_preterm37_trim1$hb, df_inf_preterm37_trim1$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm37_trim1, "iso_results/iso_preterm37_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_preterm37_trim2.rda")

#run spline model
spline_preterm37_trim2 <- knot_fun(df_inf_preterm37_trim2, "hb", "preterm37")
saveRDS(spline_preterm37_trim2, "iso_results/spline_preterm37_trim2.rds")

#run isotonic model
iso_preterm37_trim2 <- flexstepreg_glmer(df_inf_preterm37_trim2$preterm37, df_inf_preterm37_trim2$hb, df_inf_preterm37_trim2$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm37_trim2, "iso_results/iso_preterm37_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_preterm37_trim3.rda")

#run spline model
spline_preterm37_trim3 <- knot_fun(df_inf_preterm37_trim3, "hb", "preterm37")
saveRDS(spline_preterm37_trim3, "iso_results/spline_preterm37_trim3.rds")

#run isotonic model
iso_preterm37_trim3 <- flexstepreg_glmer(df_inf_preterm37_trim3$preterm37, df_inf_preterm37_trim3$hb, df_inf_preterm37_trim3$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm37_trim3, "iso_results/iso_preterm37_trim3.rds")

