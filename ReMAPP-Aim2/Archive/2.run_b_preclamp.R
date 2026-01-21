#****************************************************************************
#Preterm premature rupture of membranes
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")
#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_preclamp.rda")

#run spline model
spline_preclamp <- knot_fun(df_mat_preclamp, "hb", "preclamp")
saveRDS(spline_preclamp, "iso_results/spline_preclamp.rds")

#run isotonic model
iso_preclamp <- flexstepreg_glmer(df_mat_preclamp$preclamp, df_mat_preclamp$hb, df_mat_preclamp$SITE, 
                               covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preclamp, "iso_results/iso_preclamp.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_mat_preclamp_trim1.rda")
#run spline model
spline_preclamp_trim1 <- knot_fun(df_mat_preclamp_trim1, "hb", "preclamp")
saveRDS(spline_preclamp_trim1, "iso_results/spline_preclamp_trim1.rds")

#run isotonic model
iso_preclamp_trim1 <- flexstepreg_glmer(df_mat_preclamp_trim1$preclamp, df_mat_preclamp_trim1$hb, df_mat_preclamp_trim1$SITE, 
                                  covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preclamp_trim1, "iso_results/iso_preclamp_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_mat_preclamp_trim2.rda")

#run spline model
spline_preclamp_trim2 <- knot_fun(df_mat_preclamp_trim2, "hb", "preclamp")
saveRDS(spline_preclamp_trim2, "iso_results/spline_preclamp_trim2.rds")

#run isotonic model
iso_preclamp_trim2 <- flexstepreg_glmer(df_mat_preclamp_trim2$preclamp, df_mat_preclamp_trim2$hb, df_mat_preclamp_trim2$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preclamp_trim2, "iso_results/iso_preclamp_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_mat_preclamp_trim3.rda")

#run spline model
spline_preclamp_trim3 <- knot_fun(df_mat_preclamp_trim3, "hb", "preclamp")
saveRDS(spline_preclamp_trim3, "iso_results/spline_preclamp_trim3.rds")

#run isotonic model
iso_preclamp_trim3 <- flexstepreg_glmer(df_mat_preclamp_trim3$preclamp, df_mat_preclamp_trim3$hb, df_mat_preclamp_trim3$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preclamp_trim3, "iso_results/iso_preclamp_trim3.rds")
