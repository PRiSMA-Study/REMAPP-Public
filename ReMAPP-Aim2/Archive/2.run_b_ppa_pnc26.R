#****************************************************************************
#Maternal postpartum anemia at Pnc26
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_ppa_pnc26.rda")

#run spline model
spline_ppa_pnc26 <- knot_fun(df_mat_ppa_pnc26, "hb", "ppa_pnc26")
saveRDS(spline_ppa_pnc26, "iso_results/spline_ppa_pnc26.rds")

#run isotonic model
iso_ppa_pnc26 <- flexstepreg_glmer(df_mat_ppa_pnc26$ppa_pnc26, df_mat_ppa_pnc26$hb, df_mat_ppa_pnc26$SITE, 
                                  covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc26, "iso_results/iso_ppa_pnc26.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_mat_ppa_pnc26_trim1.rda")

#run spline model
spline_ppa_pnc26_trim1 <- knot_fun(df_mat_ppa_pnc26_trim1, "hb", "ppa_pnc26")
saveRDS(spline_ppa_pnc26_trim1, "iso_results/spline_ppa_pnc26_trim1.rds")

#run isotonic model
iso_ppa_pnc26_trim1 <- flexstepreg_glmer(df_mat_ppa_pnc26_trim1$ppa_pnc26, df_mat_ppa_pnc26_trim1$hb, df_mat_ppa_pnc26_trim1$SITE, 
                                        covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc26_trim1, "iso_results/iso_ppa_pnc26_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_mat_ppa_pnc26_trim2.rda")

#run spline model
spline_ppa_pnc26_trim2 <- knot_fun(df_mat_ppa_pnc26_trim2, "hb", "ppa_pnc26")
saveRDS(spline_ppa_pnc26_trim2, "iso_results/spline_ppa_pnc26_trim2.rds")

#run isotonic model
iso_ppa_pnc26_trim2 <- flexstepreg_glmer(df_mat_ppa_pnc26_trim2$ppa_pnc26, df_mat_ppa_pnc26_trim2$hb, df_mat_ppa_pnc26_trim2$SITE, 
                                        covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc26_trim2, "iso_results/iso_ppa_pnc26_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_mat_ppa_pnc26_trim3.rda")

#run spline model
spline_ppa_pnc26_trim3 <- knot_fun(df_mat_ppa_pnc26_trim3, "hb", "ppa_pnc26")
saveRDS(spline_ppa_pnc26_trim3, "iso_results/spline_ppa_pnc26_trim3.rds")

#run isotonic model
iso_ppa_pnc26_trim3 <- flexstepreg_glmer(df_mat_ppa_pnc26_trim3$ppa_pnc26, df_mat_ppa_pnc26_trim3$hb, df_mat_ppa_pnc26_trim3$SITE, 
                                        covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc26_trim3, "iso_results/iso_ppa_pnc26_trim3.rds")
