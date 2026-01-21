#****************************************************************************
#Maternal postpartum anemia at PNC6
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")
#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_ppa_pnc6.rda")

#run spline model
spline_ppa_pnc6 <- knot_fun(df_mat_ppa_pnc6, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6, "iso_results/spline_ppa_pnc6.rds")

#run isotonic model
iso_ppa_pnc6 <- flexstepreg_glmer(df_mat_ppa_pnc6$ppa_pnc6, df_mat_ppa_pnc6$hb, df_mat_ppa_pnc6$SITE, 
                                 covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc6, "iso_results/iso_ppa_pnc6.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_mat_ppa_pnc6_trim1.rda")

#run spline model
spline_ppa_pnc6_trim1 <- knot_fun(df_mat_ppa_pnc6_trim1, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6_trim1, "iso_results/spline_ppa_pnc6_trim1.rds")

#run isotonic model
iso_ppa_pnc6_trim1 <- flexstepreg_glmer(df_mat_ppa_pnc6_trim1$ppa_pnc6, df_mat_ppa_pnc6_trim1$hb, df_mat_ppa_pnc6_trim1$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc6_trim1, "iso_results/iso_ppa_pnc6_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_mat_ppa_pnc6_trim2.rda")

#run spline model
spline_ppa_pnc6_trim2 <- knot_fun(df_mat_ppa_pnc6_trim2, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6_trim2, "iso_results/spline_ppa_pnc6_trim2.rds")

#run isotonic model
iso_ppa_pnc6_trim2 <- flexstepreg_glmer(df_mat_ppa_pnc6_trim2$ppa_pnc6, df_mat_ppa_pnc6_trim2$hb, df_mat_ppa_pnc6_trim2$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc6_trim2, "iso_results/iso_ppa_pnc6_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_mat_ppa_pnc6_trim3.rda")

#run spline model
spline_ppa_pnc6_trim3 <- knot_fun(df_mat_ppa_pnc6_trim3, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6_trim3, "iso_results/spline_ppa_pnc6_trim3.rds")

#run isotonic model
iso_ppa_pnc6_trim3 <- flexstepreg_glmer(df_mat_ppa_pnc6_trim3$ppa_pnc6, df_mat_ppa_pnc6_trim3$hb, df_mat_ppa_pnc6_trim3$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc6_trim3, "iso_results/iso_ppa_pnc6_trim3.rds")
