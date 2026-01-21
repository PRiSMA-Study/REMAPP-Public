#****************************************************************************
#Low birth weight (<1500g)
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")


#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_lbw1500.rda")

#run spline model
spline_lbw1500 <- knot_fun(df_inf_lbw1500, "hb", "lbw1500")
saveRDS(spline_lbw1500, "iso_results/spline_lbw1500.rds")

#run isotonic model
iso_lbw1500 <- flexstepreg_glmer(df_inf_lbw1500$lbw1500, df_inf_lbw1500$hb, df_inf_lbw1500$SITE, 
                                 covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw1500, "iso_results/iso_lbw1500.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_lbw1500_trim1.rda")

#run spline model
spline_lbw1500_trim1 <- knot_fun(df_inf_lbw1500_trim1, "hb", "lbw1500")
saveRDS(spline_lbw1500_trim1, "iso_results/spline_lbw1500_trim1.rds")

#run isotonic model
iso_lbw1500_trim1 <- flexstepreg_glmer(df_inf_lbw1500_trim1$lbw1500, df_inf_lbw1500_trim1$hb, df_inf_lbw1500_trim1$SITE, 
                                 covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw1500_trim1, "iso_results/iso_lbw1500_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_lbw1500_trim2.rda")

#run spline model
spline_lbw1500_trim2 <- knot_fun(df_inf_lbw1500_trim2, "hb", "lbw1500")
saveRDS(spline_lbw1500_trim2, "iso_results/spline_lbw1500_trim2.rds")

#run isotonic model
iso_lbw1500_trim2 <- flexstepreg_glmer(df_inf_lbw1500_trim2$lbw1500, df_inf_lbw1500_trim2$hb, df_inf_lbw1500_trim2$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw1500_trim2, "iso_results/iso_lbw1500_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_lbw1500_trim3.rda")

#run spline model
spline_lbw1500_trim3 <- knot_fun(df_inf_lbw1500_trim3, "hb", "lbw1500")
saveRDS(spline_lbw1500_trim3, "iso_results/spline_lbw1500_trim3.rds")

#run isotonic model
iso_lbw1500_trim3 <- flexstepreg_glmer(df_inf_lbw1500_trim3$lbw1500, df_inf_lbw1500_trim3$hb, df_inf_lbw1500_trim3$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw1500_trim3, "iso_results/iso_lbw1500_trim3.rds")
