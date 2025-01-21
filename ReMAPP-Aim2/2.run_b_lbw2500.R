#****************************************************************************
#Low birth weight (<2500g)
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_lbw2500.rda")

#run spline model
spline_lbw2500 <- knot_fun(df_inf_lbw2500, "hb", "lbw2500")
saveRDS(spline_lbw2500, "iso_results/spline_lbw2500.rds")

#run isotonic model
iso_lbw2500 <- flexstepreg_glmer(df_inf_lbw2500$lbw2500, df_inf_lbw2500$hb, df_inf_lbw2500$SITE, 
                                 covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw2500, "iso_results/iso_lbw2500.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_lbw2500_trim1.rda")

#run spline model
spline_lbw2500_trim1 <- knot_fun(df_inf_lbw2500_trim1, "hb", "lbw2500")
saveRDS(spline_lbw2500_trim1, "iso_results/spline_lbw2500_trim1.rds")

#run isotonic model
iso_lbw2500_trim1 <- flexstepreg_glmer(df_inf_lbw2500_trim1$lbw2500, df_inf_lbw2500_trim1$hb, df_inf_lbw2500_trim1$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw2500_trim1, "iso_results/iso_lbw2500_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_lbw2500_trim2.rda")

#run spline model
spline_lbw2500_trim2 <- knot_fun(df_inf_lbw2500_trim2, "hb", "lbw2500")
saveRDS(spline_lbw2500_trim2, "iso_results/spline_lbw2500_trim2.rds")

#run isotonic model
iso_lbw2500_trim2 <- flexstepreg_glmer(df_inf_lbw2500_trim2$lbw2500, df_inf_lbw2500_trim2$hb, df_inf_lbw2500_trim2$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw2500_trim2, "iso_results/iso_lbw2500_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_lbw2500_trim3.rda")

#run spline model
spline_lbw2500_trim3 <- knot_fun(df_inf_lbw2500_trim3, "hb", "lbw2500")
saveRDS(spline_lbw2500_trim3, "iso_results/spline_lbw2500_trim3.rds")

#run isotonic model
iso_lbw2500_trim3 <- flexstepreg_glmer(df_inf_lbw2500_trim3$lbw2500, df_inf_lbw2500_trim3$hb, df_inf_lbw2500_trim3$SITE, 
                                       covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw2500_trim3, "iso_results/iso_lbw2500_trim3.rds")
