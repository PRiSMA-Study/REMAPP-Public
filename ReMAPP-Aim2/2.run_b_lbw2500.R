#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_lbw2500.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_lbw2500 <- knot_fun(df_inf_lbw2500, "hb", "lbw2500")
saveRDS(spline_lbw2500, "iso_results/spline_lbw2500.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_lbw2500 <- flexstepreg_glmer(df_inf_lbw2500$lbw2500, df_inf_lbw2500$hb, df_inf_lbw2500$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw2500, "iso_results/iso_lbw2500.rds")
