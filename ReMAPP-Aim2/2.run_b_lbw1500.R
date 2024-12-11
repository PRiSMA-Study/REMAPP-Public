#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_lbw1500.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_lbw1500 <- knot_fun(df_inf_lbw1500, "hb", "lbw1500")
saveRDS(spline_lbw1500, "iso_results/spline_lbw1500.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_lbw1500 <- flexstepreg_glmer(df_inf_lbw1500$lbw1500, df_inf_lbw1500$hb, df_inf_lbw1500$SITE, 
                                 covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_lbw1500, "iso_results/iso_lbw1500.rds")
