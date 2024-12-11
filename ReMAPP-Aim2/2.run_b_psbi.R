#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_psbi.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_psbi <- knot_fun(df_inf_psbi, "hb", "inf_psbi")
saveRDS(spline_psbi, "iso_results/spline_psbi.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_psbi <- flexstepreg_glmer(df_inf_psbi$inf_psbi, df_inf_psbi$hb, df_inf_psbi$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_psbi, "iso_results/iso_psbi.rds")
