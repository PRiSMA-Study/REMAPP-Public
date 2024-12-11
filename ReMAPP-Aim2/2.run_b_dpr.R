#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_mat_dpr.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_dpr <- knot_fun(df_mat_dpr, "hb", "dpr")
saveRDS(spline_dpr, "iso_results/spline_dpr.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_dpr <- flexstepreg_glmer(df_mat_dpr$dpr, df_mat_dpr$hb, df_mat_dpr$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_dpr, "iso_results/iso_dpr.rds")
