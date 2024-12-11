#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_mat_pprom.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_pprom <- knot_fun(df_mat_pprom, "hb", "pprom")
saveRDS(spline_pprom, "iso_results/spline_pprom.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_pprom <- flexstepreg_glmer(df_mat_pprom$pprom, df_mat_pprom$hb, df_mat_pprom$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pprom, "iso_results/iso_pprom.rds")
