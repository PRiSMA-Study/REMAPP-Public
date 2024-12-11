#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_mat_ppa_pnc6.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_ppa_pnc6 <- knot_fun(df_mat_ppa_pnc6, "hb", "ppa_pnc6")
saveRDS(spline_ppa_pnc6, "iso_results/spline_ppa_pnc6.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_ppa_pnc6 <- flexstepreg_glmer(df_mat_ppa_pnc6$ppa_pnc6, df_mat_ppa_pnc6$hb, df_mat_ppa_pnc6$SITE, 
                             covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc6, "iso_results/iso_ppa_pnc6.rds")
