#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_mat_ppa_pnc26.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_ppa_pnc26 <- knot_fun(df_mat_ppa_pnc26, "hb", "ppa_pnc26")
saveRDS(spline_ppa_pnc26, "iso_results/spline_ppa_pnc26.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_ppa_pnc26 <- flexstepreg_glmer(df_mat_ppa_pnc26$ppa_pnc26, df_mat_ppa_pnc26$hb, df_mat_ppa_pnc26$SITE, 
                                  covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_ppa_pnc26, "iso_results/iso_ppa_pnc26.rds")
