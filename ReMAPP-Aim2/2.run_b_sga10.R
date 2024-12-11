#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_sga10.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_sga10 <- knot_fun(df_inf_sga10, "hb", "sga10")
saveRDS(spline_sga10, "iso_results/spline_sga10.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_sga10 <- flexstepreg_glmer(df_inf_sga10$sga10, df_inf_sga10$hb, df_inf_sga10$SITE, 
                                 covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_sga10, "iso_results/iso_sga10.rds")
