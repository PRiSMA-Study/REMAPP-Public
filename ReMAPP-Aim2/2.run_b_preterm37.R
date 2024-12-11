#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_preterm37.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_preterm37 <- knot_fun(df_inf_preterm37, "hb", "preterm37")
saveRDS(spline_preterm37, "iso_results/spline_preterm37.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_preterm37 <- flexstepreg_glmer(df_inf_preterm37$preterm37, df_inf_preterm37$hb, df_inf_preterm37$SITE, 
                               covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm37, "iso_results/iso_preterm37.rds")
