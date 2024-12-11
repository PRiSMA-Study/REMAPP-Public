#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_preterm34.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_preterm34 <- knot_fun(df_inf_preterm34, "hb", "preterm34")
saveRDS(spline_preterm34, "iso_results/spline_preterm34.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_preterm34 <- flexstepreg_glmer(df_inf_preterm34$preterm34, df_inf_preterm34$hb, df_inf_preterm34$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm34, "iso_results/iso_preterm34.rds")
