#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_asph.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_asph <- knot_fun(df_inf_asph, "hb", "inf_asph")
saveRDS(spline_asph, "iso_results/spline_asph.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_asph <- flexstepreg_glmer(df_inf_asph$inf_asph, df_inf_asph$hb, df_inf_asph$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_asph, "iso_results/iso_asph.rds")
