#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_hyperbili.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_hyperbili <- knot_fun(df_inf_hyperbili, "hb", "hyperbili")
saveRDS(spline_hyperbili, "iso_results/spline_hyperbili.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_hyperbili <- flexstepreg_glmer(df_inf_hyperbili$hyperbili, df_inf_hyperbili$hb, df_inf_hyperbili$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_hyperbili, "iso_results/iso_hyperbili.rds")
