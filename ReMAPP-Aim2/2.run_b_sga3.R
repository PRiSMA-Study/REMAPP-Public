#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_sga3.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_sga3 <- knot_fun(df_inf_sga3, "hb", "sga3")
saveRDS(spline_sga3, "iso_results/spline_sga3.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_sga3 <- flexstepreg_glmer(df_inf_sga3$sga3, df_inf_sga3$hb, df_inf_sga3$SITE, 
                               covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_sga3, "iso_results/iso_sga3.rds")
