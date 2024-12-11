#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_stillbirth20.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_stillbirth20 <- knot_fun(df_inf_stillbirth20, "hb", "inf_stillbirth20")
saveRDS(spline_stillbirth20, "iso_results/spline_stillbirth20.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_stillbirth20 <- flexstepreg_glmer(df_inf_stillbirth20$inf_stillbirth20, df_inf_stillbirth20$hb, df_inf_stillbirth20$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_stillbirth20, "iso_results/iso_stillbirth20.rds")
