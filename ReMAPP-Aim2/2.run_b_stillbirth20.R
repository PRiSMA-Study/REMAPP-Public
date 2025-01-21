#****************************************************************************
#Stillbirth >= 20 weeks
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_stillbirth20.rda")

#run spline model
spline_stillbirth20 <- knot_fun(df_inf_stillbirth20, "hb", "inf_stillbirth20")
saveRDS(spline_stillbirth20, "iso_results/spline_stillbirth20.rds")

#run isotonic model
iso_stillbirth20 <- flexstepreg_glmer(df_inf_stillbirth20$inf_stillbirth20, df_inf_stillbirth20$hb, df_inf_stillbirth20$SITE, 
                               covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_stillbirth20, "iso_results/iso_stillbirth20.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_stillbirth20_trim1.rda")

#run spline model
spline_stillbirth20_trim1 <- knot_fun(df_inf_stillbirth20_trim1, "hb", "inf_stillbirth20")
saveRDS(spline_stillbirth20_trim1, "iso_results/spline_stillbirth20_trim1.rds")

#run isotonic model
iso_stillbirth20_trim1 <- flexstepreg_glmer(df_inf_stillbirth20_trim1$inf_stillbirth20, df_inf_stillbirth20_trim1$hb, df_inf_stillbirth20_trim1$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_stillbirth20_trim1, "iso_results/iso_stillbirth20_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_stillbirth20_trim2.rda")

#run spline model
spline_stillbirth20_trim2 <- knot_fun(df_inf_stillbirth20_trim2, "hb", "inf_stillbirth20")
saveRDS(spline_stillbirth20_trim2, "iso_results/spline_stillbirth20_trim2.rds")

#run isotonic model
iso_stillbirth20_trim2 <- flexstepreg_glmer(df_inf_stillbirth20_trim2$inf_stillbirth20, df_inf_stillbirth20_trim2$hb, df_inf_stillbirth20_trim2$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_stillbirth20_trim2, "iso_results/iso_stillbirth20_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_stillbirth20_trim3.rda")

#run spline model
spline_stillbirth20_trim3 <- knot_fun(df_inf_stillbirth20_trim3, "hb", "inf_stillbirth20")
saveRDS(spline_stillbirth20_trim3, "iso_results/spline_stillbirth20_trim3.rds")

#run isotonic model
iso_stillbirth20_trim3 <- flexstepreg_glmer(df_inf_stillbirth20_trim3$inf_stillbirth20, df_inf_stillbirth20_trim3$hb, df_inf_stillbirth20_trim3$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_stillbirth20_trim3, "iso_results/iso_stillbirth20_trim3.rds")


