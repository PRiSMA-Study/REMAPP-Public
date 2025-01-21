#****************************************************************************
#Preterm <34 weeks
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_preterm34.rda")

#run spline model
spline_preterm34 <- knot_fun(df_inf_preterm34, "hb", "preterm34")
saveRDS(spline_preterm34, "iso_results/spline_preterm34.rds")

#run isotonic model
iso_preterm34 <- flexstepreg_glmer(df_inf_preterm34$preterm34, df_inf_preterm34$hb, df_inf_preterm34$SITE, 
                               covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm34, "iso_results/iso_preterm34.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_preterm34_trim1.rda")

#run spline model
spline_preterm34_trim1 <- knot_fun(df_inf_preterm34_trim1, "hb", "preterm34")
saveRDS(spline_preterm34_trim1, "iso_results/spline_preterm34_trim1.rds")

#run isotonic model
iso_preterm34_trim1 <- flexstepreg_glmer(df_inf_preterm34_trim1$preterm34, df_inf_preterm34_trim1$hb, df_inf_preterm34_trim1$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm34_trim1, "iso_results/iso_preterm34_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_preterm34_trim2.rda")

#run spline model
spline_preterm34_trim2 <- knot_fun(df_inf_preterm34_trim2, "hb", "preterm34")
saveRDS(spline_preterm34_trim2, "iso_results/spline_preterm34_trim2.rds")

#run isotonic model
iso_preterm34_trim2 <- flexstepreg_glmer(df_inf_preterm34_trim2$preterm34, df_inf_preterm34_trim2$hb, df_inf_preterm34_trim2$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm34_trim2, "iso_results/iso_preterm34_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_preterm34_trim3.rda")

#run spline model
spline_preterm34_trim3 <- knot_fun(df_inf_preterm34_trim3, "hb", "preterm34")
saveRDS(spline_preterm34_trim3, "iso_results/spline_preterm34_trim3.rds")

#run isotonic model
iso_preterm34_trim3 <- flexstepreg_glmer(df_inf_preterm34_trim3$preterm34, df_inf_preterm34_trim3$hb, df_inf_preterm34_trim3$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_preterm34_trim3, "iso_results/iso_preterm34_trim3.rds")

