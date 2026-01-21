#****************************************************************************
#Neonatal hyperbilirubinemia
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")


#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_hyperbili.rda")

#run spline model
spline_hyperbili <- knot_fun(df_inf_hyperbili, "hb", "hyperbili")
saveRDS(spline_hyperbili, "iso_results/spline_hyperbili.rds")

#run isotonic model
iso_hyperbili <- flexstepreg_glmer(df_inf_hyperbili$hyperbili, df_inf_hyperbili$hb, df_inf_hyperbili$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_hyperbili, "iso_results/iso_hyperbili.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_hyperbili_trim1.rda")

#run spline model
spline_hyperbili_trim1 <- knot_fun(df_inf_hyperbili_trim1, "hb", "hyperbili")
saveRDS(spline_hyperbili_trim1, "iso_results/spline_hyperbili_trim1.rds")

#run isotonic model
iso_hyperbili_trim1 <- flexstepreg_glmer(df_inf_hyperbili_trim1$hyperbili, df_inf_hyperbili_trim1$hb, df_inf_hyperbili_trim1$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_hyperbili_trim1, "iso_results/iso_hyperbili_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_hyperbili_trim2.rda")

#run spline model
spline_hyperbili_trim2 <- knot_fun(df_inf_hyperbili_trim2, "hb", "hyperbili")
saveRDS(spline_hyperbili_trim2, "iso_results/spline_hyperbili_trim2.rds")

#run isotonic model
iso_hyperbili_trim2 <- flexstepreg_glmer(df_inf_hyperbili_trim2$hyperbili, df_inf_hyperbili_trim2$hb, df_inf_hyperbili_trim2$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_hyperbili_trim2, "iso_results/iso_hyperbili_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_hyperbili_trim3.rda")

#run spline model
spline_hyperbili_trim3 <- knot_fun(df_inf_hyperbili_trim3, "hb", "hyperbili")
saveRDS(spline_hyperbili_trim3, "iso_results/spline_hyperbili_trim3.rds")

#run isotonic model
iso_hyperbili_trim3 <- flexstepreg_glmer(df_inf_hyperbili_trim3$hyperbili, df_inf_hyperbili_trim3$hb, df_inf_hyperbili_trim3$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_hyperbili_trim3, "iso_results/iso_hyperbili_trim3.rds")