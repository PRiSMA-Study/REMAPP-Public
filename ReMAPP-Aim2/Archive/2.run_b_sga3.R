#****************************************************************************
#Small for gestational age (<3th)
#****************************************************************************
library(tidyverse)
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/spline_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/iso_binary.R")


#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_sga3.rda")

#run spline model
spline_sga3 <- knot_fun(df_inf_sga3, "hb", "sga3")
saveRDS(spline_sga3, "iso_results/spline_sga3.rds")

#run isotonic model
iso_sga3 <- flexstepreg_glmer(df_inf_sga3$sga3, df_inf_sga3$hb, df_inf_sga3$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_sga3, "iso_results/iso_sga3.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_sga3_trim1.rda")

#run spline model
spline_sga3_trim1 <- knot_fun(df_inf_sga3_trim1, "hb", "sga3")
saveRDS(spline_sga3_trim1, "iso_results/spline_sga3_trim1.rds")

#run isotonic model
iso_sga3_trim1 <- flexstepreg_glmer(df_inf_sga3_trim1$sga3, df_inf_sga3_trim1$hb, df_inf_sga3_trim1$SITE, 
                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_sga3_trim1, "iso_results/iso_sga3_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_sga3_trim2.rda")

#run spline model
spline_sga3_trim2 <- knot_fun(df_inf_sga3_trim2, "hb", "sga3")
saveRDS(spline_sga3_trim2, "iso_results/spline_sga3_trim2.rds")

#run isotonic model
iso_sga3_trim2 <- flexstepreg_glmer(df_inf_sga3_trim2$sga3, df_inf_sga3_trim2$hb, df_inf_sga3_trim2$SITE, 
                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_sga3_trim2, "iso_results/iso_sga3_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_sga3_trim3.rda")

#run spline model
spline_sga3_trim3 <- knot_fun(df_inf_sga3_trim3, "hb", "sga3")
saveRDS(spline_sga3_trim3, "iso_results/spline_sga3_trim3.rds")

#run isotonic model
iso_sga3_trim3 <- flexstepreg_glmer(df_inf_sga3_trim3$sga3, df_inf_sga3_trim3$hb, df_inf_sga3_trim3$SITE, 
                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_sga3_trim3, "iso_results/iso_sga3_trim3.rds")


