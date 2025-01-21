#****************************************************************************
#Composite outcome
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_compo.rda")

#run spline model
spline_compo <- knot_fun(df_inf_compo, "hb", "compo_pre_lbw_sga")
saveRDS(spline_compo, "iso_results/spline_compo.rds")

#run isotonic model
iso_compo <- flexstepreg_glmer(df_inf_compo$compo_pre_lbw_sga, df_inf_compo$hb, df_inf_compo$SITE, 
                               covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_compo, "iso_results/iso_compo.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_compo_trim1.rda")

#run spline model
spline_compo_trim1 <- knot_fun(df_inf_compo_trim1, "hb", "compo_pre_lbw_sga")
saveRDS(spline_compo_trim1, "iso_results/spline_compo_trim1.rds")

#run isotonic model
iso_compo_trim1 <- flexstepreg_glmer(df_inf_compo_trim1$compo_pre_lbw_sga, df_inf_compo_trim1$hb, df_inf_compo_trim1$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_compo_trim1, "iso_results/iso_compo_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_compo_trim2.rda")

#run spline model
spline_compo_trim2 <- knot_fun(df_inf_compo_trim2, "hb", "compo_pre_lbw_sga")
saveRDS(spline_compo_trim2, "iso_results/spline_compo_trim2.rds")

#run isotonic model
iso_compo_trim2 <- flexstepreg_glmer(df_inf_compo_trim2$compo_pre_lbw_sga, df_inf_compo_trim2$hb, df_inf_compo_trim2$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_compo_trim2, "iso_results/iso_compo_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_compo_trim3.rda")

#run spline model
spline_compo_trim3 <- knot_fun(df_inf_compo_trim3, "hb", "compo_pre_lbw_sga")
saveRDS(spline_compo_trim3, "iso_results/spline_compo_trim3.rds")

#run isotonic model
iso_compo_trim3 <- flexstepreg_glmer(df_inf_compo_trim3$compo_pre_lbw_sga, df_inf_compo_trim3$hb, df_inf_compo_trim3$SITE, 
                                     covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_compo_trim3, "iso_results/iso_compo_trim3.rds")