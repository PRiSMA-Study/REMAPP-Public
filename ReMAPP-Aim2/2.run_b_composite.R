#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_inf_compo.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_compo <- knot_fun(df_inf_compo, "hb", "compo_pre_lbw_sga")
saveRDS(spline_compo, "iso_results/spline_compo.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
# iso_compo <- flexstepreg_lmer(y, x, covar1, covar2, random_effect, alpha.adjacency = alpha) 
source("iso_code/iso_binary.R")

iso_compo <- flexstepreg_glmer(df_inf_compo$compo_pre_lbw_sga, df_inf_compo$hb, df_inf_compo$SITE, 
                                  covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_compo, "iso_results/iso_compo.rds")
