#****************************************************************************
#EPDS score
#****************************************************************************
library(tidyverse)
source("iso_code/spline_continuous.R")
source("iso_code/iso_continuous.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_dpr.rda")

#run spline model
spline_dpr_score <- knot_fun_continue(df_mat_dpr, "hb", "dpr_score")
saveRDS(spline_dpr_score, "iso_results/spline_dpr_score.rds")

#run isotonic model
iso_dpr_score <- flexstepreg_lmer(df_mat_dpr$dpr_score, df_mat_dpr$hb, df_mat_dpr$SITE, 
                                  covar2=df_mat_dpr$visit, random_effect=df_mat_dpr$MOMID, 
                                  alpha = 0.01) 
saveRDS(iso_dpr_score, "iso_results/iso_dpr_score.rds")
