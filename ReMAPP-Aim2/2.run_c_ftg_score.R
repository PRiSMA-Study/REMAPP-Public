#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_mat_ftg.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_continuous.R")

spline_ftg_score <- knot_fun_continue(df_mat_ftg, "hb", "ftg_score")
saveRDS(spline_ftg_score, "iso_results/spline_ftg_score.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
# iso_ftg_score <- flexstepreg_lmer(y, x, covar1, covar2, random_effect, alpha.adjacency = alpha) 
source("iso_code/iso_continuous.R")

iso_ftg_score <- flexstepreg_lmer(df_mat_ftg$ftg_score, df_mat_ftg$hb, df_mat_ftg$SITE, 
                                  covar2=df_mat_ftg$visit, random_effect=df_mat_ftg$MOMID, 
                                  alpha = 0.01) 
saveRDS(iso_ftg_score, "iso_results/iso_ftg_score.rds")
