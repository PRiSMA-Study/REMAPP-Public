#****************************************************************************
#Maternal Depression Likelihood
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_dpr.rda")

#run spline model
spline_dpr <- knot_fun(df_mat_dpr, "hb", "dpr")
saveRDS(spline_dpr, "iso_results/spline_dpr.rds")

#run isotonic model
iso_dpr <- flexstepreg_glmer(df_mat_dpr$dpr, df_mat_dpr$hb, df_mat_dpr$SITE, 
                             covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_dpr, "iso_results/iso_dpr.rds")

# #****************************************************************************
# #*******trim1
# #prepare data
# load("derived_data/df_mat_dpr_trim1.rda")
# 
# #run spline model
# spline_dpr_trim1 <- knot_fun(df_mat_dpr_trim1, "hb", "dpr")
# saveRDS(spline_dpr_trim1, "iso_results/spline_dpr_trim1.rds")
# 
# #run isotonic model
# iso_dpr_trim1 <- flexstepreg_glmer(df_mat_dpr_trim1$dpr, df_mat_dpr_trim1$hb, df_mat_dpr_trim1$SITE, 
#                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
# saveRDS(iso_dpr_trim1, "iso_results/iso_dpr_trim1.rds")
# 
# #****************************************************************************
# #*******trim2
# #prepare data
# load("derived_data/df_mat_dpr_trim2.rda")
# 
# #run spline model
# spline_dpr_trim2 <- knot_fun(df_mat_dpr_trim2, "hb", "dpr")
# saveRDS(spline_dpr_trim2, "iso_results/spline_dpr_trim2.rds")
# 
# #run isotonic model
# iso_dpr_trim2 <- flexstepreg_glmer(df_mat_dpr_trim2$dpr, df_mat_dpr_trim2$hb, df_mat_dpr_trim2$SITE, 
#                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
# saveRDS(iso_dpr_trim2, "iso_results/iso_dpr_trim2.rds")
# 
# #****************************************************************************
# #*******trim3
# #prepare data
# load("derived_data/df_mat_dpr_trim3.rda")
# 
# #run spline model
# spline_dpr_trim3 <- knot_fun(df_mat_dpr_trim3, "hb", "dpr")
# saveRDS(spline_dpr_trim3, "iso_results/spline_dpr_trim3.rds")
# 
# #run isotonic model
# iso_dpr_trim3 <- flexstepreg_glmer(df_mat_dpr_trim3$dpr, df_mat_dpr_trim3$hb, df_mat_dpr_trim3$SITE, 
#                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
# saveRDS(iso_dpr_trim3, "iso_results/iso_dpr_trim3.rds")
