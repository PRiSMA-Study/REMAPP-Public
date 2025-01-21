#****************************************************************************
#Postpartum hemorrahage 
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_pph.rda")

#run spline model
spline_pph <- knot_fun(df_mat_pph, "hb", "HEM_PPH")
saveRDS(spline_pph, "iso_results/spline_pph.rds")

#run isotonic model
iso_pph <- flexstepreg_glmer(df_mat_pph$HEM_PPH, df_mat_pph$hb, df_mat_pph$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pph, "iso_results/iso_pph.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_mat_pph_trim1.rda")

#run spline model
spline_pph_trim1 <- knot_fun(df_mat_pph_trim1, "hb", "HEM_PPH")
saveRDS(spline_pph_trim1, "iso_results/spline_pph_trim1.rds")

#run isotonic model
iso_pph_trim1 <- flexstepreg_glmer(df_mat_pph_trim1$HEM_PPH, df_mat_pph_trim1$hb, df_mat_pph_trim1$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pph_trim1, "iso_results/iso_pph_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_mat_pph_trim2.rda")

#run spline model
spline_pph_trim2 <- knot_fun(df_mat_pph_trim2, "hb", "HEM_PPH")
saveRDS(spline_pph_trim2, "iso_results/spline_pph_trim2.rds")

#run isotonic model
iso_pph_trim2 <- flexstepreg_glmer(df_mat_pph_trim2$HEM_PPH, df_mat_pph_trim2$hb, df_mat_pph_trim2$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pph_trim2, "iso_results/iso_pph_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_mat_pph_trim3.rda")

#run spline model
spline_pph_trim3 <- knot_fun(df_mat_pph_trim3, "hb", "HEM_PPH")
saveRDS(spline_pph_trim3, "iso_results/spline_pph_trim3.rds")

#run isotonic model
iso_pph_trim3 <- flexstepreg_glmer(df_mat_pph_trim3$HEM_PPH, df_mat_pph_trim3$hb, df_mat_pph_trim3$SITE, 
                                         covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pph_trim3, "iso_results/iso_pph_trim3.rds")
