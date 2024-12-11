#****************************************************************************
#prepare data
#****************************************************************************
library(tidyverse)

load("derived_data/df_mat_pph.rda")

#****************************************************************************
#run spline model
#****************************************************************************
source("iso_code/spline_binary.R")

spline_pph <- knot_fun(df_mat_pph, "hb", "HEM_PPH")
saveRDS(spline_pph, "iso_results/spline_pph.rds")

#****************************************************************************
#run isotonic model
#****************************************************************************
source("iso_code/iso_binary.R")

iso_pph <- flexstepreg_glmer(df_mat_pph$HEM_PPH, df_mat_pph$hb, df_mat_pph$SITE, 
                             covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pph, "iso_results/iso_pph.rds")
