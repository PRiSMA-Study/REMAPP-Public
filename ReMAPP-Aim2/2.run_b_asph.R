#****************************************************************************
#Birth asphyxia
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_inf_asph.rda")

#run spline model
spline_asph <- knot_fun(df_inf_asph, "hb", "inf_asph")
saveRDS(spline_asph, "iso_results/spline_asph.rds")

#run isotonic model
iso_asph <- flexstepreg_glmer(df_inf_asph$inf_asph, df_inf_asph$hb, df_inf_asph$SITE, 
                              covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_asph, "iso_results/iso_asph.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_inf_asph_trim1.rda")

#run spline model
spline_asph_trim1 <- knot_fun(df_inf_asph_trim1, "hb", "inf_asph")
saveRDS(spline_asph_trim1, "iso_results/spline_asph_trim1.rds")

#run isotonic model
iso_asph_trim1 <- flexstepreg_glmer(df_inf_asph_trim1$inf_asph, df_inf_asph_trim1$hb, df_inf_asph_trim1$SITE, 
                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_asph_trim1, "iso_results/iso_asph_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_inf_asph_trim2.rda")

#run spline model
spline_asph_trim2 <- knot_fun(df_inf_asph_trim2, "hb", "inf_asph")
saveRDS(spline_asph_trim2, "iso_results/spline_asph_trim2.rds")

#run isotonic model
iso_asph_trim2 <- flexstepreg_glmer(df_inf_asph_trim2$inf_asph, df_inf_asph_trim2$hb, df_inf_asph_trim2$SITE, 
                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_asph_trim2, "iso_results/iso_asph_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_inf_asph_trim3.rda")

#run spline model
spline_asph_trim3 <- knot_fun(df_inf_asph_trim3, "hb", "inf_asph")
saveRDS(spline_asph_trim3, "iso_results/spline_asph_trim3.rds")

#run isotonic model
iso_asph_trim3 <- flexstepreg_glmer(df_inf_asph_trim3$inf_asph, df_inf_asph_trim3$hb, df_inf_asph_trim3$SITE, 
                                    covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_asph_trim3, "iso_results/iso_asph_trim3.rds")