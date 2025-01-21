#****************************************************************************
#Preterm premature rupture of membranes
#****************************************************************************
library(tidyverse)
source("iso_code/spline_binary.R")
source("iso_code/iso_binary.R")

#****************************************************************************
#*******all
#prepare data
load("derived_data/df_mat_pprom.rda")

#run spline model
spline_pprom <- knot_fun(df_mat_pprom, "hb", "pprom")
saveRDS(spline_pprom, "iso_results/spline_pprom.rds")

#run isotonic model
iso_pprom <- flexstepreg_glmer(df_mat_pprom$pprom, df_mat_pprom$hb, df_mat_pprom$SITE, 
                             covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pprom, "iso_results/iso_pprom.rds")

#****************************************************************************
#*******trim1
#prepare data
load("derived_data/df_mat_pprom_trim1.rda")

#run spline model
spline_pprom_trim1 <- knot_fun(df_mat_pprom_trim1, "hb", "pprom")
saveRDS(spline_pprom_trim1, "iso_results/spline_pprom_trim1.rds")

#run isotonic model
iso_pprom_trim1 <- flexstepreg_glmer(df_mat_pprom_trim1$pprom, df_mat_pprom_trim1$hb, df_mat_pprom_trim1$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pprom_trim1, "iso_results/iso_pprom_trim1.rds")

#****************************************************************************
#*******trim2
#prepare data
load("derived_data/df_mat_pprom_trim2.rda")

#run spline model
spline_pprom_trim2 <- knot_fun(df_mat_pprom_trim2, "hb", "pprom")
saveRDS(spline_pprom_trim2, "iso_results/spline_pprom_trim2.rds")

#run isotonic model
iso_pprom_trim2 <- flexstepreg_glmer(df_mat_pprom_trim2$pprom, df_mat_pprom_trim2$hb, df_mat_pprom_trim2$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pprom_trim2, "iso_results/iso_pprom_trim2.rds")

#****************************************************************************
#*******trim3
#prepare data
load("derived_data/df_mat_pprom_trim3.rda")

#run spline model
spline_pprom_trim3 <- knot_fun(df_mat_pprom_trim3, "hb", "pprom")
saveRDS(spline_pprom_trim3, "iso_results/spline_pprom_trim3.rds")

#run isotonic model
iso_pprom_trim3 <- flexstepreg_glmer(df_mat_pprom_trim3$pprom, df_mat_pprom_trim3$hb, df_mat_pprom_trim3$SITE, 
                                   covar2=NULL, random_effect=NULL, alpha = 0.01) 
saveRDS(iso_pprom_trim3, "iso_results/iso_pprom_trim3.rds")
