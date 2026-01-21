#****************************************************************************
#read in model
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

load("derived_data/df_mat_preclamp.rda")
load("derived_data/df_mat_preclamp_trim1.rda")
load("derived_data/df_mat_preclamp_trim2.rda")
load("derived_data/df_mat_preclamp_trim3.rda")

#****************************************************************************
#Preterm premature rupture of membranes
#****************************************************************************

#*******all
#read in model
spline_preclamp <- readRDS("iso_results/spline_preclamp.rds")   
iso_preclamp <- readRDS("iso_results/iso_preclamp.rds")

#run and save plot
isoplot_preclamp <- iso_fun_glmer(df_mat_preclamp$hb, df_mat_preclamp$preclamp, "PREECLAMPSIA", spline_preclamp, iso_preclamp)

png(file = "iso_results/isoplot_preclamp.png")
isoplot_preclamp <- iso_fun_glmer(df_mat_preclamp$hb, df_mat_preclamp$preclamp, "PREECLAMPSIA", spline_preclamp, iso_preclamp)
dev.off()

#run and save output data
out_preclamp  <- outdata(iso_preclamp, "PREECLAMPSIA")
out_preclamp 
save(out_preclamp, file = "iso_results/out_preclamp.rda")

#****************************************************************************
#*******trim1
#read in model
spline_preclamp_trim1 <- readRDS("iso_results/spline_preclamp_trim1.rds")   
iso_preclamp_trim1 <- readRDS("iso_results/iso_preclamp_trim1.rds")

#run and save plot
isoplot_preclamp_trim1 <- iso_fun_glmer(df_mat_preclamp_trim1$hb, df_mat_preclamp_trim1$preclamp, "PREECLAMPSIA - trim1", spline_preclamp_trim1, iso_preclamp_trim1)

png(file = "iso_results/isoplot_preclamp_trim1.png")
isoplot_preclamp_trim1 <- iso_fun_glmer(df_mat_preclamp_trim1$hb, df_mat_preclamp_trim1$preclamp, "PREECLAMPSIA - trim1", spline_preclamp_trim1, iso_preclamp_trim1)
dev.off()

#run and save output data
out_preclamp_trim1  <- outdata(iso_preclamp_trim1, "PREECLAMPSIA")
out_preclamp_trim1 
save(out_preclamp_trim1, file = "iso_results/out_preclamp_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_preclamp_trim2 <- readRDS("iso_results/spline_preclamp_trim2.rds")   
iso_preclamp_trim2 <- readRDS("iso_results/iso_preclamp_trim2.rds")

#run and save plot
isoplot_preclamp_trim2 <- iso_fun_glmer(df_mat_preclamp_trim2$hb, df_mat_preclamp_trim2$preclamp, "PREECLAMPSIA - trim2", spline_preclamp_trim2, iso_preclamp_trim2)

png(file = "iso_results/isoplot_preclamp_trim2.png")
isoplot_preclamp_trim2 <- iso_fun_glmer(df_mat_preclamp_trim2$hb, df_mat_preclamp_trim2$preclamp, "PREECLAMPSIA - trim2", spline_preclamp_trim2, iso_preclamp_trim2)
dev.off()

#run and save output data
out_preclamp_trim2  <- outdata(iso_preclamp_trim2, "PREECLAMPSIA")
out_preclamp_trim2 
save(out_preclamp_trim2, file = "iso_results/out_preclamp_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_preclamp_trim3 <- readRDS("iso_results/spline_preclamp_trim3.rds")   
iso_preclamp_trim3 <- readRDS("iso_results/iso_preclamp_trim3.rds")

#run and save plot
isoplot_preclamp_trim3 <- iso_fun_glmer(df_mat_preclamp_trim3$hb, df_mat_preclamp_trim3$preclamp, "PREECLAMPSIA - trim3", spline_preclamp_trim3, iso_preclamp_trim3)

png(file = "iso_results/isoplot_preclamp_trim3.png")
isoplot_preclamp_trim3 <- iso_fun_glmer(df_mat_preclamp_trim3$hb, df_mat_preclamp_trim3$preclamp, "PREECLAMPSIA - trim3", spline_preclamp_trim3, iso_preclamp_trim3)
dev.off()

#run and save output data
out_preclamp_trim3  <- outdata(iso_preclamp_trim3, "PREECLAMPSIA")
out_preclamp_trim3 
save(out_preclamp_trim3, file = "iso_results/out_preclamp_trim3.rda")

