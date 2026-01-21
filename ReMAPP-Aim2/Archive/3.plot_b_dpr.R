#****************************************************************************
#Maternal Depression Likelihood
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")


load("derived_data/df_mat_dpr.rda")
load("derived_data/df_mat_dpr_trim1.rda")
load("derived_data/df_mat_dpr_trim2.rda")
load("derived_data/df_mat_dpr_trim3.rda")


#****************************************************************************
#*******all
#read in model
spline_dpr <- readRDS("iso_results/spline_dpr.rds")   
iso_dpr <- readRDS("iso_results/iso_dpr.rds")

#run and save plot
isoplot_dpr <- iso_fun_glmer(df_mat_dpr$hb, df_mat_dpr$dpr, "Depressed", spline_dpr, iso_dpr)

png(file = "iso_results/isoplot_dpr.png")
isoplot_dpr <- iso_fun_glmer(df_mat_dpr$hb, df_mat_dpr$dpr, "Depressed", spline_dpr, iso_dpr)
dev.off()

#run and save output data
out_dpr <- outdata(iso_dpr, "Depressed")
out_dpr
save(out_dpr, file = "iso_results/out_dpr.rda")

#****************************************************************************
#*******trim1
#read in model
spline_dpr_trim1 <- readRDS("iso_results/spline_dpr_trim1.rds")
iso_dpr_trim1 <- readRDS("iso_results/iso_dpr_trim1.rds")

#run and save plot
isoplot_dpr_trim1 <- iso_fun_glmer(df_mat_dpr_trim1$hb, df_mat_dpr_trim1$dpr, "Likelihood of depression - trim1", spline_dpr_trim1, iso_dpr_trim1)

png(file = "iso_results/isoplot_dpr_trim1.png")
isoplot_dpr_trim1 <- iso_fun_glmer(df_mat_dpr_trim1$hb, df_mat_dpr_trim1$dpr, "Likelihood of depression - trim1", spline_dpr_trim1, iso_dpr_trim1)
dev.off()

#run and save output data
out_dpr_trim1 <- outdata(iso_dpr_trim1, "Depress_Likel_trim1")
out_dpr_trim1
save(out_dpr_trim1, file = "iso_results/out_dpr_trim1.rda")

#****************************************************************************
#*******trim2
#read in model
spline_dpr_trim2 <- readRDS("iso_results/spline_dpr_trim2.rds")
iso_dpr_trim2 <- readRDS("iso_results/iso_dpr_trim2.rds")

#run and save plot
isoplot_dpr_trim2 <- iso_fun_glmer(df_mat_dpr_trim2$hb, df_mat_dpr_trim2$dpr, "Likelihood of depression - trim2", spline_dpr_trim2, iso_dpr_trim2)

png(file = "iso_results/isoplot_dpr_trim2.png")
isoplot_dpr_trim2 <- iso_fun_glmer(df_mat_dpr_trim2$hb, df_mat_dpr_trim2$dpr, "Likelihood of depression - trim2", spline_dpr_trim2, iso_dpr_trim2)
dev.off()

#run and save output data
out_dpr_trim2 <- outdata(iso_dpr_trim2, "Depress_Likel_trim2")
out_dpr_trim2
save(out_dpr_trim2, file = "iso_results/out_dpr_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_dpr_trim3 <- readRDS("iso_results/spline_dpr_trim3.rds")
iso_dpr_trim3 <- readRDS("iso_results/iso_dpr_trim3.rds")

#run and save plot
isoplot_dpr_trim3 <- iso_fun_glmer(df_mat_dpr_trim3$hb, df_mat_dpr_trim3$dpr, "Likelihood of depression - trim3", spline_dpr_trim3, iso_dpr_trim3)

png(file = "iso_results/isoplot_dpr_trim3.png")
isoplot_dpr_trim3 <- iso_fun_glmer(df_mat_dpr_trim3$hb, df_mat_dpr_trim3$dpr, "Likelihood of depression - trim3", spline_dpr_trim3, iso_dpr_trim3)
dev.off()

#run and save output data
out_dpr_trim3 <- outdata(iso_dpr_trim3, "Depress_Likel_trim3")
out_dpr_trim3
save(out_dpr_trim3, file = "iso_results/out_dpr_trim3.rda")