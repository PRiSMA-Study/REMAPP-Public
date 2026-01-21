#****************************************************************************
#Maternal postpartum anemia at PNC6
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

#****************************************************************************
#*******all
#*
load("derived_data/df_mat_ppa_pnc6.rda")
load("derived_data/df_mat_ppa_pnc6_trim1.rda")
load("derived_data/df_mat_ppa_pnc6_trim2.rda")
load("derived_data/df_mat_ppa_pnc6_trim3.rda")

#read in model
spline_ppa_pnc6 <- readRDS("iso_results/spline_ppa_pnc6.rds")   
iso_ppa_pnc6 <- readRDS("iso_results/iso_ppa_pnc6.rds")

#run and save plot
isoplot_ppa_pnc6 <- iso_fun_glmer(df_mat_ppa_pnc6$hb, df_mat_ppa_pnc6$ppa_pnc6, "Postpartum anemia at PNC6", spline_ppa_pnc6, iso_ppa_pnc6)

png(file = "iso_results/isoplot_ppa_pnc6.png")
isoplot_ppa_pnc6 <- iso_fun_glmer(df_mat_ppa_pnc6$hb, df_mat_ppa_pnc6$ppa_pnc6, "Postpartum anemia at PNC6", spline_ppa_pnc6, iso_ppa_pnc6)
dev.off()

#run and save output data
out_ppa_pnc6  <- outdata(iso_ppa_pnc6, "PPA-PNC6")
out_ppa_pnc6 
save(out_ppa_pnc6, file = "iso_results/out_ppa_pnc6.rda")

#****************************************************************************
#*******trim1
#read in model
spline_ppa_pnc6_trim1 <- readRDS("iso_results/spline_ppa_pnc6_trim1.rds")   
iso_ppa_pnc6_trim1 <- readRDS("iso_results/iso_ppa_pnc6_trim1.rds")

#run and save plot
isoplot_ppa_pnc6_trim1 <- iso_fun_glmer(df_mat_ppa_pnc6_trim1$hb, df_mat_ppa_pnc6_trim1$ppa_pnc6, "PPA at PNC6 - trim1", spline_ppa_pnc6_trim1, iso_ppa_pnc6_trim1)

png(file = "iso_results/isoplot_ppa_pnc6_trim1.png")
isoplot_ppa_pnc6_trim1 <- iso_fun_glmer(df_mat_ppa_pnc6_trim1$hb, df_mat_ppa_pnc6_trim1$ppa_pnc6, "PPA at PNC6 - trim1", spline_ppa_pnc6_trim1, iso_ppa_pnc6_trim1)
dev.off()

#run and save output data
out_ppa_pnc6_trim1  <- outdata(iso_ppa_pnc6_trim1, "PPA-PNC6")
out_ppa_pnc6_trim1 
save(out_ppa_pnc6_trim1, file = "iso_results/out_ppa_pnc6_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_ppa_pnc6_trim2 <- readRDS("iso_results/spline_ppa_pnc6_trim2.rds")   
iso_ppa_pnc6_trim2 <- readRDS("iso_results/iso_ppa_pnc6_trim2.rds")

#run and save plot
isoplot_ppa_pnc6_trim2 <- iso_fun_glmer(df_mat_ppa_pnc6_trim2$hb, df_mat_ppa_pnc6_trim2$ppa_pnc6, "PPA at PNC6 - trim2", spline_ppa_pnc6_trim2, iso_ppa_pnc6_trim2)

png(file = "iso_results/isoplot_ppa_pnc6_trim2.png")
isoplot_ppa_pnc6_trim2 <- iso_fun_glmer(df_mat_ppa_pnc6_trim2$hb, df_mat_ppa_pnc6_trim2$ppa_pnc6, "PPA at PNC6 - trim2", spline_ppa_pnc6_trim2, iso_ppa_pnc6_trim2)
dev.off()

#run and save output data
out_ppa_pnc6_trim2  <- outdata(iso_ppa_pnc6_trim2, "PPA-PNC6")
out_ppa_pnc6_trim2 
save(out_ppa_pnc6_trim2, file = "iso_results/out_ppa_pnc6_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_ppa_pnc6_trim3 <- readRDS("iso_results/spline_ppa_pnc6_trim3.rds")   
iso_ppa_pnc6_trim3 <- readRDS("iso_results/iso_ppa_pnc6_trim3.rds")

#run and save plot
isoplot_ppa_pnc6_trim3 <- iso_fun_glmer(df_mat_ppa_pnc6_trim3$hb, df_mat_ppa_pnc6_trim3$ppa_pnc6, "PPA at PNC6 - trim3", spline_ppa_pnc6_trim3, iso_ppa_pnc6_trim3)

png(file = "iso_results/isoplot_ppa_pnc6_trim3.png")
isoplot_ppa_pnc6_trim3 <- iso_fun_glmer(df_mat_ppa_pnc6_trim3$hb, df_mat_ppa_pnc6_trim3$ppa_pnc6, "PPA at PNC6 - trim3", spline_ppa_pnc6_trim3, iso_ppa_pnc6_trim3)
dev.off()

#run and save output data
out_ppa_pnc6_trim3  <- outdata(iso_ppa_pnc6_trim3, "PPA-PNC6")
out_ppa_pnc6_trim3 
save(out_ppa_pnc6_trim3, file = "iso_results/out_ppa_pnc6_trim3.rda")
