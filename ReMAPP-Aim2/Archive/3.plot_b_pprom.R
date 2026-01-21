#****************************************************************************
#read in model
#****************************************************************************

UploadDate = "2025-04-18"

setwd(file.path("D:/Users/williams_pj/Documents/Analysis/ReMAPP/Aim2", UploadDate))


source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

load("derived_data/df_mat_pprom.rda")
load("derived_data/df_mat_pprom_trim1.rda")
load("derived_data/df_mat_pprom_trim2.rda")
load("derived_data/df_mat_pprom_trim3.rda")

#****************************************************************************
#Preterm premature rupture of membranes
#****************************************************************************
# source("iso_code/plot_binary.R")
# source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_pprom <- readRDS("iso_results/spline_pprom.rds")   
iso_pprom <- readRDS("iso_results/iso_pprom.rds")

#run and save plot
isoplot_pprom <- iso_fun_glmer(df_mat_pprom$hb, df_mat_pprom$pprom, "PPROM", spline_pprom, iso_pprom)

png(file = "iso_results/isoplot_pprom.png")
isoplot_pprom <- iso_fun_glmer(df_mat_pprom$hb, df_mat_pprom$pprom, "PPROM", spline_pprom, iso_pprom)
dev.off()

#run and save output data
out_pprom  <- outdata(iso_pprom, "PPROM")
out_pprom 
save(out_pprom, file = "iso_results/out_pprom.rda")

#****************************************************************************
#*******trim1
#read in model
spline_pprom_trim1 <- readRDS("iso_results/spline_pprom_trim1.rds")   
iso_pprom_trim1 <- readRDS("iso_results/iso_pprom_trim1.rds")

#run and save plot
isoplot_pprom_trim1 <- iso_fun_glmer(df_mat_pprom_trim1$hb, df_mat_pprom_trim1$pprom, "PPROM - trim1", spline_pprom_trim1, iso_pprom_trim1)

png(file = "iso_results/isoplot_pprom_trim1.png")
isoplot_pprom_trim1 <- iso_fun_glmer(df_mat_pprom_trim1$hb, df_mat_pprom_trim1$pprom, "PPROM - trim1", spline_pprom_trim1, iso_pprom_trim1)
dev.off()

#run and save output data
out_pprom_trim1  <- outdata(iso_pprom_trim1, "PPROM")
out_pprom_trim1 
save(out_pprom_trim1, file = "iso_results/out_pprom_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_pprom_trim2 <- readRDS("iso_results/spline_pprom_trim2.rds")   
iso_pprom_trim2 <- readRDS("iso_results/iso_pprom_trim2.rds")

#run and save plot
isoplot_pprom_trim2 <- iso_fun_glmer(df_mat_pprom_trim2$hb, df_mat_pprom_trim2$pprom, "PPROM - trim2", spline_pprom_trim2, iso_pprom_trim2)

png(file = "iso_results/isoplot_pprom_trim2.png")
isoplot_pprom_trim2 <- iso_fun_glmer(df_mat_pprom_trim2$hb, df_mat_pprom_trim2$pprom, "PPROM - trim2", spline_pprom_trim2, iso_pprom_trim2)
dev.off()

#run and save output data
out_pprom_trim2  <- outdata(iso_pprom_trim2, "PPROM")
out_pprom_trim2 
save(out_pprom_trim2, file = "iso_results/out_pprom_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_pprom_trim3 <- readRDS("iso_results/spline_pprom_trim3.rds")   
iso_pprom_trim3 <- readRDS("iso_results/iso_pprom_trim3.rds")

#run and save plot
isoplot_pprom_trim3 <- iso_fun_glmer(df_mat_pprom_trim3$hb, df_mat_pprom_trim3$pprom, "PPROM - trim3", spline_pprom_trim3, iso_pprom_trim3)

png(file = "iso_results/isoplot_pprom_trim3.png")
isoplot_pprom_trim3 <- iso_fun_glmer(df_mat_pprom_trim3$hb, df_mat_pprom_trim3$pprom, "PPROM - trim3", spline_pprom_trim3, iso_pprom_trim3)
dev.off()

#run and save output data
out_pprom_trim3  <- outdata(iso_pprom_trim3, "PPROM")
out_pprom_trim3 
save(out_pprom_trim3, file = "iso_results/out_pprom_trim3.rda")

