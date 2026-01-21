#****************************************************************************
#Small for gestational age (<3th)
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

load("derived_data/df_inf_sga3.rda")
load("derived_data/df_inf_sga3_trim1.rda")
load("derived_data/df_inf_sga3_trim2.rda")
load("derived_data/df_inf_sga3_trim3.rda")

#****************************************************************************
#*******all
#read in model
spline_sga3 <- readRDS("iso_results/spline_sga3.rds")   
iso_sga3 <- readRDS("iso_results/iso_sga3.rds")

#run and save plot
isoplot_sga3 <- iso_fun_glmer(df_inf_sga3$hb, df_inf_sga3$sga3, "SGA <3th", spline_sga3, iso_sga3)
png(file = "iso_results/isoplot_sga3.png")
isoplot_sga3 <- iso_fun_glmer(df_inf_sga3$hb, df_inf_sga3$sga3, "SGA <3th", spline_sga3, iso_sga3)
dev.off()

#run and save output data
out_sga3  <- outdata(iso_sga3, "SGA<3th")
out_sga3 
save(out_sga3, file = "iso_results/out_sga3.rda")

#****************************************************************************
#*******trim1
#read in model
spline_sga3_trim1 <- readRDS("iso_results/spline_sga3_trim1.rds")   
iso_sga3_trim1 <- readRDS("iso_results/iso_sga3_trim1.rds")

#run and save plot
isoplot_sga3_trim1 <- iso_fun_glmer(df_inf_sga3_trim1$hb, df_inf_sga3_trim1$sga3, "SGA <3th - trim1", spline_sga3_trim1, iso_sga3_trim1)

png(file = "iso_results/isoplot_sga3_trim1.png")
isoplot_sga3_trim1 <- iso_fun_glmer(df_inf_sga3_trim1$hb, df_inf_sga3_trim1$sga3, "SGA <3th - trim1", spline_sga3_trim1, iso_sga3_trim1)
dev.off()

#run and save output data
out_sga3_trim1  <- outdata(iso_sga3_trim1, "SGA <3th - trim1")
out_sga3_trim1 
save(out_sga3_trim1, file = "iso_results/out_sga3_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_sga3_trim2 <- readRDS("iso_results/spline_sga3_trim2.rds")   
iso_sga3_trim2 <- readRDS("iso_results/iso_sga3_trim2.rds")

#run and save plot
isoplot_sga3_trim2 <- iso_fun_glmer(df_inf_sga3_trim2$hb, df_inf_sga3_trim2$sga3, "SGA <3th - trim2", spline_sga3_trim2, iso_sga3_trim2)

png(file = "iso_results/isoplot_sga3_trim2.png")
isoplot_sga3_trim2 <- iso_fun_glmer(df_inf_sga3_trim2$hb, df_inf_sga3_trim2$sga3, "SGA <3th - trim2", spline_sga3_trim2, iso_sga3_trim2)
dev.off()

#run and save output data
out_sga3_trim2  <- outdata(iso_sga3_trim2, "SGA <3th - trim2")
out_sga3_trim2 
save(out_sga3_trim2, file = "iso_results/out_sga3_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_sga3_trim3 <- readRDS("iso_results/spline_sga3_trim3.rds")   
iso_sga3_trim3 <- readRDS("iso_results/iso_sga3_trim3.rds")

#run and save plot
isoplot_sga3_trim3 <- iso_fun_glmer(df_inf_sga3_trim3$hb, df_inf_sga3_trim3$sga3, "SGA <3th - trim3", spline_sga3_trim3, iso_sga3_trim3)

png(file = "iso_results/isoplot_sga3_trim3.png")
isoplot_sga3_trim3 <- iso_fun_glmer(df_inf_sga3_trim3$hb, df_inf_sga3_trim3$sga3, "SGA <3th - trim3", spline_sga3_trim3, iso_sga3_trim3)
dev.off()

#run and save output data
out_sga3_trim3  <- outdata(iso_sga3_trim3, "SGA <3th - trim3")
out_sga3_trim3 
save(out_sga3_trim3, file = "iso_results/out_sga3_trim3.rda")
