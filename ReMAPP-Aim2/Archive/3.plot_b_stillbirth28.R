#****************************************************************************
#Small for gestational age (<10th)
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

load("derived_data/df_inf_stillbirth28.rda")
load("derived_data/df_inf_stillbirth28_trim1.rda")
load("derived_data/df_inf_stillbirth28_trim2.rda")
load("derived_data/df_inf_stillbirth28_trim3.rda")


#****************************************************************************
#*******all
#read in model
spline_stillbirth28 <- readRDS("iso_results/spline_stillbirth28.rds")   
iso_stillbirth28 <- readRDS("iso_results/iso_stillbirth28.rds")

#run and save plot
isoplot_stillbirth28 <- iso_fun_glmer(df_inf_stillbirth28$hb, df_inf_stillbirth28$inf_stillbirth28, "Stillbirth >=28 weeks", spline_stillbirth28, iso_stillbirth28)

png(file = "iso_results/isoplot_stillbirth28.png")
isoplot_stillbirth28 <- iso_fun_glmer(df_inf_stillbirth28$hb, df_inf_stillbirth28$inf_stillbirth28, "Stillbirth >=28 weeks", spline_stillbirth28, iso_stillbirth28)
dev.off()

#run and save output data
out_stillbirth28  <- outdata(iso_stillbirth28, "Stillbirth>=28weeks")
out_stillbirth28 
save(out_stillbirth28, file = "iso_results/out_stillbirth28.rda")

#****************************************************************************
#*******trim1
#read in model
spline_stillbirth28_trim1 <- readRDS("iso_results/spline_stillbirth28_trim1.rds")   
iso_stillbirth28_trim1 <- readRDS("iso_results/iso_stillbirth28_trim1.rds")

#run and save plot
isoplot_stillbirth28_trim1 <- iso_fun_glmer(df_inf_stillbirth28_trim1$hb, df_inf_stillbirth28_trim1$inf_stillbirth28, "Stillbirth >=28 weeks - trim1", spline_stillbirth28_trim1, iso_stillbirth28_trim1)

png(file = "iso_results/isoplot_stillbirth28_trim1.png")
isoplot_stillbirth28_trim1 <- iso_fun_glmer(df_inf_stillbirth28_trim1$hb, df_inf_stillbirth28_trim1$inf_stillbirth28, "Stillbirth >=28 weeks - trim1", spline_stillbirth28_trim1, iso_stillbirth28_trim1)
dev.off()

#run and save output data
out_stillbirth28_trim1  <- outdata(iso_stillbirth28_trim1, "Stillbirth >=28 weeks - trim1")
out_stillbirth28_trim1 
save(out_stillbirth28_trim1, file = "iso_results/out_stillbirth28_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_stillbirth28_trim2 <- readRDS("iso_results/spline_stillbirth28_trim2.rds")   
iso_stillbirth28_trim2 <- readRDS("iso_results/iso_stillbirth28_trim2.rds")

#run and save plot
isoplot_stillbirth28_trim2 <- iso_fun_glmer(df_inf_stillbirth28_trim2$hb, df_inf_stillbirth28_trim2$inf_stillbirth28, "Stillbirth >=28 weeks - trim2", spline_stillbirth28_trim2, iso_stillbirth28_trim2)

png(file = "iso_results/isoplot_stillbirth28_trim2.png")
isoplot_stillbirth28_trim2 <- iso_fun_glmer(df_inf_stillbirth28_trim2$hb, df_inf_stillbirth28_trim2$inf_stillbirth28, "Stillbirth >=28 weeks - trim2", spline_stillbirth28_trim2, iso_stillbirth28_trim2)
dev.off()

#run and save output data
out_stillbirth28_trim2  <- outdata(iso_stillbirth28_trim2, "Stillbirth >=28 weeks - trim2")
out_stillbirth28_trim2 
save(out_stillbirth28_trim2, file = "iso_results/out_stillbirth28_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_stillbirth28_trim3 <- readRDS("iso_results/spline_stillbirth28_trim3.rds")   
iso_stillbirth28_trim3 <- readRDS("iso_results/iso_stillbirth28_trim3.rds")

#run and save plot
isoplot_stillbirth28_trim3 <- iso_fun_glmer(df_inf_stillbirth28_trim3$hb, df_inf_stillbirth28_trim3$inf_stillbirth28, "Stillbirth >=28 weeks - trim3", spline_stillbirth28_trim3, iso_stillbirth28_trim3)

png(file = "iso_results/isoplot_stillbirth28_trim3.png")
isoplot_stillbirth28_trim3 <- iso_fun_glmer(df_inf_stillbirth28_trim3$hb, df_inf_stillbirth28_trim3$inf_stillbirth28, "Stillbirth >=28 weeks - trim3", spline_stillbirth28_trim3, iso_stillbirth28_trim3)
dev.off()

#run and save output data
out_stillbirth28_trim3  <- outdata(iso_stillbirth28_trim3, "Stillbirth >=28 weeks - trim3")
out_stillbirth28_trim3 
save(out_stillbirth28_trim3, file = "iso_results/out_stillbirth28_trim3.rda")
