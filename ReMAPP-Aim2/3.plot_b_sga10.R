#****************************************************************************
#Small for gestational age (<10th)
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_sga10 <- readRDS("iso_results/spline_sga10.rds")   
iso_sga10 <- readRDS("iso_results/iso_sga10.rds")

#run and save plot
isoplot_sga10 <- iso_fun_glmer(df_inf_sga10$hb, df_inf_sga10$sga10, "SGA <10th", spline_sga10, iso_sga10)

png(file = "iso_results/isoplot_sga10.png")
isoplot_sga10 <- iso_fun_glmer(df_inf_sga10$hb, df_inf_sga10$sga10, "SGA <10th", spline_sga10, iso_sga10)
dev.off()

#run and save output data
out_sga10  <- outdata(iso_sga10, "SGA<10th")
out_sga10 
save(out_sga10, file = "iso_results/out_sga10.rda")

#****************************************************************************
#*******trim1
#read in model
spline_sga10_trim1 <- readRDS("iso_results/spline_sga10_trim1.rds")   
iso_sga10_trim1 <- readRDS("iso_results/iso_sga10_trim1.rds")

#run and save plot
isoplot_sga10_trim1 <- iso_fun_glmer(df_inf_sga10_trim1$hb, df_inf_sga10_trim1$sga10, "SGA <10th - trim1", spline_sga10_trim1, iso_sga10_trim1)

png(file = "iso_results/isoplot_sga10_trim1.png")
isoplot_sga10_trim1 <- iso_fun_glmer(df_inf_sga10_trim1$hb, df_inf_sga10_trim1$sga10, "SGA <10th - trim1", spline_sga10_trim1, iso_sga10_trim1)
dev.off()

#run and save output data
out_sga10_trim1  <- outdata(iso_sga10_trim1, "SGA <10th - trim1")
out_sga10_trim1 
save(out_sga10_trim1, file = "iso_results/out_sga10_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_sga10_trim2 <- readRDS("iso_results/spline_sga10_trim2.rds")   
iso_sga10_trim2 <- readRDS("iso_results/iso_sga10_trim2.rds")

#run and save plot
isoplot_sga10_trim2 <- iso_fun_glmer(df_inf_sga10_trim2$hb, df_inf_sga10_trim2$sga10, "SGA <10th - trim2", spline_sga10_trim2, iso_sga10_trim2)

png(file = "iso_results/isoplot_sga10_trim2.png")
isoplot_sga10_trim2 <- iso_fun_glmer(df_inf_sga10_trim2$hb, df_inf_sga10_trim2$sga10, "SGA <10th - trim2", spline_sga10_trim2, iso_sga10_trim2)
dev.off()

#run and save output data
out_sga10_trim2  <- outdata(iso_sga10_trim2, "SGA <10th - trim2")
out_sga10_trim2 
save(out_sga10_trim2, file = "iso_results/out_sga10_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_sga10_trim3 <- readRDS("iso_results/spline_sga10_trim3.rds")   
iso_sga10_trim3 <- readRDS("iso_results/iso_sga10_trim3.rds")

#run and save plot
isoplot_sga10_trim3 <- iso_fun_glmer(df_inf_sga10_trim3$hb, df_inf_sga10_trim3$sga10, "SGA <10th - trim3", spline_sga10_trim3, iso_sga10_trim3)

png(file = "iso_results/isoplot_sga10_trim3.png")
isoplot_sga10_trim3 <- iso_fun_glmer(df_inf_sga10_trim3$hb, df_inf_sga10_trim3$sga10, "SGA <10th - trim3", spline_sga10_trim3, iso_sga10_trim3)
dev.off()

#run and save output data
out_sga10_trim3  <- outdata(iso_sga10_trim3, "SGA <10th - trim3")
out_sga10_trim3 
save(out_sga10_trim3, file = "iso_results/out_sga10_trim3.rda")

