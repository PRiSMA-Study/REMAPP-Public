#****************************************************************************
#Low birth weight (<2500g)
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_lbw2500 <- readRDS("iso_results/spline_lbw2500.rds")   
iso_lbw2500 <- readRDS("iso_results/iso_lbw2500.rds")

#run and save plot
isoplot_lbw2500 <- iso_fun_glmer(df_inf_lbw2500$hb, df_inf_lbw2500$lbw2500, "Low birth weight (<2500g)", spline_lbw2500, iso_lbw2500)

png(file = "iso_results/isoplot_lbw2500.png")
isoplot_lbw2500 <- iso_fun_glmer(df_inf_lbw2500$hb, df_inf_lbw2500$lbw2500, "Low birth weight (<2500g)", spline_lbw2500, iso_lbw2500)
dev.off()

#run and save output data
out_lbw2500  <- outdata(iso_lbw2500, "LBW<2500g")
out_lbw2500 
save(out_lbw2500, file = "iso_results/out_lbw2500.rda")

#****************************************************************************
#*******trim1
#read in model
spline_lbw2500_trim1 <- readRDS("iso_results/spline_lbw2500_trim1.rds")   
iso_lbw2500_trim1 <- readRDS("iso_results/iso_lbw2500_trim1.rds")

#run and save plot
isoplot_lbw2500_trim1 <- iso_fun_glmer(df_inf_lbw2500_trim1$hb, df_inf_lbw2500_trim1$lbw2500, "LBW (<2500g) - trim1", spline_lbw2500_trim1, iso_lbw2500_trim1)

png(file = "iso_results/isoplot_lbw2500_trim1.png")
isoplot_lbw2500_trim1 <- iso_fun_glmer(df_inf_lbw2500_trim1$hb, df_inf_lbw2500_trim1$lbw2500, "LBW (<2500g) - trim1", spline_lbw2500_trim1, iso_lbw2500_trim1)
dev.off()

#run and save output data
out_lbw2500_trim1  <- outdata(iso_lbw2500_trim1, "LBW<2500g")
out_lbw2500_trim1 
save(out_lbw2500_trim1, file = "iso_results/out_lbw2500_trim1.rda")

#****************************************************************************
#*******trim2
#read in model
spline_lbw2500_trim2 <- readRDS("iso_results/spline_lbw2500_trim2.rds")   
iso_lbw2500_trim2 <- readRDS("iso_results/iso_lbw2500_trim2.rds")

#run and save plot
isoplot_lbw2500_trim2 <- iso_fun_glmer(df_inf_lbw2500_trim2$hb, df_inf_lbw2500_trim2$lbw2500, "LBW (<2500g) - trim2", spline_lbw2500_trim2, iso_lbw2500_trim2)

png(file = "iso_results/isoplot_lbw2500_trim2.png")
isoplot_lbw2500_trim2 <- iso_fun_glmer(df_inf_lbw2500_trim2$hb, df_inf_lbw2500_trim2$lbw2500, "LBW (<2500g) - trim2", spline_lbw2500_trim2, iso_lbw2500_trim2)
dev.off()

#run and save output data 
out_lbw2500_trim2  <- outdata(iso_lbw2500_trim2, "LBW<2500g")
out_lbw2500_trim2 
save(out_lbw2500_trim2, file = "iso_results/out_lbw2500_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_lbw2500_trim3 <- readRDS("iso_results/spline_lbw2500_trim3.rds")   
iso_lbw2500_trim3 <- readRDS("iso_results/iso_lbw2500_trim3.rds")

#run and save plot
isoplot_lbw2500_trim3 <- iso_fun_glmer(df_inf_lbw2500_trim3$hb, df_inf_lbw2500_trim3$lbw2500, "LBW (<2500g) - trim3", spline_lbw2500_trim3, iso_lbw2500_trim3)

png(file = "iso_results/isoplot_lbw2500_trim3.png")
isoplot_lbw2500_trim3 <- iso_fun_glmer(df_inf_lbw2500_trim3$hb, df_inf_lbw2500_trim3$lbw2500, "LBW (<2500g) - trim3", spline_lbw2500_trim3, iso_lbw2500_trim3)
dev.off()

#run and save output data
out_lbw2500_trim3  <- outdata(iso_lbw2500_trim3, "LBW<2500g")
out_lbw2500_trim3 
save(out_lbw2500_trim3, file = "iso_results/out_lbw2500_trim3.rda")
