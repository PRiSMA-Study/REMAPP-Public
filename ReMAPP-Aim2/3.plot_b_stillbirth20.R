#****************************************************************************
#Small for gestational age (<10th)
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_stillbirth20 <- readRDS("iso_results/spline_stillbirth20.rds")   
iso_stillbirth20 <- readRDS("iso_results/iso_stillbirth20.rds")

#run and save plot
isoplot_stillbirth20 <- iso_fun_glmer(df_inf_stillbirth20$hb, df_inf_stillbirth20$inf_stillbirth20, "Stillbirth >=20 weeks", spline_stillbirth20, iso_stillbirth20)

png(file = "iso_results/isoplot_stillbirth20.png")
isoplot_stillbirth20 <- iso_fun_glmer(df_inf_stillbirth20$hb, df_inf_stillbirth20$inf_stillbirth20, "Stillbirth >=20 weeks", spline_stillbirth20, iso_stillbirth20)
dev.off()

#run and save output data
out_stillbirth20  <- outdata(iso_stillbirth20, "Stillbirth>=20weeks")
out_stillbirth20 
save(out_stillbirth20, file = "iso_results/out_stillbirth20.rda")

#****************************************************************************
#*******trim1
#read in model
spline_stillbirth20_trim1 <- readRDS("iso_results/spline_stillbirth20_trim1.rds")   
iso_stillbirth20_trim1 <- readRDS("iso_results/iso_stillbirth20_trim1.rds")

#run and save plot
isoplot_stillbirth20_trim1 <- iso_fun_glmer(df_inf_stillbirth20_trim1$hb, df_inf_stillbirth20_trim1$inf_stillbirth20, "Stillbirth >=20 weeks - trim1", spline_stillbirth20_trim1, iso_stillbirth20_trim1)

png(file = "iso_results/isoplot_stillbirth20_trim1.png")
isoplot_stillbirth20_trim1 <- iso_fun_glmer(df_inf_stillbirth20_trim1$hb, df_inf_stillbirth20_trim1$inf_stillbirth20, "Stillbirth >=20 weeks - trim1", spline_stillbirth20_trim1, iso_stillbirth20_trim1)
dev.off()

#run and save output data
out_stillbirth20_trim1  <- outdata(iso_stillbirth20_trim1, "Stillbirth >=20 weeks - trim1")
out_stillbirth20_trim1 
save(out_stillbirth20_trim1, file = "iso_results/out_stillbirth20_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_stillbirth20_trim2 <- readRDS("iso_results/spline_stillbirth20_trim2.rds")   
iso_stillbirth20_trim2 <- readRDS("iso_results/iso_stillbirth20_trim2.rds")

#run and save plot
isoplot_stillbirth20_trim2 <- iso_fun_glmer(df_inf_stillbirth20_trim2$hb, df_inf_stillbirth20_trim2$inf_stillbirth20, "Stillbirth >=20 weeks - trim2", spline_stillbirth20_trim2, iso_stillbirth20_trim2)

png(file = "iso_results/isoplot_stillbirth20_trim2.png")
isoplot_stillbirth20_trim2 <- iso_fun_glmer(df_inf_stillbirth20_trim2$hb, df_inf_stillbirth20_trim2$inf_stillbirth20, "Stillbirth >=20 weeks - trim2", spline_stillbirth20_trim2, iso_stillbirth20_trim2)
dev.off()

#run and save output data
out_stillbirth20_trim2  <- outdata(iso_stillbirth20_trim2, "Stillbirth >=20 weeks - trim2")
out_stillbirth20_trim2 
save(out_stillbirth20_trim2, file = "iso_results/out_stillbirth20_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_stillbirth20_trim3 <- readRDS("iso_results/spline_stillbirth20_trim3.rds")   
iso_stillbirth20_trim3 <- readRDS("iso_results/iso_stillbirth20_trim3.rds")

#run and save plot
isoplot_stillbirth20_trim3 <- iso_fun_glmer(df_inf_stillbirth20_trim3$hb, df_inf_stillbirth20_trim3$inf_stillbirth20, "Stillbirth >=20 weeks - trim3", spline_stillbirth20_trim3, iso_stillbirth20_trim3)

png(file = "iso_results/isoplot_stillbirth20_trim3.png")
isoplot_stillbirth20_trim3 <- iso_fun_glmer(df_inf_stillbirth20_trim3$hb, df_inf_stillbirth20_trim3$inf_stillbirth20, "Stillbirth >=20 weeks - trim3", spline_stillbirth20_trim3, iso_stillbirth20_trim3)
dev.off()

#run and save output data
out_stillbirth20_trim3  <- outdata(iso_stillbirth20_trim3, "Stillbirth >=20 weeks - trim3")
out_stillbirth20_trim3 
save(out_stillbirth20_trim3, file = "iso_results/out_stillbirth20_trim3.rda")
