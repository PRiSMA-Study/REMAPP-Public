#****************************************************************************
#Birth asphyxia
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

load("derived_data/df_inf_asph.rda")
load("derived_data/df_inf_asph_trim1.rda")
load("derived_data/df_inf_asph_trim2.rda")
load("derived_data/df_inf_asph_trim3.rda")

#****************************************************************************
#*******all
#read in model
spline_asph <- readRDS("iso_results/spline_asph.rds")   
iso_asph <- readRDS("iso_results/iso_asph.rds")

#run and save plot
isoplot_asph <- iso_fun_glmer(df_inf_asph$hb, df_inf_asph$inf_asph, "Brith asphyxia", spline_asph, iso_asph)

png(file = "iso_results/isoplot_asph.png")
isoplot_asph <- iso_fun_glmer(df_inf_asph$hb, df_inf_asph$inf_asph, "Brith asphyxia", spline_asph, iso_asph)
dev.off()

#run and save output data
out_asph <- outdata(iso_asph, "ASPH")
out_asph
save(out_asph, file = "iso_results/out_asph.rda")

#****************************************************************************
#*******trim1
#read in model
spline_asph_trim1 <- readRDS("iso_results/spline_asph_trim1.rds")   
iso_asph_trim1 <- readRDS("iso_results/iso_asph_trim1.rds")

#run and save plot
isoplot_asph_trim1 <- iso_fun_glmer(df_inf_asph_trim1$hb, df_inf_asph_trim1$inf_asph, "Brith asphyxia - trim1", spline_asph_trim1, iso_asph_trim1)

png(file = "iso_results/isoplot_asph_trim1.png")
isoplot_asph_trim1 <- iso_fun_glmer(df_inf_asph_trim1$hb, df_inf_asph_trim1$inf_asph, "Brith asphyxia - trim1", spline_asph_trim1, iso_asph_trim1)
dev.off()

#run and save output data
out_asph_trim1 <- outdata(iso_asph_trim1, "ASPH_trim1")
out_asph_trim1
save(out_asph_trim1, file = "iso_results/out_asph_trim1.rda")

#****************************************************************************
#*******trim2
#read in model
spline_asph_trim2 <- readRDS("iso_results/spline_asph_trim2.rds")   
iso_asph_trim2 <- readRDS("iso_results/iso_asph_trim2.rds")

#run and save plot
isoplot_asph_trim2 <- iso_fun_glmer(df_inf_asph_trim2$hb, df_inf_asph_trim2$inf_asph, "Brith asphyxia - trim2", spline_asph_trim2, iso_asph_trim2)

png(file = "iso_results/isoplot_asph_trim2.png")
isoplot_asph_trim2 <- iso_fun_glmer(df_inf_asph_trim2$hb, df_inf_asph_trim2$inf_asph, "Brith asphyxia - trime2", spline_asph_trim2, iso_asph_trim2)
dev.off()

#run and save output data
out_asph_trim2 <- outdata(iso_asph_trim2, "ASPH_trim2")
out_asph_trim2
save(out_asph_trim2, file = "iso_results/out_asph_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_asph_trim3 <- readRDS("iso_results/spline_asph_trim3.rds")   
iso_asph_trim3 <- readRDS("iso_results/iso_asph_trim3.rds")

#run and save plot
isoplot_asph_trim3 <- iso_fun_glmer(df_inf_asph_trim3$hb, df_inf_asph_trim3$inf_asph, "Brith asphyxia - trim3", spline_asph_trim3, iso_asph_trim3)

png(file = "iso_results/isoplot_asph_trim3.png")
isoplot_asph_trim3 <- iso_fun_glmer(df_inf_asph_trim3$hb, df_inf_asph_trim3$inf_asph, "Brith asphyxia - trim3", spline_asph_trim3, iso_asph_trim3)
dev.off()

#run and save output data
out_asph_trim3 <- outdata(iso_asph_trim3, "ASPH_trim3")
out_asph_trim3
save(out_asph_trim3, file = "iso_results/out_asph_trim3.rda")
 