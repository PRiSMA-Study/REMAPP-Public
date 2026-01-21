#****************************************************************************
#Maternal Postpartum hemorrhage
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

#****************************************************************************
#*******all

load("derived_data/df_mat_pph.rda")
load("derived_data/df_mat_pph_trim1.rda")
load("derived_data/df_mat_pph_trim2.rda")
load("derived_data/df_mat_pph_trim3.rda")

#read in model
spline_pph <- readRDS("iso_results/spline_pph.rds")   
iso_pph <- readRDS("iso_results/iso_pph.rds")

#run and save plot
isoplot_pph <- iso_fun_glmer(df_mat_pph$hb, df_mat_pph$HEM_PPH, "Postpartum hemorrhage", spline_pph, iso_pph)

png(file = "iso_results/isoplot_pph.png")
isoplot_pph <- iso_fun_glmer(df_mat_pph$hb, df_mat_pph$HEM_PPH, "Postpartum hemorrhage", spline_pph, iso_pph)
dev.off()

#run and save output data
out_pph  <- outdata(iso_pph, "PPH")
out_pph 
save(out_pph, file = "iso_results/out_pph.rda")

#****************************************************************************
#*******trim1
#read in model
spline_pph_trim1 <- readRDS("iso_results/spline_pph_trim1.rds")   
iso_pph_trim1 <- readRDS("iso_results/iso_pph_trim1.rds")

#run and save plot
isoplot_pph_trim1 <- iso_fun_glmer(df_mat_pph_trim1$hb, df_mat_pph_trim1$HEM_PPH, "PPH - trim1", spline_pph_trim1, iso_pph_trim1)

png(file = "iso_results/isoplot_pph_trim1.png")
isoplot_pph_trim1 <- iso_fun_glmer(df_mat_pph_trim1$hb, df_mat_pph_trim1$HEM_PPH, "PPH - trim1", spline_pph_trim1, iso_pph_trim1)
dev.off()

#run and save output data
out_pph_trim1  <- outdata(iso_pph_trim1, "PPH")
out_pph_trim1 
save(out_pph_trim1, file = "iso_results/out_pph_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_pph_trim2 <- readRDS("iso_results/spline_pph_trim2.rds")   
iso_pph_trim2 <- readRDS("iso_results/iso_pph_trim2.rds")

#run and save plot
isoplot_pph_trim2 <- iso_fun_glmer(df_mat_pph_trim2$hb, df_mat_pph_trim2$HEM_PPH, "PPH - trim2", spline_pph_trim2, iso_pph_trim2)

png(file = "iso_results/isoplot_pph_trim2.png")
isoplot_pph_trim2 <- iso_fun_glmer(df_mat_pph_trim2$hb, df_mat_pph_trim2$HEM_PPH, "PPH - trim2", spline_pph_trim2, iso_pph_trim2)
dev.off()

#run and save output data
out_pph_trim2  <- outdata(iso_pph_trim2, "PPH")
out_pph_trim2 
save(out_pph_trim2, file = "iso_results/out_pph_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_pph_trim3 <- readRDS("iso_results/spline_pph_trim3.rds")   
iso_pph_trim3 <- readRDS("iso_results/iso_pph_trim3.rds")

#run and save plot
isoplot_pph_trim3 <- iso_fun_glmer(df_mat_pph_trim3$hb, df_mat_pph_trim3$HEM_PPH, "PPH - trim3", spline_pph_trim3, iso_pph_trim3)

png(file = "iso_results/isoplot_pph_trim3.png")
isoplot_pph_trim3 <- iso_fun_glmer(df_mat_pph_trim3$hb, df_mat_pph_trim3$HEM_PPH, "PPH - trim3", spline_pph_trim3, iso_pph_trim3)
dev.off()

#run and save output data
out_pph_trim3  <- outdata(iso_pph_trim3, "PPH")
out_pph_trim3 
save(out_pph_trim3, file = "iso_results/out_pph_trim3.rda")
