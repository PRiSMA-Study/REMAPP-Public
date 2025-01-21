#****************************************************************************
#PSBI
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_psbi <- readRDS("iso_results/spline_psbi.rds")   
iso_psbi <- readRDS("iso_results/iso_psbi.rds")

#run and save plot
isoplot_psbi <- iso_fun_glmer(df_inf_psbi$hb, df_inf_psbi$inf_psbi, "PSBI", spline_psbi, iso_psbi)

png(file = "iso_results/isoplot_psbi.png")
isoplot_psbi <- iso_fun_glmer(df_inf_psbi$hb, df_inf_psbi$inf_psbi, "PSBI", spline_psbi, iso_psbi)
dev.off()

#run and save output data
out_psbi  <- outdata(iso_psbi, "PSBI")
out_psbi 
save(out_psbi, file = "iso_results/out_psbi.rda")

#****************************************************************************
#*******trim1
#read in model
spline_psbi_trim1 <- readRDS("iso_results/spline_psbi_trim1.rds")   
iso_psbi_trim1 <- readRDS("iso_results/iso_psbi_trim1.rds")

#run and save plot
isoplot_psbi_trim1 <- iso_fun_glmer(df_inf_psbi_trim1$hb, df_inf_psbi_trim1$inf_psbi, "PSBI - trim1", spline_psbi_trim1, iso_psbi_trim1)

png(file = "iso_results/isoplot_psbi_trim1.png")
isoplot_psbi_trim1 <- iso_fun_glmer(df_inf_psbi_trim1$hb, df_inf_psbi_trim1$inf_psbi, "PSBI - trim1", spline_psbi_trim1, iso_psbi_trim1)
dev.off()

#run and save output data
out_psbi_trim1  <- outdata(iso_psbi_trim1, "PSBI - trim1")
out_psbi_trim1 
save(out_psbi_trim1, file = "iso_results/out_psbi_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_psbi_trim2 <- readRDS("iso_results/spline_psbi_trim2.rds")   
iso_psbi_trim2 <- readRDS("iso_results/iso_psbi_trim2.rds")

#run and save plot
isoplot_psbi_trim2 <- iso_fun_glmer(df_inf_psbi_trim2$hb, df_inf_psbi_trim2$inf_psbi, "PSBI - trim2", spline_psbi_trim2, iso_psbi_trim2)

png(file = "iso_results/isoplot_psbi_trim2.png")
isoplot_psbi_trim2 <- iso_fun_glmer(df_inf_psbi_trim2$hb, df_inf_psbi_trim2$inf_psbi, "PSBI - trim2", spline_psbi_trim2, iso_psbi_trim2)
dev.off()

#run and save output data
out_psbi_trim2  <- outdata(iso_psbi_trim2, "PSBI - trim2")
out_psbi_trim2 
save(out_psbi_trim2, file = "iso_results/out_psbi_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_psbi_trim3 <- readRDS("iso_results/spline_psbi_trim3.rds")   
iso_psbi_trim3 <- readRDS("iso_results/iso_psbi_trim3.rds")

#run and save plot
isoplot_psbi_trim3 <- iso_fun_glmer(df_inf_psbi_trim3$hb, df_inf_psbi_trim3$inf_psbi, "PSBI - trim3", spline_psbi_trim3, iso_psbi_trim3)

png(file = "iso_results/isoplot_psbi_trim3.png")
isoplot_psbi_trim3 <- iso_fun_glmer(df_inf_psbi_trim3$hb, df_inf_psbi_trim3$inf_psbi, "PSBI - trim3", spline_psbi_trim3, iso_psbi_trim3)
dev.off()

#run and save output data
out_psbi_trim3  <- outdata(iso_psbi_trim3, "PSBI - trim3")
out_psbi_trim3 
save(out_psbi_trim3, file = "iso_results/out_psbi_trim3.rda")
