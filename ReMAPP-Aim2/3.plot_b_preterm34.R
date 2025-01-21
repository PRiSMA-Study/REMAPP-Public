#****************************************************************************
#Preterm <34 weeks
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_preterm34 <- readRDS("iso_results/spline_preterm34.rds")   
iso_preterm34 <- readRDS("iso_results/iso_preterm34.rds")

#run and save plot
isoplot_preterm34 <- iso_fun_glmer(df_inf_preterm34$hb, df_inf_preterm34$preterm34, "Preterm <34 weeks", spline_preterm34, iso_preterm34)

png(file = "iso_results/isoplot_preterm34.png")
isoplot_preterm34 <- iso_fun_glmer(df_inf_preterm34$hb, df_inf_preterm34$preterm34, "Preterm <34 weeks", spline_preterm34, iso_preterm34)
dev.off()

#run and save output data
out_preterm34  <- outdata(iso_preterm34, "Preterm<34weeks")
out_preterm34 
save(out_preterm34, file = "iso_results/out_preterm34.rda")

#****************************************************************************
#*******trim1
#read in model
spline_preterm34_trim1 <- readRDS("iso_results/spline_preterm34_trim1.rds")   
iso_preterm34_trim1 <- readRDS("iso_results/iso_preterm34_trim1.rds")

#run and save plot
isoplot_preterm34_trim1 <- iso_fun_glmer(df_inf_preterm34_trim1$hb, df_inf_preterm34_trim1$preterm34, "Preterm <34 weeks - trim1", spline_preterm34_trim1, iso_preterm34_trim1)

png(file = "iso_results/isoplot_preterm34_trim1.png")
isoplot_preterm34_trim1 <- iso_fun_glmer(df_inf_preterm34_trim1$hb, df_inf_preterm34_trim1$preterm34, "Preterm <34 weeks - trim1", spline_preterm34_trim1, iso_preterm34_trim1)
dev.off()

#run and save output data
out_preterm34_trim1  <- outdata(iso_preterm34_trim1, "Preterm <34 weeks - trim1")
out_preterm34_trim1 
save(out_preterm34_trim1, file = "iso_results/out_preterm34_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_preterm34_trim2 <- readRDS("iso_results/spline_preterm34_trim2.rds")   
iso_preterm34_trim2 <- readRDS("iso_results/iso_preterm34_trim2.rds")

#run and save plot
isoplot_preterm34_trim2 <- iso_fun_glmer(df_inf_preterm34_trim2$hb, df_inf_preterm34_trim2$preterm34, "Preterm <34 weeks - trim2", spline_preterm34_trim2, iso_preterm34_trim2)

png(file = "iso_results/isoplot_preterm34_trim2.png")
isoplot_preterm34_trim2 <- iso_fun_glmer(df_inf_preterm34_trim2$hb, df_inf_preterm34_trim2$preterm34, "Preterm <34 weeks - trim2", spline_preterm34_trim2, iso_preterm34_trim2)
dev.off()

#run and save output data
out_preterm34_trim2  <- outdata(iso_preterm34_trim2, "Preterm <34 weeks - trim2")
out_preterm34_trim2 
save(out_preterm34_trim2, file = "iso_results/out_preterm34_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_preterm34_trim3 <- readRDS("iso_results/spline_preterm34_trim3.rds")   
iso_preterm34_trim3 <- readRDS("iso_results/iso_preterm34_trim3.rds")

#run and save plot
isoplot_preterm34_trim3 <- iso_fun_glmer(df_inf_preterm34_trim3$hb, df_inf_preterm34_trim3$preterm34, "Preterm <34 weeks - trim3", spline_preterm34_trim3, iso_preterm34_trim3)

png(file = "iso_results/isoplot_preterm34_trim3.png")
isoplot_preterm34_trim3 <- iso_fun_glmer(df_inf_preterm34_trim3$hb, df_inf_preterm34_trim3$preterm34, "Preterm <34 weeks - trim3", spline_preterm34_trim3, iso_preterm34_trim3)
dev.off()

#run and save output data
out_preterm34_trim3  <- outdata(iso_preterm34_trim3, "Preterm <34 weeks - trim3")
out_preterm34_trim3 
save(out_preterm34_trim3, file = "iso_results/out_preterm34_trim3.rda")
