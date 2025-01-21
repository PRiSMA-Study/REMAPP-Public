#****************************************************************************
#Preterm <37 weeks
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_preterm37 <- readRDS("iso_results/spline_preterm37.rds")   
iso_preterm37 <- readRDS("iso_results/iso_preterm37.rds")

#run and save plot
isoplot_preterm37 <- iso_fun_glmer(df_inf_preterm37$hb, df_inf_preterm37$preterm37, "Preterm <37 weeks", spline_preterm37, iso_preterm37)

png(file = "iso_results/isoplot_preterm37.png")
isoplot_preterm37 <- iso_fun_glmer(df_inf_preterm37$hb, df_inf_preterm37$preterm37, "Preterm <37 weeks", spline_preterm37, iso_preterm37)
dev.off()

#run and save output data
out_preterm37  <- outdata(iso_preterm37, "Preterm<37weeks")
out_preterm37 
save(out_preterm37, file = "iso_results/out_preterm37.rda")

#****************************************************************************
#*******trim1
#read in model
spline_preterm37_trim1 <- readRDS("iso_results/spline_preterm37_trim1.rds")   
iso_preterm37_trim1 <- readRDS("iso_results/iso_preterm37_trim1.rds")

#run and save plot
isoplot_preterm37_trim1 <- iso_fun_glmer(df_inf_preterm37_trim1$hb, df_inf_preterm37_trim1$preterm37, "Preterm <37 weeks - trim1", spline_preterm37_trim1, iso_preterm37_trim1)

png(file = "iso_results/isoplot_preterm37_trim1.png")
isoplot_preterm37_trim1 <- iso_fun_glmer(df_inf_preterm37_trim1$hb, df_inf_preterm37_trim1$preterm37, "Preterm <37 weeks - trim1", spline_preterm37_trim1, iso_preterm37_trim1)
dev.off()

#run and save output data
out_preterm37_trim1  <- outdata(iso_preterm37_trim1, "Preterm <37 weeks - trim1")
out_preterm37_trim1 
save(out_preterm37_trim1, file = "iso_results/out_preterm37_trim1.rda")

#****************************************************************************
#*******trim3
#read in model
spline_preterm37_trim2 <- readRDS("iso_results/spline_preterm37_trim2.rds")   
iso_preterm37_trim2 <- readRDS("iso_results/iso_preterm37_trim2.rds")

#run and save plot
isoplot_preterm37_trim2 <- iso_fun_glmer(df_inf_preterm37_trim2$hb, df_inf_preterm37_trim2$preterm37, "Preterm <37 weeks - trim2", spline_preterm37_trim2, iso_preterm37_trim2)

png(file = "iso_results/isoplot_preterm37_trim2.png")
isoplot_preterm37_trim2 <- iso_fun_glmer(df_inf_preterm37_trim2$hb, df_inf_preterm37_trim2$preterm37, "Preterm <37 weeks - trim2", spline_preterm37_trim2, iso_preterm37_trim2)
dev.off()

#run and save output data
out_preterm37_trim2  <- outdata(iso_preterm37_trim2, "Preterm <37 weeks - trim2")
out_preterm37_trim2 
save(out_preterm37_trim2, file = "iso_results/out_preterm37_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_preterm37_trim3 <- readRDS("iso_results/spline_preterm37_trim3.rds")   
iso_preterm37_trim3 <- readRDS("iso_results/iso_preterm37_trim3.rds")

#run and save plot
isoplot_preterm37_trim3 <- iso_fun_glmer(df_inf_preterm37_trim3$hb, df_inf_preterm37_trim3$preterm37, "Preterm <37 weeks - trim3", spline_preterm37_trim3, iso_preterm37_trim3)

png(file = "iso_results/isoplot_preterm37_trim3.png")
isoplot_preterm37_trim3 <- iso_fun_glmer(df_inf_preterm37_trim3$hb, df_inf_preterm37_trim3$preterm37, "Preterm <37 weeks - trim3", spline_preterm37_trim3, iso_preterm37_trim3)
dev.off()

#run and save output data
out_preterm37_trim3  <- outdata(iso_preterm37_trim3, "Preterm <37 weeks - trim3")
out_preterm37_trim3 
save(out_preterm37_trim3, file = "iso_results/out_preterm37_trim3.rda")
