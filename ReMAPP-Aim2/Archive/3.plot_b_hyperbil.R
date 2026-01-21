#****************************************************************************
#Neonatal hyperbilirubinemia
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_binary.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")

load("derived_data/df_inf_hyperbili.rda")
load("derived_data/df_inf_hyperbili_trim1.rda")
load("derived_data/df_inf_hyperbili_trim2.rda")
load("derived_data/df_inf_hyperbili_trim3.rda")


#****************************************************************************
#*******all
#read in model
spline_hyperbili <- readRDS("iso_results/spline_hyperbili.rds")   
iso_hyperbili <- readRDS("iso_results/iso_hyperbili.rds")

#run and save plot
isoplot_hyperbili <- iso_fun_glmer(df_inf_hyperbili$hb, df_inf_hyperbili$hyperbili, "Hyperbilirubinemia", spline_hyperbili, iso_hyperbili)

png(file = "iso_results/isoplot_hyperbili .png")
isoplot_hyperbili <- iso_fun_glmer(df_inf_hyperbili$hb, df_inf_hyperbili$hyperbili, "Hyperbilirubinemia", spline_hyperbili, iso_hyperbili)
dev.off()

#run and save output data
out_hyperbili  <- outdata(iso_hyperbili, "Hyperbilirubinemia")
out_hyperbili 
save(out_hyperbili, file = "iso_results/out_hyperbili.rda")

#****************************************************************************
#*******trim1
#read in model
spline_hyperbili_trim1 <- readRDS("iso_results/spline_hyperbili_trim1.rds")   
iso_hyperbili_trim1 <- readRDS("iso_results/iso_hyperbili_trim1.rds")

#run and save plot
isoplot_hyperbili_trim1 <- iso_fun_glmer(df_inf_hyperbili_trim1$hb, df_inf_hyperbili_trim1$hyperbili, "Hyperbilirubinemia - trim1", spline_hyperbili_trim1, iso_hyperbili_trim1)

png(file = "iso_results/isoplot_hyperbili_trim1 .png")
isoplot_hyperbili_trim1 <- iso_fun_glmer(df_inf_hyperbili_trim1$hb, df_inf_hyperbili_trim1$hyperbili, "Hyperbilirubinemia - trim1", spline_hyperbili_trim1, iso_hyperbili_trim1)
dev.off()

#run and save output data
out_hyperbili_trim1  <- outdata(iso_hyperbili_trim1, "Hyperbilirubinemia")
out_hyperbili_trim1 
save(out_hyperbili_trim1, file = "iso_results/out_hyperbili_trim1.rda")

#****************************************************************************
#*******trim2
#read in model
spline_hyperbili_trim2 <- readRDS("iso_results/spline_hyperbili_trim2.rds")   
iso_hyperbili_trim2 <- readRDS("iso_results/iso_hyperbili_trim2.rds")

#run and save plot
isoplot_hyperbili_trim2 <- iso_fun_glmer(df_inf_hyperbili_trim2$hb, df_inf_hyperbili_trim2$hyperbili, "Hyperbilirubinemia - trim2", spline_hyperbili_trim2, iso_hyperbili_trim2)

png(file = "iso_results/isoplot_hyperbili_trim2 .png")
isoplot_hyperbili_trim2 <- iso_fun_glmer(df_inf_hyperbili_trim2$hb, df_inf_hyperbili_trim2$hyperbili, "Hyperbilirubinemia - trim2", spline_hyperbili_trim2, iso_hyperbili_trim2)
dev.off()

#run and save output data
out_hyperbili_trim2  <- outdata(iso_hyperbili_trim2, "Hyperbilirubinemia")
out_hyperbili_trim2 
save(out_hyperbili_trim2, file = "iso_results/out_hyperbili_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_hyperbili_trim3 <- readRDS("iso_results/spline_hyperbili_trim3.rds")   
iso_hyperbili_trim3 <- readRDS("iso_results/iso_hyperbili_trim3.rds")

#run and save plot
isoplot_hyperbili_trim3 <- iso_fun_glmer(df_inf_hyperbili_trim3$hb, df_inf_hyperbili_trim3$hyperbili, "Hyperbilirubinemia - trim3", spline_hyperbili_trim3, iso_hyperbili_trim3)

png(file = "iso_results/isoplot_hyperbili_trim3 .png")
isoplot_hyperbili_trim3 <- iso_fun_glmer(df_inf_hyperbili_trim3$hb, df_inf_hyperbili_trim3$hyperbili, "Hyperbilirubinemia - trim3", spline_hyperbili_trim3, iso_hyperbili_trim3)
dev.off()

#run and save output data
out_hyperbili_trim3  <- outdata(iso_hyperbili_trim3, "Hyperbilirubinemia")
out_hyperbili_trim3 
save(out_hyperbili_trim3, file = "iso_results/out_hyperbili_trim3.rda")
