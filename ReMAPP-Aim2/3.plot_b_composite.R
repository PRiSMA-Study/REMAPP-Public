#****************************************************************************
#Composite outcome
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_compo <- readRDS("iso_results/spline_compo.rds")   
iso_compo <- readRDS("iso_results/iso_compo.rds")

#run and save plot
isoplot_compo <- iso_fun_glmer(df_inf_compo$hb, df_inf_compo$compo_pre_lbw_sga, "Composite", spline_compo, iso_compo)

png(file = "iso_results/isoplot_compo.png")
isoplot_compo <- iso_fun_glmer(df_inf_compo$hb, df_inf_compo$compo_pre_lbw_sga, "Composite", spline_compo, iso_compo)
dev.off()

#run and save output data
out_compo <- outdata(iso_compo, "Composite outcome")
out_compo
save(out_compo, file = "iso_results/out_compo.rda")

#****************************************************************************
#*******trim1
#read in model
spline_compo_trim1 <- readRDS("iso_results/spline_compo_trim1.rds")   
iso_compo_trim1 <- readRDS("iso_results/iso_compo_trim1.rds")

#run and save plot
isoplot_compo_trim1 <- iso_fun_glmer(df_inf_compo_trim1$hb, df_inf_compo_trim1$compo_pre_lbw_sga, "Composite - Trim1", spline_compo_trim1, iso_compo_trim1)

png(file = "iso_results/isoplot_compo_trim1.png")
isoplot_compo_trim1 <- iso_fun_glmer(df_inf_compo_trim1$hb, df_inf_compo_trim1$compo_pre_lbw_sga, "Composite - Trim1", spline_compo_trim1, iso_compo_trim1)
dev.off()

#run and save output data
out_compo_trim1 <- outdata(iso_compo_trim1, "Composite_trim1")
out_compo_trim1
save(out_compo_trim1, file = "iso_results/out_compo_trim1.rda")

#****************************************************************************
#*******trim2
#read in model
spline_compo_trim2 <- readRDS("iso_results/spline_compo_trim2.rds")   
iso_compo_trim2 <- readRDS("iso_results/iso_compo_trim2.rds")

#run and save plot
isoplot_compo_trim2 <- iso_fun_glmer(df_inf_compo_trim2$hb, df_inf_compo_trim2$compo_pre_lbw_sga, "Composite - trim2", spline_compo_trim2, iso_compo_trim2)

png(file = "iso_results/isoplot_compo_trim2.png")
isoplot_compo_trim2 <- iso_fun_glmer(df_inf_compo_trim2$hb, df_inf_compo_trim2$compo_pre_lbw_sga, "Composite - trim2", spline_compo_trim2, iso_compo_trim2)
dev.off()

#run and save output data
out_compo_trim2 <- outdata(iso_compo_trim2, "Composite_trim2")
out_compo_trim2
save(out_compo_trim2, file = "iso_results/out_compo_trim2.rda")

#****************************************************************************
#*******trim3
#read in model
spline_compo_trim3 <- readRDS("iso_results/spline_compo_trim3.rds")   
iso_compo_trim3 <- readRDS("iso_results/iso_compo_trim3.rds")

#run and save plot
isoplot_compo_trim3 <- iso_fun_glmer(df_inf_compo_trim3$hb, df_inf_compo_trim3$compo_pre_lbw_sga, "Composite - trim3", spline_compo_trim3, iso_compo_trim3)

png(file = "iso_results/isoplot_compo_trim3.png")
isoplot_compo_trim3 <- iso_fun_glmer(df_inf_compo_trim3$hb, df_inf_compo_trim3$compo_pre_lbw_sga, "Composite - trim3", spline_compo_trim3, iso_compo_trim3)
dev.off()

#run and save output data
out_compo_trim3 <- outdata(iso_compo_trim3, "Composite_trim3")
out_compo_trim3
save(out_compo_trim3, file = "iso_results/out_compo_trim3.rda")
