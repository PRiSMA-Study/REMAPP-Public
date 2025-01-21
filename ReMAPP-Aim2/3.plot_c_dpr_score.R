#****************************************************************************
#EPDS score
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_dpr_score <- readRDS("iso_results/spline_dpr_score.rds")   
iso_dpr_score <- readRDS("iso_results/iso_dpr_score.rds")

#run and save plot
isoplot_dpr_score <- iso_fun_lmer(df_mat_dpr$hb, df_mat_dpr$dpr_score, "EPDS score", spline_dpr_score, iso_dpr_score)

png(file = "iso_results/isoplot_dpr_score.png")
isoplot_dpr_score <- iso_fun_lmer(df_mat_dpr$hb, df_mat_dpr$dpr_score, "EPDS score", spline_dpr_score, iso_dpr_score)
dev.off()

#run and save output data
out_dpr_score <- outdata(iso_dpr_score, "EPDS score")
out_dpr_score
save(out_dpr_score, file = "iso_results/out_dpr_score.rda")

