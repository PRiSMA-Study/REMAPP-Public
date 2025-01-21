#****************************************************************************
#Fatigue score
#****************************************************************************
source("iso_code/plot_binary.R")
source("iso_code/outdata.R")

#****************************************************************************
#*******all
#read in model
spline_ftg_score <- readRDS("iso_results/spline_ftg_score.rds")   
iso_ftg_score <- readRDS("iso_results/iso_ftg_score.rds")

#run and save plot
isoplot_ftg_score <- iso_fun_lmer(df_mat_ftg$hb, df_mat_ftg$ftg_score, "Fatigue score", spline_ftg_score, iso_ftg_score)

png(file = "iso_results/isoplot_ftg_score.png")
isoplot_ftg_score <- iso_fun_lmer(df_mat_ftg$hb, df_mat_ftg$ftg_score, "Fatigue score", spline_ftg_score, iso_ftg_score)
dev.off()

#run and save output data
out_ftg_score <- outdata(iso_ftg_score, "Fatigue score")
out_ftg_score
save(out_ftg_score, file = "iso_results/out_ftg_score.rda")


