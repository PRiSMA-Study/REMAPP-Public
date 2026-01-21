# #****************************************************************************
# #Fatigue score
# #****************************************************************************
# source("iso_code/plot_continuous.R")
# source("iso_code/outdata.R")
# 
# load("derived_data/df_mat_ftg.rda")
# 
# #****************************************************************************
# #*******all
# #read in model
# spline_ftg_score <- readRDS("iso_results/spline_ftg_score.rds")   
# iso_ftg_score <- readRDS("iso_results/iso_ftg_score.rds")
# 
# #run and save plot
# isoplot_ftg_score <- iso_fun_lmer(df_mat_ftg$hb, df_mat_ftg$ftg_score, "Fatigue score", spline_ftg_score, iso_ftg_score)
# 
# png(file = "iso_results/isoplot_ftg_score.png")
# isoplot_ftg_score <- iso_fun_lmer(df_mat_ftg$hb, df_mat_ftg$ftg_score, "Fatigue score", spline_ftg_score, iso_ftg_score)
# dev.off()
# 
# #run and save output data
# out_ftg_score <- outdata(iso_ftg_score, "Fatigue score")
# out_ftg_score
# save(out_ftg_score, file = "iso_results/out_ftg_score.rda")
# 

source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_continuous.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")


load("derived_data/df_ftg_anc.rda")
load("derived_data/df_ftg_pnc.rda")

#load("derived_data/df_mat_dpr.rda")
# load("derived_data/df_mat_dpr_trim1.rda")
# load("derived_data/df_mat_dpr_trim2.rda")
# load("derived_data/df_mat_dpr_trim3.rda")


#****************************************************************************
#*******at ANC
#read in model
spline_ftg_score <- readRDS("iso_results/spline_ftg_score.rds")   
iso_ftg_score <- readRDS("iso_results/iso_ftg_score.rds")

#run and save plot
isoplot_ftg_score <- iso_fun_lmer(df_ftg_anc$hb, df_ftg_anc$ftg_score, "Fatigue score at ANC", spline_ftg_score, iso_ftg_score)

png(file = "iso_results/isoplot_ftg_score.png")
isoplot_ftg_score <- iso_fun_lmer(df_ftg_anc$hb, df_ftg_anc$ftg_score, "Fatigue score at ANC", spline_ftg_score, iso_ftg_score)
dev.off()

#run and save output data
out_ftg_score <- outdata(iso_ftg_score, "Fatigue Score score at ANC")
out_ftg_score
save(out_ftg_score, file = "iso_results/out_ftg_score.rda")



spline_ftg_score_pnc <- readRDS("iso_results/spline_ftg_score_pnc.rds")   
iso_ftg_score_pnc <- readRDS("iso_results/iso_ftg_score_pnc.rds")

#run and save plot
isoplot_ftg_score_pnc <- iso_fun_lmer(df_ftg_pnc$hb, df_ftg_pnc$ftg_score, "Fatigue score at PNC", spline_ftg_score_pnc, iso_ftg_score_pnc)

png(file = "iso_results/isoplot_ftg_score.png")
isoplot_ftg_score_pnc <- iso_fun_lmer(df_ftg_pnc$hb, df_ftg_pnc$ftg_score, "Fatigue score at PNC", spline_ftg_score_pnc, iso_ftg_score_pnc)
dev.off()

#run and save output 
out_ftg_score_pnc <- outdata(iso_ftg_score_pnc, "Fatigue score at PNC")
out_ftg_score_pnc
save(out_ftg_score_pnc, file = "iso_results/out_ftg_score_pnc.rda")
