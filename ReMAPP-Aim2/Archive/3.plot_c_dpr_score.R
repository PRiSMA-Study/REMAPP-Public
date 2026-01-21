#****************************************************************************
#EPDS score
#****************************************************************************
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/plot_continuous.R")
source("~/REMAPP-Public/ReMAPP-Aim2/iso_code/outdata.R")


load("derived_data/df_dpr_anc.rda")
load("derived_data/df_dpr_pnc.rda")

#load("derived_data/df_mat_dpr.rda")
# load("derived_data/df_mat_dpr_trim1.rda")
# load("derived_data/df_mat_dpr_trim2.rda")
# load("derived_data/df_mat_dpr_trim3.rda")


#****************************************************************************
#*******at ANC
#read in model
spline_dpr_score <- readRDS("iso_results/spline_dpr_score.rds")   
iso_dpr_score <- readRDS("iso_results/iso_dpr_score.rds")

#run and save plot
isoplot_dpr_score <- iso_fun_lmer(df_dpr_anc$hb, df_dpr_anc$dpr_score, "EPDS score at ANC", spline_dpr_score, iso_dpr_score)

png(file = "iso_results/isoplot_dpr_score.png")
isoplot_dpr_score <- iso_fun_lmer(df_dpr_anc$hb, df_dpr_anc$dpr_score, "EPDS score at ANC", spline_dpr_score, iso_dpr_score)
dev.off()

#run and save output data
out_dpr_score <- outdata(iso_dpr_score, "EPDS score at ANC")
out_dpr_score
save(out_dpr_score, file = "iso_results/out_dpr_score.rda")



spline_dpr_score_pnc <- readRDS("iso_results/spline_dpr_score_pnc.rds")   
iso_dpr_score_pnc <- readRDS("iso_results/iso_dpr_score_pnc.rds")

#run and save plot
isoplot_dpr_score_pnc <- iso_fun_lmer(df_dpr_pnc$hb, df_dpr_pnc$dpr_score, "EPDS score at PNC", spline_dpr_score_pnc, iso_dpr_score_pnc)

png(file = "iso_results/isoplot_dpr_score.png")
isoplot_dpr_score_pnc <- iso_fun_lmer(df_dpr_pnc$hb, df_dpr_pnc$dpr_score, "EPDS score at PNC", spline_dpr_score_pnc, iso_dpr_score_pnc)
dev.off()

#run and save output data
out_dpr_score_pnc <- outdata(iso_dpr_score_pnc, "EPDS score at PNC")
out_dpr_score_pnc
save(out_dpr_score_pnc, file = "iso_results/out_dpr_score_pnc.rda")
