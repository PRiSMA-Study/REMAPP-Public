#****************************************************************************
#read in model
#****************************************************************************
spline_dpr_score <- readRDS("iso_results/spline_dpr_score.rds")   
iso_dpr_score <- readRDS("iso_results/iso_dpr_score.rds")
# load("derived_data/df_mat_dpr.rda")
#****************************************************************************
#Plot function
#****************************************************************************
iso_fun_lmer <- function(x, y, ylab, spline_model, iso_model) {
  o <- order(x)  
  plot(x, y, pch = 19, col = "#95A3A6", cex = 0.3, xlab = "Hemoglobin (g/dL)", ylab = ylab, 
       main = ylab)
  axis(side = 1, at = seq(4, 20, 1))  # Add x-axis ticks
  # Add spline
  xu <- seq(min(x), max(x), length=100)
  xu_pred <- Predict(spline_model,  hb = xu)
  lines(xu, xu_pred$yhat, col="#016795")
  # Add isotonic regression line
  lines(x[o], iso_model$estimates[o], col = "#9f2305")
  # Add legend
  legend(x = max(x) - 3, y = max(y) - 0.1, 
         # cex = 0.6, 
         lty = c(1, 1, NA), 
         pch = c(NA, NA, 20),
         col = c("#9f2305", "#016795", "#95A3A6"), 
         legend = c("isotonic", "spline", "group data"))
}

#****************************************************************************
#Run plots
#****************************************************************************
#check plots in plot window
isoplot_dpr_score <- iso_fun_lmer(df_mat_dpr$hb, df_mat_dpr$dpr_score, "EPDS score", spline_dpr_score, iso_dpr_score)

#save plots
png(file = "iso_results/isoplot_dpr_score.png")
isoplot_dpr_score <- iso_fun_lmer(df_mat_dpr$hb, df_mat_dpr$dpr_score, "EPDS score", spline_dpr_score, iso_dpr_score)
dev.off()

#****************************************************************************
#Output
#****************************************************************************
df_mean <- as.data.frame(iso_dpr_score$groups) %>% 
  bind_cols(as.data.frame(iso_dpr_score$estimates)) %>% 
  group_by(`iso_dpr_score$groups`) %>%  
  summarise(
    group_n = n(),
    mean_score = sum(`iso_dpr_score$estimates`)/n())

out_dpr_score <- as.data.frame(
  list(
  Outcome = "EPDS score",
  N_group1 = df_mean$group_n[1],
  "Risk/Score1" = df_mean$mean_score[1],
  Threshold1 = iso_dpr_score$brkPoints[1],
  N_group2 = df_mean$group_n[2],
  "Risk/Score2" = df_mean$mean_score[2],
  Threshold2 = iso_dpr_score$brkPoints[2],
  N_group3 = df_mean$group_n[3],
  "Risk/Score3" = df_mean$mean_score[3],
  Threshold3 = iso_dpr_score$brkPoints[3],
  N_group4 = df_mean$group_n[4],
  "Risk/Score4" = df_mean$mean_score[4],
  Threshold4 = iso_dpr_score$brkPoints[4],
  N_group5 = df_mean$group_n[5],
  "Risk/Score5" = df_mean$mean_score[5],
  Threshold5 = iso_dpr_score$brkPoints[5],
  N_group6 = df_mean$group_n[6],
  "Risk/Score6" = df_mean$mean_score[6],
  Threshold6 = iso_dpr_score$brkPoints[6],
  N_group7 = df_mean$group_n[7],
  "Risk/Score7" = df_mean$mean_score[7]
    ))
out_dpr_score

save(out_dpr_score, file = "iso_results/out_dpr_score.rda")