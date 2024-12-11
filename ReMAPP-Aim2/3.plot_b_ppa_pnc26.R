#****************************************************************************
#read in model
#****************************************************************************
spline_ppa_pnc26 <- readRDS("iso_results/spline_ppa_pnc26.rds")   
iso_ppa_pnc26 <- readRDS("iso_results/iso_ppa_pnc26.rds")
# load("derived_data/df_mat_ppa_pnc26.rda")
#****************************************************************************
#Plot function
#****************************************************************************
iso_fun_glmer <- function(x, y, ylab, spline_model, iso_model) {
  o <- order(x)  # Order predictor variable
  plot(x, y, pch = 19, col = "#95A3A6", cex = 0.3, xlab = "Hemoglobin (g/dL)", ylab = ylab, 
       main = ylab)
  axis(side = 1, at = seq(4, 20, 1))  # Add x-axis ticks
  # Add spline
  xu <- seq(min(x), max(x), length=30)
  xu_pred <- Predict(spline_model,  hb = xu)
  lines(xu, exp(xu_pred$yhat), col="#016795")
  # Add isotonic regression line
  lines(x[o], iso_model$estimates[o], col = "#9f2305")
  # Add legend
  legend(x = max(x) - 5, y = max(y) - 0.1, 
         # cex = 0.6, 
         lty = c(1, 1, NA), 
         pch = c(NA, NA, 20),
         col = c("#9f2305", "#016795", "#95A3A6"), 
         legend = c("Isotonic", "Spline", "Data"))
  return(iso_model)
}
#****************************************************************************
#Run plots
#****************************************************************************
#check plots in plot window
isoplot_ppa_pnc26 <- iso_fun_glmer(df_mat_ppa_pnc26$hb, df_mat_ppa_pnc26$ppa_pnc26, "Postpartum anemia at PNC26", spline_ppa_pnc26, iso_ppa_pnc26)

#save plots
png(file = "iso_results/isoplot_ppa_pnc26.png")
isoplot_ppa_pnc26 <- iso_fun_glmer(df_mat_ppa_pnc26$hb, df_mat_ppa_pnc26$ppa_pnc26, "Postpartum anemia at PNC26", spline_ppa_pnc26, iso_ppa_pnc26)
dev.off()

#****************************************************************************
#Output
#****************************************************************************
df_risk <- as.data.frame(iso_ppa_pnc26$groups) %>% 
  bind_cols(as.data.frame(iso_ppa_pnc26$estimates)) %>% 
  group_by(`iso_ppa_pnc26$groups`) %>%  
  summarise(
    group_n = n(),
    risk = sum(`iso_ppa_pnc26$estimates`)/n())

out_ppa_pnc26 <- as.data.frame(
  list(
    Outcome = "Postpartum anemia at PNC26",
    N_group1 = df_risk$group_n[1],
    "Risk/Score1" = df_risk$risk[1],
    Threshold1 = iso_ppa_pnc26$brkPoints[1],
    N_group2 = df_risk$group_n[2],
    "Risk/Score2" = df_risk$risk[2],
    Threshold2 = iso_ppa_pnc26$brkPoints[2],
    N_group3 = df_risk$group_n[3],
    "Risk/Score3" = df_risk$risk[3],
    Threshold3 = iso_ppa_pnc26$brkPoints[3],
    N_group4 = df_risk$group_n[4],
    "Risk/Score4" = df_risk$risk[4],
    Threshold4 = iso_ppa_pnc26$brkPoints[4],
    N_group5 = df_risk$group_n[5],
    "Risk/Score5" = df_risk$risk[5],
    Threshold5 = iso_ppa_pnc26$brkPoints[5],
    N_group6 = df_risk$group_n[6],
    "Risk/Score6" = df_risk$risk[6],
    Threshold6 = iso_ppa_pnc26$brkPoints[6],
    N_group7 = df_risk$group_n[7],
    "Risk/Score7" = df_risk$risk[7],
    Threshold7 = iso_ppa_pnc26$brkPoints[7],
    N_group8 = df_risk$group_n[8],
    "Risk/Score8" = df_risk$risk[8],
    Threshold8 = iso_ppa_pnc26$brkPoints[8],
    N_group9 = df_risk$group_n[9],
    "Risk/Score9" = df_risk$risk[9]
  ))
out_ppa_pnc26

save(out_ppa_pnc26, file = "iso_results/out_ppa_pnc26.rda")