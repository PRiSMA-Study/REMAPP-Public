#****************************************************************************
#read in model
#****************************************************************************
spline_pprom <- readRDS("iso_results/spline_pprom.rds")   
iso_pprom <- readRDS("iso_results/iso_pprom.rds")
# load("derived_data/df_inf_pprom.rda")
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
isoplot_pprom <- iso_fun_glmer(df_inf_pprom$hb, df_inf_pprom$inf_pprom, "PPROM", spline_pprom, iso_pprom)

#save plots
png(file = "iso_results/isoplot_pprom.png")
isoplot_pprom <- iso_fun_glmer(df_inf_pprom$hb, df_inf_pprom$inf_pprom, "PPROM", spline_pprom, iso_pprom)
dev.off()

#****************************************************************************
#Output
#****************************************************************************
df_risk <- as.data.frame(iso_pprom$groups) %>% 
  bind_cols(as.data.frame(iso_pprom$estimates)) %>% 
  group_by(`iso_pprom$groups`) %>%  
  summarise(
    group_n = n(),
    risk = sum(`iso_pprom$estimates`)/n())

out_pprom <- as.data.frame(
  list(
    Outcome = "PRPROM",
    N_group1 = df_risk$group_n[1],
    "Risk/Score1" = df_risk$risk[1],
    Threshold1 = iso_pprom$brkPoints[1],
    N_group2 = df_risk$group_n[2],
    "Risk/Score2" = df_risk$risk[2],
    Threshold2 = iso_pprom$brkPoints[2],
    N_group3 = df_risk$group_n[3],
    "Risk/Score3" = df_risk$risk[3],
    Threshold3 = iso_pprom$brkPoints[3],
    N_group4 = df_risk$group_n[4],
    "Risk/Score4"= df_risk$risk[4],
    Threshold4 = iso_pprom$brkPoints[4],
    N_group5 = df_risk$group_n[5],
    "Risk/Score5" = df_risk$risk[5],
    Threshold5 = iso_pprom$brkPoints[5],
    N_group6 = df_risk$group_n[6],
    "Risk/Score6" = df_risk$risk[6],
    Threshold6 = iso_pprom$brkPoints[6],
    N_group7 = df_risk$group_n[7],
    "Risk/Score7" = df_risk$risk[7]
  ))
out_pprom

save(out_pprom, file = "iso_results/out_pprom.rda")