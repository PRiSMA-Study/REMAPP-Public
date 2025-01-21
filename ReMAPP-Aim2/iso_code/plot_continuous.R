#*****************************************************************************
#Plots for continuous outcome
#Author: Xiaoyan Hu
#Email; xyh@gwu.edu
#*****************************************************************************
iso_fun_lmer <- function(x, y, ylab, spline_model, iso_model) {
  o <- order(x)
  plot(x, y, pch = 19, col = "#95A3A6", cex = 0.3, xlab = "Hemoglobin (g/dL)", ylab = ylab,
       main = ylab)
  axis(side = 1, at = seq(4, 20, 1))  
  # Add spline
  xu <- seq(min(x), max(x), length=100)
  xu_pred <- Predict(spline_model,  hb = xu)
  lines(xu, xu_pred$yhat, col="#016795")
  # Add isotonic regression line
  lines(x[o], iso_model$estimates[o], col = "#9f2305")
  # Add legend
  legend(x = max(x) - 4.5, y = max(y) - 0.1,
         lty = c(1, 1, NA),
         pch = c(NA, NA, 20),
         col = c("#9f2305", "#016795", "#95A3A6"),
         legend = c("isotonic", "spline", "group data"))
}

##notitle function is for panel plots - remove legend, title
iso_fun_lmer_notitle <- function(x, y, spline_model, iso_model) {
  o <- order(x)
  plot(x, y, pch = 19, col = "#95A3A6", cex = 0.3, xlab = "", ylab = "",
       main = "")
  axis(side = 1, at = seq(4, 20, 1)) 
  # Add spline
  xu <- seq(min(x), max(x), length=100)
  xu_pred <- Predict(spline_model,  hb = xu)
  lines(xu, xu_pred$yhat, col="#016795")
  # Add isotonic regression line
  lines(x[o], iso_model$estimates[o], col = "#9f2305")
}