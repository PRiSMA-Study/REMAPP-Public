#*****************************************************************************
#Plots for binary outcome
#Author: Xiaoyan Hu
#Email; xyh@gwu.edu
#*****************************************************************************
iso_fun_glmer <- function(x, y, ylab, spline_model, iso_model) {
  o <- order(x)  
  plot(x, y, pch = 18, col = "#95A3A6", cex = 0.3, xlab = "Hemoglobin (g/dL)", ylab = ylab, 
       main = ylab)
  # Add spline
  xu <- seq(min(x), max(x), length=30)
  xu_pred <- Predict(spline_model,  hb = xu)
  lines(xu, exp(xu_pred$yhat), col="#016795")
  # Add isotonic regression line
  lines(x[o], iso_model$estimates[o], col = "#9f2305")
  # Add legend
  legend(x = max(x) - 4, y = max(y) - 0.01, 
         lty = c(1, 1, NA), 
         pch = c(NA, NA, 18),
         col = c("#9f2305", "#016795", "#95A3A6"), 
         legend = c("Isotonic", "Spline", "Data"))
}

#notitle function is for panel plots 
iso_fun_glmer_notitle <- function(x, y, spline_model, iso_model) {
  o <- order(x)  
  plot(x, y, pch = 18, col = "#95A3A6", cex = 0.3, xlab = "", ylab = "", 
       main = "")
  # Add spline
  xu <- seq(min(x), max(x), length=30)
  xu_pred <- Predict(spline_model,  hb = xu)
  lines(xu, exp(xu_pred$yhat), col="#016795")
  # Add isotonic regression line
  lines(x[o], iso_model$estimates[o], col = "#9f2305")
}
