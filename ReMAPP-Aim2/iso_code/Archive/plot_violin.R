

spline_and_violin_plot <- function(x, y, xlab, ylab, spline_model, iso_model, outcome_var, title) {
  o <- order(x)
  xu <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)
  xu_pred <- Predict(spline_model, hb = xu, conf.int = 0.95)
  
  layout(matrix(1:2, ncol = 1), heights = c(1.2, 2))
  par(mar = c(0.5, 4, 3, 2))  # top plot
  
  # ───── 1. Violin Plot ─────
  outcome_levels <- sort(unique(outcome_var))
  colors <- c("#E69F00", "#56B4E9")  # extend if more classes
  
  plot(NA, xlim = range(x, na.rm = TRUE), ylim = c(0.5, length(outcome_levels) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = title)
  
  for (i in seq_along(outcome_levels)) {
    val <- outcome_levels[i]
    hb_vals <- x[outcome_var == val]
    
    if (length(hb_vals) > 1) {
      d <- density(hb_vals, na.rm = TRUE)
      scale_factor <- 0.4 / max(d$y)
      polygon(d$x, i + d$y * scale_factor, col = adjustcolor(colors[i], alpha.f = 0.4), border = NA)
      polygon(d$x, i - d$y * scale_factor, col = adjustcolor(colors[i], alpha.f = 0.4), border = NA)
      
      # Jittered raw data points
      points(jitter(hb_vals, amount = 0.05), rep(i, length(hb_vals)),
             pch = 16, col = adjustcolor("gray80", alpha.f = 0.6), cex = 0.5)
    }
  }
  
  # Violin legend with outcome labels
  legend("topright", inset = 0.02,
         legend = c(paste(ylab, "= 0"), paste(ylab, "= 1"), "Hb values"),
         fill = c(adjustcolor(colors[1], alpha.f = 0.4),
                  adjustcolor(colors[2], alpha.f = 0.4),
                  "gray80"),
         border = NA,
         pch = c(NA, NA, 16),
         col = c(NA, NA, "black"),
         pt.cex = c(NA, NA, 0.7),
         bty = "n", text.col = "black")
  
  # ───── 2. Spline Plot ─────
  par(mar = c(5, 4, 1, 2))
  y_max <- max(exp(xu_pred$upper), na.rm = TRUE)
  
  plot(NA, xlim = range(x, na.rm = TRUE), ylim = c(0, y_max),
       xlab = "Hemoglobin (g/dL)", ylab = "")
  
  polygon(c(xu, rev(xu)),
          c(exp(xu_pred$lower), rev(exp(xu_pred$upper))),
          col = adjustcolor("#016795", alpha.f = 0.2), border = NA)
  
  lines(xu, exp(xu_pred$yhat), col = "#016795", lwd = 2)
  lines(x[o], iso_model$estimates[o], col = "#9f2305", lwd = 2)
  
  legend("topright",
         legend = c("Spline fit", "Spline 95% CI", "Isotonic fit"),
         col = c("#016795", adjustcolor("#016795", alpha.f = 0.2), "#9f2305"),
         lty = c(1, NA, 1),
         lwd = c(2, NA, 2),
         pch = c(NA, 15, NA),
         pt.cex = c(NA, 2, NA),
         bty = "n")
}

