plot_boot_continuous <- function(x, y, xlab, ylab, rcs_result, iso_model, outcome_var, title) {
  o <- order(x)
  boot_ci <- rcs_result$boot_ci
  spline_model <- rcs_result$model
  
  xu <- boot_ci$x
  lower <- boot_ci$lower
  upper <- boot_ci$upper
  
  # Get predicted spline values on original scale
  pred_main <- predict(spline_model, newdata = data.frame(hb = xu))
  
  # Set up single plot
  par(mar = c(5, 4, 4, 2))
  
  # Calculate appropriate y-limits that include both data and fitted curves
  y_range_data <- range(y, na.rm = TRUE)
  y_range_fit <- range(c(lower, upper, pred_main, iso_model$estimates), na.rm = TRUE)
  y_combined_range <- range(c(y_range_data, y_range_fit), na.rm = TRUE)
  y_padding <- diff(y_combined_range) * 0.1
  y_lim <- c(y_combined_range[1] - y_padding, y_combined_range[2] + y_padding)
  
  # Create main plot area
  plot(NA, xlim = range(xu, na.rm = TRUE), ylim = y_lim,
       xlab = xlab, ylab = ylab, main = title)
  
  # Add gray background data points (semi-transparent)
  points(x, y, pch = 16, col = adjustcolor("gray70", alpha.f = 0.3), cex = 0.8)
  
  # Bootstrap CI ribbon
  polygon(c(xu, rev(xu)),
          c(lower, rev(upper)),
          col = adjustcolor("#016795", alpha.f = 0.2), border = NA)
  
  # Spline fit line
  lines(xu, pred_main, col = "#016795", lwd = 3)
  
  # Isotonic line
  lines(x[o], iso_model$estimates[o], col = "#9f2305", lwd = 3)
  
  # Enhanced legend
  legend("topright",
         legend = c("Spline fit", "Bootstrap 95% CI", "Isotonic fit", "Data points"),
         col = c("#016795", adjustcolor("#016795", alpha.f = 0.2), "#9f2305", 
                 adjustcolor("gray70", alpha.f = 0.3)),
         lty = c(1, NA, 1, NA),
         lwd = c(3, NA, 3, NA),
         pch = c(NA, 15, NA, 16),
         pt.cex = c(NA, 2, NA, 0.8),
         cex = 0.8,
         bty = "n",
         bg = adjustcolor("white", alpha.f = 0.8))
}