plot_boot_violin <- function(x, y, xlab, ylab, rcs_result, iso_model, outcome_var, title) {
  o <- order(x)
  boot_ci <- rcs_result$boot_ci
  spline_model <- rcs_result$model
  
  xu <- boot_ci$x
  lower <- plogis(boot_ci$lower)
  upper <- plogis(boot_ci$upper)
  
  # Get predicted spline values on probability scale
  pred_main <- plogis(predict(spline_model, newdata = data.frame(hb = xu)))
  
  layout(matrix(1:2, ncol = 1), heights = c(1.2, 2))
  par(mar = c(0.5, 4, 3, 2))  # top plot
  
  # ───── 1. Violin Plot ─────
  # Ensure outcome_var is binary factor with levels 0 and 1
  outcome_levels <- rev(sort(unique(outcome_var)))
  outcome_factor <- factor(outcome_var, levels = outcome_levels)
  
  # Define consistent colors for outcome levels (0 = orange, 1 = blue)
  colors <- c("#E69F00", "#56B4E9")  # Level 0 → orange, Level 1 → blue
  names(colors) <- outcome_levels
  
  plot(NA, xlim = range(xu, na.rm = TRUE), ylim = c(0.5, length(outcome_levels) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = title)
  
  for (i in seq_along(outcome_levels)) {
    val <- outcome_levels[i]
    hb_vals <- x[outcome_var == val]
    
    if (length(hb_vals) > 1) {
      d <- density(hb_vals, na.rm = TRUE)
      scale_factor <- 0.4 / max(d$y)
      polygon(d$x, i + d$y * scale_factor, col = adjustcolor(colors[i], alpha.f = 0.4), border = NA)
      polygon(d$x, i - d$y * scale_factor, col = adjustcolor(colors[i], alpha.f = 0.4), border = NA)
      
      points(jitter(hb_vals, amount = 0.05), rep(i, length(hb_vals)),
             pch = 16, col = adjustcolor("gray80", alpha.f = 0.6), cex = 0.5)
    }
  }
  
  legend("topright", inset = 0.02,
         legend = c( paste(ylab, "=", outcome_levels[2]),
                    paste(ylab, "=", outcome_levels[1]),
                    "Hb values"),
         fill = c(adjustcolor(colors[as.character(outcome_levels[2])], alpha.f = 0.4),
                  adjustcolor(colors[as.character(outcome_levels[1])], alpha.f = 0.4),
                  
                  "gray80"),
         border = NA,
         pch = c(NA, NA, 16),
         col = c(NA, NA, "black"),
         pt.cex = c(NA, NA, 0.7),
         cex = 0.7,
         bty = "n", text.col = "black")
  
  # ───── 2. Spline Plot ─────
  par(mar = c(5, 4, 1, 2))
  y_max <- max(upper, na.rm = TRUE)
  
  plot(NA, xlim = range(xu, na.rm = TRUE), ylim = c(0, 1),
       xlab = "Hemoglobin (g/dL)", ylab = "")
  
  # Bootstrap CI ribbon
  polygon(c(xu, rev(xu)),
          c(lower, rev(upper)),
          col = adjustcolor("#016795", alpha.f = 0.2), border = NA)
  
  # Spline fit line
  lines(xu, pred_main, col = "#016795", lwd = 2)
  
  # Isotonic line
  lines(x[o], iso_model$estimates[o], col = "#9f2305", lwd = 2)
  
  # Legend
  legend("topright",
         legend = c("Spline fit", "Bootstrap 95% CI", "Isotonic fit"),
         col = c("#016795", adjustcolor("#016795", alpha.f = 0.2), "#9f2305"),
         lty = c(1, NA, 1),
         lwd = c(2, NA, 2),
         pch = c(NA, 15, NA),
         pt.cex = c(NA, 2, NA),
         cex = 0.8,  # <-- smaller text
         bty = "n")
}

