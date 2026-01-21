#*******************************************************************************
# The code is adapted from the R code developed by Lai.
# Lai Y. On the adaptive partition approach to the detection of multiple change-points. 
# PLoS One. 2011 May; 6(5). https://home.gwu.edu/~ylai/research/flexstepreg_bi/Rcode.txt

# The modified code is for binary variable with/out fixed effect and/or random_effect

#The code is modified by Xiaoyan Hu, who can be reached by xyh@gwu.edu
#This code was also modified by Precious Williams, williams_pj@gwu.edu
#This change was made to use hb data from 5 to 18 and get CI using two methods 
#*******************************************************************************
library(rms)
library(boot)


knot_fun_boot <- function(data, hb_var, outcome_var) {
  
  n_boot <- 1000
  
  # Make datadist globally accessible for rms
  dd <<- datadist(data)
  options(datadist = "dd")
  
  # Helper: AIC for different numbers of knots
  getAIC <- function(k) {
    f <- as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", k, ")"))
    model <- Glm(f, family = binomial(link = "logit"), data = data)
    AIC(model)
  }
  
  # Select best number of knots based on AIC
  knot_range <- 3:5
  aic_values <- sapply(knot_range, getAIC)
  best_nk <- knot_range[which.min(aic_values)]
  
  # Fit final model with selected knots
  final_model <- Glm(
    formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", best_nk, ")")),
    family = binomial(link = "logit"),
    data = data,
    x = TRUE, y = TRUE
  )
  
  # Bootstrap function
  boot_rcs <- function(data, indices) {
    d <- data[indices, ]
    model <- Glm(
      formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", best_nk, ")")),
      family = binomial(link = "logit"),
      data = d
    )
    hb_seq <- seq(
      from = min(data[[hb_var]], na.rm = TRUE),
      to = max(data[[hb_var]], na.rm = TRUE),
      length.out = 200
    )
    pred <- predict(model, newdata = data.frame(hb = hb_seq), type = "lp")
    return(pred)
  }
  
  # Run bootstrap
  set.seed(123)
  boot_results <- boot(data = data, statistic = boot_rcs, R = n_boot)
  
  # Create same Hb sequence for confidence interval output
  hb_seq <- seq(
    from = min(data[[hb_var]], na.rm = TRUE),
    to = max(data[[hb_var]], na.rm = TRUE),
    length.out = 200
  )
  
  # Confidence intervals from bootstrap results
  boot_ci <- t(apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  
  return(list(
    model = final_model,
    boot_ci = data.frame(
      x = hb_seq,
      lower = boot_ci[, 1],
      upper = boot_ci[, 2]
    )
  ))
}
