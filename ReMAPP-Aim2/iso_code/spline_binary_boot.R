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
  
  # Create datadist and assign it globally (temporarily)
  # rms package needs to find it in the search path
  dd <- datadist(data)
  assign("dd", dd, envir = .GlobalEnv)
  options(datadist = "dd")
  
  # Helper: AIC for different numbers of knots
  getAIC <- function(k) {
    f <- as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", k, ")"))
    model <- Glm(f, family = binomial(link = "logit"), data = data)
    return(AIC(model))
  }
  
  # Select best number of knots based on AIC
  knot_range <- 3:5
  aic_values <- sapply(knot_range, getAIC)
  best_nk <- knot_range[which.min(aic_values)]
  
  # Fit final model with selected knots
  final_model <- Glm(
    formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", best_nk, ")")),
    family = binomial(link = "logit"),
    data = data
  )
  
  # Create Hb sequence for predictions
  hb_seq <- seq(
    from = min(data[[hb_var]], na.rm = TRUE),
    to = max(data[[hb_var]], na.rm = TRUE),
    length.out = 200
  )
  
  # Create prediction data frame with correct variable name
  pred_data <- data.frame(hb_seq)
  names(pred_data) <- hb_var
  
  # Bootstrap function
  boot_rcs <- function(data, indices) {
    d <- data[indices, ]
    
    # Create datadist for bootstrap sample
    dd_boot <- datadist(d)
    assign("dd_boot", dd_boot, envir = .GlobalEnv)
    old_options <- options(datadist = "dd_boot")
    on.exit({
      options(old_options)
      # Clean up temporary dd_boot
      if (exists("dd_boot", envir = .GlobalEnv)) {
        rm("dd_boot", envir = .GlobalEnv)
      }
    })
    
    # Fit model on bootstrap sample
    model <- tryCatch({
      Glm(
        formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", best_nk, ")")),
        family = binomial(link = "logit"),
        data = d
      )
    }, error = function(e) {
      return(NULL)
    })
    
    # If model fails, return NAs
    if (is.null(model) || !model$converged) {
      return(rep(NA, length(hb_seq)))
    }
    
    # Predict on the original sequence
    pred <- predict(model, newdata = pred_data, type = "lp")
    return(pred)
  }
  
  # Run bootstrap
  set.seed(123)
  boot_results <- boot(data = data, statistic = boot_rcs, R = n_boot)
  
  # Confidence intervals from bootstrap results
  boot_ci <- t(apply(boot_results$t, 2, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) {
      quantile(x_clean, probs = c(0.025, 0.975), na.rm = TRUE)
    } else {
      c(NA, NA)
    }
  }))
  
  # Clean up global dd after we're done
  if (exists("dd", envir = .GlobalEnv)) {
    rm("dd", envir = .GlobalEnv)
  }
  options(datadist = NULL)
  
  # Return same structure as original
  return(list(
    model = final_model,
    boot_ci = data.frame(
      x = hb_seq,
      lower = boot_ci[, 1],
      upper = boot_ci[, 2]
    )
  ))
}