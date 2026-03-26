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
  
  # Create datadist for the original data
  dd <- datadist(data)
  options(datadist = dd)
  
  ## Function to return AIC values of the models
  getAIC <- function(nk, model, rm_var, add_var_str){
    new_var <- sprintf(add_var_str, nk) 
    update_frml <- as.formula(sprintf(".~.-%s+%s", rm_var, new_var)) 
    res <- (update(model, update_frml))$aic 
    names(res) <- new_var
    return(res)
  }
  
  ## Start with a model having 3 knots
  base_model <- Glm(eval(parse(text=paste0(outcome_var, " ~ rcs(", hb_var, ",", 3,")"))),
                    family=gaussian(link="identity"), data = data)
  
  ## Compare the AIC of models having 3,4,5 knots
  aic_values <- sapply(3:5, getAIC, model = base_model,
                       rm_var = paste0("rcs(", hb_var, ", 3)"), 
                       add_var_str = paste0("rcs(", hb_var, ", %d)")) 
  
  # Select the best number of knots
  best_knots <- as.numeric(gsub("rcs\\(.*?, (\\d+)\\)", "\\1", 
                                names(aic_values)[which.min(aic_values)]))
  
  # Fit final model with selected knots
  final_model <- Glm(eval(parse(text=paste0(outcome_var, " ~ rcs(", hb_var, ",", best_knots,")"))),
                     family=gaussian(link="identity"), data = data)
  
  # Create prediction sequence from original data
  hb_seq <- seq(
    from = min(data[[hb_var]], na.rm = TRUE),
    to = max(data[[hb_var]], na.rm = TRUE),
    length.out = 200
  )
  
  # Bootstrap function - must be completely self-contained
  boot_rcs <- function(data, indices) {
    # Extract bootstrap sample
    d <- data[indices, ]
    
    # Create new datadist for bootstrap sample
    dd_boot <- datadist(d)
    old_options <- options(datadist = dd_boot)
    on.exit(options(old_options))
    
    # Fit model with same formula as final model
    model_boot <- Glm(
      formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", best_knots, ")")),
      family = gaussian(link = "identity"),
      data = d
    )
    
    # Predict on the original hb_seq (important!)
    pred_data <- data.frame(x = hb_seq)
    names(pred_data) <- hb_var
    
    pred <- predict(model_boot, newdata = pred_data, type = "lp")
    return(pred)
  }
  
  # Run bootstrap
  set.seed(123)
  boot_results <- boot(data = data, statistic = boot_rcs, R = n_boot)
  
  # Get predictions from the final model on the same sequence
  pred_data_final <- data.frame(x = hb_seq)
  names(pred_data_final) <- hb_var
  final_predictions <- predict(final_model, newdata = pred_data_final, type = "lp")
  
  # Confidence intervals from bootstrap results
  boot_ci <- t(apply(boot_results$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  
  # Return results
  return(list(
    model = final_model,
    predictions = data.frame(
      x = hb_seq,
      y = final_predictions
    ),
    boot_ci = data.frame(
      x = hb_seq,
      lower = boot_ci[, 1],
      upper = boot_ci[, 2]
    ),
    best_knots = best_knots,
    boot_results = boot_results
  ))
}