#*******************************************************************************
# The code is adapted from the R code developed by Lai.
# Lai Y. On the adaptive partition approach to the detection of multiple change-points. 
# PLoS One. 2011 May; 6(5). https://home.gwu.edu/~ylai/research/flexstepreg_bi/Rcode.txt

# The modified code is for binary variable with/out fixed effect and/or random_effect

#The code is modified by Xiaoyan Hu, who can be reached by xyh@gwu.edu
#This code was also modified by Precious Williams, 
#This change was made to use hb data from 5 to 18 
#*******************************************************************************

library(rms)
## Test 2: Force the data to be between 5 to 18 ----
  knot_fun <- function(data, hb_var, outcome_var) {
    
    # 3. Apply limits
    hb_min <- 5
    hb_max <- 18
    
    # 4. Trim dataset
    data_restricted <- data[data[[hb_var]] >= hb_min & data[[hb_var]] <= hb_max, ]
    
    # 5. Setup datadist
    dd <- datadist(data_restricted)
    
    options(datadist = "dd")
    
    # 6. AIC helper
    getAIC <- function(nk, model, rm_var, add_var_str) {
      new_var <- sprintf(add_var_str, nk)
      update_frml <- as.formula(sprintf(". ~ . - %s + %s", rm_var, new_var))
      res <- (update(model, update_frml))$aic
      names(res) <- new_var
      return(res)
    }
    
    # 7. Base model
    base_model <- Glm(
      formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", 3)")),
      family = binomial(link = "logit"),
      data = data_restricted
    )
    
    # 8. AIC comparison
    aic_values <- sapply(3:5, getAIC, model = base_model,
                         rm_var = paste0("rcs(", hb_var, ", 3)"),
                         add_var_str = paste0("rcs(", hb_var, ", %d)"))
    
    best_nk <- as.numeric(gsub("rcs\\(.*?, (\\d+)\\)", "\\1", names(aic_values)[which.min(aic_values)]))
    
    # 9. Final model
    final_model <- Glm(
      formula = as.formula(paste0(outcome_var, " ~ rcs(", hb_var, ", ", best_nk, ")")),
      family = binomial(link = "logit"),
      data = data_restricted,
      x = TRUE, y = TRUE
    )
    
    # Return both model and CI data
    return(final_model)
  }