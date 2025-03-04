#*******************************************************************************
# The code is adapted from the R code developed by Lai.
# Lai Y. On the adaptive partition approach to the detection of multiple change-points. 
# PLoS One. 2011 May; 6(5). https://home.gwu.edu/~ylai/research/flexstepreg_bi/Rcode.txt

# The modified code is for binary variable with/out fixed effect and/or random_effect

#The code is modified by Xiaoyan Hu, who can be reached by xyh@gwu.edu
#*******************************************************************************
library(rms)

knot_fun <- function(data, hb_var, outcome_var){
  ds <- datadist(data)
  options(datadist = ds)
  ## function to return AIC values of the models
  getAIC <- function(nk, model, rm_var, add_var_str){
    new_var <- sprintf(add_var_str, nk) 
    update_frml <- as.formula(sprintf(".~.-%s+%s", rm_var, new_var)) 
    res <- (update(model, update_frml))$aic 
    names(res) <- new_var
    return(res)
  }
  
  ## start with a model having 3 knots
  base_model <- Glm(eval(parse(text=paste0(outcome_var, " ~ rcs(", hb_var, ",", 3,")"))),
                    family=binomial(link="logit"), data = data)
  
  ## compare the AIC of models having 3,4,5 knots
  aic_values <- sapply(3:5, getAIC, model = base_model,
                       rm_var = paste0("rcs(", hb_var, ", 3)"), 
                       add_var_str = paste0("rcs(", hb_var, ", %d)")) 
  
  #select the best knots
  best_knots <- as.numeric(gsub("rcs\\(.*?, (\\d+)\\)", "\\1", names(aic_values)[which.min(aic_values)]))
  
  #final spine model by using best knots selected
  final_model <- Glm(eval(parse(text=paste0(outcome_var, " ~ rcs(", hb_var, ",", best_knots,")"))),
                     family=binomial(link="logit"), data = data,
                     x = TRUE, y = TRUE)
  
  return(final_model)
}