#*******************************************************************************
# The code is adapted from the R code developed by Lai.
# Lai Y. On the adaptive partition approach to the detection of multiple change-points. 
# PLoS One. 2011 May; 6(5). https://home.gwu.edu/~ylai/research/flexstepreg_bi/Rcode.txt

# The modified code is for binary variable with/out fixed effect and/or random_effect

#The code is modified by Xiaoyan Hu, who can be reached by xyh@gwu.edu
#This code was also modified by Precious Williams, 
#This change was made to force the knots at 10th percentile and 90th percentile 
#*******************************************************************************
library(rms)
## Test 1:Force the knots to be at the 10th and 90th percentile ----
    
    knot_fun <- function(data, hb_var, outcome_var) {
      # Get the hemoglobin values from the data
      hb_vals <- data[[hb_var]]

      # Calculate 10th and 90th percentiles (excluding NAs)
      percentiles <- quantile(hb_vals, probs = c(0.10, 0.90), na.rm = TRUE)
      manual_knots <- as.numeric(percentiles)

      # Set up datadist (important for rms functions)
      dd <- datadist(data)
      options(datadist = "dd")

      # Build formula with manual knots
      # Note: rcs() with 2 knots actually creates 1 internal knot plus boundary knots
      formula_str <- paste0(outcome_var, " ~ rcs(", hb_var, ", knots = c(",
                            paste(manual_knots, collapse = ", "), "))")

      # Fit the logistic regression model
      final_model <- lrm(as.formula(formula_str),
                         data = data,
                         x = TRUE,
                         y = TRUE)

      return(final_model)
    }


