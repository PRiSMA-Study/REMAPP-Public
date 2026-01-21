#*******************************************************************************
# The code is adapted from the R code developed by Lai.
# Lai Y. On the adaptive partition approach to the detection of multiple change-points. 
# PLoS One. 2011 May; 6(5). https://home.gwu.edu/~ylai/research/flexstepreg_bi/Rcode.txt

# The modified code is for continuous variable with/out fixed effect and/or random_effect

#*******************************************************************************
#library(lme4)
library(lmerTest)
library(dplyr)

get.numbers_lmer <- function(x) {
  n <- tapply(x, x, length)
  l <- double(length(n) + 1)
  for (i in 1:length(n)) {
    l[i + 1] <- l[i] + n[i]
  }
  return(l)
}

t2p_lmer<- function(v, n1, n2, tail, covar1=NULL, covar2=NULL, random_effect=NULL) {
  
  if (n1 < 50 || n2 < 50) {
    return(1)  
  }
  
  x <- v[1:n1]
  y <- v[(n1 + 1):(n1 + n2)]
  
  # Check if we have enough distinct values for all variables
  if (length(unique(y)) < 2) {
    return(1)
  }
  
  tryCatch({
    # Fit mixed-effects model
    if(is.null(random_effect))
    {
      if(is.null(covar2)){
        data <- data.frame(
          outcome = c(x, y),
          group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
          covar1 = factor(covar1)
        )
        if(n_distinct(data$outcome)<2 || n_distinct(data$covar1)<2){
          return(1)
        }
        model <- lmer(outcome ~ group + covar1, data = data, REML = FALSE)
      }
      else{
        data <- data.frame(
          outcome = c(x, y),
          group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
          covar1 = factor(covar1),
          covar2 = factor(covar2)
        )
        if(n_distinct(data$outcome)<2 || n_distinct(data$covar1)<2 || n_distinct(data$covar2)<2){
          return(1)
        }
        model <- lmer(outcome ~ group + covar1 + covar2, data = data, REML = FALSE)
      }
    }
    else{
      if(is.null(covar2)){
        data <- data.frame(
          outcome = c(x, y),
          group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
          covar1 = factor(covar1),
          random_effect = factor(random_effect)
        ) %>% 
          group_by(random_effect) %>% 
          mutate(n = n()) %>%
          filter(n>1) %>%
          ungroup()
        if(n_distinct(data$outcome)<2 || n_distinct(data$group)<2 || 
           n_distinct(data$covar1)<2 || n_distinct(data$random_effect)<2){
          return(1)
        }
        model <- lmer(outcome ~ group + covar1 + (1 | random_effect), data = data, REML = FALSE)
      }
      else{
        data <- data.frame(
          outcome = c(x, y),
          group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
          covar1 = factor(covar1),
          covar2 = factor(covar2),
          random_effect = factor(random_effect)
        ) %>% 
          group_by(random_effect) %>% 
          mutate(n = n()) %>%
          filter(n>1) %>%
          ungroup()
        if(n_distinct(data$outcome)<2 || n_distinct(data$group)<2 || 
           n_distinct(data$covar1)<2 || n_distinct(data$covar2)<2 || 
           n_distinct(data$random_effect)<2){
          return(1)
        }
        model <- lmer(outcome ~ group + covar1 + covar2 + (1 | random_effect), data = data, REML = FALSE)
      }
    }
    
    # Check for convergence issues
    if (is.null(random_effect)) {
      if (n_distinct(data$outcome) < 2) {
        return(1)
      } else {
        summary_coef <- summary(model)$coefficients
        row_index <- grep("^groupGroup2", rownames(summary_coef))
        if (length(row_index) > 0) {
          pvalue <- summary_coef[row_index, "Pr(>|t|)"]
        } else {
          pvalue <- 1
        }
      }
    } else if (isSingular(model)) {
      pvalue <- 1
    } else {
      summary_coef <- summary(model)$coefficients
      row_index <- grep("^groupGroup2", rownames(summary_coef))
      if (length(row_index) > 0) {
        pvalue <- summary_coef[row_index, "Pr(>|t|)"]
      } else {
        pvalue <- 1
      }
    }
    
    # Check for NaN or NA
    if (is.na(pvalue) || is.nan(pvalue)) {
      return(1)
    } else {
      return(pvalue)
    }
    
  }, error = function(e) {
    # Return p-value of 1 if any error occurs during model fitting
    return(1)
  }, warning = function(w) {
    # Return p-value of 1 if any warning occurs during model fitting
    return(1)
  })
}

flexstepreg_lmer <- function(y, x, covar1=NULL, covar2=NULL, random_effect=NULL, alpha.adjacency = 1, tail.two = "upper") {
  kk<-0
  tail.two <- match.arg(tail.two, c("upper", "lower", "two"))
  o <- order(x, decreasing = FALSE)
  
  x <- x[o]
  y <- y[o]
  covar1 <- covar1[o]
  if(!is.null(covar2)){
    covar2 <- covar2[o]
  }
  if(!is.null(random_effect)){
    random_effect <- random_effect[o]
  }
  
  mu <- mean(y)
  ss0 <- sum((y - mu)^2)
  
  # Get binned x values and group indices
  get_numbers_binned <- function(x) {
    # Step 1: Clip to [5, 18]
    x_clipped <- pmax(pmin(x, 18), 5)
    
    # Step 2: Bin to nearest center at .25 or .75 (bin width = 0.5)
    x_binned <- floor((x_clipped - 0.0) / 0.5 + 0.5) * 0.5 + 0.25
    
    # Step 3: Bin positions
    bin_counts <- rle(sort(x_binned))$lengths
    g <- c(0, cumsum(bin_counts))
    
    # Assign to global environment
    list2env(list(g = g, x_binned = x_binned), envir = .GlobalEnv)
  }
  
  
  # Call the function once and store results
  out <- get_numbers_binned(x)
  
  g <- out$g
  
  x_binned <- out$x_binned
  
  # Calculate total number of potential models for progress tracking
  total_models <- 0
  for (i in 2:(length(g) - 1)) {
    for (j in 1:(i - 1)) {
      total_models <- total_models + 1
    }
  }
  
  cat(paste("Total models to run:", total_models, "\n"))
  
  link.rank.score <- vector("list", length(g) - 1)
  
  # Initialize first element
  mu <- mean(y[(g[1] + 1):g[2]])
  score <- sum((y[(g[1] + 1):g[2]] - mu)^2)
  link.rank.score[[1]] <- list(l = 0, r = 0, s = score)
  
  for (i in 2:(length(g) - 1)) {
    mu <- mean(y[(g[1] + 1):g[i + 1]])
    score <- sum((y[(g[1] + 1):g[i + 1]] - mu)^2)
    link.rank.score[[i]] <- list(l = 0, r = 0, s = score)
    
    for (j in 1:(i - 1)) {
      mu <- mean(y[(g[j + 1] + 1):g[i + 1]])
      flag <- TRUE
      k <- 1
      
      while (flag & k <= length(link.rank.score[[j]]$l)) {
        link <- link.rank.score[[j]]$l[k]
        tempPV <- y[(g[link + 1] + 1):g[i + 1]]
        tempN1 <- g[j + 1] - g[link + 1]
        tempN2 <- g[i + 1] - g[j + 1]
        tempCovar1 <- covar1[(g[link + 1] + 1):g[i + 1]]
        
        if(!is.null(covar2)){
          tempCovar2 <- covar2[(g[link + 1] + 1):g[i + 1]]
        } else {
          tempCovar2 <- NULL
        }
        
        if(!is.null(random_effect)){
          tempRandom_effect <- random_effect[(g[link + 1] + 1):g[i + 1]]
        } else {
          tempRandom_effect <- NULL
        }
        
        # Increment and display counter
        kk <- kk + 1
        if (kk %% 100 == 0 || kk == total_models) {
          cat(paste("Running model", kk, "of", total_models, "\n"))
        }
        
        if (tempN1 < 50 || tempN2 < 50) {
          flag <- FALSE
        } else {
          p_val <- t2p_lmer(tempPV, tempN1, tempN2, tail.two, tempCovar1, tempCovar2, tempRandom_effect)
          if (p_val < alpha.adjacency) {
            flag <- FALSE
            score <- link.rank.score[[j]]$s[k] + sum((y[(g[j + 1] + 1):g[i + 1]] - mu)^2)
            link.rank.score[[i]]$l <- c(link.rank.score[[i]]$l, j)
            link.rank.score[[i]]$r <- c(link.rank.score[[i]]$r, k)
            link.rank.score[[i]]$s <- c(link.rank.score[[i]]$s, score)
          }
        }
        k <- k + 1
      }
    }
    
    # Sort by score
    o.score <- order(link.rank.score[[i]]$s, decreasing = FALSE)
    link.rank.score[[i]]$l <- link.rank.score[[i]]$l[o.score]
    link.rank.score[[i]]$r <- link.rank.score[[i]]$r[o.score]
    link.rank.score[[i]]$s <- link.rank.score[[i]]$s[o.score]
  }
  
  cat(paste("Finished! Total models run:", kk, "\n"))
  
  ss1 <- link.rank.score[[length(g) - 1]]$s[1]
  statistic <- (ss0 / ss1 - 1) * (length(y) - length(g) + 1) / (length(g) - 1 - 1)
  
  i <- length(g) - 1
  partition <- g[i + 1]
  rank <- 1
  while (rank > 0) {
    i.temp <- i
    i <- link.rank.score[[i.temp]]$l[rank]
    rank <- link.rank.score[[i.temp]]$r[rank]
    partition <- c(g[i + 1], partition)
  }
  
  estimates <- double(length(y))
  for (i in 2:length(partition)) {
    estimates[o][(partition[i - 1] + 1):partition[i]] <- mean(y[(partition[i - 1] + 1):partition[i]])
  }
  
  groups <- integer(length(y))
  for (i in 2:length(partition)) {
    groups[o][(partition[i - 1] + 1):partition[i]] <- i - 1
  }
  
  # Calculate breakpoints only 
  if (length(partition) > 1) {
    # Use pre-binned x values (already at .125/.375/.625/.875)
    brkPoints <- sapply(2:(length(partition)-1), function(i) {
      (x_binned[partition[i]] + x_binned[partition[i]+1]) / 2
    })
    
  } else {
    brkPoints <- numeric(0)
  }
  
  return(list(groups = groups, estimates = estimates, statistic = statistic, brkPoints = brkPoints))
}


