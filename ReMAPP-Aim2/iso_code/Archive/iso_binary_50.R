library(lme4)
library(lmerTest)
library(tidyverse)
library(glmmTMB)

t2p_glmer_50 <- function(v, n1, n2, tail, covar1=NULL, covar2=NULL, random_effect=NULL) {
  if (n1 < 2 || n2 < 2) {
    return
    (1) 
  }
  
  x <- v[1:n1]
  y <- v[(n1 + 1):(n1 + n2)]
  
  
  # Fit mixed-effects model 
  if(is.null(random_effect))
  {
    if(is.null(covar2)){
      data <- data.frame(
        outcome = c(x, y),
        group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
        covar1 = factor(rep(covar1))
      )
      if(n_distinct(data$outcome)<2){
        return(1)
      }
      if(n_distinct(data$covar1)<2){
        return(1)
      }
      if(sum(data$outcome==0)<5 || sum(data$outcome==1)<5){
        return(1)
      }
      if(sum(data$group=="Group1")<5 || sum(data$group=="Group2")<5){
        return(1)
      }
      model <- glm(outcome ~ group+covar1, family=binomial(link="logit"),data=data)
    }
    else{
      data <- data.frame(
        outcome = c(x, y),
        group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
        covar1 = factor(rep(covar1)),
        covar2 = factor(rep(covar2))
      )
      if(n_distinct(data$outcome)<2){
        return(1)
      }
      if(n_distinct(data$covar1)<2){
        return(1)
      }
      if(n_distinct(data$covar2)<2){
        return(1)
      }
      if(sum(data$outcome==0)<5 || sum(data$outcome==1)<5){
        return(1)
      }
      if(sum(data$group=="Group1")<5 || sum(data$group=="Group2")<5){
        return(1)
      }
      model <- glm(outcome ~ group+covar1+covar2, family=binomial(link="logit"),data=data)
    }
  }
  else{
    if(is.null(covar2)){
      data <- data.frame(
        outcome = c(x, y),
        group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
        covar1 = factor(rep(covar1)),
        random_effect = factor(rep(random_effect))
      ) %>% 
        group_by(random_effect) %>% 
        mutate(n = n()) %>%
        filter(n>1) %>%
        ungroup()
      if(n_distinct(data$outcome)<2){
        return(1)
      }
      if(n_distinct(data$group)<2){
        return(1)
      }
      if(n_distinct(data$covar1)<2){
        return(1)
      }
      if(n_distinct(data$random_effect)<2){
        return(1)
      }
      if(sum(data$outcome==0)<5 || sum(data$outcome==1)<5){
        return(1)
      }
      if(sum(data$group=="Group1")<5 || sum(data$group=="Group2")<5){
        return(1)
      }
      model <- glmmTMB(outcome ~ group+covar1+ (1|random_effect), family=binomial, data=data)
    }
    else{
      data <- data.frame(
        outcome = c(x, y),
        group = factor(rep(c("Group1", "Group2"), c(n1, n2))),
        covar1 = factor(rep(covar1)),
        covar2 = factor(rep(covar2)),
        random_effect = factor(rep(random_effect))
      ) %>% 
        group_by(random_effect) %>% 
        mutate(n = n()) %>%
        filter(n>1) %>%
        ungroup()
      if(n_distinct(data$outcome)<2){
        return(1)
      }
      if(n_distinct(data$group)<2){
        return(1)
      }
      if(n_distinct(data$covar1)<2){
        return(1)
      }
      if(n_distinct(data$covar2)<2){
        return(1)
      }
      if(n_distinct(data$random_effect)<2){
        return(1)
      }
      if(sum(data$outcome==0)<5 || sum(data$outcome==1)<5){
        return(1)
      }
      if(sum(data$group=="Group1")<5 || sum(data$group=="Group2")<5){
        return(1)
      }
      model <- glmmTMB(outcome ~ group+covar1+covar2 + (1|random_effect), family=binomial, data=data)
    }
  }
  
  if(is.null(random_effect)){
    if(n_distinct(data$outcome)<2){
      pvalue<-1
    }
    else{
      summary_coef <- summary(model)$coefficients
      row_index <- grep("^groupGroup2", rownames(summary_coef))
      # Check if row_index is not empty
      if (length(row_index) > 0) {
        tempT <- summary_coef[row_index, "z value"]
        ### "changing "t value" to "z value"
        pvalue<-2*pnorm(q=abs(tempT), lower.tail=FALSE)
      } 
      else {
        pvalue <- 1  
      }
    }
    
  }
  else{
    summary_coef <- summary(model)$coefficients
    row_index <- grep("^groupGroup2", rownames(summary_coef))
    # Check if row_index is not empty
    if (length(row_index) > 0) {
      tempT <- summary_coef[row_index, "z value"]
      ### "changing "t value" to "z value"
      pvalue<-2*pnorm(q=abs(tempT), lower.tail=FALSE)
    } 
    else {
      pvalue <- 1  # Or handle the case where the coefficient name doesn't exist
    }
  }
  
  # Check for NaN (although this is uncommon with lmer())
  if (is.na(pvalue) || is.nan(pvalue)) {
    return(1)
  } else {
    return(pvalue)
  }
}

flexstepreg_glmer <- function(y, x, covar1=NULL, covar2=NULL, random_effect=NULL, alpha.adjacency = 1, tail.two = "upper") {
  
  # y = df_inf_sga3$sga3 
  # x = df_inf_sga3$hb 
  # covar1 = df_inf_sga3$SITE  
  # covar2=NULL 
  # random_effect=NULL 
  # alpha = 0.01
  # alpha.adjacency = 1 
  # tail.two = "upper"
  
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
  
  get_numbers_binned <- function(x, bin_width = 0.5) {
    bins <- floor(x / bin_width) * bin_width  # Group into 0.50 intervals
    bin_counts <- rle(sort(bins))$lengths     # Count values in each bin
    c(0, cumsum(bin_counts))                  # Cumulative positions
  }
  
  # Usage:
  g <- get_numbers_binned(x, bin_width = 0.5)
  
  #g <- get.numbers_glmer(x)
  
  link.rank.score <- vector("list", length(g) - 1)
  
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
        tempPV <-y[(g[link + 1] + 1):g[i + 1]]
        tempN1 <- g[j + 1] - g[link + 1]
        tempN2 <- g[i + 1] - g[j + 1]
        tempCovar1 <- covar1[(g[link + 1] + 1):g[i + 1]]
        if(!is.null(covar2)){
          tempCovar2 <- covar2[(g[link + 1] + 1):g[i + 1]]
        }
        else{
          tempCovar2<-NULL
        }
        
        if(!is.null(random_effect)){
          tempRandom_effect <- random_effect[(g[link + 1] + 1):g[i + 1]]
        }   
        else{
          tempRandom_effect<-NULL
        }
        
        #temproary debugging
        kk<-kk+1
        show(kk)
        if(kk==108){
          kk=kk
        }
        
        if (tempN1 < 50 || tempN2 < 50) {
          flag <- FALSE  # Skip to the next iteration
        } else {
          if (t2p_glmer_50(tempPV, tempN1, tempN2, tail.two, tempCovar1, tempCovar2, tempRandom_effect) < alpha.adjacency) {
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
    o.score <- order(link.rank.score[[i]]$s, decreasing = FALSE)
    link.rank.score[[i]]$l <- link.rank.score[[i]]$l[o.score]
    link.rank.score[[i]]$r <- link.rank.score[[i]]$r[o.score]
    link.rank.score[[i]]$s <- link.rank.score[[i]]$s[o.score]
  }
  
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
  
  #calculate the break points
  brkPoints<-rep(0,length(partition)-2)
  for(i in 2:length(partition)-1){
    brkPoints[i-1]<-(x[partition[i]]+x[partition[i]+1])/2
  }
  
  return(list(groups = groups, estimates = estimates, statistic = statistic, brkPoints = brkPoints))
}

