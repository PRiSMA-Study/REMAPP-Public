#########################################################

# The code is adapted from the R code developed by Lai.
# Lai Y. On the adaptive partition approach to the detection of multiple change-points. 
# PLoS One. 2011 May; 6(5). https://home.gwu.edu/~ylai/research/FlexStepReg/Rcode.txt

# The modified code is for binary response variable
# (1) sse is changed to negative log-likelihood of binary response variable
# (2) two-sample t-test is changed to two-proportion z-test or Fisher exact test for event counts

#########################################################

# get number of groups of covariates
get.numbers <- function(x){
  
  # vector length of each unique x value
  n <- tapply(x, x, length)
  # an empty vector of double type
  l <- double(length(n)+1)
  # vector length accumulation
  for(i in 1:length(n)){
    l[i+1] <- l[i]+n[i]
  }
  
  return(l)
}

# two-proportion test 
t2p <- function(v, n1, n2){
  
  if(n1>1 && n2>1){
    # first group
    g1 <- v[1:n1]
    m1 <- sum(g1 == 1)
    # second group
    g2 <- v[(n1+1):(n1+n2)]
    m2 <- sum(g2 == 1)
    # chi-square or fisher exact test
    test_tb <- matrix(c(m1, n1-m1,
                        m2, n2-m2), nrow = 2, byrow = TRUE)
    chisq_res <- chisq.test(test_tb)
    # if any cell of the expected frequency table less than 5, use fisher exact test
    if (any(chisq_res$expected < 5)) {
      fisher_res <- fisher.test(test_tb)
      pvalue <- fisher_res$p.value
    } else {
      pvalue <- chisq_res$p.value
    }
    
    if( is.na(pvalue) ){
      return(1)
    }else{
       return(pvalue)
    }
    
  }else{
    return( 1 )
  }
}

# flexible isotonic regression
flexstepreg <- function(y, x, alpha.adjacency=1){
  
  o <- order(x, decreasing=FALSE)
  x <- x[o]
  y <- y[o]
  
  #sum of squares under the null
  mu <- mean(y)
  # ss0 <- sum( (y-mu)^2 )
  loglik0 <- -sum(map_dbl(y, function(i) ifelse(mu==0, NA, log(mu^i*(1-mu)^(1-i)))))
  
  #estimate under the alternative (reduced isotonic)
  g <- get.numbers(x)
  link.rank.score <- vector("list", length(g)-1)
  #print(g)
  
  # mean in the first group
  mu <- mean(y[(g[1]+1):g[2]])
  # sse of first group
  # score <- sum( (y[(g[1]+1):g[2]]-mu)^2 )
  score <- -sum(map_dbl(y[(g[1]+1):g[2]], function(i) ifelse(mu==0, NA, log(mu^i*(1-mu)^(1-i)))))
  # add frist element to the triplet
  link.rank.score[[1]] <- list(l=0,r=0,s=score)
  for(i in 2:(length(g)-1)){
    # mean of the first i groups
    mu <- mean(y[(g[1]+1):g[i+1]])
    # sse of the first i groups
    # score <- sum( (y[(g[1]+1):g[i+1]]-mu)^2 )
    score <- -sum(map_dbl(y[(g[1]+1):g[i+1]], function(i) ifelse(mu==0, NA, log(mu^i*(1-mu)^(1-i)))))
    # second element
    link.rank.score[[i]] <- list(l=0,r=0,s=score)
    for(j in 1:(i-1)){
      {
        # mean of the (j+1)-th group to i-th group
        mu <- mean(y[(g[j+1]+1):g[i+1]])
        flag <- TRUE
        k <- 1
        while(flag & k<=length(link.rank.score[[j]]$l)){
          link <- link.rank.score[[j]]$l[k]
          # if sample t test is significant
          if( t2p(y[(g[link+1]+1):g[i+1]], g[j+1]-g[link+1], g[i+1]-g[j+1]) < alpha.adjacency ){
            flag <- FALSE
            # score <- link.rank.score[[j]]$s[k] + sum( (y[(g[j+1]+1):g[i+1]]-mu)^2 )
            score <- link.rank.score[[j]]$s[k] - sum(map_dbl(y[(g[j+1]+1):g[i+1]], function(i) ifelse(mu==0, NA, log(mu^i*(1-mu)^(1-i)))))
            link.rank.score[[i]]$l <- c(link.rank.score[[i]]$l, j)
            link.rank.score[[i]]$r <- c(link.rank.score[[i]]$r, k)
            link.rank.score[[i]]$s <- c(link.rank.score[[i]]$s, score)
          }
          k <- k+1
        }
      }
    }
    # reorder by score
    o.score <- order(link.rank.score[[i]]$s, decreasing=FALSE)
    link.rank.score[[i]]$l <- link.rank.score[[i]]$l[o.score]
    link.rank.score[[i]]$r <- link.rank.score[[i]]$r[o.score]
    link.rank.score[[i]]$s <- link.rank.score[[i]]$s[o.score]
  }
  #print(link.rank.score)
  loglik1 <- link.rank.score[[length(g)-1]]$s[1]
  # statistic <- (ss0/ss1 - 1)*(length(y)-length(g)+1)/(length(g)-1-1)
  statistic <- -2*(loglik1 - loglik0)
  
  i <- length(g)-1
  partition <- g[i+1]
  rank <- 1
  while(rank>0){
    i.temp <- i
    i <- link.rank.score[[i.temp]]$l[rank]
    rank <- link.rank.score[[i.temp]]$r[rank]
    partition <- c(g[i+1], partition)
  }
  #print(partition)
  
  estimates <- double(length(y))
  for(i in 2:length(partition)){
    estimates[o][(partition[i-1]+1):partition[i]] <- mean( y[(partition[i-1]+1):partition[i]] )
  }
  
  groups <- integer(length(y))
  for(i in 2:length(partition)){
    groups[o][(partition[i-1]+1):partition[i]] <- i-1
  }
  
  return(list(groups=groups, estimates=estimates, statistic=statistic))
}

# get pvalue by permutation
pvalue_permu <- function(y, x, alpha.adjacency, stat, B){
  
  for(i in 1:B){
    stat <- c(stat, flexstepreg(sample(y), x, alpha.adjacency=alpha.adjacency)$statistic)
    # print(i)
  }
  
  # p-value
  # sum(stat >= stat[1])/length(stat)
  stat
}

# get confidence interval by bootstrap
ci_bootstrap <- function(y, x, alpha.adjacency, estimate, B){
  
  bootstrap0 <- function(idx, y, x){
    xs <- x[idx]
    ys <- y[idx]
    est <- rep(NA, length(y))
    est[i] <- flexstepreg(ys, xs, alpha.adjacency=alpha.adjacency)$estimates
    return(est)	
  }

  for(i in 1:B){
    estimate <- cbind(estimate, bootstrap0(sample(c(1:length(y)), replace=T), y, x))
    # print(i)
  }
  
  data.frame(mean = tapply(x,x,mean),
             lower = tapply(estimate[,-1], rep(x,B), quantile, 0.005, na.rm=T),
             upper = tapply(estimate[,-1], rep(x,B), quantile, 0.995, na.rm=T))
}

















