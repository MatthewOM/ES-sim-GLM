
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

essim_glm_multiple <- function(phy, trait1, trait2, a = 0.5, nsim = 1000, is) {
  
  if(missing(is)) { # If inverse splits statistics not provided, calculate it
    rootnode <- length(phy$tip.label) + 1
    is <- numeric(length(phy$tip.label))
    for (i in 1:length(is)){
      node <- i
      index <- 1
      qx <- 0
      while (node != rootnode){
        el <- phy$edge.length[phy$edge[, 2] == node]
        node <- phy$edge[, 1][phy$edge[, 2] == node]
        qx <- qx + el * (1 / (1 / a) ^ (index - 1))		
        index <- index + 1
      }
      is[i] <- 1 / qx
    }		
    names(is) <- phy$tip.label
    is <- log(is[phy$tip.label]) # log transform
  }
  
  trait1 <- trait1[phy$tip.label]
  trait2 <- trait2[phy$tip.label]
  is <- is[phy$tip.label]
  
  # Generalized Linear Models (GLM's) of correlation between splits statistic and trait
  res <- glm(is ~ trait1 + trait2, family = gaussian, method = "glm.fit")
  
  # Fit Brownian motion model to get diffusion rate and root state estimates
  vv <- vcv.phylo(as.phylo(phy))
  onev <- matrix(rep(1, length(trait1)), nrow = length(trait1), ncol = 1)
  
  root1 <- as.vector(solve(t(onev) %*% solve(vv) %*% onev) %*% (t(onev) %*% solve(vv) %*% trait1))
  rate1 <- as.vector((t(trait1 - root1) %*% solve(vv) %*% (trait1 - root1)) / length(trait1))
  
  root2 <- as.vector(solve(t(onev) %*% solve(vv) %*% onev) %*% (t(onev) %*% solve(vv) %*% trait2))
  rate2 <- as.vector((t(trait2 - root2) %*% solve(vv) %*% (trait2 - root2)) / length(trait2))
  
  # Brownian simulations 
  sims1 <- t(rmvnorm(nsim, sigma = rate1 * vv))
  sims2 <- t(rmvnorm(nsim, sigma = rate2 * vv))
  sims <- cbind(sims1, sims2)
  rownames(sims) <- rownames(vv)
  
  # GLM's correlations of simulated datasets
  sim.res <- list()
  for(i in 1:nsim){
    for(j in (nsim + 1):(nsim * 2)){
      if ((i + nsim) == j){
        sim.res[[i]] <- glm(is ~ sims[, i] + sims[, j], family = gaussian, method = "glm.fit")
      }
    }
  }
  
  sim.s1 <- list() # slope trait1
  sim.s2 <- list() # slope trait2
  for(i in 1:nsim){
    sim.s1[[i]] <- sim.res[[i]]$coefficients[2] # slope trait1
    sim.s2[[i]] <- sim.res[[i]]$coefficients[3] # slope trait2
  }
  
  intercept <- res$coefficients[1]
  slope1 <- res$coefficients[2]
  slope2 <- res$coefficients[3]
  r2 <- 1 - (res$deviance/res$null.deviance)
  aic <- res$aic[1]
  
  # Calculate the two-tailed p value (slope trait1)
  upper <- length(sim.s1[sim.s1 >= slope1]) / nsim
  lower <- length(sim.s1[sim.s1 <= slope1]) / nsim
  pval.s1 <- 2 * min(c(upper,lower)) # Remove "2" for one-tailed
  
  # Calculate the two-tailed p value (slope trait2)
  upper <- length(sim.s2[sim.s2 >= slope2]) / nsim
  lower <- length(sim.s2[sim.s2 <= slope2]) / nsim
  pval.s2 <- 2 * min(c(upper,lower)) # Remove "2" for one-tailed
  
  result <- as.vector(c(intercept, slope1, pval.s1, slope2, pval.s2, r2, aic))
  names(result) <- c("Intercept", "Slope1", "P value (Slope1)", "Slope2", "P value (Slope2)", "r^2", "AIC")
  return(result)
  
}
