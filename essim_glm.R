
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

essim_glm <- function(phy, trait, a = 0.5, nsim = 1000, is) {
  
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
  
  trait <- trait[phy$tip.label]
  is <- is[phy$tip.label]
  
  # Generalized Linear Models (GLM's) of correlation between splits statistic and trait
  res <- glm(is ~ trait, family = gaussian, method = "glm.fit")
  
  # Fit Brownian motion model to get diffusion rate and root state estimates
  vv <- vcv.phylo(as.phylo(phy))
  onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
  root <- as.vector(solve(t(onev) %*% solve(vv) %*% onev) %*% (t(onev) %*% solve(vv) %*% trait))
  rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root)) / length(trait))
  
  # Brownian simulations 
  sims <- t(rmvnorm(nsim, sigma = rate * vv))
  rownames(sims) <- rownames(vv)
  
  # GLM's correlations of simulated datasets
  sim.res <- list()
  for(i in 1:nsim){
    sim.res[[i]] <- glm(is ~ sims[, i], family = gaussian, method = "glm.fit")
  }
  
  sim.s <- list() # slope
  for(i in 1:nsim){
    sim.s[[i]] <- sim.res[[i]]$coefficients[2] # slope trait
  }
  
  intercept <- res$coefficients[1]
  slope <- res$coefficients[2]
  r2 <- 1 - (res$deviance/res$null.deviance)
  aic <- res$aic[1]
  
  # Calculate the two-tailed p value (slope)
  upper <- length(sim.s[sim.s >= slope]) / nsim
  lower <- length(sim.s[sim.s <= slope]) / nsim
  pval.s <- 2 * min(c(upper,lower)) # Remove "2" for one-tailed
  
  result <- as.vector(c(intercept, slope, pval.s, r2, aic))
  names(result) <- c("Intercept", "Slope", "P value (Slope)", "r^2", "AIC")
  return(result)
  
}
