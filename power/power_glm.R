
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

#### libraries ####

require(ape)
require(mvtnorm)
require(caper)
require(diversitree)

#### load functions ####

source("essim_glm.R")

#### Power test (50 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees50sdd <- read.tree("trees_Core/trees50sdd_0.004.txt")
states50sdd <- read.table("traits_Core/trees50sdd_0.004_states.txt", header = TRUE, row.names = 1)

essim.p.slope <- vector()

for (i in 1:length(trees50sdd)) {
  
  print(i)
  trait <- states50sdd[, i]
  names(trait) <- row.names(states50sdd)
  
  # ES-sim
  
  essim.res <- essim_glm(trees50sdd[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.50 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees50sdd)

write.csv(essim.pow.50, "power.50.glm.csv")

#### Power test (100 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees100sdd <- read.tree("trees_Core/trees100sdd_0.004.txt")
states100sdd <- read.table("traits_Core/trees100sdd_0.004_states.txt", header = TRUE, row.names = 1)

essim.p.slope <- vector()

for (i in 1:length(trees100sdd)) {
  
  print(i)
  trait <- states100sdd[, i]
  names(trait) <- row.names(states100sdd)
  
  # ES-sim
  
  essim.res <- essim_glm(trees100sdd[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.100 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees100sdd)

write.csv(essim.pow.100, "power.100.glm.csv")

#### Power test (250 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees250sdd <- read.tree("trees_Core/trees250sdd_0.004.txt")
states250sdd <- read.table("traits_Core/trees250sdd_0.004_states.txt", header = TRUE, row.names = 1)

essim.p.slope <- vector()

for (i in 1:length(trees250sdd)) {
  
  print(i)
  trait <- states250sdd[, i]
  names(trait) <- row.names(states250sdd)
  
  # ES-sim
  
  essim.res <- essim_glm(trees250sdd[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.250 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees250sdd)

write.csv(essim.pow.250, "power.250.glm.csv")

#### Power test (1250 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees1250sdd <- read.tree("trees_Core/trees1250sdd_0.004.txt")
states1250sdd <- read.table("traits_Core/trees1250sdd_0.004_states.txt", header = TRUE, row.names = 1)

essim.p.slope <- vector()

for (i in 1:length(trees1250sdd)) {
  
  print(i)
  trait <- states1250sdd[, i]
  names(trait) <- row.names(states1250sdd)
  
  # ES-sim
  
  essim.res <- essim_glm(trees1250sdd[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.1250 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees1250sdd)

write.csv(essim.pow.1250, "power.1250.glm.csv")
