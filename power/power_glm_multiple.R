
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

#### libraries ####

require(ape)
require(mvtnorm)
require(caper)
require(diversitree)

#### load functions ####

source("essim_glm_multiple.R")

#### Power test (50 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees50sdd <- read.tree("trees_Core/trees50sdd_multiple_0.004.txt")
states50sdd <- read.table("traits_Core/trees50sdd_states_multiple_0.004.txt", header = TRUE, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:length(trees50sdd)) {
  
  print(i)
  trait1 <- states50sdd[, ((i*2)-1)]
  names(trait1) <- row.names(states50sdd)
  trait2 <- states50sdd[, (i*2)]
  names(trait2) <- row.names(states50sdd)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees50sdd[[i]], trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec <- vec[vec=="TRUE"]

essim.pow.50.slope <- length(vec)/length(trees50sdd)

write.csv(essim.pow.50, "power/power.50.glm.multiple_0.004.csv")

#### Power test (100 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees100sdd <- read.tree("trees_Core/trees100sdd_multiple_0.004.txt")
states100sdd <- read.table("traits_Core/trees100sdd_states_multiple_0.004.txt", header = TRUE, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:length(trees100sdd)) {
  
  print(i)
  trait1 <- states100sdd[, ((i*2)-1)]
  names(trait1) <- row.names(states100sdd)
  trait2 <- states100sdd[, (i*2)]
  names(trait2) <- row.names(states100sdd)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees100sdd[[i]], trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec <- vec[vec=="TRUE"]

essim.pow.100.slope <- length(vec)/length(trees100sdd)

write.csv(essim.pow.100, "power/power.100.glm.multiple_0.004.csv")

#### Power test (250 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees250sdd <- read.tree("trees_Core/trees250sdd_multiple_0.004.txt")
states250sdd <- read.table("traits_Core/trees250sdd_states_multiple_0.004.txt", header = TRUE, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:length(trees250sdd)) {
  
  print(i)
  trait1 <- states250sdd[, ((i*2)-1)]
  names(trait1) <- row.names(states250sdd)
  trait2 <- states250sdd[, (i*2)]
  names(trait2) <- row.names(states250sdd)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees250sdd[[i]], trait1, trait2, nsim = 2500)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec <- vec[vec=="TRUE"]

essim.pow.250.slope <- length(vec)/length(trees250sdd)

write.csv(essim.pow.250, "power/power.250.glm.multiple_0.004.csv")

#### Power test (1250 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees1250sdd <- read.tree("trees_Core/trees1250sdd_multiple_0.004.txt")
states1250sdd <- read.table("traits_Core/trees1250sdd_states_multiple_0.004.txt", header = TRUE, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:length(trees1250sdd)) {
  
  print(i)
  trait1 <- states1250sdd[, ((i*2)-1)]
  names(trait1) <- row.names(states1250sdd)
  trait2 <- states1250sdd[, (i*2)]
  names(trait2) <- row.names(states1250sdd)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees1250sdd[[i]], trait1, trait2, nsim = 12500)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec <- vec[vec=="TRUE"]

essim.pow.1250.slope <- length(vec)/length(trees1250sdd)

write.csv(essim.pow.1250, "power/power.1250.glm.multiple_0.004.csv")
