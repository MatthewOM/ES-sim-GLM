
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

trees1250sdd <- read.tree("trees_Core/trees1250sdd_0.004.txt")
states1250sdd <- read.table("traits_Core/trees1250sdd_0.004_states.txt", header = TRUE, row.names = 1)

## create subset data from the full tree to see how that impacts the power analysis with essim.R and essim_pgls.R

trees.sample.50 <- vector("list", 100) # Vector to place trees
class(trees.sample.50) <- "multiPhylo"
states.sample.50 <- vector("list", 100)

for (i in 1:100) {
  dtips <- sample(trees1250sdd[[i]]$tip.label, 1200)
  trees.sample.50[[i]] <- drop.tip(trees1250sdd[[i]], dtips)
  states.sample.50[[i]] <- as.matrix(states1250sdd[ !(rownames(states1250sdd) %in% dtips), i])
  rownames(states.sample.50[[i]]) <- rownames(states1250sdd[ !(rownames(states1250sdd) %in% dtips),])
}

essim.p.slope <- vector()

for (i in 1:length(trees.sample.50)) {
  
  print(i)
  trait <- states.sample.50[[i]]
  names(trait) <- row.names(states.sample.50[[i]])
  
  # ES-sim
  
  essim.res <- essim_glm(trees.sample.50[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.50 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees.sample.50)

write.csv(essim.pow.50, "power.50.glm.subset.csv")

#### Power test (100 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees1250sdd <- read.tree("trees_Core/trees1250sdd_0.004.txt")
states1250sdd <- read.table("traits_Core/trees1250sdd_0.004_states.txt", header = TRUE, row.names = 1)

## create subset data from the full tree to see how that impacts the power analysis with essim.R and essim_pgls.R

trees.sample.100 <- vector("list", 100) # Vector to place trees
class(trees.sample.100) <- "multiPhylo"
states.sample.100 <- vector("list", 100)

for (i in 1:100) {
  dtips <- sample(trees1250sdd[[i]]$tip.label, 1200)
  trees.sample.100[[i]] <- drop.tip(trees1250sdd[[i]], dtips)
  states.sample.100[[i]] <- as.matrix(states1250sdd[ !(rownames(states1250sdd) %in% dtips), i])
  rownames(states.sample.100[[i]]) <- rownames(states1250sdd[ !(rownames(states1250sdd) %in% dtips),])
}

essim.p.slope <- vector()

for (i in 1:length(trees.sample.100)) {
  
  print(i)
  trait <- states.sample.100[[i]]
  names(trait) <- row.names(states.sample.100[[i]])
  
  # ES-sim
  
  essim.res <- essim_glm(trees.sample.100[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.100 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees.sample.100)

write.csv(essim.pow.100, "power.100.glm.subset.csv")

#### Power test (250 tips) ####

## power will be tested on the simulations with trait-dependence, SSD (state-dependent diversification)

trees1250sdd <- read.tree("trees_Core/trees1250sdd_0.004.txt")
states1250sdd <- read.table("traits_Core/trees1250sdd_0.004_states.txt", header = TRUE, row.names = 1)

## create subset data from the full tree to see how that impacts the power analysis with essim.R and essim_pgls.R

trees.sample.250 <- vector("list", 100) # Vector to place trees
class(trees.sample.250) <- "multiPhylo"
states.sample.250 <- vector("list", 100)

for (i in 1:100) {
  dtips <- sample(trees1250sdd[[i]]$tip.label, 1000)
  trees.sample.250[[i]] <- drop.tip(trees1250sdd[[i]], dtips)
  states.sample.250[[i]] <- as.matrix(states1250sdd[ !(rownames(states1250sdd) %in% dtips), i])
  rownames(states.sample.250[[i]]) <- rownames(states1250sdd[ !(rownames(states1250sdd) %in% dtips),])
}

essim.p.slope <- vector()

for (i in 1:length(trees.sample.250)) {
  
  print(i)
  trait <- states.sample.250[[i]]
  names(trait) <- row.names(states.sample.250[[i]])
  
  # ES-sim
  
  essim.res <- essim_glm(trees.sample.250[[i]], trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.pow.250 <- length(essim.p.slope[essim.p.slope < 0.05])/length(trees.sample.250)

write.csv(essim.pow.250, "power.250.glm.subset.csv")
