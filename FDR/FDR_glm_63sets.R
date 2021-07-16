
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

#### libraries ####

require(ape)
require(mvtnorm)
require(caper)
require(diversitree)

#### load functions ####

source("essim_glm.R")

#### FDR test (bisse) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.bisse <- read.tree("trees_FDR/trees_bisse.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.bisse <- list()

for(i in 1:50){
  # load .txt frequency tables
  states.bisse[[i]] <- read.table(states.list[[i]], header = T, row.names = 1)
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.1cladefixed, "FDR/fdr.bisse.1cladefixed.csv")

## 1cladenosignal ##

states.bisse <- list()

for(i in 51:100){
  for (j in 1:50) {
    if((i - 50) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.1cladenosignal, "FDR/fdr.bisse.1cladenosignal.csv")

## BM ##

states.bisse <- list()

for(i in 101:150){
  for (j in 1:50) {
    if((i - 100) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.BM <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.BM, "FDR/fdr.bisse.BM.csv")

## BMjump ##

states.bisse <- list()

for(i in 151:200){
  for (j in 1:50) {
    if((i - 150) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.BMjump, "FDR/fdr.bisse.BMjump.csv")

## BMmultirate ##

states.bisse <- list()

for(i in 201:250){
  for (j in 1:50) {
    if((i - 200) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.BMmultirate, "FDR/fdr.bisse.BMmultirate.csv")

## discretetrait ##

states.bisse <- list()

for(i in 251:300){
  for (j in 1:50) {
    if((i - 250) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.discretetrait, "FDR/fdr.bisse.discretetrait.csv")

## nosignal ##

states.bisse <- list()

for(i in 301:350){
  for (j in 1:50) {
    if((i - 300) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.nosignal, "FDR/fdr.bisse.nosignal.csv")

## OUstrong ##

states.bisse <- list()

for(i in 351:400){
  for (j in 1:50) {
    if((i - 350) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.OUstrong, "FDR/fdr.bisse.OUstrong.csv")

## OUweak ##

states.bisse <- list()

for(i in 401:450){
  for (j in 1:50) {
    if((i - 400) == j) {
      # load .txt frequency tables
      states.bisse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.bisse)) {
  
  print(i)
  states <- states.bisse[[i]]
  
  for (j in 1:length(trees.bisse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.bisse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.bisse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.bisse.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.bisse.OUweak, "FDR/fdr.bisse.OUweak.csv")

#### FDR test (carnivore) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.carnivore <- read.tree("trees_FDR/trees_carnivore.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[451]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.1cladefixed, "FDR/fdr.carnivore.1cladefixed.csv")

## 1cladenosignal ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[452]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.1cladenosignal, "FDR/fdr.carnivore.1cladenosignal.csv")

## BM ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[453]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.BM <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.BM, "FDR/fdr.carnivore.BM.csv")

## BMjump ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[454]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.BMjump, "FDR/fdr.carnivore.BMjump.csv")

## BMmultirate ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[455]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.BMmultirate, "FDR/fdr.carnivore.BMmultirate.csv")

## discretetrait ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[456]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.discretetrait, "FDR/fdr.carnivore.discretetrait.csv")

## nosignal ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[457]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.nosignal, "FDR/fdr.carnivore.nosignal.csv")

## OUstrong ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[458]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.OUstrong, "FDR/fdr.carnivore.OUstrong.csv")

## OUweak ##

states.carnivore <- list()

# load .txt frequency tables
states.carnivore <- read.table(states.list[[459]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.carnivore)) {
  
  print(i)
  trait <- states.carnivore[, i]
  names(trait) <- row.names(states.carnivore)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.carnivore, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.carnivore.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.carnivore)

write.csv(essim.fdr.carnivore.OUweak, "FDR/fdr.carnivore.OUweak.csv")

#### FDR test (constant) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.constant <- read.tree("trees_FDR/trees_constant.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.constant <- list()

for(i in 460:509){
  for (j in 1:50) {
    if((i - 459) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.1cladefixed, "FDR/fdr.constant.1cladefixed.csv")

## 1cladenosignal ##

states.constant <- list()

for(i in 510:559){
  for (j in 1:50) {
    if((i - 509) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.1cladenosignal, "FDR/fdr.constant.1cladenosignal.csv")

## BM ##

states.constant <- list()

for(i in 560:609){
  for (j in 1:50) {
    if((i - 559) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.BM <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.BM, "FDR/fdr.constant.BM.csv")

## BMjump ##

states.constant <- list()

for(i in 610:659){
  for (j in 1:50) {
    if((i - 609) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.BMjump, "FDR/fdr.constant.BMjump.csv")

## BMmultirate ##

states.constant <- list()

for(i in 660:709){
  for (j in 1:50) {
    if((i - 659) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.BMmultirate, "FDR/fdr.constant.BMmultirate.csv")

## discretetrait ##

states.constant <- list()

for(i in 710:759){
  for (j in 1:50) {
    if((i - 709) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.discretetrait, "FDR/fdr.constant.discretetrait.csv")

## nosignal ##

states.constant <- list()

for(i in 760:809){
  for (j in 1:50) {
    if((i - 759) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.nosignal, "FDR/fdr.constant.nosignal.csv")

## OUstrong ##

states.constant <- list()

for(i in 810:859){
  for (j in 1:50) {
    if((i - 809) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.OUstrong, "FDR/fdr.constant.OUstrong.csv")

## OUweak ##

states.constant <- list()

for(i in 860:909){
  for (j in 1:50) {
    if((i - 859) == j) {
      # load .txt frequency tables
      states.constant[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.constant)) {
  
  print(i)
  states <- states.constant[[i]]
  
  for (j in 1:length(trees.constant)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.constant[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.constant[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.constant.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.constant.OUweak, "FDR/fdr.constant.OUweak.csv")

#### FDR test (coral) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.coral <- read.tree("trees_FDR/trees_coral.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.coral <- list()

for(i in 910:959){
  for (j in 1:50) {
    if((i - 909) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.1cladefixed, "FDR/fdr.coral.1cladefixed.csv")

## 1cladenosignal ##

states.coral <- list()

for(i in 960:1009){
  for (j in 1:50) {
    if((i - 959) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.1cladenosignal, "FDR/fdr.coral.1cladenosignal.csv")

## BM ##

states.coral <- list()

for(i in 1010:1059){
  for (j in 1:50) {
    if((i - 1009) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.BM <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.BM, "FDR/fdr.coral.BM.csv")

## BMjump ##

states.coral <- list()

for(i in 1060:1109){
  for (j in 1:50) {
    if((i - 1059) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.BMjump, "FDR/fdr.coral.BMjump.csv")

## BMmultirate ##

states.coral <- list()

for(i in 1110:1159){
  for (j in 1:50) {
    if((i - 1109) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.BMmultirate, "FDR/fdr.coral.BMmultirate.csv")

## discretetrait ##

states.coral <- list()

for(i in 1160:1209){
  for (j in 1:50) {
    if((i - 1159) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.discretetrait, "FDR/fdr.coral.discretetrait.csv")

## nosignal ##

states.coral <- list()

for(i in 1210:1259){
  for (j in 1:50) {
    if((i - 1209) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.nosignal, "FDR/fdr.coral.nosignal.csv")

## OUstrong ##

states.coral <- list()

for(i in 1260:1309){
  for (j in 1:50) {
    if((i - 1259) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.OUstrong, "FDR/fdr.coral.OUstrong.csv")

## OUweak ##

states.coral <- list()

for(i in 1310:1359){
  for (j in 1:50) {
    if((i - 1309) == j) {
      # load .txt frequency tables
      states.coral[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.coral)) {
  
  print(i)
  states <- states.coral[[i]]
  
  for (j in 1:length(trees.coral)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.coral[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.coral[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.coral.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.coral.OUweak, "FDR/fdr.coral.OUweak.csv")

#### FDR test (ddk4) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.ddk4 <- read.tree("trees_FDR/trees_ddk4.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.ddk4 <- list()

for(i in 1360:1409){
  for (j in 1:50) {
    if((i - 1359) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.1cladefixed, "FDR/fdr.ddk4.1cladefixed.csv")

## 1cladenosignal ##

states.ddk4 <- list()

for(i in 1410:1459){
  for (j in 1:50) {
    if((i - 1409) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.1cladenosignal, "FDR/fdr.ddk4.1cladenosignal.csv")

## BM ##

states.ddk4 <- list()

for(i in 1460:1509){
  for (j in 1:50) {
    if((i - 1459) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.BM <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.BM, "FDR/fdr.ddk4.BM.csv")

## BMjump ##

states.ddk4 <- list()

for(i in 1510:1559){
  for (j in 1:50) {
    if((i - 1509) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.BMjump, "FDR/fdr.ddk4.BMjump.csv")

## BMmultirate ##

states.ddk4 <- list()

for(i in 1560:1609){
  for (j in 1:50) {
    if((i - 1559) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.BMmultirate, "FDR/fdr.ddk4.BMmultirate.csv")

## discretetrait ##

states.ddk4 <- list()

for(i in 1610:1659){
  for (j in 1:50) {
    if((i - 1609) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.discretetrait, "FDR/fdr.ddk4.discretetrait.csv")

## nosignal ##

states.ddk4 <- list()

for(i in 1660:1709){
  for (j in 1:50) {
    if((i - 1659) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.nosignal, "FDR/fdr.ddk4.nosignal.csv")

## OUstrong ##

states.ddk4 <- list()

for(i in 1710:1759){
  for (j in 1:50) {
    if((i - 1709) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.OUstrong, "FDR/fdr.ddk4.OUstrong.csv")

## OUweak ##

states.ddk4 <- list()

for(i in 1760:1809){
  for (j in 1:50) {
    if((i - 1759) == j) {
      # load .txt frequency tables
      states.ddk4[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.ddk4)) {
  
  print(i)
  states <- states.ddk4[[i]]
  
  for (j in 1:length(trees.ddk4)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.ddk4[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.ddk4[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.ddk4.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.ddk4.OUweak, "FDR/fdr.ddk4.OUweak.csv")

#### FDR test (quasse) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.quasse <- read.tree("trees_FDR/trees_quasse.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.quasse <- list()

for(i in 1810:1859){
  for (j in 1:50) {
    if((i - 1809) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.1cladefixed, "FDR/fdr.quasse.1cladefixed.csv")

## 1cladenosignal ##

states.quasse <- list()

for(i in 1860:1909){
  for (j in 1:50) {
    if((i - 1859) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.1cladenosignal, "FDR/fdr.quasse.1cladenosignal.csv")

## BM ##

states.quasse <- list()

for(i in 1910:1959){
  for (j in 1:50) {
    if((i - 1909) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.BM <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.BM, "FDR/fdr.quasse.BM.csv")

## BMjump ##

states.quasse <- list()

for(i in 1960:2009){
  for (j in 1:50) {
    if((i - 1959) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.BMjump, "FDR/fdr.quasse.BMjump.csv")

## BMmultirate ##

states.quasse <- list()

for(i in 2010:2059){
  for (j in 1:50) {
    if((i - 2009) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.BMmultirate, "FDR/fdr.quasse.BMmultirate.csv")

## discretetrait ##

states.quasse <- list()

for(i in 2060:2109){
  for (j in 1:50) {
    if((i - 2059) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.discretetrait, "FDR/fdr.quasse.discretetrait.csv")

## nosignal ##

states.quasse <- list()

for(i in 2110:2159){
  for (j in 1:50) {
    if((i - 2109) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.nosignal, "FDR/fdr.quasse.nosignal.csv")

## OUstrong ##

states.quasse <- list()

for(i in 2160:2209){
  for (j in 1:50) {
    if((i - 2159) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.OUstrong, "FDR/fdr.quasse.OUstrong.csv")

## OUweak ##

states.quasse <- list()

for(i in 2210:2259){
  for (j in 1:50) {
    if((i - 2209) == j) {
      # load .txt frequency tables
      states.quasse[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.quasse)) {
  
  print(i)
  states <- states.quasse[[i]]
  
  for (j in 1:length(trees.quasse)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.quasse[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.quasse[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.quasse.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.quasse.OUweak, "FDR/fdr.quasse.OUweak.csv")

#### FDR test (slowdown) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.slowdown <- read.tree("trees_FDR/trees_slowdown.tre")
states.list <- list.files(paste("~/traits_FDR",
                                sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.slowdown <- list()

for(i in 2260:2309){
  for (j in 1:50) {
    if((i - 2259) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.1cladefixed, "FDR/fdr.slowdown.1cladefixed.csv")

## 1cladenosignal ##

states.slowdown <- list()

for(i in 2310:2359){
  for (j in 1:50) {
    if((i - 2309) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.1cladenosignal, "FDR/fdr.slowdown.1cladenosignal.csv")

## BM ##

states.slowdown <- list()

for(i in 2360:2409){
  for (j in 1:50) {
    if((i - 2359) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.BM <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.BM, "FDR/fdr.slowdown.BM.csv")

## BMjump ##

states.slowdown <- list()

for(i in 2410:2459){
  for (j in 1:50) {
    if((i - 2409) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.BMjump, "FDR/fdr.slowdown.BMjump.csv")

## BMmultirate ##

states.slowdown <- list()

for(i in 2460:2509){
  for (j in 1:50) {
    if((i - 2459) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.BMmultirate, "FDR/fdr.slowdown.BMmultirate.csv")

## discretetrait ##

states.slowdown <- list()

for(i in 2510:2559){
  for (j in 1:50) {
    if((i - 2509) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.discretetrait, "FDR/fdr.slowdown.discretetrait.csv")

## nosignal ##

states.slowdown <- list()

for(i in 2560:2609){
  for (j in 1:50) {
    if((i - 2559) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.nosignal, "FDR/fdr.slowdown.nosignal.csv")

## OUstrong ##

states.slowdown <- list()

for(i in 2610:2659){
  for (j in 1:50) {
    if((i - 2609) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.OUstrong, "FDR/fdr.slowdown.OUstrong.csv")

## OUweak ##

states.slowdown <- list()

for(i in 2660:2709){
  for (j in 1:50) {
    if((i - 2659) == j) {
      # load .txt frequency tables
      states.slowdown[[j]] <- read.table(states.list[[i]], header = T, row.names = 1)
    }
  }
}

essim.p.slope <- vector()

for (i in 1:length(trees.slowdown)) {
  
  print(i)
  states <- states.slowdown[[i]]
  
  for (j in 1:length(trees.slowdown)) {
    
    print(j)
    trait <- states[, j]
    names(trait) <- row.names(states.slowdown[[i]])
    
    # ES-sim
    
    essim.res <- essim_glm(trees.slowdown[[i]], trait, nsim = 1000)
    
    essim.p.slope <- c(essim.p.slope, essim.res[4])
    
  }
}

essim.fdr.slowdown.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/2500

write.csv(essim.fdr.slowdown.OUweak, "FDR/fdr.slowdown.OUweak.csv")
