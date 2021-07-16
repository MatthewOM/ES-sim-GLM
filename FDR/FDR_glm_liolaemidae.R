
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

#### libraries ####

require(ape)
require(mvtnorm)
require(caper)
require(diversitree)

#### load functions ####

source("essim_glm.R")
source("essim_glm_multiple.R")

#### FDR test univariate (liolaemidae) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.liolaemidae <- read.nexus("Dataset S2.tre")
states.list <- list.files(paste("traits_FDR", sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[1]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.1cladefixed <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.1cladefixed, "FDR/fdr.liolaemidae.1cladefixed.csv")

## 1cladenosignal ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[3]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.1cladenosignal <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.1cladenosignal, "FDR/fdr.liolaemidae.1cladenosignal.csv")

## BM ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[5]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.BM <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.BM, "FDR/fdr.liolaemidae.BM.csv")

## BMjump ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[7]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.BMjump <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.BMjump, "FDR/fdr.liolaemidae.BMjump.csv")

## BMmultirate ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[9]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.BMmultirate <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.BMmultirate, "FDR/fdr.liolaemidae.BMmultirate.csv")

## discretetrait ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[11]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.discretetrait <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.discretetrait, "FDR/fdr.liolaemidae.discretetrait.csv")

## nosignal ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[13]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.nosignal <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.nosignal, "FDR/fdr.liolaemidae.nosignal.csv")

## OUstrong ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[15]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])
  
}

essim.fdr.liolaemidae.OUstrong <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.OUstrong, "FDR/fdr.liolaemidae.OUstrong.csv")

## OUweak ##

states.liolaemidae <- list()

# load .txt frequency tables
states.liolaemidae <- read.table(states.list[[17]], header = T, row.names = 1)

essim.p.slope <- vector()

for (i in 1:ncol(states.liolaemidae)) {
  
  print(i)
  trait <- states.liolaemidae[, i]
  names(trait) <- row.names(states.liolaemidae)
  
  # ES-sim
  
  essim.res <- essim_glm(trees.liolaemidae, trait, nsim = 1000)
  
  essim.p.slope <- c(essim.p.slope, essim.res[4])

}

essim.fdr.liolaemidae.OUweak <- length(essim.p.slope[essim.p.slope < 0.05])/ncol(states.liolaemidae)

write.csv(essim.fdr.liolaemidae.OUweak, "FDR/fdr.liolaemidae.OUweak.csv")

#### FDR test multivariate (liolaemidae) ####

## fdr will be tested on the simulations without trait-dependence, non-SSD (state-dependent diversification)

trees.liolaemidae <- read.nexus("~/U.A/3_PhD/Thesis/essim_glm/literature/Esquerre_et_al_2019_Data/LiolaemidaeDec17Combined_greenMCC.tre")
states.list <- list.files(paste("traits_FDR", sep = "/"), pattern = ".txt", full.name = T, recursive = F)

## 1cladefixed ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[1]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[2]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.1cladefixed.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.1cladefixed.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.1cladefixed <- rbind(liolaemidae.1cladefixed.slope.0.05, liolaemidae.1cladefixed.slope.0.1)

write.csv(essim.fdr.liolaemidae.1cladefixed, "FDR/fdr.liolaemidae.1cladefixed.multiple.csv")

## 1cladenosignal ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[3]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[4]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])

}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.1cladenosignal.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.1cladenosignal.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.1cladenosignal <- rbind(liolaemidae.1cladenosignal.slope.0.05, liolaemidae.1cladenosignal.slope.0.1)

write.csv(essim.fdr.liolaemidae.1cladenosignal, "FDR/fdr.liolaemidae.1cladenosignal.multiple.csv")

## BM ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[5]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[6]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.BM.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.BM.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.BM <- rbind(liolaemidae.BM.slope.0.05, liolaemidae.BM.slope.0.1)

write.csv(essim.fdr.liolaemidae.BM, "FDR/fdr.liolaemidae.BM.multiple.csv")

## BMjump ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[7]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[8]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.BMjump.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.BMjump.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.BMjump <- rbind(liolaemidae.BMjump.slope.0.05, liolaemidae.BMjump.slope.0.1)

write.csv(essim.fdr.liolaemidae.BMjump, "FDR/fdr.liolaemidae.BMjump.multiple.csv")

## BMmultirate ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[9]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[10]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.BMmultirate.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.BMmultirate.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.BMmultirate <- rbind(liolaemidae.BMmultirate.slope.0.05, liolaemidae.BMmultirate.slope.0.1)

write.csv(essim.fdr.liolaemidae.BMmultirate, "FDR/fdr.liolaemidae.BMmultirate.multiple.csv")

## discretetrait ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[11]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[12]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.discretetrait.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.discretetrait.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.discretetrait <- rbind(liolaemidae.discretetrait.slope.0.05, liolaemidae.discretetrait.slope.0.1)

write.csv(essim.fdr.liolaemidae.discretetrait, "FDR/fdr.liolaemidae.discretetrait.multiple.csv")

## nosignal ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[13]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[14]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.nosignal.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.nosignal.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.nosignal <- rbind(liolaemidae.nosignal.slope.0.05, liolaemidae.nosignal.slope.0.1)

write.csv(essim.fdr.liolaemidae.nosignal, "FDR/fdr.liolaemidae.nosignal.multiple.csv")

## OUstrong ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[15]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[16]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])
  
}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.OUstrong.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.OUstrong.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.OUstrong <- rbind(liolaemidae.OUstrong.slope.0.05, liolaemidae.OUstrong.slope.0.1)

write.csv(essim.fdr.liolaemidae.OUstrong, "FDR/fdr.liolaemidae.OUstrong.multiple.csv")

## OUweak ##

states.liolaemidae1 <- list()
states.liolaemidae2 <- list()

# load .txt frequency tables
states.liolaemidae1 <- read.table(states.list[[17]], header = T, row.names = 1)
states.liolaemidae2 <- read.table(states.list[[18]], header = T, row.names = 1)

essim.p.slope1 <- vector()
essim.p.slope2 <- vector()

for (i in 1:ncol(states.liolaemidae1)) {
  
  print(i)
  trait1 <- states.liolaemidae1[, i]
  names(trait1) <- row.names(states.liolaemidae1)
  trait2 <- states.liolaemidae2[, i]
  names(trait2) <- row.names(states.liolaemidae2)
  
  # ES-sim
  
  essim.res <- essim_glm_multiple(trees.liolaemidae, trait1, trait2, nsim = 1000)
  
  essim.p.slope1 <- c(essim.p.slope1, essim.res[4])
  essim.p.slope2 <- c(essim.p.slope2, essim.res[6])

}

vec1 <- essim.p.slope1 < 0.05 & essim.p.slope2 < 0.05; vec1 <- vec1[vec1=="TRUE"]
vec2 <- essim.p.slope1 <  0.1 & essim.p.slope2 <  0.1; vec2 <- vec2[vec2=="TRUE"]

liolaemidae.OUweak.slope.0.05 <- length(vec1)/ncol(states.liolaemidae1)
liolaemidae.OUweak.slope.0.1  <- length(vec2)/ncol(states.liolaemidae1)

essim.fdr.liolaemidae.OUweak <- rbind(liolaemidae.OUweak.slope.0.05, liolaemidae.OUweak.slope.0.1)

write.csv(essim.fdr.liolaemidae.OUweak, "FDR/fdr.liolaemidae.OUweak.multiple.csv")
