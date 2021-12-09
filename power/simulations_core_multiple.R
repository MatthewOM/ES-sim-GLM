# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

# A script to simulate two traits on trees.

#### libraries ####

require(ape)
require(mvtnorm)
require(caper)
require(diversitree)
require(mvMORPH)

source("essim.R")

#### Simulate two traits (50 tips) ####

trees50sdd <- vector("list", 100) # Vector to place trees
class(trees50sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {	
  repeat{
    print(i)
    m <- 0.004 # strength of relationship between trait and speciation
    ntax <- 50 # number of tips in tree
    linear.x <- make.linear.x(-5, 5)
    lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
    mu <- function(x) constant.x(x, 0)
    char <- make.brownian.with.drift(0, 0.06)
    phy <- NULL
    while(is.null(phy)) { # to avoid non-coalescence and keep the loop running
      try(
        phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=FALSE)
        )
      } 
    trees50sdd[[i]] <- phy
    count <- vector()
    repeat{ # get 2nd trait
      trait <- mvSIM(phy, model = "BM1", param = list(sigma = 0.06, theta = 0))
      rrr <- essim(phy, trait[,1], nsim = 1000)
      cor <- cor.test(phy$tip.state, trait[,1])
      count <- c(count, cor$parameter)
      if(rrr[2]<0.001 || length(count)==100){
        # check if trait from mvSIM is significantly related with speciation using essim with a significant correlation 
        # allow repeat function to run only 1000 times using count vector
        break
      }
    }
    if(rrr[2]<0.001){
      print(c(rrr, cor$estimate))
      print(c(abs(rrr[1])>0.5, rrr[2]<0.001, abs(cor$estimate)<0.5))
      break
    }
  }
  states <- cbind(states, phy$tip.state)
  states <- cbind(states, trait)
}

write.tree(trees50sdd, "trees_Core/trees50sdd_multiple_0.004.txt")
write.table(as.matrix(states), "traits_Core/trees50sdd_states_multiple_0.004.txt")

#### Simulate two traits (100 tips) ####

trees100sdd <- vector("list", 100) # Vector to place trees
class(trees100sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {	
  repeat{
    print(i)
    m <- 0.004 # strength of relationship between trait and speciation
    ntax <- 100 # number of tips in tree
    linear.x <- make.linear.x(-5, 5)
    lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
    mu <- function(x) constant.x(x, 0)
    char <- make.brownian.with.drift(0, 0.06)
    phy <- NULL
    while(is.null(phy)) { # to avoid non-coalescence and keep the loop running
      try(
        phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=FALSE)
      )
    } 
    trees100sdd[[i]] <- phy
    count <- vector()
    repeat{ # get 2nd trait
      trait <- mvSIM(phy, model = "BM1", param = list(sigma = 0.06, theta = 0))
      rrr <- essim(phy, trait[,1], nsim = 1000)
      cor <- cor.test(phy$tip.state, trait[,1])
      count <- c(count, cor$parameter)
      if(rrr[2]<0.001 || length(count)==100){
        # check if trait from mvSIM is significantly related with speciation using essim with a significant correlation 
        # allow repeat function to run only 1000 times using count vector
        break
      }
    }
    if(rrr[2]<0.001){
      print(c(rrr, cor$estimate))
      print(c(abs(rrr[1])>0.5, rrr[2]<0.001, abs(cor$estimate)<0.5))
      break
    }
  }
  states <- cbind(states, phy$tip.state)
  states <- cbind(states, trait)
}

write.tree(trees100sdd, "trees_Core/trees100sdd_multiple_0.004.txt")
write.table(as.matrix(states), "traits_Core/trees100sdd_states_multiple_0.004.txt")

#### Simulate two traits (250 tips) ####

trees250sdd <- vector("list", 100) # Vector to place trees
class(trees250sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {	
  repeat{
    print(i)
    m <- 0.004 # strength of relationship between trait and speciation
    ntax <- 250 # number of tips in tree
    linear.x <- make.linear.x(-5, 5)
    lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
    mu <- function(x) constant.x(x, 0)
    char <- make.brownian.with.drift(0, 0.06)
    phy <- NULL
    while(is.null(phy)) { # to avoid non-coalescence and keep the loop running
      try(
        phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=FALSE)
      )
    } 
    trees250sdd[[i]] <- phy
    count <- vector()
    repeat{ # get 2nd trait
      trait <- mvSIM(phy, model = "BM1", param = list(sigma = 0.06, theta = 0))
      rrr <- essim(phy, trait[,1], nsim = 1000)
      cor <- cor.test(phy$tip.state, trait[,1])
      count <- c(count, cor$parameter)
      if(rrr[2]<0.001 || length(count)==100){
        # check if trait from mvSIM is significantly related with speciation using essim with a significant correlation 
        # allow repeat function to run only 1000 times using count vector
        break
      }
    }
    if(rrr[2]<0.001){
      print(c(rrr, cor$estimate))
      print(c(abs(rrr[1])>0.5, rrr[2]<0.001, abs(cor$estimate)<0.5))
      break
    }
  }
  states <- cbind(states, phy$tip.state)
  states <- cbind(states, trait)
}

write.tree(trees250sdd, "trees_Core/trees250sdd_multiple_0.004.txt")
write.table(as.matrix(states), "traits_Core/trees250sdd_states_multiple_0.004.txt")

#### Simulate two traits (1250 tips) ####

trees1250sdd <- vector("list", 100) # Vector to place trees
class(trees1250sdd) <- "multiPhylo"
states <- vector()
for (i in 1:100) {	
  repeat{
    print(i)
    m <- 0.004 # strength of relationship between trait and speciation
    ntax <- 1250 # number of tips in tree
    linear.x <- make.linear.x(-5, 5)
    lambda <- function(x) linear.x(x, 0.025, m) # x, c, m (y=mx+c)
    mu <- function(x) constant.x(x, 0)
    char <- make.brownian.with.drift(0, 0.06)
    phy <- NULL
    while(is.null(phy)) { # to avoid non-coalescence and keep the loop running
      try(
        phy <- tree.quasse(c(lambda, mu, char), max.taxa=ntax, x0=0, single.lineage=FALSE, verbose=FALSE)
      )
    } 
    trees1250sdd[[i]] <- phy
    count <- vector()
    repeat{ # get 2nd trait
      trait <- mvSIM(phy, model = "BM1", param = list(sigma = 0.06, theta = 0))
      rrr <- essim(phy, trait[,1], nsim = 1000)
      cor <- cor.test(phy$tip.state, trait[,1])
      count <- c(count, cor$parameter)
      if(rrr[2]<0.001 || length(count)==100){
        # check if trait from mvSIM is significantly related with speciation using essim with a significant correlation 
        # allow repeat function to run only 1000 times using count vector
        break
      }
    }
    if(rrr[2]<0.001){
      print(c(rrr, cor$estimate))
      print(c(abs(rrr[1])>0.5, rrr[2]<0.001, abs(cor$estimate)<0.5))
      break
    }
  }
  states <- cbind(states, phy$tip.state)
  states <- cbind(states, trait)
}

write.tree(trees1250sdd, "trees_Core/trees1250sdd_multiple_0.004.txt")
write.table(as.matrix(states), "traits_Core/trees1250sdd_states_multiple_0.004.txt")
