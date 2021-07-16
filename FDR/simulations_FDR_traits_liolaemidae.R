
# Copyright (c) 2017, Michael G. Harvey
# All rights reserved.

# A script to simulate traits on trees.

library(ape)
library(caper)
library(phytools)
library(diversitree)

nsim <- 1000 # Number of traits to simulate per tree per treeset

### Simulate traits ###

# BM (uncorrelated)

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		vv <- vcv.phylo(tree)
		sims <- t(rmvnorm(nsim, sigma=0.06*vv))
		rownames(sims) <- rownames(vv)		
		
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_BM_2.txt"))
	}	
}

# BM (multirate)

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
				
		# Isolate clade of > 10% of tips
		subtrees <- subtrees(tree)
		large.subtrees <- vector()
		for (k in 1:length(subtrees)) {
			if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
				large.subtrees <- c(large.subtrees, subtrees[k])				
			}
		}
		subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]

		# Increase branch lengths for that clade
		branches <- which.edge(tree, subclade$tip.label)
		tree$edge.length[branches] <- tree$edge.length[branches]*5
				
		# Trait simulations 
		vv <- vcv.phylo(tree)
		sims <- t(rmvnorm(nsim, sigma=0.06*vv))
		rownames(sims) <- rownames(vv)	
		
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_BMmultirate_2.txt"))
	}	
}

# BM - jump in mean

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
				
		# Trait simulations 
		vv <- vcv.phylo(tree)
		sims <- t(rmvnorm(nsim, sigma=0.06*vv))
		rownames(sims) <- rownames(vv)	
		
		# Isolate clade of > 10% of tips
		subtrees <- subtrees(tree)
		large.subtrees <- vector()
		for (k in 1:length(subtrees)) {
			if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
				large.subtrees <- c(large.subtrees, subtrees[k])				
			}
		}
		subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]

		# Jump in mean for that clade
		sims[rownames(sims) %in% subclade$tip.label,] <- sims[rownames(sims) %in% subclade$tip.label,]+0.3
						
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_BMjump_2.txt"))
	}	
}

# No phylogenetic signal

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		sim.vals <- rnorm(length(tree$tip.label)*1000, mean = 0, sd = max(branching.times(tree)))
		sims <- matrix(unlist(sim.vals), ncol=1000)
		rownames(sims) <- tree$tip.label	
		
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_nosignal_2.txt"))
	}	
}

# 1 clade with no phylogenetic signal

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		vv <- vcv.phylo(tree)
		sims <- t(rmvnorm(nsim, sigma=0.06*vv))
		rownames(sims) <- rownames(vv)		
		
		# Isolate clade of > 10% of tips
		subtrees <- subtrees(tree)
		large.subtrees <- vector()
		for (k in 1:length(subtrees)) {
			if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
				large.subtrees <- c(large.subtrees, subtrees[k])				
			}
		}
		subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]

		# No signal in that one clade
		norm.vals <- rnorm(length(subclade$tip.label)*1000, mean = 0, sd = max(branching.times(subclade)))
		sims[rownames(sims) %in% subclade$tip.label,] <- norm.vals

		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_1cladenosignal_2.txt"))
	}	
}

# 1 clade fixed

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		vv <- vcv.phylo(tree)
		sims <- t(rmvnorm(nsim, sigma=0.06*vv))
		rownames(sims) <- rownames(vv)		
		
		# Isolate clade of > 10% of tips
		subtrees <- subtrees(tree)
		large.subtrees <- vector()
		for (k in 1:length(subtrees)) {
			if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
				large.subtrees <- c(large.subtrees, subtrees[k])				
			}
		}
		subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]

		# Fix that clade to the value of one of its tips
		val <- sims[sample(1:length(subclade$tip.label), 1),]
		rep.val <- rep(val, each=nrow(sims[rownames(sims) %in% subclade$tip.label,]))
		sims[rownames(sims) %in% subclade$tip.label,] <- rep.val

		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_1cladefixed_2.txt"))
	}	
}

# Discrete trait (2 normal distributions)

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
				
		# Trait simulations 		
		statecols <- vector()
		for(k in 1:nsim) {
			states <- sim.character(tree, c(0.1,0.1), x0=0, model="mk2")
			states[states == 0] <- rnorm(length(states[states == 0]), mean = 0, sd = 2) # Normal dist. for state 0
			states[states == 1] <- rnorm(length(states[states == 1]), mean = 5, sd = 4) # Normal dist. for state 1
			statecols <- c(statecols, states[tree$tip.label])	
		}
		
		sims <- matrix(unlist(statecols), ncol=1000)		
		rownames(sims) <- tree$tip.label	
								
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_discretetrait_2.txt"))
	}	
}

# OU process (strong)

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.2, theta = 0, nsim = nsim) 		
		
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_OUstrong_2.txt"))
	}	
}

# OU process (weak)

treefiles <- list.files("",pattern="Dataset S2.tre", full.names=T)
for(i in 1:length(treefiles)){
	trees <- read.nexus(treefiles[i])
	if(length(trees) < 50) {
		trees <- list(trees)
	}
	for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.002, theta = 0, nsim = nsim) 		
		
		write.table(as.matrix(sims), paste0("traits_FDR/traits_", "Liolaemidae", "_OUweak_2.txt"))
	}	
}
