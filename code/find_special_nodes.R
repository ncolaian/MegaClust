#This script will take a rooted tree and return a list of nodes at a RED value of 0.35 that contain at least 10% of the genomes
#These nodes will then be used to reroot the trees and average RED values for other nodes in a later script

# to calculate Relative Evolutionary Divergence (RED)
# This is an operational approximation of relative time with extant taxa existing
# in the present (RED=1) or the the last common ancestor some fixed time in the past (RED=0)
# The internal nodes are linearly interpolated between these values according to lineage specific rates of evolution
# red is calculated by:

#### p + (d/u) x (1-p)

### p = Red of the parent node
### d = branch length to the parent node
### u = avg branch length from the parent node to all extant taxa descendant from n

#I also want to keep it consistent to the tree we are using in the figure.

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(ape)
library(phytools)
library(castor)

#This file will get the tips of all the strains associated with an MCA

params = matrix(c(
  "tree", "t", 1, "character",
  "out_dir", "o", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

get_nodes <- function(tree) {
  red_vals <- (length(tree$tip.label)+1):(length(tree$edge.length)+1)
  red_vals <- c(red_vals, get_reds(tree) )
  red_vals <- matrix(red_vals, ncol=2)
  
  #get the nodes that reach the cutoff threshold
  over_cutoff <- red_vals[which(red_vals[,2] > 0.35),]
  over_cutoff <- over_cutoff[order(over_cutoff[,2]),]
  tips_accounted_for <- c()
  nodelist <- c()
  for ( i in 1:nrow(over_cutoff) ) {
    tips <- extract.clade(tree,over_cutoff[i,1])$tip.label
    if ( sum(tips %in% tips_accounted_for) == length(tips) ) {
      next
    }
  if ( length(tips) > (length(tree$tip.label) / 10) ) {
      nodelist <- append(nodelist, over_cutoff[i,1])
  }
    tips_accounted_for <- append(tips_accounted_for, tips)
  }
  return(nodelist)
}


#basically I will be trying to identify the nodes in which meet this criteria

tree <- read.tree(opt$tree)
nodes <- get_nodes(tree)
write.table(nodes, paste(opt$out_dir, "/nodes_to_root.tsv", sep=""), sep="\t", row.names = F, quote = F, col.names = F)





