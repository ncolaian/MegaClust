#this script will take a rooted tree and calculate the red scores for each of the nodes on the tree
#it will then return the clusters that meet some RED criteria

#The original paper randomly assigns the root and averages the nodes. I will not do that.
#Since they are defining species, I think that the root should not play a factor when trying
#to identify closely related species in the tree

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
  "percent_ident", "p", 1, "character",
  "out_dir", "o", 1, "character",
  "nodelist", "n", 2, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

grab_group_with_max_redval <- function(tree, red_cutoff, nodelist) {
  red_vals <- (length(tree$tip.label)+1):(length(tree$edge.length)+1)
  red_vals <- c(red_vals, get_reds(tree) )
  red_vals <- matrix(red_vals, ncol=2)
  if ( length(nodelist) != 0 ) {
    for (node in nodelist){
      red_vals[1, 2] <- red_vals[1, 2] + red_vals[which(red_vals[,1] == node), 2]
    }
    for (node in nodelist) {
      tree1 <- root( tree, node = as.integer(node) )
      reds <- tree1$node.label
      for(i in 2:length(reds)){
        reds[i] = as.integer(as.character(reds[i]))
      }
      reds <- c(reds, get_reds(tree1))
      reds <- matrix(reds, ncol=2)
      for (j in 2:(nrow(red_vals)-1)){
        red_vals[j,2] = red_vals[j,2] + as.numeric(reds[which(reds[,1] == red_vals[j,1]), 2])
        
      }
    }
  }
  test <- red_vals[order(red_vals[,2]),]
  red_vals[,2] <- ( red_vals[,2] / (1+length(nodelist)) )
  #get the nodes that reach the cutoff threshold
  over_cutoff <- red_vals[which(red_vals[,2] > red_cutoff),]
  over_cutoff <- over_cutoff[order(over_cutoff[,2]),]
  tips_accounted_for <- c()
  groups <- c()
  for ( i in 1:nrow(over_cutoff) ) {
    tips <- extract.clade(tree,over_cutoff[i,1])$tip.label
    if ( sum(tips %in% tips_accounted_for) == length(tips) ) {
      next
    }
    tips_accounted_for <- append(tips_accounted_for, tips)
    groups <- append(groups, c(over_cutoff[i,1], paste(tips, collapse = ",")))
  }
  return(groups)
}


#basically I will be trying to identify the nodes in which meet this criteria

tree <- read.tree(opt$tree)
tree$node.label = (length(tree$tip.label)+1):(length(tree$edge.length)+1)
nodes <- readLines(opt$nodelist)
out_groups <- grab_group_with_max_redval(tree,opt$percent_ident, nodes)
out_groups <- matrix(out_groups, ncol=2, byrow = T)
write.table(out_groups, paste(opt$out_dir, "/red_groups.tsv", sep=""), sep="\t", row.names = F, quote = F, col.names = F)
