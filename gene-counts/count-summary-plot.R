#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

args = c("NewData/mt-barcodes-manual.csv", "Prelims/mt-tree-manual.phy", "NewData/mt-stats-means-manual.csv", "NewData/mt-simple-manual-indices.csv", "NewData/pt-barcodes-manual.csv", "Prelims/pt-tree-manual.phy", "NewData/pt-stats-means-manual.csv", "NewData/pt-simple-manual-indices.csv", "Plots/")

if(length(args) < 9) {
  stop("Need 9 argumnets -- MT barcodes, tree, stats, indices; PT barcodes, tree, stats, indices; output directory")
}

message("Loading libraries...")

library(ggplot2)
#library(ggtree)
#library(ggtreeExtra)
library(gridExtra)
library(phytools)
#library(igraph)
#library(cowplot)

mtbarcodefilename = args[1]
mttreefilename = args[2]
mtgenestatsfilename = args[3]
mtindexfilename = args[4]

ptbarcodefilename = args[5]
pttreefilename = args[6]
ptgenestatsfilename = args[7]
ptindexfilename = args[8]

out.dir = args[9]

sf = 3

#################### MT

#### get barcode data

message("Getting data...")

# read barcode data
mt.df <- read.csv(mtbarcodefilename, header=T)
mt.ngenes = ncol(mt.df)-1

mt.species.list = trimws(mt.df$Species)
mt.df$GeneCount = rowSums( mt.df[,2:ncol(mt.df)] ) 

# read taxonomy tree
mt.treeString = tolower(paste(readLines(mttreefilename), collapse=""))
mt.treeString = gsub(" ", "_", mt.treeString)
mt.tree <- read.tree(text=mt.treeString)
mt.tree$tip.label = gsub("_", " ", mt.tree$tip.label)
mt.to.drop = setdiff(mt.tree$tip.label, mt.species.list)
mt.tree = drop.tip(mt.tree, mt.to.drop)

# list of all node labels
mt.all.labels = c(mt.tree$tip.label, mt.tree$node.label)

# numbers of leaves and internal nodes
mt.n.tips = length(mt.tree$tip.label)
mt.n.nodes = length(mt.tree$node.label)

# identify tree root
mt.luca.node = grep("eukaryota", mt.all.labels)

# identify those clades that are directly descended from root...
mt.clades.refs = mt.tree$edge[which(mt.tree$edge[,1]==mt.luca.node),2]

# ... and which are not leaves
mt.clades.refs = mt.clades.refs[mt.clades.refs > mt.n.tips]

# get and sort labels for these clades
mt.clades = mt.all.labels[mt.clades.refs]
mt.clades = mt.clades[order(mt.clades)]

message("Building phylogenetic info...")

# build list of lists of organisms in each clade
mt.clades.list = list()
mt.counts.list = list()
mt.counts.df = data.frame()
for(i in 1:length(mt.clades)) {
  # identify root of clade
  mt.head.node = which(mt.all.labels == mt.clades[i])
  # get all descendants
  mt.head.des = getDescendants(mt.tree, mt.head.node)
  # get labels of descendants that are leaves
  mt.species = gsub("_", " ", mt.all.labels[mt.head.des[mt.head.des <= mt.n.tips]])
  # append these labels to this clade's list
  mt.clades.list = append(mt.clades.list, list(mt.species))
  counts = c()
  for(species in mt.species) {
    ref = which(mt.df$Species == species)
    if(length(ref) > 0) {
      counts = append(counts, mt.df$GeneCount[ref])
    }
  }
  mt.counts.list = append(mt.counts.list, list(counts))
  mt.counts.df = rbind(mt.counts.df, data.frame(Clade=mt.clades[i], Counts = counts))
}


mt.plot = ggplot(mt.counts.df, aes(x=Clade, y=Counts)) + geom_jitter() + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("mtDNA protein-coding gene count")

#################### MT

#### get barcode data

message("Getting data...")

# read barcode data
pt.df <- read.csv(ptbarcodefilename, header=T)
pt.ngenes = ncol(pt.df)-1

pt.species.list = trimws(pt.df$Species)
pt.df$GeneCount = rowSums( pt.df[,2:ncol(pt.df)] ) 

# read taxonomy tree
pt.treeString = tolower(paste(readLines(pttreefilename), collapse=""))
pt.treeString = gsub(" ", "_", pt.treeString)
pt.tree <- read.tree(text=pt.treeString)
pt.tree$tip.label = gsub("_", " ", pt.tree$tip.label)
pt.to.drop = setdiff(pt.tree$tip.label, pt.species.list)
pt.tree = drop.tip(pt.tree, pt.to.drop)

# list of all node labels
pt.all.labels = c(pt.tree$tip.label, pt.tree$node.label)

# numbers of leaves and internal nodes
pt.n.tips = length(pt.tree$tip.label)
pt.n.nodes = length(pt.tree$node.label)

# identify tree root
pt.luca.node = grep("eukaryota", pt.all.labels)

# identify those clades that are directly descended from root...
pt.clades.refs = pt.tree$edge[which(pt.tree$edge[,1]==pt.luca.node),2]

# ... and which are not leaves
pt.clades.refs = pt.clades.refs[pt.clades.refs > pt.n.tips]

# get and sort labels for these clades
pt.clades = pt.all.labels[pt.clades.refs]
pt.clades = pt.clades[order(pt.clades)]

message("Building phylogenetic info...")

# build list of lists of organisms in each clade
pt.clades.list = list()
pt.counts.list = list()
pt.counts.df = data.frame()
for(i in 1:length(pt.clades)) {
  # identify root of clade
  pt.head.node = which(pt.all.labels == pt.clades[i])
  # get all descendants
  pt.head.des = getDescendants(pt.tree, pt.head.node)
  # get labels of descendants that are leaves
  pt.species = gsub("_", " ", pt.all.labels[pt.head.des[pt.head.des <= pt.n.tips]])
  # append these labels to this clade's list
  pt.clades.list = append(pt.clades.list, list(pt.species))
  counts = c()
  for(species in pt.species) {
    ref = which(pt.df$Species == species)
    if(length(ref) > 0) {
      counts = append(counts, pt.df$GeneCount[ref])
    }
  }
  pt.counts.list = append(pt.counts.list, list(counts))
  pt.counts.df = rbind(pt.counts.df, data.frame(Clade=pt.clades[i], Counts = counts))
}


pt.plot = ggplot(pt.counts.df, aes(x=Clade, y=Counts)) + geom_jitter() + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("ptDNA protein-coding gene count")

png("count-summary-scatters.png", width=800*sf, height=300*sf, res=72*sf)
grid.arrange(mt.plot, pt.plot, nrow=1)
dev.off()