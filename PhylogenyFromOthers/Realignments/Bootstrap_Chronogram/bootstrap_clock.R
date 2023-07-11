library(ape)
library(phytools)
phy <- read.nexus("bootstrap_tree_classified.treefile")
outgroups <- c("Corticium_vaceleti", "Leuclathrina_translucida")
root_descendant_node <- ape::getMRCA(phy, outgroups)
root_descendant_edge_length <- phy$edge.length[which(phy$edge[,2]==root_descendant_node)]
phy <- phytools::reroot(phy, root_descendant_node, position=0.5*root_descendant_edge_length)
pdf(file="tree.pdf", width=10, height=10); plot(phy, type="fan", use.edge.length = FALSE, cex=0.2); dev.off(); system("open tree.pdf")

write.tree(phy, file="intree.tre")
cat("treefile = intree.tre
smooth = 100
numsites = 11462
mrca = ALL Corticium_vaceleti Hexadella_indica
min = ALL 535
max = ALL 535
outfile = intree.dated.tre
thorough
prime
randomcv
", file="treeplRun.txt")
system("treepl treeplRun.txt")

#change lines 33-35 based on output of last code block
cat("treefile = intree.tre
smooth = 100
numsites = 11462
mrca = ALL Corticium_vaceleti Hexadella_indica
min = ALL 535
max = ALL 535
outfile = intree.dated.final.tre
thorough
opt = 5
optad = 2
optcvad = 3
", file="treeplRun2.txt")
system("treepl treeplRun2.txt")

phy_clock <- read.tree("intree.dated.final.tre")
pdf(file="clock.pdf", width=10, height=10); plot(phy, type="fan", use.edge.length = FALSE, cex=0.2); dev.off(); system("open clock.pdf")

#the part of the output used for the second code block:
#PLACE THE LINES BELOW IN THE CONFIGURATION FILE
#opt = 5
#optad = 2
#optcvad = 3