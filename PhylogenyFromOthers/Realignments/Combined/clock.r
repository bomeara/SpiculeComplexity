library(ape)
library(phytools)
phy <-  read.nexus("classifiedtree.treefile")
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
min = ALL 545
max = ALL 545
outfile = intree.dated.tre
thorough
prime
randomcv
", file="treeplRun.txt")
system("treepl treeplRun.txt")

cat("treefile = intree.tre
smooth = 100
numsites = 11462
mrca = ALL Corticium_vaceleti Hexadella_indica
min = ALL 545
max = ALL 545
outfile = intree.dated.final.tre
thorough
opt = 5
optad = 5
optcvad = 0
", file="treeplRun2.txt")
system("treepl treeplRun2.txt")

phy_clock <- read.tree("intree.dated.final.tre")
pdf(file="clock.pdf", width=10, height=10); plot(phy, type="fan", use.edge.length = FALSE, cex=0.2); dev.off(); system("open clock.pdf")