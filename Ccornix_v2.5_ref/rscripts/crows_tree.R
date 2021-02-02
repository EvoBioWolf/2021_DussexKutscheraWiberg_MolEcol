###
# Plotting and modifying the phylogenetic tree
# Last Modified: 29.04.2017
###
# Clear all
rm(list = ls(all=TRUE))

library(ggplot2)
source("~/RData/RScripts/ggplot_theme.R")

#install.packages("devtools")
#library(devtools)
#install.packages("ape")
library(ape)
#install.packages("ade4")
library(ade4)
#install_github("igraph/rigraph")
library(igraph)
#install.packages("phangorn")
library(igraph)
#install.packages("phytools")
library(phytools)
#source("https://bioconductor.org/biocLite.R")
##?BiocUpgrade
#biocLite()
#biocLite("BiocUpgrade")
#biocLite("ggtree")
library(ggtree)
# Set colours
or <-"#E69F00"
bl <-"#0072B2"

# Load the trees
crow_tree_af0.5<-read.tree("~/PhD/crows_project/results/darren_phylogeny/RAxML_bipartitionsBranchLabels.crows_tree_Ccornix-ref-af0.5_GTRGAMMA_rc4_boot1k")  
crow_tree_mskVars<-read.raxml("~/PhD/crows_project/results/darren_phylogeny/RAxML_bipartitionsBranchLabels.crows_tree_Ccornix-ref-mskVars_GTRGAMMA_rc4_boot1k")
crow_tree_af0.5<-as.phylo(crow_tree_af0.5)
crow_tree_mskVars<-as.phylo(crow_tree_mskVars)
str(crow_tree_af0.5)
str(crow_tree_mskVars)
# Make a multiphylo object
trees<-as.multiPhylo(c(crow_tree_af0.5,crow_tree_mskVars))
str(trees)
# Calculate the Robinson-Fould distance between the two trees. 
multiRF(trees)
crow_tree<-crow_tree_af0.5
is.ultrametric(crow_tree)
str(crow_tree)


crow_tree_af0.5<-read.raxml("~/PhD/crows_project/results/darren_phylogeny/RAxML_bipartitionsBranchLabels.crows_tree_Ccornix-ref-mskVars_GTRGAMMA_rc4_boot1k")
crow_tree<-crow_tree_af0.5
str(crow_tree)
# Plot with ggtree
treeplot<- ggtree(crow_tree) + geom_tiplab(size=3)
treeplot + geom_treescale(x=0.06,width = 0.01) + ggplot2::xlim(0,0.08)
crow_tree@phylo$tip.label
# Make nice labels for the tree
species_names<-c(
  'paste(bolditalic("C. corone "),bold(Group))',
  'paste(bolditalic("C. brachyrhynchos"))',
  'paste(bolditalic("C. woodfordi"))',
  'paste(bolditalic("C. tasmanicus"))',
  'paste(bolditalic("C. moneduloides"))',
  'paste(bolditalic("C. kubaryi"))',
  'paste(bolditalic("C. hawaiiensis"))',
  'paste(bolditalic("C. corax"))',
  'paste(bolditalic("C. splendens"))',
  'paste(bolditalic("C. frugilegus"))',
  'paste(bolditalic("C. dauuricus"))',
  'paste(bolditalic("C. monedula"))',
  'paste(bolditalic("T. guttata"))')
crow_tree@phylo$tip.label<-species_names
tree_data<-data.frame(species=as.factor(crow_tree@phylo$tip.label),
                      tool=as.factor(c("n","n","n","n","y","n",
                                       "y","n","n","n","n","n","n")))
rownames(tree_data)<-crow_tree@phylo$tip.label

treeplot<- ggtree(crow_tree) %<+% tree_data + 
  geom_tiplab(size=4,align=TRUE,aes(color=tool),
              linesize=0.5,parse=TRUE)+
  scale_colour_manual("",values=c("grey50",or),guide=FALSE)
treeplot<-treeplot+geom_treescale(x = 0.07,offset = -0.5,width = 0.01)
str(crow_tree)
crow_tree@data
treeplot+ggplot2::xlim(0.04,0.08)+geom_nodepoint(aes(subset=bootstrap > 95),
                                             colour="red",size=3.5,alpha=1/2)
treeplot+ggplot2::xlim(0.04,0.08)+geom_nodelab(aes(label=bootstrap))

# For PAML:
# Load tree again
crow_tree<-read.tree("~/PhD/crows_project/results/darren_phylogeny/RAxML_bipartitionsBranchLabels.crows_tree_Ccornix-ref-mskVars_GTRGAMMA_rc4_boot1k")
crow_tree$tip.label
plot(crow_tree)

# Prune the tree
# 5 species set
species5 <- c("3sp_Ccorone","Cdau","Cfru","Csple",
              "Cmon","Tgut")
drop_species5 <- crow_tree$tip.label[which(!(crow_tree$tip.label %in% species5))]
species5_tree <- drop.tip(crow_tree,drop_species5)
plot(as.phylo(unroot(species5_tree)))
write.tree(as.phylo(unroot(species5_tree)),
            "~/PhD/crows_project/results/phylogeny/5species_crows_Nulltree_ur.nex")
# 7 species set
species7 <- c("3sp_Ccorone","Cdau","Cfru","Csple",
              "Cmon","Ctas","Ccorx","Tgut")
drop_species7 <- crow_tree$tip.label[which(!(crow_tree$tip.label %in% species7))]
species7_tree <- drop.tip(crow_tree,drop_species7)
plot(as.phylo(unroot(species7_tree)))
write.tree(as.phylo(unroot(species7_tree)),
           "~/PhD/crows_project/results/phylogeny/7species_crows_Nulltree_ur.nex")
# 8 species set
species8 <- c("3sp_Ccorone","Cdau","Cfru","Csple",
              "Cmon","Ctas","Ccorx","Chaw","Tgut")
drop_species8 <- crow_tree$tip.label[which(!(crow_tree$tip.label %in% species8))]
species8_tree <- drop.tip(crow_tree,drop_species8)
plot(as.phylo(unroot(species8_tree)))
write.tree(as.phylo(unroot(species8_tree)),
           "~/PhD/crows_project/results/phylogeny/8species_crows_Nulltree_ur.nex")


# For CAFE:
# Make dataset of dates for some nodes (C. cornix v C. moneduloides node)
dates<-data.frame(node=c(14,19),age.min=c(36,10),age.max=c(50,11),soft.bounds=c(NA))
# Make the tree ultrametric
crow_tree_um<-chronos(crow_tree,lambda = 0.5,calibration = dates)
plot(crow_tree_um)
add.scale.bar()
crow_tree_um$edge.length<-round(crow_tree_um$edge.length,0)
tips<-c("Tgut","3sp_Ccornix","Cmon")
crow_tree_um <- drop.tip(crow_tree_um,
                         tip = crow_tree_um$tip.label[
                           which(!(crow_tree_um$tip.label %in% tips))])

str(crow_tree_um)
ggtree(crow_tree_um)+scale_x_ggtree()
add.scale.bar(x = 1,y=1.2,lwd = 2)
text("(Million Years)",y=1,x=4,cex = 0.6)

str(crow_tree_um)
# Write the ultrametric pruned tree to a file.
write.tree(crow_tree_um,
           "~/PhD/crows_project/results/darren_phylogeny/crows_CAFEtree.new")


