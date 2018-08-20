##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                             WGCNA                                                                      #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# QC PIPELINE VERSION: v02
# DATE: 30/10/2017

#
# NOTES - 
# Run WGCNA on processed data
# seperate by tissue
#
#
#

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

Tissue="Entorhinal_Cortex"
# Tissue="Frontal_Cortex"
# Tissue="Temporal_Cortex"
# Tissue="Cerebellum"

##### LOAD LIBRARIES ####

library(WGCNA)
library(dynamicTreeCut)
library(flashClust)
library(Hmisc)
library(qvalue)
library(org.Hs.eg.db)

#install.packages("flashClust")

#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")

# enable multi threading

enableWGCNAThreads()

##### SET DIRECTORIES #####

temp_work_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/6.WGCNA/"
data_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/1.Data/2.Process_data/clean_data"

setwd(temp_work_dir)

##### CREATE DIRECTORIES #####

# directory with tissue

dir.create(paste(temp_work_dir, Tissue, sep="/"))
work_dir=paste(temp_work_dir,Tissue, sep="/")

# directory for correlation plots
dir.create(paste(work_dir,"Correlation_plots", sep="/"))
Correlation_plots_dir=paste(work_dir,"Correlation_plots", sep="/")

# directory for network preservation plots
dir.create(paste(work_dir,"Network_preservation_plots", sep="/"))
Network_preservation_plots_dir=paste(work_dir,"Network_preservation_plots", sep="/")

# directory for module overlap
dir.create(paste(work_dir,"Module_overlap_plots", sep="/"))
Module_overlap_plots_dir=paste(work_dir,"Module_overlap_plots", sep="/")

# directory for gene list in each module
dir.create(paste(work_dir,"Modules", sep="/"))
dir.create(paste(work_dir,"Modules/CO", sep="/"))
dir.create(paste(work_dir,"Modules/CO_AD", sep="/"))
dir.create(paste(work_dir,"Modules/AD", sep="/"))

CO_module_dir=paste(work_dir,"Modules/CO", sep="/")
CO_AD_module_dir=paste(work_dir,"Modules/CO_AD", sep="/")
AD_module_dir=paste(work_dir,"Modules/AD", sep="/")

# directory for module relationship

dir.create(paste(work_dir,"Module_relationship_plots", sep="/"))
Module_relationship_plots_dir=paste(work_dir,"Module_relationship_plots", sep="/")

# directory for kME

dir.create(paste(work_dir,"kME", sep="/"))
dir.create(paste(work_dir,"kME/CO", sep="/"))
dir.create(paste(work_dir,"kME/CO_AD", sep="/"))
dir.create(paste(work_dir,"kME/AD", sep="/"))

CO_kME_dir=paste(work_dir,"kME/CO", sep="/")
CO_AD_kME_dir=paste(work_dir,"kME/CO_AD", sep="/")
AD_kME_dir=paste(work_dir,"kME/AD", sep="/")

# directory for hub genes
dir.create(paste(work_dir, "Hub_genes", sep="/"))
hub_gene_dir<-paste(work_dir, "Hub_genes", sep="/")

# directory for Knonwn AD genes
dir.create(paste(work_dir, "Known_AD_genes", sep="/"))
Known_AD_genes_dir<-paste(work_dir, "Known_AD_genes", sep="/")

##### ADDITIONAL SCRIPT REQUIRED #####

setwd(temp_work_dir)

# download from https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/JMiller/
# Other required files (.zip)

source("tutorialFunctions.R")

setwd(work_dir)

##### LOAD DATA ####

setwd(data_dir)

expression_data<-read.table("MRC-LBB_expression_table_QCd.txt", head=T, as.is=T, sep="\t", check.names = F)
phenotype_data<-read.table("MRC-LBB_phenotype.txt", head=T, as.is=T, sep="\t")

head(expression_data)[1:5]
head(phenotype_data)

dim(expression_data)
dim(phenotype_data)

##### EXTARCT TISSUE OF INTEREST #####

#split data by tissue + diagnosis - require columns for samples and rows for probes

CO_exprs<-as.data.frame(t(subset(expression_data, rownames(expression_data) %in% rownames(subset(phenotype_data, PHENOTYPE=="CO" & TISSUE==Tissue)))))

CO_AD_exprs<-EC_CO_AD<-as.data.frame(t(subset(expression_data, rownames(expression_data) %in% rownames(subset(phenotype_data, PHENOTYPE=="CO_AD" & TISSUE==Tissue)))))

AD_exprs<-EC_AD<-as.data.frame(t(subset(expression_data, rownames(expression_data) %in% rownames(subset(phenotype_data, PHENOTYPE=="AD" & TISSUE==Tissue)))))

dim(CO_exprs)
dim(CO_AD_exprs)
dim(AD_exprs)

##### PICK SOFT POWER #####

# soft threshold is a value used to power correaltion of the genes.

# create function

pickSoftThreshold_function<-function(datExpr){
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
}

#apply function and set softpower threshold
pickSoftThreshold_function(as.data.frame(t(AD_exprs)))
AD_softPower=8

pickSoftThreshold_function(as.data.frame(t(CO_AD_exprs)))
CO_AD_softPower=5

pickSoftThreshold_function(as.data.frame(t(CO_exprs)))
CO_softPower=6

##### CORRELATING NETWROK PROPERTIES #####

#checking comaprability of datsets.

setwd(Correlation_plots_dir)

# within same brain region - but across AD vs AD_CO + AD_CO vs CO

# create function 

plot_correlation<- function(dataset1, dataset1_SF, dataset2, dataset2_SF, dataset3, dataset3_SF){
  
  rankExprA1= rank(rowMeans(dataset1)) 
  rankExprA2= rank(rowMeans(dataset2)) 
  
  rankConnA1= rank(softConnectivity(t(dataset1),type="signed",power=dataset1_SF)) 
  rankConnA2= rank(softConnectivity(t(dataset2),type="signed",power=dataset2_SF)) 
  rankExprB1= rank(rowMeans(dataset2)) 
  rankExprB2= rank(rowMeans(dataset3)) 
  
  rankConnB1= rank(softConnectivity(t(dataset2),type="signed",power=dataset2_SF)) 
  rankConnB2= rank(softConnectivity(t(dataset3),type="signed",power=dataset3_SF)) 
  
  par(mfrow=c(2,2)) 
  verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (A1)",  
                     ylab="Ranked Expression (A2)") 
  verboseScatterplot(rankConnA1,rankConnA2, xlab="Ranked Connectivity (A1)",  
                     ylab="Ranked Connectivity (A2)") 
  verboseScatterplot(rankExprB1,rankExprB2, xlab="Ranked Expression (B1)",  
                     ylab="Ranked Expression (B2)") 
  verboseScatterplot(rankConnB1,rankConnB2, xlab="Ranked Connectivity (B1)",  
                     ylab="Ranked Connectivity (B2)") 
}

#apply function
plot_correlation(AD_exprs, AD_softPower, CO_AD_exprs, CO_AD_softPower, CO_exprs, CO_softPower)

#plot to pdf
pdf(paste(Tissue, "correlation_plot.PDF", sep="_"))
plot_correlation(AD_exprs, AD_softPower, CO_AD_exprs, CO_AD_softPower, CO_exprs, CO_softPower)
dev.off()

setwd(work_dir)

##### CALCULATE WGCNA VARIABLES #####

# calculate WGCNA variables for AD

# calculate correlation network adjacency
AD_adjacency = adjacency(t(AD_exprs),power=AD_softPower,type="signed"); 
#set diagnal of matrix to zero
diag(AD_adjacency)=0 
# convert adjacency matrix to topology matrix- then dissimilarity (1-TOM)
AD_dissTOM   = 1-TOMsimilarity(AD_adjacency, TOMType="signed") 
# hierarchial cluster tree
AD_geneTree  = flashClust(as.dist(AD_dissTOM), method="average") 

### calculate WGCNA variables for CO_AD
CO_AD_adjacency = adjacency(t(CO_AD_exprs),power=CO_AD_softPower,type="signed"); 
diag(CO_AD_adjacency)=0 
CO_AD_dissTOM   = 1-TOMsimilarity(CO_AD_adjacency, TOMType="signed") 
CO_AD_geneTree  = flashClust(as.dist(CO_AD_dissTOM), method="average") 

### calculate WGCNA variables for CO
CO_adjacency = adjacency(t(CO_exprs),power=CO_softPower,type="signed"); 
diag(CO_adjacency)=0 
CO_dissTOM   = 1-TOMsimilarity(CO_adjacency, TOMType="signed") 
CO_geneTree  = flashClust(as.dist(CO_dissTOM), method="average") 

##### BASE MODULE COLOURS #####

#colours - CO to be baseline
#CO
CO_mColorh=NULL 

# detect clusters, assign colours
for (ds in 0:3){ 
  tree = cutreeHybrid(dendro = CO_geneTree, pamStage=FALSE, 
                      minClusterSize = (30-3*ds), cutHeight = 0.99,  
                      deepSplit = ds, distM = CO_dissTOM) 
  CO_mColorh=cbind(CO_mColorh,labels2colors(tree$labels)); 
} 

# use 1st assigned colours
CO_modules =  CO_mColorh[,1] 

##### SAME MODULE COLOURS FOR CO_AD AND AD AS CO #####

# based on CO. same probe in CO will have same colour in CO_AD + AD

#module colours for CO_AD

CO_AD_mColorh<-NULL

# detect clusters, assign colours
for (ds in 0:3){ 
  CO_AD_tree = cutreeHybrid(dendro = CO_AD_geneTree, pamStage=FALSE, 
                            minClusterSize = (30-3*ds), cutHeight = 0.99,  
                            deepSplit = ds, distM = CO_AD_dissTOM) 
  CO_AD_mColorh=cbind(CO_AD_mColorh,labels2colors(CO_AD_tree$labels)); 
} 

# use 1st assigned colours
CO_AD_modules =  CO_AD_mColorh[,1] 

#modules colours for AD

AD_mColorh<-NULL

# detect clusters, assign colours
for (ds in 0:3){ 
  AD_tree = cutreeHybrid(dendro = AD_geneTree, pamStage=FALSE, 
                         minClusterSize = (30-3*ds), cutHeight = 0.99,  
                         deepSplit = ds, distM = AD_dissTOM) 
  AD_mColorh=cbind(AD_mColorh,labels2colors(AD_tree$labels)); 
} 

# use 1st assigned colours
AD_modules =  AD_mColorh[,1] 

# assign colours based on CO

##### PRESERVED MODULE PLOT #####

setwd(Network_preservation_plots_dir)

plotDendroAndColors(CO_geneTree, CO_modules, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors CO")  

plotDendroAndColors(CO_AD_geneTree, CO_modules, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors CO_AD")  

plotDendroAndColors(AD_geneTree, CO_modules, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors AD")  

#plot preservation
pdf(paste(Tissue, "Dendogram_and_Module_preservation_CO.pdf", sep="_"))
plotDendroAndColors(CO_geneTree, CO_modules, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors CO")  
dev.off()

pdf(paste(Tissue, "Dendogram_and_Module_preservation_CO_AD.pdf", sep="_"))
plotDendroAndColors(CO_AD_geneTree, CO_modules, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors CO_AD")  
dev.off()

pdf(paste(Tissue, "Dendogram_and_Module_preservation_AD.pdf", sep="_"))
plotDendroAndColors(AD_geneTree, CO_modules, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors AD")  
dev.off()

##### PRESERVED MODULE STATS #####

### calaculate how well a module is preserved
# all data expression into 1 list
multiExpr  = list(CO=list(data=t(CO_exprs)),CO_AD=list(data=t(CO_AD_exprs)), AD=list(data=t(AD_exprs))) 
# coloours to use
multiColor = list(CO = CO_modules) 
# calculate module preservation statisitcs
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed", 
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400) 

# extract median Rank + Zsummary stats
CO_AD_mp_stats<-cbind(mp$preservation$observed$ref.CO$inColumnsAlsoPresentIn.CO_AD[c(1,2)], mp$preservation$Z$ref.CO$inColumnsAlsoPresentIn.CO_AD[2])
CO_AD_mp_stats[order(-CO_AD_mp_stats$Zsummary.pres),]

AD_mp_stats<-cbind(mp$preservation$observed$ref.CO$inColumnsAlsoPresentIn.AD[c(1,2)], mp$preservation$Z$ref.CO$inColumnsAlsoPresentIn.AD[2])
AD_mp_stats[order(-AD_mp_stats$Zsummary.pres),]

##### PLOT MODULE PRESERVED STATS #####

#create function
plot_preservation_stats<-function(stats, stats2, label_pos, label_pos2, label_pos3, label_pos4){
  par(mfrow = c(2,2))
  
  plot(stats[,1], stats[, 2],  col = 1, bg = rownames(stats), pch = 21,  
       main = "AsymAD preservation median rank",
       cex = 2.4,
       ylab = "Preservation median rank", 
       xlab ="Module size" , 
       ylim = c(max(stats[,2])+2, 0),
       xlim = c(10, (max(stats[,1])+200))) 
  labelPoints(stats[, 1], stats[, 2], rownames(stats), offs = label_pos, xpd=T, protectEdges=T, jiggle=T)
  abline(lm(stats[, 2]~stats[, 1]), col="red")
  
  plot(stats[,1], stats[, 3],  col = 1, bg = rownames(stats), pch = 21,  
       main = "AsymAD preservation Zsummary",
       cex = 2.4,
       ylab = "Preservation Zsummary", 
       xlab ="Module size" , 
       ylim = c(0, max(stats[,3])+5),
       xlim = c(10, (max(stats[,1])+200)))
  labelPoints(stats[, 1], stats[, 3], rownames(stats), offs = label_pos2, xpd=T, protectEdges=T, jiggle=T)
  abline(lm(stats[, 3]~stats[, 1]), col="red")
  abline(h=2, col ="blue", lty = 2)
  abline(h=10,col ="darkgreen", lty = 2)

  plot(stats2[,1], stats2[, 2],  col = 1, bg = rownames(stats2), pch = 21,  
       main = "AD preservation median rank",
       cex = 2.4,
       ylab = "Preservation median rank", 
       xlab ="Module size" , 
       ylim = c(max(stats2[,2])+2, 0),
       xlim = c(10, (max(stats2[,1])+200)))
  labelPoints(stats2[, 1], stats2[, 2], rownames(stats2), offs = label_pos3, xpd=T, protectEdges=T, jiggle=T)
  abline(lm(stats2[, 2]~stats2[, 1]), col="red")
  
  plot(stats2[,1], stats2[, 3],  col = 1, bg = rownames(stats2), pch = 21,  
       main = "AD preservation Zsummary",
       cex = 2.4,
       ylab = "Preservation Zsummary", 
       xlab ="Module size" , 
       ylim = c(0, max(stats2[,3])+5),
       xlim = c(10, (max(stats2[,1])+200)))
  labelPoints(stats2[, 1], stats2[, 3], rownames(stats2), offs = label_pos4, xpd=T, protectEdges=T, jiggle=T)
  abline(lm(stats2[, 3]~stats2[, 1]), col="red")
  abline(h=2, col ="blue", lty = 2)
  abline(h=10,col ="darkgreen", lty = 2)
  
  }

plot_preservation_stats(CO_AD_mp_stats, AD_mp_stats, 0, 0.1, 0, 0.1)


# plot to PDF

setwd(Network_preservation_plots_dir)

tiff("Module_preservation_stats.tiff", width = 8, height = 8, unit="in",  res=300)
plot_preservation_stats(CO_AD_mp_stats, AD_mp_stats, 0, 0.1, 0, 0.1)
dev.off()

##### MODULE COLOURS BASED ON MATCHING ######

# CO_AD and AD colours matched to CO. 
# e.g control group has a blue module with 100 genes. 
# If a module in CO_AD overlaps most with the control blue module, 
# it is assigned the colour blue as well, 
# if the blue module has already been assigned by another module, 
# it is assigned a new colour.


# function taken from WGCNA tutorial

#apply function

CO_AD_modules_new = matchModules(rownames(CO_exprs), CO_modules, rownames(CO_exprs), CO_AD_modules)
AD_modules_new = matchModules(rownames(CO_exprs), CO_modules, rownames(CO_exprs), AD_modules)

##### MODULE OVERLAP #####

# across all 3 phenotype - mapping colours to CO

plotDendroAndColors(CO_geneTree, cbind(CO_modules, CO_AD_modules_new, AD_modules_new), c("CO", "CO_AD", "AD"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="dendrogram and module colors mapped to CO")  

# plot to  PDF

setwd(Module_overlap_plots_dir)

pdf(paste(Tissue, "Module_overlap.pdf", sep="_"))
plotDendroAndColors(CO_geneTree, cbind(CO_modules, CO_AD_modules_new, AD_modules_new), c("CO", "Pre-MCI", "AD"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="dendrogram and module colors mapped to CO")  
dev.off()

tiff(paste(Tissue, "Module_overlap.tiff", sep="_"), width = 8, height = 8, unit="in",  res=300)
plotDendroAndColors(CO_geneTree, cbind(CO_modules, CO_AD_modules_new, AD_modules_new), c("CO", "AsymAD", "AD"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Entorhinal Cortex Module Preservation")  
dev.off()



##### EXPORT GENE LIST IN EACH MODULE FOR CO #####

setwd(CO_module_dir)

CO_module_colour_dataframe<-cbind(rownames(CO_exprs), as.data.frame(CO_modules))
colnames(CO_module_colour_dataframe)<-c("Entrez_ID", "Module_colour")
head(CO_module_colour_dataframe)
table(CO_module_colour_dataframe$Module_colour)

for (x in 1:length(unique(CO_module_colour_dataframe$Module_colour))) {
  # number of colours
  colours<-unique(CO_module_colour_dataframe$Module_colour)[x]
  #extract colour
  module<-subset(CO_module_colour_dataframe, Module_colour==colours)
  # write out entrez gene ID and save as colour
  write.table(module[1], file=paste(colours, ".txt", sep=""), row.names=F, col.names=F, quote=F)
}

write.table(as.data.frame(rownames(CO_exprs)), file="background_list.txt", row.names=F, col.names=F, quote=F)

##### EXPORT GENE LIST IN EACH MODULE FOR CO_AD #####

setwd(CO_AD_module_dir)

CO_AD_module_colour_dataframe<-cbind(rownames(CO_exprs), as.data.frame(CO_AD_modules_new))
colnames(CO_AD_module_colour_dataframe)<-c("Entrez_ID", "Module_colour")
head(CO_AD_module_colour_dataframe)
table(CO_AD_module_colour_dataframe$Module_colour)

for (x in 1:length(unique(CO_AD_module_colour_dataframe$Module_colour))) {
  # number of colours
  colours<-unique(CO_AD_module_colour_dataframe$Module_colour)[x]
  #extract colour
  module<-subset(CO_AD_module_colour_dataframe, Module_colour==colours)
  # write out entrez gene ID and save as colour
  write.table(module[1], file=paste(colours, ".txt", sep=""), row.names=F, col.names=F, quote=F)
}

write.table(as.data.frame(rownames(CO_exprs)), file="background_list.txt", row.names=F, col.names=F, quote=F)

##### EXPORT GENE LIST IN EACH MODULE FOR AD #####

setwd(AD_module_dir)

AD_module_colour_dataframe<-cbind(rownames(CO_exprs), as.data.frame(AD_modules_new))
colnames(AD_module_colour_dataframe)<-c("Entrez_ID", "Module_colour")
head(AD_module_colour_dataframe)
table(AD_module_colour_dataframe$Module_colour)

for (x in 1:length(unique(AD_module_colour_dataframe$Module_colour))) {
  # number of colours
  colours<-unique(AD_module_colour_dataframe$Module_colour)[x]
  #extract colour
  module<-subset(AD_module_colour_dataframe, Module_colour==colours)
  # write out entrez gene ID and save as colour
  write.table(module[1], file=paste(colours, ".txt", sep=""), row.names=F, col.names=F, quote=F)
}

write.table(as.data.frame(rownames(CO_exprs)), file="background_list.txt", row.names=F, col.names=F, quote=F)

##### HEATMAP FUNCTION #####

### function taken from WGCNA

#---------------------------------------------------------------------------------------------------------
# HeatmapWithTextLabels.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# ReverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


ReverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

ReverseVector = function(Vector)
{
  ind = seq(from=length(Vector), to=1, by=-1);
  Vector[ind];
  #Vector
}

#--------------------------------------------------------------------------
#
# HeatmapWithTextLabels = function ( Matrix, xLabels, yLabels, ... ) { 
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix 
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent 
# with the dimensions in.
# If ColorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

HeatmapWithTextLabels = function ( Matrix, xLabels, yLabels = NULL, xSymbols = NULL, ySymbols = NULL, 
                                   InvertColors=FALSE, ColorLabels = NULL, xColorLabels = FALSE,
                                   yColorLabels = FALSE,
                                   SetMargins = TRUE,
                                   colors = NULL, NumMatrix = NULL, cex.Num = NULL, cex.lab = NULL, 
                                   plotLegend = TRUE, ... ) 
{
  if (!is.null(ColorLabels)) {xColorLabels = ColorLabels; yColorLabels = ColorLabels; }
  
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(Matrix)[1]==dim(Matrix)[2])) 
    yLabels = xLabels; 
  if (SetMargins)
  {
    if (ColorLabels)
    {
      par(mar=c(2,2,3,5)+0.2);
    } else {
      par(mar = c(7,7,3,5)+0.2);
    }
  }
  if (is.null(colors)) 
    #if (IncludeSign)
    #{
    #colors = GreenRedWhite(50);
    #} else {
    colors = heat.colors(30);
  #}
  if (InvertColors) colors = ReverseVector(colors);
  if (plotLegend)
  {
    image.plot(t(ReverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  } else {
    image(z = t(ReverseRows(Matrix)), xaxt = "n", xlab="", yaxt="n", ylab="", col=colors, ...);
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
  nxlabels = length(xLabels)
  #   plot x axis labels using:
  #   par("usr")[3] - 0.25 as the vertical placement
  #   srt = 45 as text rotation angle
  #   adj = 1 to place right end of text at tick mark
  #   pd = TRUE to allow for text outside the plot region
  plotbox = par("usr");
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;
  
  # print(paste("plotbox:", plotbox[1], plotbox[2], plotbox[3], plotbox[4]));
  
  nylabels = length(yLabels)
  axis(2, labels = FALSE, tick = FALSE)
  xspacing = 1/(nxlabels-1); yspacing = 1/(nylabels-1)
  # print(paste("nxlabels:", nxlabels));
  if (!xColorLabels)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(((1:nxlabels)-1)*xspacing , ymin - 0.02, srt = 45, 
         adj = 1, labels = xLabels, xpd = TRUE, cex = cex.lab)
  } else {
    rect(((1:nxlabels)-1)*xspacing - xspacing/2, ymin-xspacing*1.2,
         ((1:nxlabels)-1)*xspacing + xspacing/2, ymin-xspacing*0.2,
         density = -1,  col = substring(xLabels, 3), border = substring(xLabels, 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( ((1:nxlabels)-1)*xspacing, ymin-xspacing*1.3, xSymbols, adj = c(0.5, 0), xpd = TRUE);
  }
  if (!yColorLabels)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - 0.01*xrange, ((1:nylabels)-1)/(nylabels-1), srt = 0, 
         adj = 1, labels = ReverseVector(yLabels), xpd = TRUE, cex = cex.lab )
  } else {
    rect(xmin-yspacing*1.2, ((nylabels:1)-1)*yspacing - yspacing/2,
         xmin-yspacing*0.2, ((nylabels:1)-1)*yspacing + yspacing/2, 
         density = -1,  col = substring(yLabels, 3), border = substring(yLabels, 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin-yspacing*1.2, ((nylabels:1)-1)*yspacing, ySymbols, adj = c(1, 0.5), xpd = TRUE);
  }
  
  if (!is.null(NumMatrix))
  {
    if (is.null(cex.Num)) cex.Num = par("cex");
    #if (dim(NumMatrix)!=dim(Matrix))
    #  stop("HeatmapWithTextLabels: NumMatrix was given, but has dimensions incompatible with Matrix.");
    for (rw in 1:dim(Matrix)[1])
      for (cl in 1:dim(Matrix)[2])
      {
        text((cl-1)*xspacing, (dim(Matrix)[1]-rw)*yspacing, 
             as.character(NumMatrix[rw,cl]), xpd = TRUE, cex = cex.Num, adj = c(0.5, 0.5));
      }
  }
  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
}

##### MODULE RELATIONSHIP CO vs CO_AD #####

# We calculate the overlaps
# of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
# assign a p-value to each of the pairwise overlaps.

# Isolate the module labels in the order they appear in module eigengenes
CO_ModuleLabels = unique(CO_modules)
CO_AD_ModuleLabels = unique(CO_AD_modules_new)

#number of modules
nCO_Mod = length(CO_ModuleLabels)
nCO_AD_Mod =  length(CO_AD_ModuleLabels)

# Initialize tables of p-values and of the corresponding counts
CO_vs_CO_AD_pTable = matrix(0, nrow = nCO_Mod, ncol = nCO_AD_Mod);
CO_vs_CO_AD_CountTbl = matrix(0, nrow = nCO_Mod, ncol = nCO_AD_Mod)

# Execute all pairwaise comparisons
for (CO in 1:nCO_Mod)
  for (CO_AD in 1:nCO_AD_Mod)
  {
    CO_Members = (CO_modules == CO_ModuleLabels[CO]);
    CO_AD_Members = (CO_AD_modules_new == CO_AD_ModuleLabels[CO_AD]);
    CO_vs_CO_AD_pTable[CO, CO_AD] = -log10(fisher.test(CO_Members, CO_AD_Members, alternative = "greater")$p.value);
    CO_vs_CO_AD_CountTbl[CO, CO_AD] = sum(CO_modules == CO_ModuleLabels[CO] & CO_AD_modules_new ==
                                            CO_AD_ModuleLabels[CO_AD])
  }


# Truncate p values smaller than 10^{-50} to 10^{-50}
CO_vs_CO_AD_pTable[is.infinite(CO_vs_CO_AD_pTable)] = 1.3*max(CO_vs_CO_AD_pTable[is.finite(CO_vs_CO_AD_pTable)]);
CO_vs_CO_AD_pTable[CO_vs_CO_AD_pTable>50 ] = 50 ;

# Marginal counts (really module sizes)
CO_ModTotals = apply(CO_vs_CO_AD_CountTbl, 1, sum)
CO_AD_ModTotals = apply(CO_vs_CO_AD_CountTbl, 2, sum)

# Actual plotting
sizeGrWindow(10,7 );

par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_CO_AD_pTable,
               xLabels = paste(" ", CO_AD_ModuleLabels),
               yLabels = paste(" ", CO_ModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("AsymAD ", CO_AD_ModuleLabels, ": ", CO_AD_ModTotals, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels, ": ", CO_ModTotals, sep=""),
               textMatrix = CO_vs_CO_AD_CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO and AsymAD modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);


setwd(Module_relationship_plots_dir)

tiff(paste(Tissue, "correspondence of CO and CO_AD modules.tiff", sep="_"), width = 10, height = 8, unit="in",  res=300)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(9, 10, 5, 2));
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_CO_AD_pTable,
               xLabels = paste(" ", CO_AD_ModuleLabels),
               yLabels = paste(" ", CO_ModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("AsymAD ", CO_AD_ModuleLabels, ": ", CO_AD_ModTotals, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels, ": ", CO_ModTotals, sep=""),
               textMatrix = CO_vs_CO_AD_CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of CO and AsymAD modules in the EC",
               cex.text = 1, 
               cex.lab = 1, 
               setStdMargins = FALSE);
dev.off()

##### MODULE RELATIONSHIP CO_AD vs AD #####

# We calculate the overlaps
# of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
# assign a p-value to each of the pairwise overlaps.

# Isolate the module labels in the order they appear in module eigengenes
AD_ModuleLabels = unique(AD_modules_new)


#number of modules
nAD_Mod =  length(AD_ModuleLabels)

# Initialize tables of p-values and of the corresponding counts
CO_AD_vs_AD_pTable = matrix(0, nrow = nCO_AD_Mod, ncol = nAD_Mod);
CO_AD_vs_AD_CountTbl = matrix(0, nrow = nCO_AD_Mod, ncol = nAD_Mod)

# Execute all pairwaise comparisons
for (CO_AD in 1:nCO_AD_Mod)
  for (AD in 1:nAD_Mod)
  {
    CO_AD_Members = (CO_AD_modules_new == CO_AD_ModuleLabels[CO_AD]);
    AD_Members = (AD_modules_new == AD_ModuleLabels[AD]);
    CO_AD_vs_AD_pTable[CO_AD, AD] = -log10(fisher.test(CO_AD_Members, AD_Members, alternative = "greater")$p.value);
    CO_AD_vs_AD_CountTbl[CO_AD, AD] = sum(CO_AD_modules_new == CO_AD_ModuleLabels[CO_AD] & AD_modules_new ==
                                            AD_ModuleLabels[AD])
  }

# Truncate p values smaller than 10^{-50} to 10^{-50}
CO_AD_vs_AD_pTable[is.infinite(CO_AD_vs_AD_pTable)] = 1.3*max(CO_AD_vs_AD_pTable[is.finite(CO_AD_vs_AD_pTable)]);
CO_AD_vs_AD_pTable[CO_AD_vs_AD_pTable>50 ] = 50 ;

# Marginal counts (really module sizes)
CO_AD_ModTotals = apply(CO_AD_vs_AD_CountTbl, 1, sum)
AD_ModTotals = apply(CO_AD_vs_AD_CountTbl, 2, sum)

# Actual plotting
sizeGrWindow(10,7 );

par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_AD_vs_AD_pTable,
               xLabels = paste(" ", AD_ModuleLabels),
               yLabels = paste(" ", CO_AD_ModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels, ": ", AD_ModTotals, sep=""),
               ySymbols = paste("CO_AD ", CO_AD_ModuleLabels, ": ", CO_AD_ModTotals, sep=""),
               textMatrix = CO_AD_vs_AD_CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO_AD and AD modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 0.9, 
               setStdMargins = FALSE);


setwd(Module_relationship_plots_dir)

tiff(paste(Tissue, "correspondence of CO_AD and AD modules.tiff", sep="_"), width = 10, height = 8, unit="in",  res=300)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(9, 10, 5, 2));
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_AD_vs_AD_pTable,
               xLabels = paste(" ", AD_ModuleLabels),
               yLabels = paste(" ", CO_AD_ModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels, ": ", AD_ModTotals, sep=""),
               ySymbols = paste("AsymAD ", CO_AD_ModuleLabels, ": ", CO_AD_ModTotals, sep=""),
               textMatrix = CO_AD_vs_AD_CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of AsymAD and AD modules in the EC", 
               cex.text = 1, 
               cex.lab = 0.88, 
               setStdMargins = FALSE);
dev.off()

##### MODULE RELATIONSHIP CO vs AD #####

# We calculate the overlaps
# of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
# assign a p-value to each of the pairwise overlaps.

# Initialize tables of p-values and of the corresponding counts
CO_vs_AD_pTable = matrix(0, nrow = nCO_Mod, ncol = nAD_Mod);
CO_vs_AD_CountTbl = matrix(0, nrow = nCO_Mod, ncol = nAD_Mod)

# Execute all pairwaise comparisons
for (CO in 1:nCO_Mod)
  for (AD in 1:nAD_Mod)
  {
    CO_Members = (CO_modules == CO_ModuleLabels[CO]);
    AD_Members = (AD_modules_new == AD_ModuleLabels[AD]);
    CO_vs_AD_pTable[CO, AD] = -log10(fisher.test(CO_Members, AD_Members, alternative = "greater")$p.value);
    CO_vs_AD_CountTbl[CO, AD] = sum(CO_modules == CO_ModuleLabels[CO] & AD_modules_new ==
                                      AD_ModuleLabels[AD])
  }


# Truncate p values smaller than 10^{-50} to 10^{-50}
CO_vs_AD_pTable[is.infinite(CO_vs_AD_pTable)] = 1.3*max(CO_vs_AD_pTable[is.finite(CO_vs_AD_pTable)]);
CO_vs_AD_pTable[CO_vs_AD_pTable>50 ] = 50 ;

# Marginal counts (really module sizes)
CO_ModTotals = apply(CO_vs_AD_CountTbl, 1, sum)
AD_ModTotals = apply(CO_vs_AD_CountTbl, 2, sum)

# Actual plotting
sizeGrWindow(10,7 );

par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_AD_pTable,
               xLabels = paste(" ", AD_ModuleLabels),
               yLabels = paste(" ", CO_ModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels, ": ", AD_ModTotals, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels, ": ", CO_ModTotals, sep=""),
               textMatrix = CO_vs_AD_CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO and AD modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);


setwd(Module_relationship_plots_dir)

tiff(paste(Tissue, "correspondence of CO and AD modules.tiff", sep="_"), width = 10, height = 8, unit="in",  res=300)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(9, 10, 5, 2));
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_AD_pTable,
               xLabels = paste(" ", AD_ModuleLabels),
               yLabels = paste(" ", CO_ModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels, ": ", AD_ModTotals, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels, ": ", CO_ModTotals, sep=""),
               textMatrix = CO_vs_AD_CountTbl,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of CO and AD modules in the EC",
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);
dev.off()

##### MODULE MEMBERSHIP (kME) CO #####

setwd(CO_kME_dir)

# calculate module eigengenes
CO_PCs    = moduleEigengenes(t(CO_exprs),  colors=CO_modules)  
CO_ME    = CO_PCs$eigengenes 
CO_colors = names(table(CO_modules)) 

# calculate kME
CO_geneModuleMembership = signedKME(t(CO_exprs), CO_ME) 
colnames(CO_geneModuleMembership)=paste("PC", CO_colors ,".cor",sep="")

# calculate p value for gene kME to module
CO_MMPvalue=corPvalueStudent(as.matrix(CO_geneModuleMembership),dim(CO_exprs)[[2]])
colnames(CO_MMPvalue)=paste("PC",CO_colors,".pval",sep="")

#create kME table with cor + pvalue for each module
CO_Gene       = rownames(CO_exprs) 
CO_kMEtable  = cbind(CO_Gene,CO_Gene,CO_modules) 
for (i in 1:length(CO_colors)) 
  CO_kMEtable = cbind(CO_kMEtable, CO_geneModuleMembership[,i], CO_MMPvalue[,i]) 
colnames(CO_kMEtable)=c("PSID","Gene","Module",sort(c(colnames(CO_geneModuleMembership), 
                                                      colnames(CO_MMPvalue)))) 

# write out kME table

write.csv(CO_kMEtable[,2:dim(CO_kMEtable)[2]],"CO_kMEtable.csv",row.names=FALSE)

##### MODULE MEMBERSHIP (kME) CO_AD #####

##CO_AD
setwd(CO_AD_kME_dir)

# First calculate MEs for CO_AD
CO_AD_PCs = moduleEigengenes(t(CO_AD_exprs),  colors=CO_AD_modules_new)  
CO_AD_ME = CO_AD_PCs$eigengenes 
CO_AD_colors = names(table(CO_AD_modules_new)) 

# calculate kME for CO_AD
CO_AD_geneModuleMembership = signedKME(t(CO_AD_exprs), CO_AD_ME) 
colnames(CO_AD_geneModuleMembership)=paste("PC",CO_AD_colors,".cor",sep="");  

# calculate P value
CO_AD_MMPvalue=corPvalueStudent(as.matrix(CO_AD_geneModuleMembership),dim(CO_AD_exprs)[[2]]);  
colnames(CO_AD_MMPvalue)=paste("PC",CO_AD_colors,".pval",sep=""); 

#create kME table with cor + pvalue for each module
CO_AD_Gene = rownames(CO_AD_exprs) 
CO_AD_kMEtable  = cbind(CO_AD_Gene,CO_AD_Gene,CO_AD_modules_new) 

for (i in 1:length(CO_AD_colors)) 
  CO_AD_kMEtable = cbind(CO_AD_kMEtable, CO_AD_geneModuleMembership[,i], CO_AD_MMPvalue[,i]) 
colnames(CO_AD_kMEtable)=c("PSID","Gene","Module",sort(c(colnames(CO_AD_geneModuleMembership), 
                                                         colnames(CO_AD_MMPvalue)))) 

write.csv(CO_AD_kMEtable[,2:dim(CO_AD_kMEtable)[2]],"CO_AD_kMEtable.csv",row.names=FALSE)

##### MODULE MEMBERSHIP (kME) AD #####

setwd(AD_kME_dir)

# First calculate MEs for AD
AD_PCs = moduleEigengenes(t(AD_exprs),  colors=AD_modules_new)  
AD_ME = AD_PCs$eigengenes 
AD_colors = names(table(AD_modules_new)) 

# calculate kME for AD
AD_geneModuleMembership = signedKME(t(AD_exprs), AD_ME) 
colnames(AD_geneModuleMembership)=paste("PC",AD_colors,".cor",sep="");  

# calculate P value
AD_MMPvalue=corPvalueStudent(as.matrix(AD_geneModuleMembership),dim(AD_exprs)[[2]]);  
colnames(AD_MMPvalue)=paste("PC",AD_colors,".pval",sep=""); 

#create kME table with cor + pvalue for each module
AD_Gene = rownames(AD_exprs) 
AD_kMEtable  = cbind(AD_Gene,AD_Gene,AD_modules_new) 

for (i in 1:length(AD_colors)) 
  AD_kMEtable = cbind(AD_kMEtable, AD_geneModuleMembership[,i], AD_MMPvalue[,i]) 
colnames(AD_kMEtable)=c("PSID","Gene","Module",sort(c(colnames(AD_geneModuleMembership), 
                                                      colnames(AD_MMPvalue)))) 

write.csv(AD_kMEtable[,2:dim(AD_kMEtable)[2]],"AD_kMEtable.csv",row.names=FALSE)

##### CONVERT ENTREZ ID TO GENE SYMBOL #####

# aconvert entrez gene id to gene symbol
entrez_gene_symbol <- org.Hs.egSYMBOL

# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(entrez_gene_symbol)

# Convert to a list
gene_symbol_lookup_table <- as.data.frame(entrez_gene_symbol[mapped_genes])

# add X to entrez id
gene_symbol_lookup_table$gene_id<-paste("X", gene_symbol_lookup_table$gene_id, sep="")

head(gene_symbol_lookup_table)

# create convert_probe_id_to_entrez_id function 

convert_probe_id_to_entrez_id <- function(expression_dataset, probe_mapping_file){
  # transform dataset # - removed this step
  expression_dataset_t<-as.data.frame(t(expression_dataset))
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset_t[colnames(expression_dataset_t)%in%probe_mapping_file$gene_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$symbol[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$gene_id)]
  # transform
  data_frame_in_probe_mapper<-as.data.frame(t(data_frame_in_probe_mapper))
  return(data_frame_in_probe_mapper)
}

# using illuminaHumanv4.db_mapping_entrez_id_unique instead of Illumina_HT_12_v4_probe_list_entrez_id_unique from this point onwards

CO_geneModuleMembership_gene_symbol<-convert_probe_id_to_entrez_id(CO_geneModuleMembership, gene_symbol_lookup_table)


head(CO_geneModuleMembership_gene_symbol)[1:5]

dim(CO_geneModuleMembership_gene_symbol)
head(CO_geneModuleMembership)


CO_geneModuleMembership

##### MODULE HUB GENE FUNCTION #####

# create function to extract top hub genes per module
extract_hub_gene<-function(geneModuleMembership){
  # create empty vector
  top_hub_genes<-NULL
  # convert dataframe to genesymbol
  geneModuleMembership<-convert_probe_id_to_entrez_id(geneModuleMembership, gene_symbol_lookup_table)
  # iterate through all modules
  for (gene in 1:ncol(geneModuleMembership)){
    # extract module
    module_connectivity<-geneModuleMembership[gene]
    # rank by abs(cor)
    module_connectivity$rank<-rank(abs(geneModuleMembership[,gene]))
    # order by highest abs(cor)
    module_connectivity<-module_connectivity[order(module_connectivity$rank),]
    # extract top ten genes
    #top_hub_genes<-cbind(top_hub_genes, head(rownames(module_connectivity),10))
    top_hub_genes<-cbind(top_hub_genes, rownames(module_connectivity))
    # assign module colour
    colnames(top_hub_genes)[gene]<-colnames(geneModuleMembership)[gene]
  }
  # create rank column
  rank=as.data.frame(1:nrow(geneModuleMembership))
  colnames(rank)<-"Rank"
  top_hub_genes<-cbind(rank, top_hub_genes)
  return(top_hub_genes)
}

##### MODULE TOP 10 HUB GENES - CO #####

CO_top_10_module_hub_genes<-extract_hub_gene(CO_geneModuleMembership)
head(CO_top_10_module_hub_genes)

##### MODULE TOP 10 HUB GENES - CO_AD #####

CO_AD_top_10_module_hub_genes<-extract_hub_gene(CO_AD_geneModuleMembership)
head(CO_AD_top_10_module_hub_genes)

##### MODULE TOP 10 HUB GENES - AD #####

AD_top_10_module_hub_genes<-extract_hub_gene(AD_geneModuleMembership)
head(AD_top_10_module_hub_genes)

##### EXPORT HUB GENE LIST #####

setwd(hub_gene_dir)

write.table(CO_top_10_module_hub_genes, file="CO_module_hub_genes.txt", quote=F, row.names=F)
write.table(CO_AD_top_10_module_hub_genes, file="CO_AD_module_hub_genes.txt", quote=F, row.names=F)
write.table(AD_top_10_module_hub_genes, file="AD_module_hub_genes.txt", quote=F, row.names=F)

##### KNOWN AD GENES #####

knonwn_AD_genes<-read.table("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/8.Known_AD_genes/Alzgene_list.txt")
head(knonwn_AD_genes)
dim(knonwn_AD_genes)

# using previous created gene symbol/entrez gene mapping file gene_symbol_lookup_table

# Convert to a list
gene_symbol_lookup_table2 <- as.data.frame(entrez_gene_symbol[mapped_genes])

# keep 

knonwn_AD_genes_entrez_ID<-gene_symbol_lookup_table2[gene_symbol_lookup_table2$symbol %in% knonwn_AD_genes$V1,]
head(knonwn_AD_genes_entrez_ID)
dim(knonwn_AD_genes_entrez_ID)

# top 10
Well_known_AD_genes<-c("ABCA7",
                          "ADAM10",
                          "APOE",
                          "APP",
                          "BACE1",
                          "BIN1",
                          "C9ORF72",
                          "CD33",
                          "CLU",
                          "CR1",
                          "FUS",
                          "GBA",
                          "LRRK2",
                          "PICALM", 
                          "PSEN1",
                          "PSEN2",
                          "GRN",
                          "SNCA",
                          "SORL1",
                          "MAPT",
                          "TDP-43",
                          "TREM2")

Well_known_AD_genes_entrez_ID<-subset(knonwn_AD_genes_entrez_ID, knonwn_AD_genes_entrez_ID$symbol %in% Well_known_AD_genes)
Well_known_AD_genes_entrez_ID


# table of which module

#check files

head(CO_modules)
head(CO_AD_modules_new)
head(AD_modules_new)

head(CO_module_colour_dataframe)
head(CO_AD_module_colour_dataframe)
head(AD_module_colour_dataframe)

head(knonwn_AD_genes_entrez_ID)
Well_known_AD_genes_entrez_ID

#check top 10 known AD genes

Well_known_AD_genes_entrez_ID

CO_module_colour_dataframe[CO_module_colour_dataframe$Entrez_ID %in% Well_known_AD_genes_entrez_ID$gene_id,]

CO_AD_module_colour_dataframe[CO_AD_module_colour_dataframe$Entrez_ID %in% Well_known_AD_genes_entrez_ID$gene_id,]

AD_module_colour_dataframe[AD_module_colour_dataframe$Entrez_ID %in% Well_known_AD_genes_entrez_ID$gene_id,]

#extract well known AD gene changes
well_known_AD_gene_changes<-cbind(CO_module_colour_dataframe[CO_module_colour_dataframe$Entrez_ID %in% Well_known_AD_genes_entrez_ID$gene_id,],
      CO_AD_module_colour_dataframe[CO_AD_module_colour_dataframe$Entrez_ID %in% Well_known_AD_genes_entrez_ID$gene_id,][2],
      AD_module_colour_dataframe[AD_module_colour_dataframe$Entrez_ID %in% Well_known_AD_genes_entrez_ID$gene_id,][2])

# merge gene symbol and change colnames
colnames(well_known_AD_gene_changes)<-c("gene_id", "CO", "CO_AD", "AD")
well_known_AD_gene_changes<-merge(well_known_AD_gene_changes, knonwn_AD_genes_entrez_ID, by="gene_id")

well_known_AD_gene_changes<-well_known_AD_gene_changes[c(1,5,2,3,4)]

well_known_AD_gene_changes

# exttract all AD known genes

CO_modules_known_AD<-CO_module_colour_dataframe[CO_module_colour_dataframe$Entrez_ID %in% knonwn_AD_genes_entrez_ID$gene_id,]

CO_AD_modules_known_AD<-CO_AD_module_colour_dataframe[CO_AD_module_colour_dataframe$Entrez_ID %in% knonwn_AD_genes_entrez_ID$gene_id,]

AD_modules_known_AD<-AD_module_colour_dataframe[AD_module_colour_dataframe$Entrez_ID %in% knonwn_AD_genes_entrez_ID$gene_id,]

head(CO_modules_known_AD)
head(CO_AD_modules_known_AD)
head(AD_modules_known_AD)

# check same order - should be TRUE
all(CO_modules_known_AD$Entrez_ID==CO_AD_modules_known_AD$Entrez_ID)==TRUE
all(CO_modules_known_AD$Entrez_ID==AD_modules_known_AD$Entrez_ID)==TRUE

# keep only colours

CO_modules_known_AD<-as.character(CO_modules_known_AD$Module_colour)
CO_AD_modules_known_AD<-as.character(CO_AD_modules_known_AD$Module_colour)
AD_modules_known_AD<-as.character(AD_modules_known_AD$Module_colour)

length(CO_modules_known_AD)
length(CO_AD_modules_known_AD)
length(AD_modules_known_AD)

##### MODULE RELATIONSHIP CO vs CO_AD - KNOWN AD GENES #####

# We calculate the overlaps
# of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
# assign a p-value to each of the pairwise overlaps.

# Isolate the module labels in the order they appear in module eigengenes
CO_ModuleLabels_known_AD = unique(CO_modules_known_AD)
CO_AD_ModuleLabels_known_AD = unique(CO_AD_modules_known_AD)

#number of modules
nCO_Mod_known_AD = length(CO_ModuleLabels_known_AD)
nCO_AD_Mod_known_AD =  length(CO_AD_ModuleLabels_known_AD)

# Initialize tables of p-values and of the corresponding counts
CO_vs_CO_AD_pTable_known_AD = matrix(0, nrow = nCO_Mod_known_AD, ncol = nCO_AD_Mod_known_AD);
CO_vs_CO_AD_CountTbl_known_AD = matrix(0, nrow = nCO_Mod_known_AD, ncol = nCO_AD_Mod_known_AD)

# Execute all pairwaise comparisons
for (CO in 1:nCO_Mod_known_AD)
  for (CO_AD in 1:nCO_AD_Mod_known_AD)
  {
    CO_Members_known_AD = (CO_modules_known_AD == CO_ModuleLabels_known_AD[CO]);
    CO_AD_Members_known_AD = (CO_AD_modules_known_AD == CO_AD_ModuleLabels_known_AD[CO_AD]);
    CO_vs_CO_AD_pTable_known_AD[CO, CO_AD] = -log10(fisher.test(CO_Members_known_AD, CO_AD_Members_known_AD, alternative = "greater")$p.value);
    CO_vs_CO_AD_CountTbl_known_AD[CO, CO_AD] = sum(CO_modules_known_AD == CO_ModuleLabels_known_AD[CO] & CO_AD_modules_known_AD ==
                                                     CO_AD_ModuleLabels_known_AD[CO_AD])
  }


# Truncate p values smaller than 10^{-50} to 10^{-50}
CO_vs_CO_AD_pTable_known_AD[is.infinite(CO_vs_CO_AD_pTable_known_AD)] = 1.3*max(CO_vs_CO_AD_pTable_known_AD[is.finite(CO_vs_CO_AD_pTable_known_AD)]);
CO_vs_CO_AD_pTable_known_AD[CO_vs_CO_AD_pTable_known_AD>50 ] = 50 ;

# Marginal counts (really module sizes)
CO_ModTotals_known_AD = apply(CO_vs_CO_AD_CountTbl_known_AD, 1, sum)
CO_AD_ModTotals_known_AD = apply(CO_vs_CO_AD_CountTbl_known_AD, 2, sum)

# Actual plotting
sizeGrWindow(10,7 );

par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_CO_AD_pTable_known_AD,
               xLabels = paste(" ", CO_AD_ModuleLabels_known_AD),
               yLabels = paste(" ", CO_ModuleLabels_known_AD),
               colorLabels = TRUE,
               xSymbols = paste("CO_AD ", CO_AD_ModuleLabels_known_AD, ": ", CO_AD_ModTotals_known_AD, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels_known_AD, ": ", CO_ModTotals_known_AD, sep=""),
               textMatrix = CO_vs_CO_AD_CountTbl_known_AD,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO and CO_AD known AD gene modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);


setwd(Known_AD_genes_dir)

pdf(paste(Tissue, "correspondence of CO and Pre-MCI known AD gene modules.PDF", sep="_"))
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(9, 10, 5, 2));
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_CO_AD_pTable_known_AD,
               xLabels = paste(" ", CO_AD_ModuleLabels_known_AD),
               yLabels = paste(" ", CO_ModuleLabels_known_AD),
               colorLabels = TRUE,
               xSymbols = paste("CO_AD ", CO_AD_ModuleLabels_known_AD, ": ", CO_AD_ModTotals_known_AD, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels_known_AD, ": ", CO_ModTotals_known_AD, sep=""),
               textMatrix = CO_vs_CO_AD_CountTbl_known_AD,
               colors = blueWhiteRed(100)[50:100],
               main = paste("Correspondence of known AD genes in CO and pre-MCI modules", sep=" "),
               cex.text = 0.8, 
               cex.lab = 0.8, 
               setStdMargins = FALSE);
dev.off()

##### MODULE RELATIONSHIP CO_AD vs AD - KNOWN AD GENES #####

# We calculate the overlaps
# of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
# assign a p-value to each of the pairwise overlaps.

# Isolate the module labels in the order they appear in module eigengenes
AD_ModuleLabels_known_AD = unique(AD_modules_known_AD)


#number of modules
nAD_Mod_known_AD =  length(AD_ModuleLabels_known_AD)

# Initialize tables of p-values and of the corresponding counts
CO_AD_vs_AD_pTable_known_AD = matrix(0, nrow = nCO_AD_Mod_known_AD, ncol = nAD_Mod_known_AD);
CO_AD_vs_AD_CountTbl_known_AD = matrix(0, nrow = nCO_AD_Mod_known_AD, ncol = nAD_Mod_known_AD)

# Execute all pairwaise comparisons
for (CO_AD in 1:nCO_AD_Mod_known_AD)
  for (AD in 1:nAD_Mod_known_AD)
  {
    CO_AD_Members_known_AD = (CO_AD_modules_known_AD == CO_AD_ModuleLabels_known_AD[CO_AD]);
    AD_Members_known_AD = (AD_modules_known_AD == AD_ModuleLabels_known_AD[AD]);
    CO_AD_vs_AD_pTable_known_AD[CO_AD, AD] = -log10(fisher.test(CO_AD_Members_known_AD, AD_Members_known_AD, alternative = "greater")$p.value);
    CO_AD_vs_AD_CountTbl_known_AD[CO_AD, AD] = sum(CO_AD_modules_known_AD == CO_AD_ModuleLabels_known_AD[CO_AD] & AD_modules_known_AD ==
                                                     AD_ModuleLabels_known_AD[AD])
  }

# Truncate p values smaller than 10^{-50} to 10^{-50}
CO_AD_vs_AD_pTable_known_AD[is.infinite(CO_AD_vs_AD_pTable_known_AD)] = 1.3*max(CO_AD_vs_AD_pTable_known_AD[is.finite(CO_AD_vs_AD_pTable_known_AD)]);
CO_AD_vs_AD_pTable_known_AD[CO_AD_vs_AD_pTable_known_AD>50 ] = 50 ;

# Marginal counts (really module sizes)
CO_AD_ModTotals_known_AD = apply(CO_AD_vs_AD_CountTbl_known_AD, 1, sum)
AD_ModTotals_known_AD = apply(CO_AD_vs_AD_CountTbl_known_AD, 2, sum)

# Actual plotting
sizeGrWindow(10,7 );

par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_AD_vs_AD_pTable_known_AD,
               xLabels = paste(" ", AD_ModuleLabels_known_AD),
               yLabels = paste(" ", CO_AD_ModuleLabels_known_AD),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels_known_AD, ": ", AD_ModTotals_known_AD, sep=""),
               ySymbols = paste("CO_AD ", CO_AD_ModuleLabels_known_AD, ": ", CO_AD_ModTotals_known_AD, sep=""),
               textMatrix = CO_AD_vs_AD_CountTbl_known_AD,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO_AD and AD known AD gene modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);


setwd(Known_AD_genes_dir)

pdf(paste(Tissue, "correspondence of CO_AD and AD known AD gene modules.PDF", sep="_"))
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_AD_vs_AD_pTable_known_AD,
               xLabels = paste(" ", AD_ModuleLabels_known_AD),
               yLabels = paste(" ", CO_AD_ModuleLabels_known_AD),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels_known_AD, ": ", AD_ModTotals_known_AD, sep=""),
               ySymbols = paste("CO_AD ", CO_AD_ModuleLabels_known_AD, ": ", CO_AD_ModTotals_known_AD, sep=""),
               textMatrix = CO_AD_vs_AD_CountTbl_known_AD,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO_AD and AD known AD gene modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);
dev.off()

##### MODULE RELATIONSHIP CO vs AD - KNOWN AD GENES #####

# We calculate the overlaps
# of each pair of female-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
# assign a p-value to each of the pairwise overlaps.

# Initialize tables of p-values and of the corresponding counts
CO_vs_AD_pTable_known_AD = matrix(0, nrow = nCO_Mod_known_AD, ncol = nAD_Mod_known_AD);
CO_vs_AD_CountTbl_known_AD = matrix(0, nrow = nCO_Mod_known_AD, ncol = nAD_Mod_known_AD)

# Execute all pairwaise comparisons
for (CO in 1:nCO_Mod_known_AD)
  for (AD in 1:nAD_Mod_known_AD)
  {
    CO_Members_known_AD = (CO_modules_known_AD == CO_ModuleLabels_known_AD[CO]);
    AD_Members_known_AD = (AD_modules_known_AD == AD_ModuleLabels_known_AD[AD]);
    CO_vs_AD_pTable_known_AD[CO, AD] = -log10(fisher.test(CO_Members_known_AD, AD_Members_known_AD, alternative = "greater")$p.value);
    CO_vs_AD_CountTbl_known_AD[CO, AD] = sum(CO_modules_known_AD == CO_ModuleLabels_known_AD[CO] & AD_modules_known_AD ==
                                               AD_ModuleLabels_known_AD[AD])
  }


# Truncate p values smaller than 10^{-50} to 10^{-50}
CO_vs_AD_pTable_known_AD[is.infinite(CO_vs_AD_pTable_known_AD)] = 1.3*max(CO_vs_AD_pTable_known_AD[is.finite(CO_vs_AD_pTable_known_AD)]);
CO_vs_AD_pTable_known_AD[CO_vs_AD_pTable_known_AD>50 ] = 50 ;

# Marginal counts (really module sizes)
CO_ModTotals_known_AD = apply(CO_vs_AD_CountTbl_known_AD, 1, sum)
AD_ModTotals_known_AD = apply(CO_vs_AD_CountTbl_known_AD, 2, sum)

# Actual plotting
sizeGrWindow(10,7 );

par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_AD_pTable_known_AD,
               xLabels = paste(" ", AD_ModuleLabels_known_AD),
               yLabels = paste(" ", CO_ModuleLabels_known_AD),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels_known_AD, ": ", AD_ModTotals_known_AD, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels_known_AD, ": ", CO_ModTotals_known_AD, sep=""),
               textMatrix = CO_vs_AD_CountTbl_known_AD,
               colors = blueWhiteRed(100)[50:100],
               main = paste(Tissue, "correspondence of CO and AD  known AD gene modules", sep=" "),
               cex.text = 1.0, 
               cex.lab = 1.0, 
               cex.main=2,
               setStdMargins = FALSE);


setwd(Known_AD_genes_dir)

tiff(paste(Tissue, "correspondence of CO and AD known AD gene modules.tiff", sep="_"), width = 10, height = 8, unit="in",  res=1200)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings

labeledHeatmap(Matrix = CO_vs_AD_pTable_known_AD,
               xLabels = paste(" ", AD_ModuleLabels_known_AD),
               yLabels = paste(" ", CO_ModuleLabels_known_AD),
               colorLabels = TRUE,
               xSymbols = paste("AD ", AD_ModuleLabels_known_AD, ": ", AD_ModTotals_known_AD, sep=""),
               ySymbols = paste("CO ", CO_ModuleLabels_known_AD, ": ", CO_ModTotals_known_AD, sep=""),
               textMatrix = CO_vs_AD_CountTbl_known_AD,
               colors = blueWhiteRed(100)[50:100],
               main = "Correspondence of CO and AD modules containing known AD genes in the EC",
               cex.text = 1.0, 
               cex.lab = 1.0, 
               setStdMargins = FALSE);
dev.off()

##### SAVE IMAGE #####

setwd(work_dir)

save.image(paste(Tissue, "WGCNA.Rdata", sep="_"))
