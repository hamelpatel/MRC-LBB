########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                                                RUN DE                                                #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 02/10/2017

##### DESCRIPTION OF ANALYSIS ####
## Run DE using limma
##
## 3 groups - 
## AD vs CO
## AD vs CO_AD 
## CO_AD vs AD
##
## 4 tissues - 
## Entorhinal cortex
## Temporal cortex
## Frontal cortex
## Cerebellum
## 
#####

##### LOAD LIBRARIES ######

library(limma)
library(calibrate)

##### SET DIRECTORIES

work_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/3.DE_analysis/1.Results/"
data_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/1.Data/2.Process_data/clean_data"

setwd(work_dir)

# create directories - plots
dir.create(paste(work_dir,"volcano_plots", sep="/"))
volcano_plots_dir=paste(work_dir,"volcano_plots", sep="/")

# pathway analysis input
dir.create("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis")
dir.create("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis/1.Create_Input/")
pathway_input_dir=paste("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis/1.Create_Input")

# PPI analysis
dir.create("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/5.PPI_Network")
dir.create("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/5.PPI_Network/1.Create_input/")
PPI_input_dir=paste("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/5.PPI_Network/1.Create_input/")

##### LOAD DATA #####

setwd(data_dir)

expression_data<-read.table("MRC-LBB_expression_table_QCd.txt", head=T, as.is=T, sep="\t", check.names = F)
phenotype_data<-read.table("MRC-LBB_phenotype.txt", head=T, as.is=T, sep="\t")


head(expression_data)[1:5]
head(phenotype_data)

dim(expression_data)
dim(phenotype_data)

##### DE #####

#merge pheno+exprs data together
exprs_pheno<-merge(phenotype_data, expression_data, by="row.names")
rownames(exprs_pheno)<-exprs_pheno$Row.names
exprs_pheno$Row.names<-NULL

#create function 

head(exprs_pheno)[1:5]

run_diff_exprs_analysis<-function(phenotype_1, phenotype_2, Tissue, x){
  # subset expression table to disease group + tissue
  dataset<-subset(exprs_pheno, PHENOTYPE==phenotype_1 & TISSUE==Tissue | PHENOTYPE==phenotype_2 & TISSUE==Tissue)
  # change coding of disease to 1 + 0
  dataset[dataset$PHENOTYPE==phenotype_1,1]<-0
  dataset[dataset$PHENOTYPE==phenotype_2,1]<-1
  # sort by disease - want disease 1st. i.e 0
  dataset_sorted<-dataset[order(dataset$PHENOTYPE),]
  #split dataset into expres and diagnosis + gender
  dataset_exprs<-dataset_sorted[5:dim(dataset_sorted)[2]]
  dataset_pheno<-dataset_sorted[1:2]
  #replace sample name to numbers
  #rownames(dataset_exprs)<-c(1:dim(dataset)[1])
  # setup experimental desgin
  design <- model.matrix(~0+PHENOTYPE+Clinical_Gender, data=dataset_pheno)
  # colnames
  colnames(design)<-c("case", "control", "gender")
  #run diff expression
  dataset_exprs_fit <- lmFit(t(dataset_exprs), design, method="robust", maxit=x)
  # make contrast
  dataset_exprs_contrast_matrix<- makeContrasts(case-control, levels=design)
  # extract linear model fit
  dataset_exprs_contrast_fit <-contrasts.fit(dataset_exprs_fit, dataset_exprs_contrast_matrix)
  dataset_exprs_contrast_ebayes <- eBayes(dataset_exprs_contrast_fit, robust=T)
  # DE results
  dataset_exprs_top_genes <- topTable(dataset_exprs_contrast_ebayes, number=(dim(dataset_exprs)[2]), coef=1, adjust.method="fdr", confint=TRUE) 
  #return table
  return(dataset_exprs_top_genes)
}

# apply function

table(exprs_pheno$PHENOTYPE, exprs_pheno$TISSUE)

# AD vs Control - Temporal_cortex

Temporal_cortex_AD_v_CO<-run_diff_exprs_analysis("AD", "CO", "Temporal_Cortex", 100)
head(Temporal_cortex_AD_v_CO)
dim(Temporal_cortex_AD_v_CO)
dim(subset(Temporal_cortex_AD_v_CO, adj.P.Val<=0.05))

# AD vs Control_AD - Temporal_cortex

Temporal_cortex_AD_v_CO_AD<-run_diff_exprs_analysis("AD", "CO_AD", "Temporal_Cortex", 100)
head(Temporal_cortex_AD_v_CO_AD)
dim(Temporal_cortex_AD_v_CO_AD)
dim(subset(Temporal_cortex_AD_v_CO_AD, adj.P.Val<=0.05))

# Control_AD vs AD - Temporal_cortex

Temporal_cortex_CO_AD_v_CO<-run_diff_exprs_analysis("CO_AD", "CO", "Temporal_Cortex", 100)
head(Temporal_cortex_CO_AD_v_CO)
dim(Temporal_cortex_CO_AD_v_CO)
dim(subset(Temporal_cortex_CO_AD_v_CO, adj.P.Val<=0.05))


# AD vs Control - Frontal_Cortex

Frontal_Cortex_AD_v_CO<-run_diff_exprs_analysis("AD", "CO", "Frontal_Cortex", 100)
head(Frontal_Cortex_AD_v_CO)
dim(Frontal_Cortex_AD_v_CO)
dim(subset(Frontal_Cortex_AD_v_CO, adj.P.Val<=0.05))

# AD vs Control_AD - Frontal_Cortex

Frontal_Cortex_AD_v_CO_AD<-run_diff_exprs_analysis("AD", "CO_AD", "Frontal_Cortex", 100)
head(Frontal_Cortex_AD_v_CO_AD)
dim(Frontal_Cortex_AD_v_CO_AD)
dim(subset(Frontal_Cortex_AD_v_CO_AD, adj.P.Val<=0.05))

# Control_AD vs AD - Frontal_Cortex

Frontal_Cortex_CO_AD_v_CO<-run_diff_exprs_analysis("CO_AD", "CO", "Frontal_Cortex", 100)
head(Frontal_Cortex_CO_AD_v_CO)
dim(Frontal_Cortex_CO_AD_v_CO)
dim(subset(Frontal_Cortex_CO_AD_v_CO, adj.P.Val<=0.05))

# AD vs Control - Entorhinal_Cortex

Entorhinal_Cortex_AD_v_CO<-run_diff_exprs_analysis("AD", "CO", "Entorhinal_Cortex", 100)
head(Entorhinal_Cortex_AD_v_CO)
dim(Entorhinal_Cortex_AD_v_CO)
dim(subset(Entorhinal_Cortex_AD_v_CO, adj.P.Val<=0.05))

# AD vs Control_AD - Entorhinal_Cortex

Entorhinal_Cortex_AD_v_CO_AD<-run_diff_exprs_analysis("AD", "CO_AD", "Entorhinal_Cortex", 100)
head(Entorhinal_Cortex_AD_v_CO_AD)
dim(Entorhinal_Cortex_AD_v_CO_AD)
dim(subset(Entorhinal_Cortex_AD_v_CO_AD, adj.P.Val<=0.05))

# Control_AD vs AD - Entorhinal_Cortex

Entorhinal_Cortex_CO_AD_v_CO<-run_diff_exprs_analysis("CO_AD", "CO", "Entorhinal_Cortex", 100)
head(Entorhinal_Cortex_CO_AD_v_CO)
dim(Entorhinal_Cortex_CO_AD_v_CO)
dim(subset(Entorhinal_Cortex_CO_AD_v_CO, adj.P.Val<=0.05))

# AD vs Control - Cerebellum

Cerebellum_AD_v_CO<-run_diff_exprs_analysis("AD", "CO", "Cerebellum", 100)
head(Cerebellum_AD_v_CO)
dim(Cerebellum_AD_v_CO)
dim(subset(Cerebellum_AD_v_CO, adj.P.Val<=0.05))

# AD vs Control_AD - Cerebellum

Cerebellum_AD_v_CO_AD<-run_diff_exprs_analysis("AD", "CO_AD", "Cerebellum", 100)
head(Cerebellum_AD_v_CO_AD)
dim(Cerebellum_AD_v_CO_AD)
dim(subset(Cerebellum_AD_v_CO_AD, adj.P.Val<=0.05))

# Control_AD vs AD - Cerebellum

Cerebellum_CO_AD_v_CO<-run_diff_exprs_analysis("CO_AD", "CO", "Cerebellum", 100)
head(Cerebellum_CO_AD_v_CO)
dim(Cerebellum_CO_AD_v_CO)
dim(subset(Cerebellum_CO_AD_v_CO, adj.P.Val<=0.05))

##### CREATE VOLCANO PLOTS #####

setwd(volcano_plots_dir)

#create function for volcano plot

create_volcano_plot<-function(DE_results, label){
  # Make a basic volcano plot
  with(DE_results, plot(logFC, -log10(adj.P.Val), pch=20, main=label, xlim=c(-5,5), ylim=c(-0.5, 30)))
  # Add colored points: red if adj.P.Val<0.05, orange of log2FC>1, green if both)
  with(subset(DE_results, adj.P.Val<=.05 ), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  #with(subset(DE_results, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="red"))
  with(subset(DE_results, adj.P.Val<=.05 & abs(logFC)>=1), points(logFC, -log10(adj.P.Val), pch=20, col="green"))
  abline(h=-log10(0.05), v=c(1, -1))
  # Label points with the textxy function from the calibrate plot
  with(subset(DE_results, adj.P.Val<=.05 & abs(logFC)>=1), textxy(logFC, -log10(adj.P.Val), labs=rownames(subset(DE_results, adj.P.Val<=.05 & abs(logFC)>=1)), cex=1))
}

#apply function to plot volcano plot

create_volcano_plot(Temporal_cortex_AD_v_CO, "Temporal_cortex_AD_v_CO")
create_volcano_plot(Temporal_cortex_AD_v_CO_AD, "Temporal_cortex_AD_v_CO_AD")
create_volcano_plot(Temporal_cortex_CO_AD_v_CO, "Temporal_cortex_CO_AD_v_CO")

create_volcano_plot(Frontal_Cortex_AD_v_CO, "Frontal_Cortex_AD_v_CO")
create_volcano_plot(Frontal_Cortex_AD_v_CO_AD, "Frontal_Cortex_AD_v_CO_AD")
create_volcano_plot(Frontal_Cortex_CO_AD_v_CO, "Frontal_Cortex_CO_AD_v_CO")

create_volcano_plot(Entorhinal_Cortex_AD_v_CO, "Entorhinal_Cortex_AD_v_CO")
create_volcano_plot(Entorhinal_Cortex_AD_v_CO_AD, "Entorhinal_Cortex_AD_v_CO_AD")
create_volcano_plot(Entorhinal_Cortex_CO_AD_v_CO, "Entorhinal_Cortex_CO_AD_v_CO")

create_volcano_plot(Cerebellum_AD_v_CO, "Cerebellum_AD_v_CO")
create_volcano_plot(Cerebellum_AD_v_CO_AD, "Cerebellum_AD_v_CO_AD")
create_volcano_plot(Cerebellum_CO_AD_v_CO, "Cerebellum_CO_AD_v_CO")

#plot to pdf

pdf("Temporal_cortex_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Temporal_cortex_AD_v_CO, "Temporal_cortex_AD_v_CO")
dev.off()

pdf("Temporal_cortex_AD_v_CO_AD_volcano_plot.pdf")
create_volcano_plot(Temporal_cortex_AD_v_CO_AD, "Temporal_cortex_AD_v_CO_AD")
dev.off()

pdf("Temporal_cortex_CO_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Temporal_cortex_CO_AD_v_CO, "Temporal_cortex_CO_AD_v_CO")
dev.off()

pdf("Frontal_Cortex_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Frontal_Cortex_AD_v_CO, "Frontal_Cortex_AD_v_CO")
dev.off()

pdf("Frontal_Cortex_AD_v_CO_AD_volcano_plot.pdf")
create_volcano_plot(Frontal_Cortex_AD_v_CO_AD, "Frontal_Cortex_AD_v_CO_AD")
dev.off()

pdf("Frontal_Cortex_CO_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Frontal_Cortex_CO_AD_v_CO, "Frontal_Cortex_CO_AD_v_CO")
dev.off()

pdf("Entorhinal_Cortex_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Entorhinal_Cortex_AD_v_CO, "Entorhinal_Cortex_AD_v_CO")
dev.off()

pdf("Entorhinal_Cortex_AD_v_CO_AD_volcano_plot.pdf")
create_volcano_plot(Entorhinal_Cortex_AD_v_CO_AD, "Entorhinal_Cortex_AD_v_CO_AD")
dev.off()

pdf("Entorhinal_Cortex_CO_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Entorhinal_Cortex_CO_AD_v_CO, "Entorhinal_Cortex_CO_AD_v_CO")
dev.off()

pdf("Cerebellum_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Cerebellum_AD_v_CO, "Cerebellum_AD_v_CO")
dev.off()

pdf("Cerebellum_AD_v_CO_AD_volcano_plot.pdf")
create_volcano_plot(Cerebellum_AD_v_CO_AD, "Cerebellum_AD_v_CO_AD")
dev.off()

pdf("Cerebellum_CO_AD_v_CO_volcano_plot.pdf")
create_volcano_plot(Cerebellum_CO_AD_v_CO, "Cerebellum_CO_AD_v_CO")
dev.off()

##### OVERLAP OF DEG ACROSS ALL FOUR BRAIN BRAIN REGIONS #####
# create function to extract consensus genes across all four brain regions

consensus_DEG_across_4<-function(dataset1, dataset2, dataset3, dataset4){
  # AD vs CO
  #merge all sig results
  case_vs_control_combined<-c(rownames(subset(dataset1, adj.P.Val<=0.05)), 
                              rownames(subset(dataset2, adj.P.Val<=0.05)),
                              rownames(subset(dataset3, adj.P.Val<=0.05)),
                              rownames(subset(dataset4, adj.P.Val<=0.05)))
  
  #count number of duplicates
  case_vs_control_combined_count<-as.data.frame(table(case_vs_control_combined))
  # order
  case_vs_control_combined_count<-case_vs_control_combined_count[order(-case_vs_control_combined_count$Freq),]
  # check 
  head(case_vs_control_combined_count)
  # count same gene 
  dim(subset(case_vs_control_combined_count, Freq==4))
  
  print(c("Number of DEG across all four brain regions: ", dim(subset(case_vs_control_combined_count, Freq==4))[1]))
  
  # of which are consensus in logFC direction
  
  dataset1_subset<-(subset(dataset1, rownames(dataset1) %in% subset(case_vs_control_combined_count, Freq==4)[,1]))[1]
  dataset2_subset<-(subset(dataset2, rownames(dataset2) %in% subset(case_vs_control_combined_count, Freq==4)[,1]))[1]
  dataset3_subset<-(subset(dataset3, rownames(dataset3) %in% subset(case_vs_control_combined_count, Freq==4)[,1]))[1]
  dataset4_subset<-(subset(dataset4, rownames(dataset4) %in% subset(case_vs_control_combined_count, Freq==4)[,1]))[1]
  
  #merge 
  
  temp<-merge(dataset1_subset, dataset2_subset, by="row.names")
  temp2<-merge(dataset3_subset, dataset4_subset, by="row.names")
  merged<-merge(temp, temp2, by="Row.names")
  head(merged)
  
  rownames(merged)<-merged$Row.names
  merged$Row.names<-NULL
  
  head(merged)
  colnames(merged)<-c("A", "B", "C", "D")
  
  subset(merged, A>=0 & B>=0 & C>=0 & D>=0)
  subset(merged, A<=0 & B<=0 & C<=0 & D<=0)
  
  results<-list()
  results[[1]]<-subset(merged, A>=0 & B>=0 & C>=0 & D>=0)
  results[[2]]<-subset(merged, A<=0 & B<=0 & C<=0 & D<=0)
  
  print(c("Number consensus up-regulated DEG across all four brain regions: ", dim(results[[1]])[1]))
  print(c("Number consensus down-regulated DEG across all four brain regions: ", dim(results[[2]])[1]))
  
  return(results)
}



AD_vs_CO_across_four_brain_regions<-consensus_DEG_across_4(Temporal_cortex_AD_v_CO,
                                                           Frontal_Cortex_AD_v_CO,
                                                           Entorhinal_Cortex_AD_v_CO,
                                                           Cerebellum_AD_v_CO)

AD_vs_CO_across_four_brain_regions[[1]]
AD_vs_CO_across_four_brain_regions[[2]]

# AD vs CO_AD

AD_vs_CO_AD_across_four_brain_regions<-consensus_DEG_across_4(Temporal_cortex_AD_v_CO_AD,
                                                           Frontal_Cortex_AD_v_CO_AD,
                                                           Entorhinal_Cortex_AD_v_CO_AD,
                                                           Cerebellum_AD_v_CO_AD)

AD_vs_CO_AD_across_four_brain_regions[[1]]
AD_vs_CO_AD_across_four_brain_regions[[2]]

# CO_AD vs CO

CO_AD_vs_CO_across_four_brain_regions<-consensus_DEG_across_4(Temporal_cortex_CO_AD_v_CO,
                                                           Frontal_Cortex_CO_AD_v_CO,
                                                           Entorhinal_Cortex_CO_AD_v_CO,
                                                           Cerebellum_CO_AD_v_CO)

CO_AD_vs_CO_across_four_brain_regions[[1]]
CO_AD_vs_CO_across_four_brain_regions[[2]]

##### WRITE DE RESULTS ####

setwd(work_dir)

write.csv(Temporal_cortex_AD_v_CO, "Temporal_cortex_AD_v_CO_DE_results.txt")
write.csv(Temporal_cortex_AD_v_CO_AD, "Temporal_cortex_AD_v_CO_AD_DE_results.txt")
write.csv(Temporal_cortex_CO_AD_v_CO, "Temporal_cortex_CO_AD_v_CO_DE_results.txt")

write.csv(Frontal_Cortex_AD_v_CO, "Frontal_Cortex_AD_v_CO_DE_results.txt")
write.csv(Frontal_Cortex_AD_v_CO_AD, "Frontal_Cortex_AD_v_CO_AD_DE_results.txt")
write.csv(Frontal_Cortex_CO_AD_v_CO, "Frontal_Cortex_CO_AD_v_CO_DE_results.txt")

write.csv(Entorhinal_Cortex_AD_v_CO, "Entorhinal_Cortex_AD_v_CO_DE_results.txt")
write.csv(Entorhinal_Cortex_AD_v_CO_AD, "Entorhinal_Cortex_AD_v_CO_AD_DE_results.txt")
write.csv(Entorhinal_Cortex_CO_AD_v_CO, "Entorhinal_Cortex_CO_AD_v_CO_DE_results.txt")

write.csv(Cerebellum_AD_v_CO, "Cerebellum_AD_v_CO_DE_results.txt")
write.csv(Cerebellum_AD_v_CO_AD, "Cerebellum_AD_v_CO_AD_DE_results.txt")
write.csv(Cerebellum_CO_AD_v_CO, "Cerebellum_CO_AD_v_CO_DE_results.txt")

##### WRITE SIG DE RESULTS FOR PATHWAY ANALYSIS #####

setwd(pathway_input_dir)

write.table(rownames(subset(Temporal_cortex_AD_v_CO, adj.P.Val<=0.05)), "Temporal_cortex_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Temporal_cortex_AD_v_CO_AD, adj.P.Val<=0.05)), "Temporal_cortex_AD_v_CO_AD_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Temporal_cortex_CO_AD_v_CO, adj.P.Val<=0.05)), "Temporal_cortex_CO_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)

write.table(rownames(subset(Frontal_Cortex_AD_v_CO, adj.P.Val<=0.05)), "Frontal_Cortex_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Frontal_Cortex_AD_v_CO_AD, adj.P.Val<=0.05)), "Frontal_Cortex_AD_v_CO_AD_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Frontal_Cortex_CO_AD_v_CO, adj.P.Val<=0.05)), "Frontal_Cortex_CO_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)

write.table(rownames(subset(Entorhinal_Cortex_AD_v_CO, adj.P.Val<=0.05)), "Entorhinal_Cortex_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Entorhinal_Cortex_AD_v_CO_AD, adj.P.Val<=0.05)), "Entorhinal_Cortex_AD_v_CO_AD_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Entorhinal_Cortex_CO_AD_v_CO, adj.P.Val<=0.05)), "Entorhinal_Cortex_CO_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)

write.table(rownames(subset(Cerebellum_AD_v_CO, adj.P.Val<=0.05)), "Cerebellum_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Cerebellum_AD_v_CO_AD, adj.P.Val<=0.05)), "Cerebellum_AD_v_CO_AD_pathway_input.txt", row.names=F, quote=F, col.names=F)
write.table(rownames(subset(Cerebellum_CO_AD_v_CO, adj.P.Val<=0.05)), "Cerebellum_CO_AD_v_CO_pathway_input.txt", row.names=F, quote=F, col.names=F)

##### WRITE SIG DE RESULTS FOR NETWORK ANALYSIS #####

setwd(PPI_input_dir)

write.table(subset(Temporal_cortex_AD_v_CO, adj.P.Val<=0.05)[1], "Temporal_cortex_AD_v_CO_PPI_input.txt", quote=F, col.names=F)
write.table(subset(Temporal_cortex_AD_v_CO_AD, adj.P.Val<=0.05)[1], "Temporal_cortex_AD_v_CO_AD_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Temporal_cortex_CO_AD_v_CO, adj.P.Val<=0.05)[1], "Temporal_cortex_CO_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)

write.table(subset(Frontal_Cortex_AD_v_CO, adj.P.Val<=0.05)[1], "Frontal_Cortex_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Frontal_Cortex_AD_v_CO_AD, adj.P.Val<=0.05)[1], "Frontal_Cortex_AD_v_CO_AD_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Frontal_Cortex_CO_AD_v_CO, adj.P.Val<=0.05)[1], "Frontal_Cortex_CO_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)

write.table(subset(Entorhinal_Cortex_AD_v_CO, adj.P.Val<=0.05)[1], "Entorhinal_Cortex_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Entorhinal_Cortex_AD_v_CO_AD, adj.P.Val<=0.05)[1], "Entorhinal_Cortex_AD_v_CO_AD_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Entorhinal_Cortex_CO_AD_v_CO, adj.P.Val<=0.05)[1], "Entorhinal_Cortex_CO_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)

write.table(subset(Cerebellum_AD_v_CO, adj.P.Val<=0.05)[1], "Cerebellum_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Cerebellum_AD_v_CO_AD, adj.P.Val<=0.05)[1], "Cerebellum_AD_v_CO_AD_PPI_input.txt",  quote=F, col.names=F)
write.table(subset(Cerebellum_CO_AD_v_CO, adj.P.Val<=0.05)[1], "Cerebellum_CO_AD_v_CO_PPI_input.txt",  quote=F, col.names=F)

##### MERGE RESULTS #####

merge_results<-function(CO_AD_vs_CO, AD_vs_CO_AD, AD_vs_CO, brain_region){
  #subset to sig genes only + keep only logFC + adj p value
  CO_AD_vs_CO<-CO_AD_vs_CO[c(1,7)]
  AD_vs_CO_AD<-AD_vs_CO_AD[c(1,7)]
  AD_vs_CO<-AD_vs_CO[c(1,7)]
  # rename column names
  colnames(CO_AD_vs_CO)<-c(paste(brain_region, "CO_AD_vs_CO_logFC", sep="_"), paste(brain_region, "CO_AD_vs_CO_adj.P.val", sep="_"))
  colnames(AD_vs_CO_AD)<-c(paste(brain_region, "AD_vs_CO_AD_logFC", sep="_"), paste(brain_region, "AD_vs_CO_AD_adj.P.val", sep="_"))
  colnames(AD_vs_CO)<-c(paste(brain_region, "AD_vs_CO_logFC", sep="_"), paste(brain_region, "AD_vs_CO_adj.P.val", sep="_"))
  #merge by rownames
  merged_results<-merge(CO_AD_vs_CO, AD_vs_CO_AD, by="row.names")
  rownames(merged_results)<-merged_results$Row.names
  merged_results$Row.names<-NULL
  #merge again
  merged_results<-merge(merged_results, AD_vs_CO, by="row.names")
  rownames(merged_results)<-merged_results$Row.names
  merged_results$Row.names<-NULL
  return(merged_results)
}

#merge TC results
TC_results<-merge_results(Temporal_cortex_CO_AD_v_CO, Temporal_cortex_AD_v_CO_AD, Temporal_cortex_AD_v_CO, "TC")
head(TC_results)

#merge EC results
EC_results<-merge_results(Entorhinal_Cortex_CO_AD_v_CO, Entorhinal_Cortex_AD_v_CO_AD, Entorhinal_Cortex_AD_v_CO, "EC")
head(EC_results)

#merge FC results
FC_results<-merge_results(Frontal_Cortex_CO_AD_v_CO, Frontal_Cortex_AD_v_CO_AD, Frontal_Cortex_AD_v_CO, "FC")
head(FC_results)

#merge CB results
CB_results<-merge_results(Cerebellum_CO_AD_v_CO, Cerebellum_AD_v_CO_AD, Cerebellum_AD_v_CO, "CB")
head(CB_results)

# merge all brain results

MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names")
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}

DE_results_table <- Reduce(MyMerge, list(TC_results, EC_results, FC_results, CB_results))
head(DE_results_table)

##### WRITE OUT MERGED RESULTS #####

setwd(work_dir)

write.csv(DE_results_table, "All_AD_DE_results_table.txt", quote=F)

##### SAVE IMAGE #####

setwd(work_dir)

save.image("DE_analysis.Rdata")
