########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                                       DE RESULT EXPLORATION                                          #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 13/11/2017

##### DESCRIPTION OF ANALYSIS ####
## DE analysis completed on all brain regions. Looking into unique, common + known AD genes
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

rm(list=ls())

##### LOAD LIBRARIES ######

library(org.Hs.eg.db)

library("qqman")

library(ggplot2)

library(RColorBrewer)
library(extrafont)

##### SET DIRECTORIES #####

work_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/3.DE_analysis/2.Result_exploration"

setwd(work_dir)

# create dir
dir.create(paste(work_dir,"manhattan_plot", sep="/"))
manhattan_plot_dir=paste(work_dir,"manhattan_plot", sep="/")

# create dir
dir.create(paste(work_dir,"boxplot_plot", sep="/"))
boxplot_plot_dir=paste(work_dir,"boxplot_plot", sep="/")

##### READ IN DE RESULTS #####

DE_results<-read.csv("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/3.DE_analysis/1.Results/All_AD_DE_results_table.txt", head=T, as.is=T, row.names=1)

head(DE_results)
dim(DE_results)

##### READ IN KNOWN AD GENES #####

knonwn_AD_genes<-read.table("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/8.Known_AD_genes/Alzgene_list.txt")
head(knonwn_AD_genes)
dim(knonwn_AD_genes)

#convert gene symbol to entrez ID

# convert entrez gene id to gene symbol
entrez_gene_symbol <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(entrez_gene_symbol)
# Convert to a list
gene_symbol_lookup_table <- as.data.frame(entrez_gene_symbol[mapped_genes])

# keep 

knonwn_AD_genes_entrez_ID<-gene_symbol_lookup_table[gene_symbol_lookup_table$symbol %in% knonwn_AD_genes$V1,]
head(knonwn_AD_genes_entrez_ID)
dim(knonwn_AD_genes_entrez_ID)

# top 10
Top_10_knonwn_AD_genes<-c("APOE",
                          "BIN1",
                          "CLU",
                          "ABCA7", 
                          "CR1",
                          "PICALM", 
                          "MS4A6A",
                          "CD33",   
                          "MS4A4E", 
                          "CD2AP")

Top_10_knonwn_AD_genes_entrez_ID<-subset(knonwn_AD_genes_entrez_ID, knonwn_AD_genes_entrez_ID$symbol %in% Top_10_knonwn_AD_genes)
Top_10_knonwn_AD_genes_entrez_ID

##### DE RESULTS FOR TC BREAKDOWN #####

# n.o DE results for

#1. CO vs AD - general
dim(subset(DE_results, TC_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_adj.P.val<=0.05))))

#2. AD vs CO_AD - late changes
dim(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05))))

#3. CO_AD vs CO - early changes
dim(subset(DE_results, TC_CO_AD_vs_CO_adj.P.val<=0.05))

# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_CO_AD_vs_CO_adj.P.val<=0.05))))

#4. sig in both CO_AD vs AD + AD vs CO - consistent changes in both
dim(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05 & TC_CO_AD_vs_CO_adj.P.val<=0.05))

# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05 & TC_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05 & TC_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05 & TC_CO_AD_vs_CO_adj.P.val<=0.05))))

subset(DE_results, rownames(DE_results) %in% (subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, TC_AD_vs_CO_AD_adj.P.val<=0.05 & TC_CO_AD_vs_CO_adj.P.val<=0.05))))$gene_id)[1:6]

##### DE RESULTS FOR EC BREAKDOWN #####

# n.o DE results for

#1. CO vs AD - general
dim(subset(DE_results, EC_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_adj.P.val<=0.05))))
subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_adj.P.val<=0.05)))

#2. AD vs CO_AD - late changes
dim(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05))))

#3. CO_AD vs CO - early changes
dim(subset(DE_results, EC_CO_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_CO_AD_vs_CO_adj.P.val<=0.05))))

#4. sig in both CO_AD vs AD + AD vs CO - consistent changes in both
dim(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05 & EC_CO_AD_vs_CO_adj.P.val<=0.05))


# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05 & EC_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05 & EC_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05 & EC_CO_AD_vs_CO_adj.P.val<=0.05))))


# unique known AD genes across all 

length(unique(rbind(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_adj.P.val<=0.05))),
                    subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_AD_vs_CO_AD_adj.P.val<=0.05))),
                    subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, EC_CO_AD_vs_CO_adj.P.val<=0.05)))))$gene_id)

##### DE RESULTS FOR FC BREAKDOWN #####

# n.o DE results for

#1. CO vs AD - general
dim(subset(DE_results, FC_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_adj.P.val<=0.05))))
subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_adj.P.val<=0.05)))

#2. AD vs CO_AD - late changes
dim(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05))))

#3. CO_AD vs CO - early changes
dim(subset(DE_results, FC_CO_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_CO_AD_vs_CO_adj.P.val<=0.05))))

#4. sig in both CO_AD vs AD + AD vs CO - consistent changes in both
dim(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05 & FC_CO_AD_vs_CO_adj.P.val<=0.05))

# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05 & FC_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05 & FC_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05 & FC_CO_AD_vs_CO_adj.P.val<=0.05))))


# unique known AD genes across all 

length(unique(rbind(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_adj.P.val<=0.05))),
                    subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_AD_vs_CO_AD_adj.P.val<=0.05))),
                    subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, FC_CO_AD_vs_CO_adj.P.val<=0.05)))))$gene_id)

##### DE RESULTS FOR CB BREAKDOWN #####

# n.o DE results for

#1. CO vs AD - general
dim(subset(DE_results, CB_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_adj.P.val<=0.05))))
subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_adj.P.val<=0.05)))

#2. AD vs CO_AD - late changes
dim(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05))))

#3. CO_AD vs CO - early changes
dim(subset(DE_results, CB_CO_AD_vs_CO_adj.P.val<=0.05))
# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_CO_AD_vs_CO_adj.P.val<=0.05))))

#4. sig in both CO_AD vs AD + AD vs CO - consistent changes in both
dim(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05 & CB_CO_AD_vs_CO_adj.P.val<=0.05))

# known AD
subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05 & CB_CO_AD_vs_CO_adj.P.val<=0.05)))
dim(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05 & CB_CO_AD_vs_CO_adj.P.val<=0.05))))
dim(subset(Top_10_knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05 & CB_CO_AD_vs_CO_adj.P.val<=0.05))))

# unique known AD genes across all 

length(unique(rbind(subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_adj.P.val<=0.05))),
                    subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_AD_vs_CO_AD_adj.P.val<=0.05))),
                    subset(knonwn_AD_genes_entrez_ID, gene_id %in% rownames(subset(DE_results, CB_CO_AD_vs_CO_adj.P.val<=0.05)))))$gene_id)

##### BOXPLOT OF LOGFC #####

# boxplot of log FC
par(mar=c(6, 4.1, 4.1, 2.1))

# x-axis labels
logFC_labels<-colnames(DE_results[c(seq(1,23,2))])

#remove "_logFC"
logFC_labels<-unlist(strsplit(logFC_labels, "_logFC"))

# boxplot
boxplot(DE_results[c(seq(1,23,2))],
        xaxt= "n",
        main="Boxplot of DE logFC results across brain", 
        ylab="LogFC")

# x axis with ticks but without labels
axis(1, at=seq(1, length(logFC_labels), 1), labels = FALSE)

text(x =  seq_along(logFC_labels), 
     y = par("usr")[3]-0.1, 
     srt = 45, 
     adj = 1,
     labels = logFC_labels, 
     xpd = TRUE)

setwd(work_dir)

pdf("boxplot_of_logFC_across_phenotype.pdf")
# boxplot of log FC
par(mar=c(6, 4.1, 4.1, 2.1))

# boxplot
boxplot(DE_results[c(seq(1,23,2))],
        xaxt= "n",
        main="Boxplot of DE logFC results across brain", 
        ylab="LogFC")

# x axis with ticks but without labels
axis(1, at=seq(1, length(logFC_labels), 1), labels = FALSE)

text(x =  seq_along(logFC_labels), 
     y = par("usr")[3]-0.1, 
     srt = 45, 
     adj = 1,
     labels = logFC_labels, 
     xpd = TRUE)

dev.off()

##### BOXPLOT OF P VAL ######

# boxplot of adj.p.value
par(mar=c(6, 4.1, 4.1, 2.1))

# x-axis labels
P_val_labels<-colnames(DE_results[c(seq(2,24,2))])

#remove "_P_val"
P_val_labels<-unlist(strsplit(P_val_labels, "_adj.P.val"))

# boxplot
boxplot(-log10(DE_results[c(seq(1,23,2))]),
        xaxt= "n",
        main="Boxplot of DE adj.P.value results across brain", 
        ylab="-log10")

# x axis with ticks but without labels
axis(1, at=seq(1, length(P_val_labels), 1), labels = FALSE)

text(x =  seq_along(P_val_labels), 
     y = par("usr")[3]-0.1, 
     srt = 45, 
     adj = 1,
     labels = P_val_labels, 
     xpd = TRUE)

abline(h=-log10(0.05), col="red")


#plot

setwd(work_dir)

pdf("boxplot_of_adj.p.val_across_phenotype.pdf")
# boxplot of adj.p.value
par(mar=c(6, 4.1, 4.1, 2.1))

# boxplot
boxplot(-log10(DE_results[c(seq(1,23,2))]),
        xaxt= "n",
        main="Boxplot of DE adj.P.value results across brain", 
        ylab="-log10")

# x axis with ticks but without labels
axis(1, at=seq(1, length(P_val_labels), 1), labels = FALSE)

text(x =  seq_along(P_val_labels), 
     y = par("usr")[3]-0.1, 
     srt = 45, 
     adj = 1,
     labels = P_val_labels, 
     xpd = TRUE)

abline(h=-log10(0.05), col="red")

dev.off()

##### MANHATTAN PLOT #####

head(DE_results)

# names of genes
Entrez_ID<-rownames(DE_results)

# Get the entrez gene identifiers that are mapped to chromosome locations
entrez_chromo <- as.list(org.Hs.egCHRLOC[mappedkeys(org.Hs.egCHRLOC)])

head(entrez_chromo)

#empty dataframe
entrez_chromo_df<-data.frame(matrix(nrow=length(entrez_chromo), ncol=3))
colnames(entrez_chromo_df)<-c("Entrez_ID", "Chromosome", "position")

# convert to dataframe
for (x in 1:length(entrez_chromo)) {
  #entrez ID
  entrez_chromo_df[x,1]<-names(entrez_chromo)[[x]]
  #chromosoe
  entrez_chromo_df[x,2]<-names(entrez_chromo[[x]])[1]
  #position
  entrez_chromo_df[x,3]<-as.numeric(entrez_chromo[[x]][1])
}

#check
head(entrez_chromo_df)

# subset to ones in data
entrez_chromo_df<-entrez_chromo_df[entrez_chromo_df$Entrez_ID %in% Entrez_ID,]
dim(entrez_chromo_df)

length(Entrez_ID)

#check any duplicated
anyDuplicated(entrez_chromo_df$Entrez_ID)

#change X to 22 + Y to 23
table(entrez_chromo_df$Chromosome)

entrez_chromo_df[entrez_chromo_df$Chromosome=="X",2]<-22
entrez_chromo_df[entrez_chromo_df$Chromosome=="Y",2]<-23

#keep only chro 1:23
entrez_chromo_df<-entrez_chromo_df[entrez_chromo_df$Chromosome %in% seq(1,23,1),]

#check
table(entrez_chromo_df$Chromosome)

dim(entrez_chromo_df)
length(Entrez_ID)

#merge with DE
rownames(entrez_chromo_df)<-entrez_chromo_df$Entrez_ID
entrez_chromo_df$Entrez_ID<-NULL

DE_results_with_chromo<-merge(entrez_chromo_df, DE_results, by="row.names")
colnames(DE_results_with_chromo)[1]<-"Entrez_ID"

head(DE_results_with_chromo)
dim(DE_results_with_chromo)

# add gene symbol

head(gene_symbol_lookup_table)
colnames(gene_symbol_lookup_table)[1]<-"Entrez_ID"

#merge
DE_results_with_chromo<-merge(gene_symbol_lookup_table, DE_results_with_chromo, by="Entrez_ID")
dim(DE_results_with_chromo)

# order by chromo then pos
DE_results_with_chromo<-DE_results_with_chromo[order(DE_results_with_chromo$Chromosome, DE_results_with_chromo$position),]

DE_results_with_chromo$position<-abs(DE_results_with_chromo$position)


#plot

plot_manhattan<-function(dataset, column_to_extract, title){
  #column to extract
  column<-grep(column_to_extract, colnames(dataset))
  # extract entrez id, chr, position , and p value and rename
  dataset_subset<-dataset[c(2,3,4,column)]
  colnames(dataset_subset)<-c("SNP", "CHR", "BP", "P")
  #make numeric
  dataset_subset$CHR<-as.numeric(dataset_subset$CHR)
  #plot
  manhattan(dataset_subset, 
            main=title,
            annotatePval = 0.005, 
            annotateTop = TRUE, 
            suggestiveline = 
              ,
            genomewideline = -log10(0.05),
            col = c("blue4", "orange3"))
}


x11()
plot_manhattan(DE_results_with_chromo, "TC_CO_AD_vs_CO_adj.P.val", "Temporal Cortex CO_AD vs CO")

plot_manhattan(DE_results_with_chromo, "TC_AD_vs_CO_AD_adj.P.val", "Temporal CortexAD vs CO_AD")

plot_manhattan(DE_results_with_chromo, "TC_AD_vs_CO_adj.P.val", "Temporal Cortex AD vs CO")

plot_manhattan(DE_results_with_chromo, "EC_CO_AD_vs_CO_adj.P.val", "Entorhinal Cortex CO_AD vs CO")

plot_manhattan(DE_results_with_chromo, "EC_AD_vs_CO_AD_adj.P.val", "Entorhinal CortexAD vs CO_AD")

plot_manhattan(DE_results_with_chromo, "EC_AD_vs_CO_adj.P.val", "Entorhinal Cortex AD vs CO")

plot_manhattan(DE_results_with_chromo, "FC_CO_AD_vs_CO_adj.P.val", "Frontal Cortex CO_AD vs CO")

plot_manhattan(DE_results_with_chromo, "FC_AD_vs_CO_AD_adj.P.val", "Frontal CortexAD vs CO_AD")

plot_manhattan(DE_results_with_chromo, "FC_AD_vs_CO_adj.P.val", "Frontal Cortex AD vs CO")

plot_manhattan(DE_results_with_chromo, "CB_CO_AD_vs_CO_adj.P.val", "Cerebellum CO_AD vs CO")

plot_manhattan(DE_results_with_chromo, "CB_AD_vs_CO_AD_adj.P.val", "CerebellumAD vs CO_AD")

plot_manhattan(DE_results_with_chromo, "CB_AD_vs_CO_adj.P.val", "Cerebellum AD vs CO")

#plot to pdf
setwd(manhattan_plot_dir)


pdf("Temporal_Cortex_CO_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "TC_CO_AD_vs_CO_adj.P.val", "Temporal Cortex CO_AD vs CO")
dev.off()

pdf("Temporal_Cortex_AD_vs_CO_AD_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "TC_AD_vs_CO_AD_adj.P.val", "Temporal Corte xAD vs CO_AD")
dev.off()

pdf("Temporal_Cortex_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "TC_AD_vs_CO_adj.P.val", "Temporal Cortex AD vs CO")
dev.off()

pdf("Entorhinal_Cortex_CO_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "EC_CO_AD_vs_CO_adj.P.val", "Early stage AD")
dev.off()

pdf("Entorhinal_Cortex_AD_vs_CO_AD_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "EC_AD_vs_CO_AD_adj.P.val", "Late stage AD")
dev.off()

pdf("Entorhinal_Cortex_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "EC_AD_vs_CO_adj.P.val", "Entorhinal Cortex AD vs CO")
dev.off()

pdf("Frontal_Cortex_CO_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "FC_CO_AD_vs_CO_adj.P.val", "Frontal Cortex CO_AD vs CO")
dev.off()

pdf("Frontal_Cortex_AD_vs_CO_AD_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "FC_AD_vs_CO_AD_adj.P.val", "Frontal Cortex AD vs CO_AD")
dev.off()

pdf("Frontal_Cortex_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "FC_AD_vs_CO_adj.P.val", "Frontal Cortex AD vs CO")
dev.off()

pdf("Cerebellum_CO_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "CB_CO_AD_vs_CO_adj.P.val", "Cerebellum CO_AD vs CO")
dev.off()

pdf("Cerebellum_AD_vs_CO_AD_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "CB_AD_vs_CO_AD_adj.P.val", "Cerebellum AD vs CO_AD")
dev.off()

pdf("Cerebellum_AD_vs_CO_manhattan.pdf")
plot_manhattan(DE_results_with_chromo, "CB_AD_vs_CO_adj.P.val", "Cerebellum AD vs CO")
dev.off()

##### BOXPLOT OF INTERESTING GENES ######

# read in QC'd data

exprs<-read.table("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/1.Data/2.Process_data/clean_data/MRC-LBB_expression_table_QCd.txt", as.is=T, check.names = F)
pheno<-read.table("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/1.Data/2.Process_data/clean_data/MRC-LBB_phenotype.txt", as.is=T)

head(exprs)[1:5]
head(pheno)
dim(exprs)
dim(pheno)

exprs_with_ph<-merge(pheno[c(1,4)], exprs, by="row.names")

rownames(exprs_with_ph)<-exprs_with_ph$Row.names
exprs_with_ph$Row.names<-NULL

head(exprs_with_ph)[1:5]

# change brain name

exprs_with_ph[exprs_with_ph$TISSUE=="Temporal_Cortex", 2]<-"TC"
exprs_with_ph[exprs_with_ph$TISSUE=="Entorhinal_Cortex", 2]<-"EC"
exprs_with_ph[exprs_with_ph$TISSUE=="Frontal_Cortex", 2]<-"FC"
exprs_with_ph[exprs_with_ph$TISSUE=="Cerebellum", 2]<-"CB"



# boxplot of MOSPD3

# create boxplot function

create_boxplot<-function(gene_to_extract, gene_name, label){
  #extract gene of interest
  Gene_of_interest<-exprs_with_ph[c(1,2, grep(gene_to_extract, names(exprs_with_ph)))]
  #change colnames
  colnames(Gene_of_interest)[3]<-"Gene"
  colnames(Gene_of_interest)[1]<-"Phenotype"
  # change phenotype CO_AD to Pre-MCI
  Gene_of_interest[Gene_of_interest$Phenotype=="CO_AD", 1]<-"AsymAD"
  #change order
  Gene_of_interest$Phenotype<-factor(Gene_of_interest$Phenotype, levels= c("CO", "AsymAD", "AD"))
  Gene_of_interest$TISSUE<-factor(Gene_of_interest$TISSUE, levels= c("EC", "TC", "FC", "CB"))
  # boxplot
  ggplot(Gene_of_interest, aes(x = TISSUE, y= Gene, fill = Phenotype)) + 
    stat_boxplot(geom ='errorbar') +
    geom_boxplot() +
    labs(x="Brain Region", y = "Log2(expression)") +
    #ggtitle(paste(gene_name)) +
    labs(title = gene_name, tag = label) +
              theme(plot.title = element_text(size= 14, hjust = 0.5, face="bold")) +
    scale_fill_manual(values=rep(c("palegreen", "paleturquoise", "slateblue1"), 4))
}

create_boxplot("64598", "MOSPD3", "A)")
create_boxplot("10577", "NPC2")
create_boxplot("351", "APP")
create_boxplot("1843", "DUSP1")
create_boxplot("2824", "GPM6B")
create_boxplot("54407", "SLC38A2")
create_boxplot("63926", "ANKEF1")
create_boxplot("348", "APOE")
create_boxplot("274", "BIN1")
create_boxplot("9865", "TRIL")

setwd(boxplot_plot_dir)

# pdf("MOSPD3.pdf")
# create_boxplot("64598", "MOSPD3")
# dev.off()
# 
# pdf("NPC2.pdf")
# create_boxplot("10577", "NPC2")
# dev.off()
# 
# pdf("APP.pdf")
# create_boxplot("351", "APP")
# dev.off()
# 
# pdf("DUSP1.pdf")
# create_boxplot("1843", "DUSP1")
# dev.off()
# 
# pdf("GPM6B.pdf")
# create_boxplot("2824", "GPM6B")
# dev.off()
# 
# pdf("SLC38A2.pdf")
# create_boxplot("54407", "SLC38A2")
# dev.off()
# 
# pdf("ANKEF1.pdf")
# create_boxplot("63926", "ANKEF1")
# dev.off()
# 
# pdf("APOE.pdf")
# create_boxplot("348", "APOE")
# dev.off()
# 
# pdf("BIN1.pdf")
# create_boxplot("274", "BIN1")
# dev.off()
# 
# pdf("TRIL.pdf")
# create_boxplot("9865", "TRIL")
# dev.off()



# TIFF format

tiff("MOSPD3.tiff", width = 8, height = 10, unit="in",  res=300)
create_boxplot("64598", "MOSPD3")
dev.off()

tiff("NPC2.tiff", width = 8, height = 10, unit="in",  res=300)
create_boxplot("10577", "NPC2")
dev.off()

tiff("TRIL.tiff", width = 8, height = 10, unit="in",  res=300)
create_boxplot("9865", "TRIL")
dev.off()

tiff("TRIL.tiff", width = 8, height = 10, unit="in",  res=300)


CairoTIFF("APP.tiff")
create_boxplot("351", "APP")
dev.off()

CairoTIFF("DUSP1.tiff")
create_boxplot("1843", "DUSP1")
dev.off()

CairoTIFF("GPM6B.tiff")
create_boxplot("2824", "GPM6B")
dev.off()

CairoTIFF("SLC38A2.tiff")
create_boxplot("54407", "SLC38A2")
dev.off()

CairoTIFF("ANKEF1.tiff")
create_boxplot("63926", "ANKEF1")
dev.off()

CairoTIFF("APOE.tiff")
create_boxplot("348", "APOE")
dev.off()

##### MULTIPLOT #####

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


setwd(boxplot_plot_dir)


tiff("Gene_boxplots.tiff", width = 20, height = 10, unit="in",  res=300)
multiplot(create_boxplot("9865", "TRIL", "A)"),
          create_boxplot("64598", "MOSPD3", "B)"), 
          create_boxplot("10577", "NPC2", "C)"),
          cols=3)
dev.off()





##### SAVE IMAGE #####

setwd(work_dir)

save.image("DE_result_exploration.Rdata")

