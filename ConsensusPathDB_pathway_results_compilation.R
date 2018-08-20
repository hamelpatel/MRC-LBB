########################################################################################################
####                                                                                                ####
###                                                                                                  ###
##                                                                                                    ##
#                             CONSENSUSPATHDB PATHWAYS RESULTS COMPILATION                             #
##                                                                                                    ##
###                                                                                                  ###
####                                                                                                ####
########################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 03/10/2017

##### DESCRIPTION OF ANALYSIS ####
## each brain region sig DEG list queried in ConsensusPathDB with background gene list
## 
## 
## 
## NOTE: 
## 
## 
#####

##### LIBRARIES #####

library(gplots)
library(d3heatmap)
library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)

##### SET PARAMETERS #####

rm(list=ls())

##### DIRECTORIES #####

work_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis/2.ConsensusPathDB_results/"

setwd(work_dir)

##### READ DATA  - ONLY PATHWAY ######

setwd(work_dir)

# files in directory
file_names<-list.files()

# read in files
for (i in 1:length(file_names)) 
  assign(file_names[i], 
         read.table(file_names[i], head=T, sep="\t", quote="", fill=T))

ls()

#general check
head(Cerebellum_AD_v_CO_AD_ORA_results.tab)[1:6]

##### SUBSET TO SIG ONLY #####

# subset to significant data to where q value<=0.05 - (FDR adjusted p value)

for (x in 1:length(file_names)) {
  # subset to FDR sig only and assign a new object name
  assign(paste(gsub(".tab", "_sig_only", file_names[x])), subset(get(file_names[x]), q.value<=0.05))
}

dim(Cerebellum_AD_v_CO_AD_ORA_results_sig_only)
dim(Cerebellum_AD_v_CO_AD_ORA_results.tab)

#dim each sig file
sig_file_names<-gsub(".tab", "_sig_only", file_names)

for (x in 1:length(sig_file_names)) {
  print(sig_file_names[x])
  print(dim(get(sig_file_names[x])))
}

##### ANALYSIS 1: AD vs CO #####

# top 10 sig results

head(Temporal_cortex_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]
head(Frontal_Cortex_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]
head(Entorhinal_Cortex_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]
head(Cerebellum_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]

# common pathways

# merge results - @ 0.05
Analysis_1_all_DEG_across_brain_regions<-as.data.frame(table(rbind(Cerebellum_AD_v_CO_ORA_results_sig_only[3],
                                                                   Temporal_cortex_AD_v_CO_ORA_results_sig_only[3],
                                                                   Entorhinal_Cortex_AD_v_CO_ORA_results_sig_only[3],
                                                                   Frontal_Cortex_AD_v_CO_ORA_results_sig_only[3])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions<-Analysis_1_all_DEG_across_brain_regions[order(-Analysis_1_all_DEG_across_brain_regions$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions)<-1:nrow(Analysis_1_all_DEG_across_brain_regions)

head(Analysis_1_all_DEG_across_brain_regions)

subset(Analysis_1_all_DEG_across_brain_regions, Freq==4)

## repeat with only 3 brain regions - removing cerebellum

# merge results
Analysis_1_all_DEG_across_brain_regions_no_cerebellum<-as.data.frame(table(rbind(Temporal_cortex_AD_v_CO_ORA_results_sig_only[3],
                                                                                 Entorhinal_Cortex_AD_v_CO_ORA_results_sig_only[3],
                                                                                 Frontal_Cortex_AD_v_CO_ORA_results_sig_only[3])))

#arrange by count
Analysis_1_all_DEG_across_brain_regions_no_cerebellum<-Analysis_1_all_DEG_across_brain_regions_no_cerebellum[order(-Analysis_1_all_DEG_across_brain_regions_no_cerebellum$Freq),]
rownames(Analysis_1_all_DEG_across_brain_regions_no_cerebellum)<-1:nrow(Analysis_1_all_DEG_across_brain_regions_no_cerebellum)

head(Analysis_1_all_DEG_across_brain_regions_no_cerebellum)

subset(Analysis_1_all_DEG_across_brain_regions_no_cerebellum, Freq==3)
dim(subset(Analysis_1_all_DEG_across_brain_regions_no_cerebellum, Freq==3))

# check known AD pathways:

#grepping th efollowing pathways from results to see if they exist:

# immune
# tranlsation/transcription
# calcium
# RA
# synapse
# metabloism
# neurotransmitters
# wnt

Cerebellum_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Cerebellum_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Temporal_cortex_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Temporal_cortex_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Entorhinal_Cortex_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Entorhinal_Cortex_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Frontal_Cortex_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Frontal_Cortex_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]

##### ANALYSIS 2: CO_AD vs CO #####

# top 10 sig results

head(Temporal_cortex_CO_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]
head(Frontal_Cortex_CO_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]
head(Entorhinal_Cortex_CO_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]
head(Cerebellum_CO_AD_v_CO_ORA_results_sig_only,10)[c(3,4,2)]

# common pathways

# merge results - @ 0.05
Analysis_2_all_DEG_across_brain_regions<-as.data.frame(table(rbind(#Cerebellum_CO_AD_v_CO_ORA_results_sig_only[3],
                                                                   Temporal_cortex_CO_AD_v_CO_ORA_results_sig_only[3],
                                                                   Entorhinal_Cortex_CO_AD_v_CO_ORA_results_sig_only[3],
                                                                   Frontal_Cortex_CO_AD_v_CO_ORA_results_sig_only[3])))

#arrange by count
Analysis_2_all_DEG_across_brain_regions<-Analysis_2_all_DEG_across_brain_regions[order(-Analysis_2_all_DEG_across_brain_regions$Freq),]
rownames(Analysis_2_all_DEG_across_brain_regions)<-1:nrow(Analysis_2_all_DEG_across_brain_regions)

head(Analysis_2_all_DEG_across_brain_regions)

subset(Analysis_2_all_DEG_across_brain_regions, Freq==4)

# check known AD pathways:

#grepping th efollowing pathways from results to see if they exist:

# immune
# tranlsation/transcription
# calcium
# RA
# synapse
# metabloism
# neurotransmitters
# wnt

# Cerebellum_CO_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Cerebellum_CO_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
# Temporal_cortex_CO_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Temporal_cortex_CO_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
# Entorhinal_Cortex_CO_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Entorhinal_Cortex_CO_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
# Frontal_Cortex_CO_AD_v_CO_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Frontal_Cortex_CO_AD_v_CO_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]

##### ANALYSIS 3: AD vs CO_AD #####

# top 10 sig results

head(Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only,10)[c(3,4,2)]
head(Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only,10)[c(3,4,2)]
head(Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only,10)[c(3,4,2)]
head(Cerebellum_AD_v_CO_AD_ORA_results_sig_only,10)[c(3,4,2)]

# common pathways

# merge results - @ 0.05
Analysis_3_all_DEG_across_brain_regions<-as.data.frame(table(rbind(Cerebellum_AD_v_CO_AD_ORA_results_sig_only[3],
                                                                   Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only[3],
                                                                   Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only[3],
                                                                   Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only[3])))

#arrange by count
Analysis_3_all_DEG_across_brain_regions<-Analysis_3_all_DEG_across_brain_regions[order(-Analysis_3_all_DEG_across_brain_regions$Freq),]
rownames(Analysis_3_all_DEG_across_brain_regions)<-1:nrow(Analysis_3_all_DEG_across_brain_regions)

head(Analysis_3_all_DEG_across_brain_regions)

subset(Analysis_3_all_DEG_across_brain_regions, Freq==4)

## repeat with only 3 brain regions - removing cerebellum

# merge results
Analysis_3_all_DEG_across_brain_regions_no_cerebellum<-as.data.frame(table(rbind(Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only[3],
                                                                                 Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only[3],
                                                                                 Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only[3])))

#arrange by count
Analysis_3_all_DEG_across_brain_regions_no_cerebellum<-Analysis_3_all_DEG_across_brain_regions_no_cerebellum[order(-Analysis_3_all_DEG_across_brain_regions_no_cerebellum$Freq),]
rownames(Analysis_3_all_DEG_across_brain_regions_no_cerebellum)<-1:nrow(Analysis_3_all_DEG_across_brain_regions_no_cerebellum)

head(Analysis_3_all_DEG_across_brain_regions_no_cerebellum)

subset(Analysis_3_all_DEG_across_brain_regions_no_cerebellum, Freq==3)


## repeat using q value of 0.01

Cerebellum_AD_v_CO_AD_ORA_results_sig_only_001<-subset(Cerebellum_AD_v_CO_AD_ORA_results_sig_only, q.value<=0.01)

Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only_001<-subset(Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only, q.value<=0.01)
Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only_001<-subset(Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only, q.value<=0.01)
Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only_001<-subset(Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only, q.value<=0.01)

# merge results
Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001<-as.data.frame(table(rbind(Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only_001[3],
                                                                                     Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only_001[3],
                                                                                     Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only_001[3])))

#arrange by count
Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001<-Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001[order(-Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001$Freq),]
rownames(Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001)<-1:nrow(Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001)

head(Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001)

subset(Analysis_3_all_DEG_across_brain_regions_no_cerebellum_001, Freq==3)

# check known AD pathways:

#grepping th efollowing pathways from results to see if they exist:

# immune
# tranlsation/transcription
# calcium
# RA
# synapse
# metabloism
# neurotransmitters
# wnt

Cerebellum_AD_v_CO_AD_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Cerebellum_AD_v_CO_AD_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Temporal_cortex_AD_v_CO_AD_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Entorhinal_Cortex_AD_v_CO_AD_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]
Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only[grep("\\translation\\b|\\transcription\\b|\\immune\\b|\\calcium\\b|\\RA\\b|\\rheumatoid\\b|\\chemical\\b|\\synapse\\b|\\neuro\\b|\\wnt\\b|\\meta\\b", Frontal_Cortex_AD_v_CO_AD_ORA_results_sig_only[,3],ignore.case=TRUE),][c(3,2)]

##### SAVE IMAGE #####

setwd(work_dir)

save.image("CONSENSUSPATHDB_ORA_RESULTS_COMPILATION.Rdata")
