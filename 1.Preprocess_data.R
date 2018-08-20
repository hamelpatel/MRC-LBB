
##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                               PREPROCESS DATA - PIEPLINE V02                                                           #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

## AUTHOR: HAMEL PATEL
## EMAIL: Hamel.patel@kcl.ac.uk
## DATE: 29/09/2017


##### DESCRIPTION OF ANALYSIS ####
# REPROCESSING DATA USING NEW PIPELINE
# Added removal of duplicate samples with low RIN number
#
#
#
#####

rm(list=ls())

##### LOAD LIBRARIES #####

library(affy)
library(WGCNA)
library(illuminaHumanv4.db) 
library(sva)
library(reshape)
library(massiR)
library(ggplot2)
library(gridExtra)

##### SET DIRECTORY #####

work_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/1.Data/2.Process_data/"
data_dir="/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/1.Data/1.Raw_data"

setwd(work_dir)

#create directories

dir.create(paste(work_dir,"pca_plots", sep="/"))
pca_dir=paste(work_dir,"pca_plots", sep="/")

dir.create(paste(work_dir,"clean_data", sep="/"))
clean_data_dir=paste(work_dir,"clean_data", sep="/")

dir.create(paste(work_dir,"sample_network_plots", sep="/"))
sample_network_dir=paste(work_dir,"sample_network_plots", sep="/")

dir.create(paste(work_dir,"Initial_box_and_density_plots", sep="/"))
initial_plots=paste(work_dir,"Initial_box_and_density_plots", sep="/")

dir.create(paste(work_dir,"Gender_specific_gene_plots", sep="/"))
Gender_plots_dir=paste(work_dir,"Gender_specific_gene_plots", sep="/")

dir.create(paste(work_dir,"Probe_detection_threshold_test", sep="/"))
Probe_det_plots_dir=paste(work_dir,"Probe_detection_threshold_test", sep="/")

##### LOAD DATA #####

setwd(data_dir)

load("mrc_brains.eset_bg_log2_rsn_combat_adj.RData")

ls()

setwd(work_dir)

# data is in eset_bg_log2_rsn_combat_adj

brain_data_normalised_as_data_frame<-as.data.frame(t(exprs(eset_bg_log2_rsn_combat_adj)))
phenotype_data<-pData(eset_bg_log2_rsn_combat_adj)

head(brain_data_normalised_as_data_frame)[1:5]
head(phenotype_data)
names(phenotype_data)

dim(brain_data_normalised_as_data_frame)
dim(phenotype_data)

##### INITIAL PLOTS ####

setwd(initial_plots)

boxplot(t(brain_data_normalised_as_data_frame))
pdf(file="raw_data_boxplot.pdf")
boxplot(t(brain_data_normalised_as_data_frame))
dev.off()

plotDensity(t(brain_data_normalised_as_data_frame))
pdf(file="raw_data_density_plot.pdf")
plotDensity(t(brain_data_normalised_as_data_frame))
dev.off()

setwd(work_dir)

##### CHECK FOR DUPLICATE SAMPLES IDs #####

#check for duplicate samples
anyDuplicated(rownames(brain_data_normalised_as_data_frame))

#check for duplicate probe ID
anyDuplicated(colnames(brain_data_normalised_as_data_frame))

##### IDENTIFY DUPLIACTE PATIENTS #####

# table of study ID vs brain region = any 2 are dupliactes
duplicate_study_ID_table<-as.data.frame(table(phenotype_data$Study_ID, phenotype_data$Brain_Region))

# extract duplicates
duplicate_study_ID_table<-subset(duplicate_study_ID_table, Freq>=2)

# create empty list to hold bad sample idS
duplicate_samples_to_remove<-character()

# create loop to extract sample ID, and keep only one where RIN number is highest

for (x in 1:dim(duplicate_study_ID_table)[1]) {
  # extract 1st sample
  duplicate<-subset(phenotype_data, Brain_Region==duplicate_study_ID_table[x,2] & Study_ID==duplicate_study_ID_table[x,1])
  #susbet to RIN number + probe detected
  duplicate<-duplicate[c(6,43)]
  # sort by RIN and then probe detected
  duplicate<-duplicate[order(-duplicate[1], -duplicate[2] ),]
  #write out bad sample ID
  duplicate_samples_to_remove<-c(duplicate_samples_to_remove, rownames(duplicate)[2])
}

duplicate_samples_to_remove

subset(phenotype_data, Brain_Region=="Cerebellum" & Study_ID=="A007_10")

# new exprs and pheno table without these IDs

phenotype_data_no_dup<-phenotype_data[!(rownames(phenotype_data) %in% duplicate_samples_to_remove),]
brain_data_normalised_as_data_frame_no_dup<-brain_data_normalised_as_data_frame[!(rownames(brain_data_normalised_as_data_frame) %in% duplicate_samples_to_remove),]

dim(brain_data_normalised_as_data_frame_no_dup)

##### GENDER CHECK ##### 

# Tried to convert probe ID from massiR to entrez and then back to affy U133a chip - less probes hit.
# Using directly affy_hg_u133_plus_2 probe list
#
#

head(phenotype_data_no_dup)
names(phenotype_data_no_dup)

# gender_information
gender_info<-phenotype_data_no_dup[64]
head(gender_info)
table(gender_info)

#subset male/female
female_samples<-subset(gender_info, clinical_gender=="FEMALE")

male_samples<-subset(gender_info, clinical_gender=="MALE")

head(female_samples)
head(male_samples)

head(brain_data_normalised_as_data_frame_no_dup)[1:5]
dim(brain_data_normalised_as_data_frame_no_dup)

#convert nuid back to illumina id

head(anno_probes)
dim(anno_probes)

nuID_ilmID<-anno_probes[c(1,2)]

head(nuID_ilmID)
nuID_ilmID<-unique(nuID_ilmID)
rownames(nuID_ilmID)<- nuID_ilmID$nuID
nuID_ilmID$nuID<-NULL

brain_data_normalised_as_data_frame_no_dup_illID<-merge(nuID_ilmID, t(brain_data_normalised_as_data_frame_no_dup), by="row.names")

dim(brain_data_normalised_as_data_frame_no_dup)
dim(brain_data_normalised_as_data_frame_no_dup_illID)

head(brain_data_normalised_as_data_frame_no_dup_illID)[1:5]

rownames(brain_data_normalised_as_data_frame_no_dup_illID)<-brain_data_normalised_as_data_frame_no_dup_illID$PROBE_ID
brain_data_normalised_as_data_frame_no_dup_illID$Row.names<-NULL
brain_data_normalised_as_data_frame_no_dup_illID$PROBE_ID<-NULL

head(brain_data_normalised_as_data_frame_no_dup_illID)[1:5]

# get Y choromosome genes
data(y.probes)
names(y.probes)

y_chromo_probes <- data.frame(y.probes["illumina_humanht_12"])

# split data by tissue

Cerebellum_exprs<-brain_data_normalised_as_data_frame_no_dup_illID[,colnames(brain_data_normalised_as_data_frame_no_dup_illID) %in% rownames(subset(phenotype_data_no_dup, Brain_Region=="Cerebellum"))]
Entorhinal_Cortex<-brain_data_normalised_as_data_frame_no_dup_illID[,colnames(brain_data_normalised_as_data_frame_no_dup_illID) %in% rownames(subset(phenotype_data_no_dup, Brain_Region=="Entorhinal_Cortex"))]
Frontal_Cortex<-brain_data_normalised_as_data_frame_no_dup_illID[,colnames(brain_data_normalised_as_data_frame_no_dup_illID) %in% rownames(subset(phenotype_data_no_dup, Brain_Region=="Frontal_Cortex"))]
Temporal_Cortex<-brain_data_normalised_as_data_frame_no_dup_illID[,colnames(brain_data_normalised_as_data_frame_no_dup_illID) %in% rownames(subset(phenotype_data_no_dup, Brain_Region=="Temporal_Cortex"))]

dim(Cerebellum_exprs)
dim(Entorhinal_Cortex)
dim(Frontal_Cortex)
dim(Temporal_Cortex)

# function to predict

predict_gender<-function(dataset){
  # extract Y chromosome genes from dataset
  eset.select.out <- massi_select(dataset, y_chromo_probes)
  #dip test
  dip.result <- massi_dip(eset.select.out)
  # run gender predict
  eset.results <- massi_cluster(eset.select.out)
  #extract gender prediction
  predicted_gender<-(eset.results$massi.results)[c(1,5)]
  rownames(predicted_gender)<-predicted_gender$ID
  predicted_gender$ID<-NULL
  colnames(predicted_gender)<-"Predicted_Gender"
  head(predicted_gender)
  return(predicted_gender)
}

Cerebellum_predicted_gender<-predict_gender(Cerebellum_exprs)
Entorhinal_Cortex_predicted_gender<-predict_gender(Entorhinal_Cortex)
Frontal_Cortex_predicted_gender<-predict_gender(Frontal_Cortex)
Temporal_Cortex_predicted_gender<-predict_gender(Temporal_Cortex)

# merge predicted genders
predicted_gender<-rbind(Cerebellum_predicted_gender,
                        Entorhinal_Cortex_predicted_gender,
                        Frontal_Cortex_predicted_gender,
                        Temporal_Cortex_predicted_gender)

#compare to clinical Gender
#standardise

gender_info$clinical_gender<-as.character(gender_info$clinical_gender)
gender_info[gender_info$clinical_gender=="MALE",]<-"male"
gender_info[gender_info$clinical_gender=="FEMALE",]<-"female"

#merge
gender_comparison<-merge(gender_info, predicted_gender, by="row.names")
rownames(gender_comparison)<-gender_comparison$Row.names
gender_comparison$Row.names<-NULL
colnames(gender_comparison)<-c("Clinical_Gender", "Predicted_Gender")
head(gender_comparison)

#differences

Gender_missmatch<-gender_comparison[which(gender_comparison$Clinical_Gender!=gender_comparison$Predicted_Gender),]
Gender_missmatch
dim(Gender_missmatch)

# look at missmatch 

table(phenotype_data[25][rownames(Gender_missmatch),])

gender_check_results<-eset.results$massi.results

gender_check_results[gender_check_results$ID %in% rownames(Gender_missmatch),]

# XIST gene of gender missmatches

phenotype_data_no_dup[rownames(Gender_missmatch),]

# remove gender missmatches if missmatch across all 

#brain_data_gender_missmatch_removed<-brain_data_normalised_as_data_frame_no_dup_illID[!colnames(brain_data_normalised_as_data_frame_no_dup_illID) %in% rownames(Gender_missmatch)]

# extract missmatch study ID + count

missmatch1<-as.data.frame(sort(table(phenotype_data_no_dup[rownames(Gender_missmatch),]$Study_ID)))

# extract number of smaples from study ID with missmatch
missmatch2<-as.data.frame(sort(table(subset(phenotype_data_no_dup, phenotype_data_no_dup$Study_ID %in% phenotype_data_no_dup[rownames(Gender_missmatch),]$Study_ID)$Study_ID)))

# any missmtach in gender found across all samples for same brain region?
Study_ID_with_gender_issues<-merge(missmatch1, missmatch2, by="row.names")
rownames(Study_ID_with_gender_issues)<-Study_ID_with_gender_issues$Row.names
Study_ID_with_gender_issues$Row.names<-NULL

colnames(Study_ID_with_gender_issues)<-c("V1", "V2")

Study_ID_gender_issues_to_remove<-rownames(Study_ID_with_gender_issues[Study_ID_with_gender_issues$V1==Study_ID_with_gender_issues$V2,])

Study_ID_gender_issues_to_remove

#get sample IDs

gendermissmatch_to_remove<-rownames(phenotype_data_no_dup[phenotype_data_no_dup$Study_ID %in% Study_ID_gender_issues_to_remove,])

# new exprs and pheno table without these IDs

phenotype_data_good_samples<-phenotype_data_no_dup[!(rownames(phenotype_data_no_dup) %in% gendermissmatch_to_remove),]
brain_data_normalised_as_data_frame_good_samples<-brain_data_normalised_as_data_frame_no_dup_illID[,!(colnames(brain_data_normalised_as_data_frame_no_dup_illID) %in% gendermissmatch_to_remove),]

dim(brain_data_normalised_as_data_frame_good_samples)
dim(phenotype_data_good_samples)

##### DIAGNOSIS OF SAMPLES #####

Diagnosis<-phenotype_data_good_samples[25]
table(Diagnosis)

Tissue<-phenotype_data_good_samples[24]
table(Tissue)

Diagnosis_case<-subset(Diagnosis, PHENOTYPE=="AD")

Diagnosis_control<-subset(Diagnosis, PHENOTYPE=="CO")

Diagnosis_control_case<-subset(Diagnosis, PHENOTYPE=="CO_AD")

Entorhnal_Cortex<-subset(Tissue, TISSUE=="Entorhinal_Cortex")
Frontal_Cortex<-subset(Tissue, TISSUE=="Frontal_Cortex")
Temporal_Cortex<-subset(Tissue, TISSUE=="Temporal_Cortex")
Cerebellum<-subset(Tissue, TISSUE=="Cerebellum")

Entorhnal_Cortex
Frontal_Cortex
Temporal_Cortex
Cerebellum

##### PROBE ID DETECTION #####

# using dataframe before probe removal

# separate case control 

case_ID<-rownames(Diagnosis_case)

control_ID<-rownames(Diagnosis_control)

control_case_ID<-rownames(Diagnosis_control_case)

#brain_data_normalised_as_data_frame_good_samples<-brain_data

head(brain_data_normalised_as_data_frame_good_samples)[1:5]
dim(brain_data_normalised_as_data_frame_good_samples)

case_exprs<-as.data.frame(t(brain_data_normalised_as_data_frame_good_samples[,colnames(brain_data_normalised_as_data_frame_good_samples)%in%case_ID]))
control_exprs<-as.data.frame(t(brain_data_normalised_as_data_frame_good_samples[,colnames(brain_data_normalised_as_data_frame_good_samples)%in%control_ID]))
control_case_exprs<-as.data.frame(t(brain_data_normalised_as_data_frame_good_samples[,colnames(brain_data_normalised_as_data_frame_good_samples)%in%control_case_ID]))

# check
head(case_exprs[1:5])
dim(case_exprs)
head(control_exprs[1:5])
dim(control_exprs)
head(control_case_exprs[1:5])
dim(control_case_exprs)

# separate by brain region and disorder

Entorhnal_Cortex_case_exprs<-case_exprs[rownames(case_exprs)%in%rownames(Entorhnal_Cortex),]
Frontal_Cortex_case_exprs<-case_exprs[rownames(case_exprs)%in%rownames(Frontal_Cortex),]
Temporal_Cortex_case_exprs<-case_exprs[rownames(case_exprs)%in%rownames(Temporal_Cortex),]
Cerebellum_case_exprs<-case_exprs[rownames(case_exprs)%in%rownames(Cerebellum),]

Entorhnal_Cortex_control_exprs<-control_exprs[rownames(control_exprs)%in%rownames(Entorhnal_Cortex),]
Frontal_Cortex_control_exprs<-control_exprs[rownames(control_exprs)%in%rownames(Frontal_Cortex),]
Temporal_Cortex_control_exprs<-control_exprs[rownames(control_exprs)%in%rownames(Temporal_Cortex),]
Cerebellum_control_exprs<-control_exprs[rownames(control_exprs)%in%rownames(Cerebellum),]

Entorhnal_Cortex_control_case_exprs<-control_case_exprs[rownames(control_case_exprs)%in%rownames(Entorhnal_Cortex),]
Frontal_Cortex_control_case_exprs<-control_case_exprs[rownames(control_case_exprs)%in%rownames(Frontal_Cortex),]
Temporal_Cortex_control_case_exprs<-control_case_exprs[rownames(control_case_exprs)%in%rownames(Temporal_Cortex),]
Cerebellum_control_case_exprs<-control_case_exprs[rownames(control_case_exprs)%in%rownames(Cerebellum),]

# check dataframe

dim(Entorhnal_Cortex_case_exprs)
dim(Frontal_Cortex_case_exprs)
dim(Temporal_Cortex_case_exprs)
dim(Cerebellum_case_exprs)

head(Entorhnal_Cortex_case_exprs)[1:5]
head(Frontal_Cortex_case_exprs)[1:5]
head(Temporal_Cortex_case_exprs)[1:5]
head(Cerebellum_case_exprs)[1:5]

dim(Entorhnal_Cortex_control_exprs)
dim(Frontal_Cortex_control_exprs)
dim(Temporal_Cortex_control_exprs)
dim(Cerebellum_control_exprs)

head(Entorhnal_Cortex_control_exprs)[1:5]
head(Frontal_Cortex_control_exprs)[1:5]
head(Temporal_Cortex_control_exprs)[1:5]
head(Cerebellum_control_exprs)[1:5]

dim(Entorhnal_Cortex_control_case_exprs)
dim(Frontal_Cortex_control_case_exprs)
dim(Temporal_Cortex_control_case_exprs)
dim(Cerebellum_control_case_exprs)

head(Entorhnal_Cortex_control_case_exprs)[1:5]
head(Frontal_Cortex_control_case_exprs)[1:5]
head(Temporal_Cortex_control_case_exprs)[1:5]
head(Cerebellum_control_case_exprs)[1:5]

# separate by gender

Entorhnal_Cortex_case_exprs_F<-Entorhnal_Cortex_case_exprs[rownames(Entorhnal_Cortex_case_exprs)%in%rownames(female_samples),]
Frontal_Cortex_case_exprs_F<-Frontal_Cortex_case_exprs[rownames(Frontal_Cortex_case_exprs)%in%rownames(female_samples),]
Temporal_Cortex_case_exprs_F<-Temporal_Cortex_case_exprs[rownames(Temporal_Cortex_case_exprs)%in%rownames(female_samples),]
Cerebellum_case_exprs_F<-Cerebellum_case_exprs[rownames(Cerebellum_case_exprs)%in%rownames(female_samples),]

Entorhnal_Cortex_control_exprs_F<-Entorhnal_Cortex_control_exprs[rownames(Entorhnal_Cortex_control_exprs)%in%rownames(female_samples),]
Frontal_Cortex_control_exprs_F<-Frontal_Cortex_control_exprs[rownames(Frontal_Cortex_control_exprs)%in%rownames(female_samples),]
Temporal_Cortex_control_exprs_F<-Temporal_Cortex_control_exprs[rownames(Temporal_Cortex_control_exprs)%in%rownames(female_samples),]
Cerebellum_control_exprs_F<-Cerebellum_control_exprs[rownames(Cerebellum_control_exprs)%in%rownames(female_samples),]

Entorhnal_Cortex_case_exprs_M<-Entorhnal_Cortex_case_exprs[rownames(Entorhnal_Cortex_case_exprs)%in%rownames(male_samples),]
Frontal_Cortex_case_exprs_M<-Frontal_Cortex_case_exprs[rownames(Frontal_Cortex_case_exprs)%in%rownames(male_samples),]
Temporal_Cortex_case_exprs_M<-Temporal_Cortex_case_exprs[rownames(Temporal_Cortex_case_exprs)%in%rownames(male_samples),]
Cerebellum_case_exprs_M<-Cerebellum_case_exprs[rownames(Cerebellum_case_exprs)%in%rownames(male_samples),]

Entorhnal_Cortex_control_exprs_M<-Entorhnal_Cortex_control_exprs[rownames(Entorhnal_Cortex_control_exprs)%in%rownames(male_samples),]
Frontal_Cortex_control_exprs_M<-Frontal_Cortex_control_exprs[rownames(Frontal_Cortex_control_exprs)%in%rownames(male_samples),]
Temporal_Cortex_control_exprs_M<-Temporal_Cortex_control_exprs[rownames(Temporal_Cortex_control_exprs)%in%rownames(male_samples),]
Cerebellum_control_exprs_M<-Cerebellum_control_exprs[rownames(Cerebellum_control_exprs)%in%rownames(male_samples),]

Entorhnal_Cortex_control_case_exprs_F<-Entorhnal_Cortex_control_case_exprs[rownames(Entorhnal_Cortex_control_case_exprs)%in%rownames(female_samples),]
Frontal_Cortex_control_case_exprs_F<-Frontal_Cortex_control_case_exprs[rownames(Frontal_Cortex_control_case_exprs)%in%rownames(female_samples),]
Temporal_Cortex_control_case_exprs_F<-Temporal_Cortex_control_case_exprs[rownames(Temporal_Cortex_control_case_exprs)%in%rownames(female_samples),]
Cerebellum_control_case_exprs_F<-Cerebellum_control_case_exprs[rownames(Cerebellum_control_case_exprs)%in%rownames(female_samples),]

Entorhnal_Cortex_control_case_exprs_M<-Entorhnal_Cortex_control_case_exprs[rownames(Entorhnal_Cortex_control_case_exprs)%in%rownames(male_samples),]
Frontal_Cortex_control_case_exprs_M<-Frontal_Cortex_control_case_exprs[rownames(Frontal_Cortex_control_case_exprs)%in%rownames(male_samples),]
Temporal_Cortex_control_case_exprs_M<-Temporal_Cortex_control_case_exprs[rownames(Temporal_Cortex_control_case_exprs)%in%rownames(male_samples),]
Cerebellum_control_case_exprs_M<-Cerebellum_control_case_exprs[rownames(Cerebellum_control_case_exprs)%in%rownames(male_samples),]


# calculate 90th percentile for each sample in each group

extract_good_probe_list<-function(dataset, probe_percentile_threshold, sample_threshold) {
  # dataset - expression dataset as dataframe
  # probe_percentile_threshold - percentile at which to use as cut-off for detected probes = 90th (compared to negative bead controls) - 0.9 = 90th percentile
  # number of samples in which probe must be expressed in - i.e 0.8 = 80% of samples
  # calculate quantile threshold for each sample
  dataset<-as.data.frame(t(dataset))
  sample_quantiles<-apply(dataset, 2, quantile, probs=c(probe_percentile_threshold))
  # count length of quantile - will be number of samples
  number_of_samples<-length(sample_quantiles)
  #boxplot of quanties and print indivdual
  boxplot(as.data.frame(sample_quantiles))
  print(as.data.frame(sample_quantiles))
  # convert probes values to NA in for each sample if probe expression value below sample percentile cut-off
  for (x in 1:number_of_samples) {
    is.na(dataset[x]) <- dataset[x] >= sample_quantiles[x]
  }
  # convert to dataframe
  dataset_dataframe<-as.data.frame(dataset)
  # count number of NA
  dataset_count<-as.data.frame(rowSums(is.na(dataset_dataframe)))
  colnames(dataset_count)<-"count"
  # subset good probes
  good_probes<-rownames(subset(dataset_count, dataset_count$count >= (number_of_samples*sample_threshold)))
  # return good probes
  return(good_probes)
}

#apply function

Entorhnal_Cortex_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Entorhnal_Cortex_case_exprs_F, 0.9, 0.8)
length(Entorhnal_Cortex_case_exprs_F_expressed_probes_list)

Frontal_Cortex_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Frontal_Cortex_case_exprs_F, 0.9, 0.8)
length(Frontal_Cortex_case_exprs_F_expressed_probes_list)

Temporal_Cortex_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Temporal_Cortex_case_exprs_F, 0.9, 0.8)
length(Temporal_Cortex_case_exprs_F_expressed_probes_list)

Cerebellum_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Cerebellum_case_exprs_F, 0.9, 0.8)
length(Cerebellum_case_exprs_F_expressed_probes_list)

Entorhnal_Cortex_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Entorhnal_Cortex_control_exprs_F, 0.9, 0.8)
length(Entorhnal_Cortex_control_exprs_F_expressed_probes_list)

Frontal_Cortex_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Frontal_Cortex_control_exprs_F, 0.9, 0.8)
length(Frontal_Cortex_control_exprs_F_expressed_probes_list)

Temporal_Cortex_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Temporal_Cortex_control_exprs_F, 0.9, 0.8)
length(Temporal_Cortex_control_exprs_F_expressed_probes_list)

Cerebellum_control_exprs_F_expressed_probes_list<-extract_good_probe_list(Cerebellum_control_exprs_F, 0.9, 0.8)
length(Cerebellum_control_exprs_F_expressed_probes_list)

Entorhnal_Cortex_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Entorhnal_Cortex_case_exprs_M, 0.9, 0.8)
length(Entorhnal_Cortex_case_exprs_M_expressed_probes_list)

Frontal_Cortex_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Frontal_Cortex_case_exprs_M, 0.9, 0.8)
length(Frontal_Cortex_case_exprs_M_expressed_probes_list)

Temporal_Cortex_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Temporal_Cortex_case_exprs_M, 0.9, 0.8)
length(Temporal_Cortex_case_exprs_M_expressed_probes_list)

Cerebellum_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Cerebellum_case_exprs_M, 0.9, 0.8)
length(Cerebellum_case_exprs_M_expressed_probes_list)

Entorhnal_Cortex_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Entorhnal_Cortex_control_exprs_M, 0.9, 0.8)
length(Entorhnal_Cortex_control_exprs_M_expressed_probes_list)

Frontal_Cortex_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Frontal_Cortex_control_exprs_M, 0.9, 0.8)
length(Frontal_Cortex_control_exprs_M_expressed_probes_list)

Temporal_Cortex_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Temporal_Cortex_control_exprs_M, 0.9, 0.8)
length(Temporal_Cortex_control_exprs_M_expressed_probes_list)

Cerebellum_control_exprs_M_expressed_probes_list<-extract_good_probe_list(Cerebellum_control_exprs_M, 0.9, 0.8)
length(Cerebellum_control_exprs_M_expressed_probes_list)

Entorhnal_Cortex_control_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Entorhnal_Cortex_control_case_exprs_F, 0.9, 0.8)
length(Entorhnal_Cortex_control_case_exprs_F_expressed_probes_list)

Frontal_Cortex_control_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Frontal_Cortex_control_case_exprs_F, 0.9, 0.8)
length(Frontal_Cortex_control_case_exprs_F_expressed_probes_list)

Temporal_Cortex_control_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Temporal_Cortex_control_case_exprs_F, 0.9, 0.8)
length(Temporal_Cortex_control_case_exprs_F_expressed_probes_list)

Cerebellum_control_case_exprs_F_expressed_probes_list<-extract_good_probe_list(Cerebellum_control_case_exprs_F, 0.9, 0.8)
length(Cerebellum_control_case_exprs_F_expressed_probes_list)

Entorhnal_Cortex_control_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Entorhnal_Cortex_control_case_exprs_M, 0.9, 0.8)
length(Entorhnal_Cortex_control_case_exprs_M_expressed_probes_list)

Frontal_Cortex_control_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Frontal_Cortex_control_case_exprs_M, 0.9, 0.8)
length(Frontal_Cortex_control_case_exprs_M_expressed_probes_list)

Temporal_Cortex_control_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Temporal_Cortex_control_case_exprs_M, 0.9, 0.8)
length(Temporal_Cortex_control_case_exprs_M_expressed_probes_list)

Cerebellum_control_case_exprs_M_expressed_probes_list<-extract_good_probe_list(Cerebellum_control_case_exprs_M, 0.9, 0.8)
length(Cerebellum_control_case_exprs_M_expressed_probes_list)

# merge list of good probes from both case + control + and all brain regions, sort and keep unique values

good_probe_list<-unique(sort(c(Entorhnal_Cortex_case_exprs_F_expressed_probes_list,
                               Frontal_Cortex_case_exprs_F_expressed_probes_list,
                               Temporal_Cortex_case_exprs_F_expressed_probes_list,
                               Cerebellum_case_exprs_F_expressed_probes_list,
                               Entorhnal_Cortex_control_exprs_F_expressed_probes_list,
                               Frontal_Cortex_control_exprs_F_expressed_probes_list,
                               Temporal_Cortex_control_exprs_F_expressed_probes_list,
                               Cerebellum_control_exprs_F_expressed_probes_list,
                               Entorhnal_Cortex_case_exprs_M_expressed_probes_list,
                               Frontal_Cortex_case_exprs_M_expressed_probes_list,
                               Temporal_Cortex_case_exprs_M_expressed_probes_list,
                               Cerebellum_case_exprs_M_expressed_probes_list,
                               Entorhnal_Cortex_control_exprs_M_expressed_probes_list,
                               Frontal_Cortex_control_exprs_M_expressed_probes_list,
                               Temporal_Cortex_control_exprs_M_expressed_probes_list,
                               Cerebellum_control_exprs_M_expressed_probes_list,
                               Entorhnal_Cortex_control_case_exprs_F_expressed_probes_list,
                               Frontal_Cortex_control_case_exprs_F_expressed_probes_list,
                               Temporal_Cortex_control_case_exprs_F_expressed_probes_list,
                               Cerebellum_control_case_exprs_F_expressed_probes_list,
                               Entorhnal_Cortex_control_case_exprs_M_expressed_probes_list,
                               Frontal_Cortex_control_case_exprs_M_expressed_probes_list,
                               Temporal_Cortex_control_case_exprs_M_expressed_probes_list,
                               Cerebellum_control_case_exprs_M_expressed_probes_list)))

length(good_probe_list)
dim(brain_data_normalised_as_data_frame_good_samples)

# extract good probes from dataset

brain_exprs_good_probes<-brain_data_normalised_as_data_frame_good_samples[,colnames(brain_data_normalised_as_data_frame_good_samples)%in%good_probe_list]

Entorhnal_Cortex_case_exprs_good_probes<-Entorhnal_Cortex_case_exprs[,colnames(Entorhnal_Cortex_case_exprs)%in%good_probe_list]

Frontal_Cortex_case_exprs_good_probes<-Frontal_Cortex_case_exprs[,colnames(Frontal_Cortex_case_exprs)%in%good_probe_list]

Temporal_Cortex_case_exprs_good_probes<-Temporal_Cortex_case_exprs[,colnames(Temporal_Cortex_case_exprs)%in%good_probe_list]

Cerebellum_case_exprs_good_probes<-Cerebellum_case_exprs[,colnames(Cerebellum_case_exprs)%in%good_probe_list]

Entorhnal_Cortex_control_exprs_good_probes<-Entorhnal_Cortex_control_exprs[,colnames(Entorhnal_Cortex_control_exprs)%in%good_probe_list]

Frontal_Cortex_control_exprs_good_probes<-Frontal_Cortex_control_exprs[,colnames(Frontal_Cortex_control_exprs)%in%good_probe_list]

Temporal_Cortex_control_exprs_good_probes<-Temporal_Cortex_control_exprs[,colnames(Temporal_Cortex_control_exprs)%in%good_probe_list]

Cerebellum_control_exprs_good_probes<-Cerebellum_control_exprs[,colnames(Cerebellum_control_exprs)%in%good_probe_list]

Entorhnal_Cortex_control_case_exprs_good_probes<-Entorhnal_Cortex_control_case_exprs[,colnames(Entorhnal_Cortex_control_case_exprs)%in%good_probe_list]

Frontal_Cortex_control_case_exprs_good_probes<-Frontal_Cortex_control_case_exprs[,colnames(Frontal_Cortex_control_case_exprs)%in%good_probe_list]

Temporal_Cortex_control_case_exprs_good_probes<-Temporal_Cortex_control_case_exprs[,colnames(Temporal_Cortex_control_case_exprs)%in%good_probe_list]

Cerebellum_control_case_exprs_good_probes<-Cerebellum_control_case_exprs[,colnames(Cerebellum_control_case_exprs)%in%good_probe_list]

# check dataframe - all should be 5376

head(Entorhnal_Cortex_case_exprs_good_probes)[1:5]
dim(Entorhnal_Cortex_case_exprs_good_probes)

head(Frontal_Cortex_case_exprs_good_probes)[1:5]
dim(Frontal_Cortex_case_exprs_good_probes)

head(Temporal_Cortex_case_exprs_good_probes)[1:5]
dim(Temporal_Cortex_case_exprs_good_probes)

head(Cerebellum_case_exprs_good_probes)[1:5]
dim(Cerebellum_case_exprs_good_probes)

head(Entorhnal_Cortex_control_exprs_good_probes)[1:5]
dim(Entorhnal_Cortex_control_exprs_good_probes)

head(Frontal_Cortex_control_exprs_good_probes)[1:5]
dim(Frontal_Cortex_control_exprs_good_probes)

head(Temporal_Cortex_control_exprs_good_probes)[1:5]
dim(Temporal_Cortex_control_exprs_good_probes)

head(Cerebellum_control_exprs_good_probes)[1:5]
dim(Cerebellum_control_exprs_good_probes)

head(Entorhnal_Cortex_control_case_exprs_good_probes)[1:5]
dim(Entorhnal_Cortex_control_case_exprs_good_probes)

head(Frontal_Cortex_control_case_exprs_good_probes)[1:5]
dim(Frontal_Cortex_control_case_exprs_good_probes)

head(Temporal_Cortex_control_case_exprs_good_probes)[1:5]
dim(Temporal_Cortex_control_case_exprs_good_probes)

head(Cerebellum_control_case_exprs_good_probes)[1:5]
dim(Cerebellum_control_case_exprs_good_probes)

# create dataframe with good probes - contain all brain regions, case and control

brain_data_normalised_exprs_good_probes<-brain_data_normalised_as_data_frame_good_samples[rownames(brain_data_normalised_as_data_frame_good_samples)%in%good_probe_list,]
dim(brain_data_normalised_exprs_good_probes)
head(brain_data_normalised_exprs_good_probes)[1:5]

##### PROBE DETECTION TEST ######

setwd(Probe_det_plots_dir)

##### GENDER SPECIFIC PROBE PLOTS #####

Gene_symbols_probes <- mappedkeys(illuminaHumanv4SYMBOL)

# Convert to a list
Gene_symbols <- as.data.frame(illuminaHumanv4SYMBOL[Gene_symbols_probes])

head(Gene_symbols)
dim(Gene_symbols)

#Expressed in females only
XIST_probe_ID<-subset(Gene_symbols, symbol=="XIST")
XIST_probe_ID

#Expressed in Males only
PRKY_probe_ID<-subset(Gene_symbols, symbol=="PRKY")
PRKY_probe_ID

RPS4Y1_probe_ID<-subset(Gene_symbols, symbol=="RPS4Y1")
RPS4Y1_probe_ID

KDM5D_probe_ID<-subset(Gene_symbols, symbol=="KDM5D")
KDM5D_probe_ID

USP9Y_probe_ID<-subset(Gene_symbols, symbol=="USP9Y")
USP9Y_probe_ID

UTY_probe_ID<-subset(Gene_symbols, symbol=="UTY")
UTY_probe_ID

# HK genes expressed in all cells + males + females
MKRN1_probe_ID<-subset(Gene_symbols, symbol=="MKRN1")
MKRN1_probe_ID

ADIPOR1_probe_ID<-subset(Gene_symbols, symbol=="ADIPOR1")
ADIPOR1_probe_ID

BNIP3L_probe_ID<-subset(Gene_symbols, symbol=="BNIP3L")
BNIP3L_probe_ID

#RNF10_probe_ID<-subset(Gene_symbols, symbol=="RNF10")
#RNF10_probe_ID


# merge all genes 
gene_list<-rbind(XIST_probe_ID,
                 PRKY_probe_ID,
                 RPS4Y1_probe_ID,
                 KDM5D_probe_ID,
                 USP9Y_probe_ID,
                 UTY_probe_ID,
                 MKRN1_probe_ID,
                 ADIPOR1_probe_ID,
                 BNIP3L_probe_ID)
#                 RNF10_probe_ID)

gene_list

# create table of genes and state if KH or gender specific
gene_table<-read.table(text =
                      "Gene Expressed_in
                       ADIPOR1 All
                       BNIP3L All
                       KDM5D Males
                       MKRN1 All
                       PRKY Males
                       RPS4Y1 Males
                       USP9Y Males
                       UTY Males
                       XIST Females" , header=T)


gene_table

plot_gender_specific_genes<-function(Expression_table, Gender, genes_to_extract, threshold, boxplot_title){
  #extract gene of interest
  Expression_table_gene_check<-as.data.frame(t(Expression_table[rownames(Expression_table)%in% genes_to_extract$probe_id,]))
  # change colnames TO GENE SYMBOL using genes to extract file
  for (x in 1:dim(Expression_table_gene_check)[2]){
    colnames(Expression_table_gene_check)[x]<-gene_list[genes_to_extract$probe_id==colnames(Expression_table_gene_check)[x],2]
  }
  # add in gender information
  Expression_table_gene_check_gender<-merge(Gender, Expression_table_gene_check, by="row.names")
  rownames(Expression_table_gene_check_gender)<-Expression_table_gene_check_gender$Row.names
  Expression_table_gene_check_gender$Row.names<-NULL
  #melt dataframe for plot
  Expression_table_gene_check_gender_melt<-melt(Expression_table_gene_check_gender, by=Gender)
  # change variable colun from factor to character
  Expression_table_gene_check_gender_melt$variable<-as.character(Expression_table_gene_check_gender_melt$variable)
  # order dataframe by variable
  Expression_table_gene_check_gender_melt<-Expression_table_gene_check_gender_melt[order(Expression_table_gene_check_gender_melt$variable),]
  # calculate user defined percentie threshold
  sample_quantiles<-apply(Expression_table, 2, quantile, probs=threshold)
  # mean of used defined threshold across samples
  mean_threshold=mean(sample_quantiles)
  #plot
  plot1<-qplot(variable, value, colour=get(colnames(Gender)), data = Expression_table_gene_check_gender_melt, geom = c("boxplot", "jitter")) + 
    geom_hline(yintercept = mean_threshold) +
    ggtitle(boxplot_title) +
    theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="Gene",y="Expression", colour = colnames(Gender)) 
  # 2nd legend
  plot2<-tableGrob(gene_table, rows=NULL)
  # plot
  grid.arrange(plot1, plot2,
               nrow=2,
               heights=c(3,1))
}


#plot


plot_gender_specific_genes(brain_data_normalised_as_data_frame_no_dup_illID, gender_comparison[1], gene_list, 0.9, "FULL Dataset")
plot_gender_specific_genes(brain_data_normalised_as_data_frame_no_dup_illID, gender_comparison[2], gene_list, 0.9, "FULL Dataset")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_case_exprs)), gender_comparison[2], gene_list, 0.9, "Entorhnal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_control_exprs)), gender_comparison[2], gene_list, 0.9, "Entorhnal_Cortex_control")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Entorhnal_Cortex_control_case")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_case_exprs)), gender_comparison[2], gene_list, 0.9, "Temporal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_control_exprs)), gender_comparison[2], gene_list, 0.9, "Temporal_Cortex_control")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Temporal_Cortex_control_case")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_case_exprs)), gender_comparison[2], gene_list, 0.9, "Cerebellum_case")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_control_exprs)), gender_comparison[2], gene_list, 0.9, "Cerebellum_control")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Cerebellum_control_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_case_exprs)), gender_comparison[2], gene_list, 0.9, "Frontal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_control_exprs)), gender_comparison[2], gene_list, 0.9, "Frontal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Frontal_Cortex_control_case")

#plot to pdf

setwd(Gender_plots_dir)
pdf("Gender_specific_gene_plot_and_detectable_cut_off_threshold_used_before_probe_removal_predicted_gender.pdf")
plot_gender_specific_genes(brain_data_normalised_as_data_frame_no_dup_illID, gender_comparison[2], gene_list, 0.9, "FULL Dataset")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_case_exprs)), gender_comparison[2], gene_list, 0.9, "Entorhnal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_control_exprs)), gender_comparison[2], gene_list, 0.9, "Entorhnal_Cortex_control")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Entorhnal_Cortex_control_case")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_case_exprs)), gender_comparison[2], gene_list, 0.9, "Temporal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_control_exprs)), gender_comparison[2], gene_list, 0.9, "Temporal_Cortex_control")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Temporal_Cortex_control_case")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_case_exprs)), gender_comparison[2], gene_list, 0.9, "Cerebellum_case")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_control_exprs)), gender_comparison[2], gene_list, 0.9, "Cerebellum_control")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Cerebellum_control_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_case_exprs)), gender_comparison[2], gene_list, 0.9, "Frontal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_control_exprs)), gender_comparison[2], gene_list, 0.9, "Frontal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_control_case_exprs)), gender_comparison[2], gene_list, 0.9, "Frontal_Cortex_control_case")
dev.off()

pdf("Gender_specific_gene_plot_and_detectable_cut_off_threshold_used_before_probe_removal_clinical_gender.pdf")
plot_gender_specific_genes(brain_data_normalised_as_data_frame_no_dup_illID, gender_comparison[1], gene_list, 0.9, "FULL Dataset")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_case_exprs)), gender_comparison[1], gene_list, 0.9, "Entorhnal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_control_exprs)), gender_comparison[1], gene_list, 0.9, "Entorhnal_Cortex_control")
plot_gender_specific_genes(as.data.frame(t(Entorhnal_Cortex_control_case_exprs)), gender_comparison[1], gene_list, 0.9, "Entorhnal_Cortex_control_case")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_case_exprs)), gender_comparison[1], gene_list, 0.9, "Temporal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_control_exprs)), gender_comparison[1], gene_list, 0.9, "Temporal_Cortex_control")
plot_gender_specific_genes(as.data.frame(t(Temporal_Cortex_control_case_exprs)), gender_comparison[1], gene_list, 0.9, "Temporal_Cortex_control_case")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_case_exprs)), gender_comparison[1], gene_list, 0.9, "Cerebellum_case")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_control_exprs)), gender_comparison[1], gene_list, 0.9, "Cerebellum_control")
plot_gender_specific_genes(as.data.frame(t(Cerebellum_control_case_exprs)), gender_comparison[1], gene_list, 0.9, "Cerebellum_control_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_case_exprs)), gender_comparison[1], gene_list, 0.9, "Frontal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_control_exprs)), gender_comparison[1], gene_list, 0.9, "Frontal_Cortex_case")
plot_gender_specific_genes(as.data.frame(t(Frontal_Cortex_control_case_exprs)), gender_comparison[1], gene_list, 0.9, "Frontal_Cortex_control_case")
dev.off()

##### PCA ####

pca_data<-function(data, legend_position){
  #run PCA
  pca<-prcomp(data)
  # order of samples in expression data
  sample_order<-colnames(data)
  # merge pheno info together
  pheno<-merge(Diagnosis, gender_comparison, by="row.names")
  rownames(pheno)<-pheno$Row.names
  pheno$Row.names<-NULL
  pheno<-merge(pheno, Tissue, by="row.names")
  rownames(pheno)<-pheno$Row.names
  pheno$Row.names<-NULL
  # subset pca_pheno to match expression data
  pca_pheno<-subset(pheno, rownames(pheno) %in% colnames(data))
  # match order
  ordered_pca_pheno<-pca_pheno[match(sample_order, rownames(pca_pheno)),]
  Diagnosis_pca_colour<-labels2colors(as.character(ordered_pca_pheno$PHENOTYPE))
  Gender_pca_colour<-labels2colors(as.character(ordered_pca_pheno$Predicted_Gender))
  clinical_Gender_pca_colour<-labels2colors(as.character(ordered_pca_pheno$Clinical_Gender))
  Ethnicity_pca_colour<-labels2colors(as.character(ordered_pca_pheno$TISSUE))
  # pca plot - Diagnosis
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Diagnosis",col="black", pch=21,bg=Diagnosis_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$PHENOTYPE), fill=unique(Diagnosis_pca_colour), title="Diagnosis")
  # pca plot - Gender
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Predicted Gender",col="black", pch=21,bg=Gender_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$Predicted_Gender), fill=unique(Gender_pca_colour), title="Predicted Gender")
  # pca plot - Gender clinical
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Clinical Gender",col="black", pch=21,bg=clinical_Gender_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$Clinical_Gender), fill=unique(clinical_Gender_pca_colour), title="Clinical Gender")
  # pca plot - Ethnicity
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Tissue",col="black", pch=21,bg=Ethnicity_pca_colour)
  legend(legend_position, legend=unique(ordered_pca_pheno$TISSUE), fill=unique(Ethnicity_pca_colour), title="Ethnicity")
}

#apply function
pca_data(brain_data_normalised_exprs_good_probes, 'bottomright')

#plot to pdf
setwd(pca_dir)

pdf("1.PCA_plot_before_SVA.pdf")
pca_data(brain_data_normalised_exprs_good_probes, 'bottomright')
dev.off()

##### SVA PER BRAIN REGION PER CASE/CONTROL #####

# add diagnosis column in - CASE CONTROL

Entorhnal_Cortex_case_exprs_good_probes<-cbind(Diagnosis = "AD", Entorhnal_Cortex_case_exprs_good_probes)
Frontal_Cortex_case_exprs_good_probes<-cbind(Diagnosis = "AD", Frontal_Cortex_case_exprs_good_probes)
Temporal_Cortex_case_exprs_good_probes<-cbind(Diagnosis = "AD", Temporal_Cortex_case_exprs_good_probes)
Cerebellum_case_exprs_good_probes<-cbind(Diagnosis = "AD", Cerebellum_case_exprs_good_probes)
Entorhnal_Cortex_control_exprs_good_probes<-cbind(Diagnosis = "CONTROL", Entorhnal_Cortex_control_exprs_good_probes)
Frontal_Cortex_control_exprs_good_probes<-cbind(Diagnosis = "CONTROL", Frontal_Cortex_control_exprs_good_probes)
Temporal_Cortex_control_exprs_good_probes<-cbind(Diagnosis = "CONTROL", Temporal_Cortex_control_exprs_good_probes)
Cerebellum_control_exprs_good_probes<-cbind(Diagnosis = "CONTROL", Cerebellum_control_exprs_good_probes)
Entorhnal_Cortex_control_case_exprs_good_probes<-cbind(Diagnosis = "CONTROL_CASE", Entorhnal_Cortex_control_case_exprs_good_probes)
Frontal_Cortex_control_case_exprs_good_probes<-cbind(Diagnosis = "CONTROL_CASE", Frontal_Cortex_control_case_exprs_good_probes)
Temporal_Cortex_control_case_exprs_good_probes<-cbind(Diagnosis = "CONTROL_CASE", Temporal_Cortex_control_case_exprs_good_probes)
Cerebellum_control_case_exprs_good_probes<-cbind(Diagnosis = "CONTROL_CASE", Cerebellum_control_case_exprs_good_probes)

# check order of gene in each dataframe

any(colnames(Entorhnal_Cortex_case_exprs_good_probes)==colnames(Entorhnal_Cortex_control_exprs_good_probes))==F
any(colnames(Frontal_Cortex_case_exprs_good_probes)==colnames(Frontal_Cortex_control_exprs_good_probes))==F
any(colnames(Temporal_Cortex_case_exprs_good_probes)==colnames(Temporal_Cortex_control_exprs_good_probes))==F
any(colnames(Cerebellum_case_exprs_good_probes)==colnames(Cerebellum_control_exprs_good_probes))==F

# merge case control in to one dataset per brain region

Entorhnal_Cortex_good_probes<-rbind(Entorhnal_Cortex_case_exprs_good_probes, Entorhnal_Cortex_control_exprs_good_probes, Entorhnal_Cortex_control_case_exprs_good_probes)
Frontal_Cortex_good_probes<-rbind(Frontal_Cortex_case_exprs_good_probes, Frontal_Cortex_control_exprs_good_probes, Frontal_Cortex_control_case_exprs_good_probes)
Temporal_Cortex_good_probes<-rbind(Temporal_Cortex_case_exprs_good_probes, Temporal_Cortex_control_exprs_good_probes, Temporal_Cortex_control_case_exprs_good_probes)
Cerebellum_good_probes<-rbind(Cerebellum_case_exprs_good_probes, Cerebellum_control_exprs_good_probes, Cerebellum_control_case_exprs_good_probes)

#add predicted gender in 

Entorhnal_Cortex_good_probes<-merge(gender_comparison[1], Entorhnal_Cortex_good_probes, by="row.names")
rownames(Entorhnal_Cortex_good_probes)<-Entorhnal_Cortex_good_probes$Row.names
Entorhnal_Cortex_good_probes$Row.names<-NULL

Frontal_Cortex_good_probes<-merge(gender_comparison[1], Frontal_Cortex_good_probes, by="row.names")
rownames(Frontal_Cortex_good_probes)<-Frontal_Cortex_good_probes$Row.names
Frontal_Cortex_good_probes$Row.names<-NULL

Temporal_Cortex_good_probes<-merge(gender_comparison[1], Temporal_Cortex_good_probes, by="row.names")
rownames(Temporal_Cortex_good_probes)<-Temporal_Cortex_good_probes$Row.names
Temporal_Cortex_good_probes$Row.names<-NULL

Cerebellum_good_probes<-merge(gender_comparison[1], Cerebellum_good_probes, by="row.names")
rownames(Cerebellum_good_probes)<-Cerebellum_good_probes$Row.names
Cerebellum_good_probes$Row.names<-NULL

head(Entorhnal_Cortex_good_probes)[1:5]
head(Frontal_Cortex_good_probes)[1:5]
head(Temporal_Cortex_good_probes)[1:5]
head(Cerebellum_good_probes)[1:5]

# create sva function

check_SV_in_data<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_pheno<-sorted_by_diagnosis[c(1,2)]
  dataset_exprs<-t(sorted_by_diagnosis[3:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis + Gender
  mod = model.matrix(~Diagnosis+Clinical_Gender, data=dataset_pheno)
  #mod = model.matrix(~Diagnosis+Predicted_Gender, data=dataset_pheno)
  # check number of SV in data
  num.sv(dataset_exprs, mod, method="leek")
}

# check sv

check_SV_in_data(Entorhnal_Cortex_good_probes)
check_SV_in_data(Frontal_Cortex_good_probes)
check_SV_in_data(Temporal_Cortex_good_probes)
check_SV_in_data(Cerebellum_good_probes)

# create function to adjust for SV

adjust_for_sva<-function(dataset){
  # create sva compatable matrix - sample in columns, probes in rows - pheno info seperate - sort by diagnosis 1st to keep AD top
  sorted_by_diagnosis<-dataset[order(dataset$Diagnosis),]
  # separate expresion and pheno
  dataset_sva_pheno<-sorted_by_diagnosis[c(1,2)]
  dataset_sva_exprs<-t(sorted_by_diagnosis[3:dim(sorted_by_diagnosis)[2]])
  #full model matrix for Diagnosis
  mod = model.matrix(~Diagnosis+Clinical_Gender, data=dataset_sva_pheno)
  mod0 = model.matrix(~1, data=dataset_sva_pheno)
  # number of SV
  num.sv(dataset_sva_exprs, mod, method="leek")
  n.sv=num.sv(dataset_sva_exprs, mod, method="leek")
  # exit if n.sv=0
  if(n.sv==0){stop("No Significant Variable found, exiting....")}
  # apply sva - removed n.sv
  svobj = sva(dataset_sva_exprs, mod, mod0, n.sv=n.sv, method="two-step")
  # adjust for sva
  X = cbind(mod, svobj$sv)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(dataset_sva_exprs))
  P = ncol(mod)
  clean_data<-dataset_sva_exprs - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  # merge clean data with pheno
  clean_data_with_pheno<-merge(dataset_sva_pheno, as.data.frame(t(clean_data)), by="row.names")
  rownames(clean_data_with_pheno)<-clean_data_with_pheno$Row.names
  clean_data_with_pheno$Row.names<-NULL
  # check SVA on adjusted data
  cat("\n")
  cat("number of surrogate variables after adjustment:")
  cat("\n")
  print(num.sv(clean_data, mod, method="leek"))
  # return clean data with pheno
  return(clean_data_with_pheno)
}

# run sva

Temporal_Cortex_good_probes_sva_adjusted<-adjust_for_sva(Temporal_Cortex_good_probes)

#Temporal cortex no adjustment Cortex no SVA adjusment
#Temporal_Cortex_good_probes_sva_adjusted<-Temporal_Cortex_good_probes
Entorhinal_Cortex_good_probes_sva_adjusted<-Entorhnal_Cortex_good_probes
Frontal_Cortex_good_probes_sva_adjusted<-Frontal_Cortex_good_probes
Cerebellum_good_probes_sva_adjusted<-Cerebellum_good_probes

# separate case and control

Frontal_Cortex_good_probes_sva_adjusted_case<-Frontal_Cortex_good_probes_sva_adjusted[Frontal_Cortex_good_probes_sva_adjusted$Diagnosis=="AD",]
Temporal_Cortex_good_probes_sva_adjusted_case<-Temporal_Cortex_good_probes_sva_adjusted[Temporal_Cortex_good_probes_sva_adjusted$Diagnosis=="AD",]
Cerebellum_good_probes_sva_adjusted_case<-Cerebellum_good_probes_sva_adjusted[Cerebellum_good_probes_sva_adjusted$Diagnosis=="AD",]
Entorhinal_Cortex_good_probes_sva_adjusted_case<-Entorhinal_Cortex_good_probes_sva_adjusted[Entorhinal_Cortex_good_probes_sva_adjusted$Diagnosis=="AD",]

Frontal_Cortex_good_probes_sva_adjusted_control<-Frontal_Cortex_good_probes_sva_adjusted[Frontal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Temporal_Cortex_good_probes_sva_adjusted_control<-Temporal_Cortex_good_probes_sva_adjusted[Temporal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Cerebellum_good_probes_sva_adjusted_control<-Cerebellum_good_probes_sva_adjusted[Cerebellum_good_probes_sva_adjusted$Diagnosis=="CONTROL",]
Entorhinal_Cortex_good_probes_sva_adjusted_control<-Entorhinal_Cortex_good_probes_sva_adjusted[Entorhinal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL",]

Frontal_Cortex_good_probes_sva_adjusted_control_case<-Frontal_Cortex_good_probes_sva_adjusted[Frontal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL_CASE",]
Temporal_Cortex_good_probes_sva_adjusted_control_case<-Temporal_Cortex_good_probes_sva_adjusted[Temporal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL_CASE",]
Cerebellum_good_probes_sva_adjusted_control_case<-Cerebellum_good_probes_sva_adjusted[Cerebellum_good_probes_sva_adjusted$Diagnosis=="CONTROL_CASE",]
Entorhinal_Cortex_good_probes_sva_adjusted_control_case<-Entorhinal_Cortex_good_probes_sva_adjusted[Entorhinal_Cortex_good_probes_sva_adjusted$Diagnosis=="CONTROL_CASE",]
#check sample numbers

dim(Frontal_Cortex_good_probes_sva_adjusted)
dim(Frontal_Cortex_good_probes_sva_adjusted_case)
dim(Frontal_Cortex_good_probes_sva_adjusted_control)
dim(Frontal_Cortex_good_probes_sva_adjusted_control_case)

dim(Temporal_Cortex_good_probes_sva_adjusted)
dim(Temporal_Cortex_good_probes_sva_adjusted_case)
dim(Temporal_Cortex_good_probes_sva_adjusted_control)
dim(Temporal_Cortex_good_probes_sva_adjusted_control_case)

dim(Cerebellum_good_probes_sva_adjusted)
dim(Cerebellum_good_probes_sva_adjusted_case)
dim(Cerebellum_good_probes_sva_adjusted_control)
dim(Cerebellum_good_probes_sva_adjusted_control_case)

dim(Entorhinal_Cortex_good_probes_sva_adjusted)
dim(Entorhinal_Cortex_good_probes_sva_adjusted_case)
dim(Entorhinal_Cortex_good_probes_sva_adjusted_control)
dim(Entorhinal_Cortex_good_probes_sva_adjusted_control_case)

##### SAMPLE NETWORK PLOT ######

# sample plot function - taken from steve expression pipeline

sampleNetwork_plot <- function(dataset) {
  datExprs<-t(dataset[3:dim(dataset)[2]])
  diagnosis<-dataset[2]
  gp_col <- "group"
  cat(" setting up data for qc plots","\r","\n")
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- colnames(datExprs)
  IAC=cor(datExprs, method="p",use="p")
  diag(IAC)=0
  A.IAC=((1+IAC)/2)^2 ## ADJACENCY MATRIX
  cat(" fundamentalNetworkConcepts","\r","\n")
  FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
  K2=FNC$ScaledConnectivity
  Z.K=(K2-mean(K2))/sd(K2)
  Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
  rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
  rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
  # set colours
  cat(" colorvec [",paste(gp_col),"]","\r","\n")
  if(gp_col=="chip") { colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode)) }
  if(gp_col=="group") { colorvec <- labels2colors(diagnosis[1]) }
  mean_IAC <- mean(IAC[upper.tri(IAC)])
  ## samplenetwork
  local(
    {colLab <<- function(n,treeorder) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colorvec[treeorder][i], lab.font = i%%3))
      }
      n
    }
    i <- 0
    })
  cat(" begin SampleNetwork plots","\r","\n")
  group_colours<-unique(cbind(colorvec, diagnosis))
  ## Cluster for pics
  cluster1 <- hclust(as.dist(1-A.IAC),method="average")
  cluster1order <- cluster1$order
  cluster2 <- as.dendrogram(cluster1,hang=0.1)
  cluster3 <- dendrapply(cluster2,colLab,cluster1order)
  ## PLOTS
  ## cluster IAC
  par(mfrow=c(2,2))
  par(mar=c(5,6,4,2))
  plot(cluster3,nodePar=list(lab.cex=1,pch=NA),
       main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),
       xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
  mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
  ## Connectivity
  par(mar=c(5,5,4,2))
  plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
  text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
  abline(h=-2)
  abline(h=-3)
  par(mar=c(5,5,4,2))
  plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec ,cex.main=1.8,cex.lab=1.4)
  abline(lm(Z.C~Z.K),col="black",lwd=2)
  mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
  abline(v=-2,lty=2,col="grey")
  abline(h=-2,lty=2,col="grey")
  ##blank plot for legend
  par(mar=c(5,5,4,2))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend(0.6, 1.4, unique(diagnosis[,1]), fill=unique(colorvec))
} #taken from steves expression pipeline

# create functio to ID outliers

names_of_outliers<-function(dataset, threshold){
  datExprs<-t(dataset[3:dim(dataset)[2]])
  IAC = cor(datExprs, method = "p", use = "p")
  diag(IAC) = 0
  A.IAC = ((1 + IAC)/2)^2  ## ADJACENCY MATRIX
  # fundamentalNetworkConcepts
  FNC = fundamentalNetworkConcepts(A.IAC)  ## WGCNA
  K2 = FNC$ScaledConnectivity
  Z.K = round((K2 - mean(K2))/sd(K2), 3)
  Z.K_outliers <- Z.K < threshold
  Z.K_outliers <- names(Z.K_outliers[Z.K_outliers == TRUE])
  n_outliers <- length(Z.K_outliers)
  return(Z.K_outliers)
}

# create function to run network analysis on each expression dataset, plot and remove bad samples 

run_sample_network_plot<-function(dataset, threshold){
  #sample network plot
  sampleNetwork_plot(dataset)
  #identify sample below Z.K -2  
  dataset_removal_1<-names_of_outliers(dataset, threshold)
  #remove samples with ZK below -2
  dataset_QC<-dataset[!(rownames(dataset)%in%dataset_removal_1),]
  #sample network plot
  sampleNetwork_plot(dataset_QC)
  #create empty count list to record samples removed
  count<-dataset_removal_1
  # reiterate above till no samples fall below threshold
  while (length(dataset_removal_1)>0) {
    # remove bad samples - 1st iteration removes none
    dataset_QC<-dataset_QC[!(rownames(dataset_QC)%in%dataset_removal_1),]
    #identify sample below Z.K -2  
    dataset_removal_1<-names_of_outliers(dataset_QC, threshold)
    #record samples removed
    count<-c(count, dataset_removal_1)
  }
  #final network plot
  sampleNetwork_plot(dataset_QC)
  # print to screen number of samples removed
  cat("\n")
  print(c("Total number of samples removed...", length(count)))
  # return clean expression set
  return(dataset_QC)
}


# # run sample network on entorhinal Cortex

Entorhinal_Cortex_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_case, -3)
Entorhinal_Cortex_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_control, -3)
Entorhinal_Cortex_good_probes_sva_adjusted_control_case_QC<-run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_control_case, -3)

# run sample network on Frontal_Cortex

Frontal_Cortex_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_case, -3)
Frontal_Cortex_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_control, -3)
Frontal_Cortex_good_probes_sva_adjusted_control_case_QC<-run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_control_case, -3)

# run sample network on Temporal_Cortex_

Temporal_Cortex_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Temporal_Cortex_good_probes_sva_adjusted_case, -3)
Temporal_Cortex_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Temporal_Cortex_good_probes_sva_adjusted_control, -3)
Temporal_Cortex_good_probes_sva_adjusted_control_case_QC<-run_sample_network_plot(Temporal_Cortex_good_probes_sva_adjusted_control_case, -3)

# run sample network on Cerebellum_

Cerebellum_good_probes_sva_adjusted_case_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_case, -3)
Cerebellum_good_probes_sva_adjusted_control_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control, -3)
Cerebellum_good_probes_sva_adjusted_control_case_QC<-run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control_case, -3)

##### PLOT SAMPLE NETWORK ANALYSIS TO PDF #####

setwd(sample_network_dir)

pdf("Entorhnal_Cortex_case_sample_network_analysis.pdf")
run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Entorhnal_Cortex_control_sample_network_analysis.pdf")
run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Entorhnal_Cortex_control_case_sample_network_analysis.pdf")
run_sample_network_plot(Entorhinal_Cortex_good_probes_sva_adjusted_control_case, -3)
dev.off()

pdf("Frontal_Cortex_case_sample_network_analysis.pdf")
run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Frontal_Cortex_control_sample_network_analysis.pdf")
run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Frontal_Cortex_control_case_sample_network_analysis.pdf")
run_sample_network_plot(Frontal_Cortex_good_probes_sva_adjusted_control_case, -3)
dev.off()

pdf("Temporal_Cortex_case_sample_network_analysis.pdf")
run_sample_network_plot(Temporal_Cortex_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Temporal_Cortex_control_sample_network_analysis.pdf")
run_sample_network_plot(Temporal_Cortex_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Temporal_Cortex_control_case_sample_network_analysis.pdf")
run_sample_network_plot(Temporal_Cortex_good_probes_sva_adjusted_control_case, -3)
dev.off()

pdf("Cerebellum_case_sample_network_analysis.pdf")
run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_case, -3)
dev.off()

pdf("Cerebellum_control_sample_network_analysis.pdf")
run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control, -3)
dev.off()

pdf("Cerebellum_control_case_sample_network_analysis.pdf")
run_sample_network_plot(Cerebellum_good_probes_sva_adjusted_control_case, -3)
dev.off()

##### CREATE QC'd DATASET #####

# extract sample ID's from QC'd sample network file

# check colnames same in all dataframes
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Frontal_Cortex_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Temporal_Cortex_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Temporal_Cortex_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Cerebellum_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Cerebellum_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Entorhinal_Cortex_good_probes_sva_adjusted_case_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Entorhinal_Cortex_good_probes_sva_adjusted_control_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Frontal_Cortex_good_probes_sva_adjusted_control_case_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Temporal_Cortex_good_probes_sva_adjusted_control_case_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Cerebellum_good_probes_sva_adjusted_control_case_QC))==F
any(colnames(Frontal_Cortex_good_probes_sva_adjusted_case_QC)==colnames(Entorhinal_Cortex_good_probes_sva_adjusted_control_case_QC))==F

# cbind all dataframes

brain_data_QCd<-rbind(Frontal_Cortex_good_probes_sva_adjusted_case_QC,
                      Frontal_Cortex_good_probes_sva_adjusted_control_QC,
                      Temporal_Cortex_good_probes_sva_adjusted_case_QC,
                      Temporal_Cortex_good_probes_sva_adjusted_control_QC,
                      Cerebellum_good_probes_sva_adjusted_case_QC,
                      Cerebellum_good_probes_sva_adjusted_control_QC,
                      Entorhinal_Cortex_good_probes_sva_adjusted_case_QC,
                      Entorhinal_Cortex_good_probes_sva_adjusted_control_QC,
                      Frontal_Cortex_good_probes_sva_adjusted_control_case_QC,
                      Temporal_Cortex_good_probes_sva_adjusted_control_case_QC,
                      Cerebellum_good_probes_sva_adjusted_control_case_QC,
                      Entorhinal_Cortex_good_probes_sva_adjusted_control_case_QC)

dim(brain_data_QCd)
head(brain_data_QCd)[1:5]

##### PCA ON CLEAN DATA ####

setwd(pca_dir)

pca_data(as.data.frame(t(brain_data_QCd[3:dim(brain_data_QCd)[2]])), 'bottomright')

pdf("2.PCA_plot_after_QC.pdf")
pca_data(as.data.frame(t(brain_data_QCd[3:dim(brain_data_QCd)[2]])), 'bottomright')
dev.off()

##### CONVERT PROBE ID TO ENTREZ ID #####
#- changed 26/04/16 - convert from illumina manifest, and then missing values from biomart
# biomart converted most DE gene illumina probe ID to multiple Entrez gene IDS - so was droppped
# further investigation found this is fault of biomart.

# create convert_probe_id_to_entrez_id function 

convert_probe_id_to_entrez_id <- function(expression_dataset, probe_mapping_file){
  # transform dataset # - removed this step
  expression_dataset_t<-as.data.frame(expression_dataset)
  # keep only probes which appear in probe_mapping_file
  data_frame_in_probe_mapper<-expression_dataset_t[colnames(expression_dataset_t)%in%probe_mapping_file$probe_id]
  # match probe id in data_frame_in_probe_mapper to that in probe_mapping_file and convert to entrez id
  colnames(data_frame_in_probe_mapper)<-probe_mapping_file$entrezgene[match(colnames(data_frame_in_probe_mapper), probe_mapping_file$probe_id)]
  return(data_frame_in_probe_mapper)
}

# Get the probe identifiers that are mapped to an ENTREZ Gene ID using illuminaHumanv4.db
mapped_probes <- mappedkeys(illuminaHumanv4ENTREZID)

# Convert to a list
illuminaHumanv4.db_mapping <- as.data.frame(illuminaHumanv4ENTREZID[mapped_probes])
# arrange order of column by entrezgene probe_id
illuminaHumanv4.db_mapping<-illuminaHumanv4.db_mapping[c(2,1)]
colnames(illuminaHumanv4.db_mapping)[1]<-"entrezgene"

head(illuminaHumanv4.db_mapping)
dim(illuminaHumanv4.db_mapping)

#check any duplicated probe IDs
anyDuplicated(illuminaHumanv4.db_mapping$probe_id)

#check any dupliacted entrezgene IDs
anyDuplicated(illuminaHumanv4.db_mapping$entrezgene)

illuminaHumanv4.db_mapping_entrez_id_unique<-illuminaHumanv4.db_mapping[!(duplicated(illuminaHumanv4.db_mapping$probe_id) | duplicated(illuminaHumanv4.db_mapping$probe_id, fromLast= TRUE)),]

head(illuminaHumanv4.db_mapping_entrez_id_unique)
dim(illuminaHumanv4.db_mapping_entrez_id_unique)

# using illuminaHumanv4.db_mapping_entrez_id_unique instead of Illumina_HT_12_v4_probe_list_entrez_id_unique from this point onwards

brain_data_QCd_entrez_id<-convert_probe_id_to_entrez_id(brain_data_QCd, illuminaHumanv4.db_mapping_entrez_id_unique)
dim(brain_data_QCd)
dim(brain_data_QCd_entrez_id)
length(which(duplicated(colnames(brain_data_QCd_entrez_id))))

head(brain_data_QCd_entrez_id[100:110])

##### COLLAPSE MULTIPLE ENTREZ ID BY SELECTING ONE WITH HIGHEST AVERAGE EXPRESSION ACROSS SAMPLES ######

## create function

select_duplicate_probe_by_top_expr <- function(x) {
  # transpose data frame - keep as dataframe
  x_t<-as.data.frame(t(x))
  # calculate mean expression per probe across samples - create new column - probe mean column
  x_t$probe_mean_expression<-rowMeans(x_t)
  #copy rownames (probe id) to column and truncate
  x_t$trunc_entrez_id<-trunc(as.numeric(as.character(rownames(x_t))))
  # order data frame by truncated probe id and then expression level
  x_t<-x_t[order(x_t$trunc_entrez_id, -x_t$probe_mean_expression), ]
  # remove all duplicate probe id - keep one with highest mean expression
  x_t_unique<-x_t[!duplicated(x_t$trunc_entrez_id),]
  #unique entrez column back to row name
  rownames(x_t_unique)<-x_t_unique$trunc_entrez_id
  #remove unwanted column
  x_t_unique$trunc_entrez_id<-NULL
  #remove unwanted column
  x_t_unique$probe_mean_expression<-NULL
  #transpose dataframe back
  x_unique<-as.data.frame(t(x_t_unique))
  return(x_unique)
}

## apply function to dataframes - check number of probes - check duplicates

brain_data_QCd_entrez_id_unique<-select_duplicate_probe_by_top_expr(brain_data_QCd_entrez_id)
dim(brain_data_QCd_entrez_id_unique)
length(which(duplicated(colnames(brain_data_QCd_entrez_id_unique))))

head(brain_data_QCd_entrez_id_unique[1:5])

##### MERGE EXPRESSION DATA AND PHENOTYPE DATA #####

# create pheno file
pheno<-merge(Diagnosis, gender_comparison, by="row.names")
rownames(pheno)<-pheno$Row.names
pheno$Row.names<-NULL
pheno<-merge(pheno, Tissue, by="row.names")
rownames(pheno)<-pheno$Row.names
pheno$Row.names<-NULL
head(pheno)

#subset pheno to data
pheno<-subset(pheno, rownames(pheno) %in% rownames(brain_data_QCd_entrez_id_unique))
dim(pheno)
dim(brain_data_QCd_entrez_id_unique)

##### WRITE #####

setwd(clean_data_dir)

write.table(brain_data_QCd_entrez_id_unique, file ="MRC-LBB_expression_table_QCd.txt", sep="\t", quote=F)
write.table(pheno, file ="MRC-LBB_phenotype.txt", sep="\t", quote=F)

##### FULL GENE LIST IN ENTREZ GENE ID #####

# added 26/09/17 - need as background list for pathway analysis

dir.create("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis")
dir.create("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis/1.Create_Input")
setwd("/media/hamel/Workspace/Dropbox/Projects/MRC-LBB/4.Pathway_Analysis/1.Create_Input")

write.table(unique(illuminaHumanv4.db_mapping[,1]), row.names=F, col.names=F, quote=F, file= "background_gene_list.txt")

##### SAVE #####

setwd(work_dir)

save.image("MRC-LBB_processing.Rdata")


