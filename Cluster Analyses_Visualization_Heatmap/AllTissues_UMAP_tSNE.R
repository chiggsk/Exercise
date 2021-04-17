rm(list=ls())

library(dplyr)
library(gplots)
library(forcats)
library(tibble)
library(tidyverse)
library(RColorBrewer)
library(varhandle)
library(ggforce)
library(M3C)

####Gastrocnemius####

dataFolder <- setwd("/Users/Chigoziri/MoTrPAC")
phenotypeDataFile <- "phenotype/merged_dmaqc_data2019-10-15.txt"
gastroCountsFile <- "transcriptomics/T55_gastrocnemius/rna-seq/results/MoTrPAC_rsem_genes_count_gastrocnemius_v1.txt"
phenotypeData <- read.table(file=phenotypeDataFile, header = TRUE, sep="\t", quote="", comment.char="")
gastroCounts <- floor(read.csv(file=gastroCountsFile, header = TRUE, row.names = 1, sep="\t", quote="", comment.char=""))

saveFolder <- setwd("/Users/Chigoziri/MoTrPAC/gastroCtrl")

phenotypeData[,1] <- as.character(phenotypeData[,1])

animalGroup <- phenotypeData$animal.key.anirandgroup
animalGroupFactor <- fct_infreq(factor(animalGroup))

names(gastroCounts) <- substring(names(gastroCounts), 2)
sampleIDsGastro <- colnames(gastroCounts)
vialIDs <- as.character(phenotypeData$vial_label)

gastroPhenotypes <- subset(phenotypeData,phenotypeData$vial_label %in% sampleIDsGastro)
ControlIPEPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Control - IPE")
Control7hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Control - 7 hr")
ExerciseIPEPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE")
Exercise0.5hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr")
Exercise1hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr")
Exercise4hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr")
Exercise7hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr")
Exercise24hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr")
Exercise48hrPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondPhenotypes <- subset(gastroPhenotypes,gastroPhenotypes$animal.key.anirandgroup %in% "Control - IPE" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr" |
                    gastroPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondSex <- data.frame(VialLabel=allCondPhenotypes$vial_label,Sex=allCondPhenotypes$animal.registration.sex)
allCondSex[,2][allCondSex[,2]=="1"] <- "Female"
allCondSex[,2][allCondSex[,2]=="2"] <- "Male"


####Split data into experimental groups
setwd("/Users/Chigoziri/MoTrPAC/gastroCtrl")
load("gastro_LRTData.RData")

allCondMatch <- match(coldata$sampleIDs,allCondSex[,1])
allCondArrange <- allCondSex[allCondMatch,]

coldata$Sex <- allCondArrange[,2]
sexFactor <- factor(coldata$Sex,levels=c("Female","Male"))

gastroData <- rld_mat
gastroDataJoin <- rownames_to_column(as.data.frame(gastroData))

gastroColdata <- coldata

plotLabels <- factor(c(rep("Ctrl",sum(coldata$Condition=="Control_IPE")),
                rep("0 hr",sum(coldata$Condition=="Exercise_IPE")),
                rep("0.5 hr",sum(coldata$Condition=="Exercise0.5hr")),
                rep("1 hr",sum(coldata$Condition=="Exercise01hr")),
                rep("4 hr",sum(coldata$Condition=="Exercise04hr")),
                rep("7 hr",sum(coldata$Condition=="Exercise07hr")),
                rep("24 hr",sum(coldata$Condition=="Exercise24hr")),
                rep("48 hr",sum(coldata$Condition=="Exercise48hr"))),
                levels=c("Ctrl","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))

tiff("umap_gastro_sex.tiff", units="in", width=9, height=6, res=300)
umap(mydata=rld_mat, labels = plotLabels,
     axistextsize=20,legendtextsize=16,
     dotsize=3,legendtitle="Time") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  #geom_mark_hull(aes(fill = plotLabels),expand = unit(4, "mm"),concavity = 3) +
  geom_mark_hull(aes(fill = sexFactor, labels = sexFactor),expand = unit(4, "mm"),concavity = 3) +
  labs(fill = "Sex") +
  #guides(fill = FALSE) +
  theme(legend.title = element_text(size = 20),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()


######Heart######

dataFolder <- setwd("/Users/Chigoziri/MoTrPAC")
phenotypeDataFile <- "phenotype/merged_dmaqc_data2019-10-15.txt"
heartCountsFile <- "transcriptomics/T58_heart/rna-seq/results/MoTrPAC_rsem_genes_count_heart_v1.txt"
phenotypeData <- read.table(file=phenotypeDataFile, header = TRUE, sep="\t", quote="", comment.char="")
heartCounts <- floor(read.csv(file=heartCountsFile, header = TRUE, row.names = 1, sep="\t", quote="", comment.char=""))

saveFolder <- setwd("/Users/Chigoziri/MoTrPAC/heartCtrl")

phenotypeData[,1] <- as.character(phenotypeData[,1])

animalGroup <- phenotypeData$animal.key.anirandgroup
animalGroupFactor <- fct_infreq(factor(animalGroup))

names(heartCounts) <- substring(names(heartCounts), 2)
sampleIDsHeart <- colnames(heartCounts)
vialIDs <- as.character(phenotypeData$vial_label)

heartPhenotypes <- subset(phenotypeData,phenotypeData$vial_label %in% sampleIDsHeart)
ControlIPEPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Control - IPE")
Control7hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Control - 7 hr")
ExerciseIPEPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE")
Exercise0.5hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr")
Exercise1hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr")
Exercise4hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr")
Exercise7hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr")
Exercise24hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr")
Exercise48hrPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")


allCondPhenotypes <- subset(heartPhenotypes,heartPhenotypes$animal.key.anirandgroup %in% "Control - IPE" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr" |
                              heartPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondSex <- data.frame(VialLabel=allCondPhenotypes$vial_label,Sex=allCondPhenotypes$animal.registration.sex)
allCondSex[,2][allCondSex[,2]=="1"] <- "Female"
allCondSex[,2][allCondSex[,2]=="2"] <- "Male"

####Split data into experimental groups

setwd("/Users/Chigoziri/MoTrPAC/heartCtrl")
load("heart_LRTData.RData")

allCondMatch <- match(coldata$sampleIDs,allCondSex[,1])
allCondArrange <- allCondSex[allCondMatch,]

coldata$Sex <- allCondArrange[,2]
sexFactor <- factor(coldata$Sex,levels=c("Female","Male"))

heartData <- rld_mat
heartDataJoin <- rownames_to_column(as.data.frame(heartData))

heartColdata <- coldata

plotLabels <- factor(c(rep("Ctrl",sum(coldata$Condition=="Control_IPE")),
                       rep("0 hr",sum(coldata$Condition=="Exercise_IPE")),
                       rep("0.5 hr",sum(coldata$Condition=="Exercise0.5hr")),
                       rep("1 hr",sum(coldata$Condition=="Exercise01hr")),
                       rep("4 hr",sum(coldata$Condition=="Exercise04hr")),
                       rep("7 hr",sum(coldata$Condition=="Exercise07hr")),
                       rep("24 hr",sum(coldata$Condition=="Exercise24hr")),
                       rep("48 hr",sum(coldata$Condition=="Exercise48hr"))),
                     levels=c("Ctrl","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))

tiff("umap_heart.tiff", units="in", width=9, height=6, res=300)
umap(mydata=rld_mat, labels = plotLabels,
     axistextsize=20,legendtextsize=16,
     dotsize=3,legendtitle="Time") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  geom_mark_hull(aes(fill = plotLabels),expand = unit(4, "mm"),concavity = 3) +
  #geom_mark_hull(aes(fill = sexFactor, labels = sexFactor),expand = unit(4, "mm"),concavity = 3) +
  labs(fill = "Time") +
  #guides(fill = FALSE) +
  theme(legend.title = element_text(size = 20),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()


####Liver####

dataFolder <- setwd("/Users/Chigoziri/MoTrPAC")
phenotypeDataFile <- "phenotype/merged_dmaqc_data2019-10-15.txt"
liverCountsFile <- "transcriptomics/T68_liver/rna-seq/results/MoTrPAC_rsem_genes_count_liver_v1.txt"
phenotypeData <- read.table(file=phenotypeDataFile, header = TRUE, sep="\t", quote="", comment.char="")
liverCounts <- floor(read.csv(file=liverCountsFile, header = TRUE, row.names = 1, sep="\t", quote="", comment.char=""))

saveFolder <- setwd("/Users/Chigoziri/MoTrPAC/liverCtrl")

phenotypeData[,1] <- as.character(phenotypeData[,1])

animalGroup <- phenotypeData$animal.key.anirandgroup
animalGroupFactor <- fct_infreq(factor(animalGroup))

names(liverCounts) <- substring(names(liverCounts), 2)
sampleIDsLiver <- colnames(liverCounts)
vialIDs <- as.character(phenotypeData$vial_label)

liverPhenotypes <- subset(phenotypeData,phenotypeData$vial_label %in% sampleIDsLiver)
ControlIPEPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Control - IPE")
Control7hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Control - 7 hr")
ExerciseIPEPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE")
Exercise0.5hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr")
Exercise1hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr")
Exercise4hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr")
Exercise7hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr")
Exercise24hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr")
Exercise48hrPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondPhenotypes <- subset(liverPhenotypes,liverPhenotypes$animal.key.anirandgroup %in% "Control - IPE" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr" |
                              liverPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondSex <- data.frame(VialLabel=allCondPhenotypes$vial_label,Sex=allCondPhenotypes$animal.registration.sex)
allCondSex[,2][allCondSex[,2]=="1"] <- "Female"
allCondSex[,2][allCondSex[,2]=="2"] <- "Male"


####Split data into experimental groups
setwd("/Users/Chigoziri/MoTrPAC/liverCtrl")
load("liver_LRTData.RData")

allCondMatch <- match(coldata$sampleIDs,allCondSex[,1])
allCondArrange <- allCondSex[allCondMatch,]

coldata$Sex <- allCondArrange[,2]
sexFactor <- factor(coldata$Sex,levels=c("Female","Male"))

liverData <- rld_mat
liverDataJoin <- rownames_to_column(as.data.frame(liverData))

liverColdata <- coldata


plotLabels <- factor(c(rep("Ctrl",sum(coldata$Condition=="Control_IPE")),
                       rep("0 hr",sum(coldata$Condition=="Exercise_IPE")),
                       rep("0.5 hr",sum(coldata$Condition=="Exercise0.5hr")),
                       rep("1 hr",sum(coldata$Condition=="Exercise01hr")),
                       rep("4 hr",sum(coldata$Condition=="Exercise04hr")),
                       rep("7 hr",sum(coldata$Condition=="Exercise07hr")),
                       rep("24 hr",sum(coldata$Condition=="Exercise24hr")),
                       rep("48 hr",sum(coldata$Condition=="Exercise48hr"))),
                     levels=c("Ctrl","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))

tiff("umap_liver.tiff", units="in", width=9, height=6, res=300)
umap(mydata=rld_mat, labels = plotLabels,
     axistextsize=20,legendtextsize=16,
     dotsize=3,legendtitle="Time") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  geom_mark_hull(aes(fill = plotLabels),expand = unit(4, "mm"),concavity = 3) +
  #geom_mark_hull(aes(fill = sexFactor, labels = sexFactor),expand = unit(4, "mm"),concavity = 3) +
  labs(fill = "Time") +
  #guides(fill = FALSE) +
  theme(legend.title = element_text(size = 20),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()


####Adipose####

dataFolder <- setwd("/Users/Chigoziri/MoTrPAC")
phenotypeDataFile <- "phenotype/merged_dmaqc_data2019-10-15.txt"
adiposeCountsFile <- "transcriptomics/T70_white_adipose/rna-seq/results/MoTrPAC_rsem_genes_count_white_adipose_v1.txt"
phenotypeData <- read.table(file=phenotypeDataFile, header = TRUE, sep="\t", quote="", comment.char="")
adiposeCounts <- floor(read.csv(file=adiposeCountsFile, header = TRUE, row.names = 1, sep="\t", quote="", comment.char=""))

saveFolder <- setwd("/Users/Chigoziri/MoTrPAC/adiposeCtrl")

phenotypeData[,1] <- as.character(phenotypeData[,1])

animalGroup <- phenotypeData$animal.key.anirandgroup
animalGroupFactor <- fct_infreq(factor(animalGroup))

names(adiposeCounts) <- substring(names(adiposeCounts), 2)
sampleIDsAdipose <- colnames(adiposeCounts)
vialIDs <- as.character(phenotypeData$vial_label)

adiposePhenotypes <- subset(phenotypeData,phenotypeData$vial_label %in% sampleIDsAdipose)
ControlIPEPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Control - IPE")
Control7hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Control - 7 hr")
ExerciseIPEPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - IPE")
Exercise0.5hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr")
Exercise1hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr")
Exercise4hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr")
Exercise7hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr")
Exercise24hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr")
Exercise48hrPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

adiposeGroups <- adiposePhenotypes$animal.key.anirandgroup
adiposeGroupsFactor <- fct_infreq(factor(adiposeGroups))

allCondPhenotypes <- subset(adiposePhenotypes,adiposePhenotypes$animal.key.anirandgroup %in% "Control - IPE" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - IPE" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr" |
                              adiposePhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondSex <- data.frame(VialLabel=allCondPhenotypes$vial_label,Sex=allCondPhenotypes$animal.registration.sex)
allCondSex[,2][allCondSex[,2]=="1"] <- "Female"
allCondSex[,2][allCondSex[,2]=="2"] <- "Male"

setwd("/Users/Chigoziri/MoTrPAC/adiposeCtrl")
load("adipose_LRTData.RData")


allCondMatch <- match(coldata$sampleIDs,allCondSex[,1])
allCondArrange <- allCondSex[allCondMatch,]

coldata$Sex <- allCondArrange[,2]
sexFactor <- factor(coldata$Sex,levels=c("Female","Male"))

adiposeData <- rld_mat
adiposeDataJoin <- rownames_to_column(as.data.frame(adiposeData))

adiposeColdata <- coldata

plotLabels <- factor(c(rep("Ctrl",sum(coldata$Condition=="Control_IPE")),
                       rep("0 hr",sum(coldata$Condition=="Exercise_IPE")),
                       rep("0.5 hr",sum(coldata$Condition=="Exercise0.5hr")),
                       rep("1 hr",sum(coldata$Condition=="Exercise01hr")),
                       rep("4 hr",sum(coldata$Condition=="Exercise04hr")),
                       rep("7 hr",sum(coldata$Condition=="Exercise07hr")),
                       rep("24 hr",sum(coldata$Condition=="Exercise24hr")),
                       rep("48 hr",sum(coldata$Condition=="Exercise48hr"))),
                     levels=c("Ctrl","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))

tiff("umap_adipose_sex.tiff", units="in", width=9, height=6, res=300)
umap(mydata=rld_mat, labels = plotLabels,
     axistextsize=20,legendtextsize=16,
     dotsize=3,legendtitle="Time") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  #geom_mark_hull(aes(fill = plotLabels),expand = unit(4, "mm"),concavity = 3) +
  geom_mark_hull(aes(fill = sexFactor, labels = sexFactor),expand = unit(4, "mm"),concavity = 3) +
  labs(fill = "Sex") +
  #guides(fill = FALSE) +
  theme(legend.title = element_text(size = 20),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()


####Blood####

dataFolder <- setwd("/Users/Chigoziri/MoTrPAC")
phenotypeDataFile <- "phenotype/merged_dmaqc_data2019-10-15.txt"
bloodCountsFile <- "transcriptomics/T30_blood_rna/rna-seq/results/MoTrPAC_rsem_genes_count_PaxGeneRNA_v1.txt"
phenotypeData <- read.table(file=phenotypeDataFile, header = TRUE, sep="\t", quote="", comment.char="")
bloodCounts <- floor(read.csv(file=bloodCountsFile, header = TRUE, row.names = 1, sep="\t", quote="", comment.char=""))

saveFolder <- setwd("/Users/Chigoziri/MoTrPAC/bloodCtrl")

phenotypeData[,1] <- as.character(phenotypeData[,1])

animalGroup <- phenotypeData$animal.key.anirandgroup
animalGroupFactor <- fct_infreq(factor(animalGroup))

names(bloodCounts) <- substring(names(bloodCounts), 2)
sampleIDsBlood <- colnames(bloodCounts)
vialIDs <- as.character(phenotypeData$vial_label)

bloodPhenotypes <- subset(phenotypeData,phenotypeData$vial_label %in% sampleIDsBlood)
ControlIPEPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Control - IPE")
Control7hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Control - 7 hr")
ExerciseIPEPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE")
Exercise0.5hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr")
Exercise1hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr")
Exercise4hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr")
Exercise7hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr")
Exercise24hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr")
Exercise48hrPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

bloodGroups <- bloodPhenotypes$animal.key.anirandgroup
bloodGroupsFactor <- fct_infreq(factor(bloodGroups))

allCondPhenotypes <- subset(bloodPhenotypes,bloodPhenotypes$animal.key.anirandgroup %in% "Control - IPE" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - IPE" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 0.5 hr" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 1 hr" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 4 hr" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 7 hr" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 24 hr" |
                              bloodPhenotypes$animal.key.anirandgroup %in% "Exercise - 48 hr")

allCondSex <- data.frame(VialLabel=allCondPhenotypes$vial_label,Sex=allCondPhenotypes$animal.registration.sex)
allCondSex[,2][allCondSex[,2]=="1"] <- "Female"
allCondSex[,2][allCondSex[,2]=="2"] <- "Male"

setwd("/Users/Chigoziri/MoTrPAC/bloodCtrl")
load("blood_LRTData.RData")


allCondMatch <- match(coldata$sampleIDs,allCondSex[,1])
allCondArrange <- allCondSex[allCondMatch,]

coldata$Sex <- allCondArrange[,2]
sexFactor <- factor(coldata$Sex,levels=c("Female","Male"))

bloodData <- rld_mat
bloodDataJoin <- rownames_to_column(as.data.frame(bloodData))

bloodColdata <- coldata

plotLabels <- factor(c(rep("Ctrl",sum(coldata$Condition=="Control_IPE")),
                       rep("0 hr",sum(coldata$Condition=="Exercise_IPE")),
                       rep("0.5 hr",sum(coldata$Condition=="Exercise0.5hr")),
                       rep("1 hr",sum(coldata$Condition=="Exercise01hr")),
                       rep("4 hr",sum(coldata$Condition=="Exercise04hr")),
                       rep("7 hr",sum(coldata$Condition=="Exercise07hr")),
                       rep("24 hr",sum(coldata$Condition=="Exercise24hr")),
                       rep("48 hr",sum(coldata$Condition=="Exercise48hr"))),
                     levels=c("Ctrl","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))

tiff("umap_blood_sex.tiff", units="in", width=9, height=6, res=300)
umap(mydata=rld_mat, labels = plotLabels,
     axistextsize=20,legendtextsize=16,
     dotsize=3,legendtitle="Time") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  #geom_mark_hull(aes(fill = plotLabels),expand = unit(4, "mm"),concavity = 3) +
  geom_mark_hull(aes(fill = sexFactor, labels = sexFactor),expand = unit(4, "mm"),concavity = 3) +
  labs(fill = "Sex") +
  #guides(fill = FALSE) +
  theme(legend.title = element_text(size = 20),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()


####All Tissues####
data1 <- inner_join(gastroDataJoin,heartDataJoin)
data2 <- inner_join(data1,liverDataJoin)
data3 <- inner_join(data2,adiposeDataJoin)
data4 <- inner_join(data3,bloodDataJoin)
data4 <- column_to_rownames(data4,var="rowname")

plotLabels <- factor(c(rep("Gastrocnemius - Ctrl",sum(gastroColdata$Condition=="Control_IPE")),
                       rep("Gastrocnemius - 0 hr",sum(gastroColdata$Condition=="Exercise_IPE")),
                       rep("Gastrocnemius - 0.5 hr",sum(gastroColdata$Condition=="Exercise0.5hr")),
                       rep("Gastrocnemius - 1 hr",sum(gastroColdata$Condition=="Exercise01hr")),
                       rep("Gastrocnemius - 4 hr",sum(gastroColdata$Condition=="Exercise04hr")),
                       rep("Gastrocnemius - 7 hr",sum(gastroColdata$Condition=="Exercise07hr")),
                       rep("Gastrocnemius - 24 hr",sum(gastroColdata$Condition=="Exercise24hr")),
                       rep("Gastrocnemius - 48 hr",sum(gastroColdata$Condition=="Exercise48hr")),
                       rep("Heart - Ctrl",sum(heartColdata$Condition=="Control_IPE")),
                       rep("Heart - 0 hr",sum(heartColdata$Condition=="Exercise_IPE")),
                       rep("Heart - 0.5 hr",sum(heartColdata$Condition=="Exercise0.5hr")),
                       rep("Heart - 1 hr",sum(heartColdata$Condition=="Exercise01hr")),
                       rep("Heart - 4 hr",sum(heartColdata$Condition=="Exercise04hr")),
                       rep("Heart - 7 hr",sum(heartColdata$Condition=="Exercise07hr")),
                       rep("Heart - 24 hr",sum(heartColdata$Condition=="Exercise24hr")),
                       rep("Heart - 48 hr",sum(heartColdata$Condition=="Exercise48hr")),
                       rep("Liver - Ctrl",sum(liverColdata$Condition=="Control_IPE")),
                       rep("Liver - 0 hr",sum(liverColdata$Condition=="Exercise_IPE")),
                       rep("Liver - 0.5 hr",sum(liverColdata$Condition=="Exercise0.5hr")),
                       rep("Liver - 1 hr",sum(liverColdata$Condition=="Exercise01hr")),
                       rep("Liver - 4 hr",sum(liverColdata$Condition=="Exercise04hr")),
                       rep("Liver - 7 hr",sum(liverColdata$Condition=="Exercise07hr")),
                       rep("Liver - 24 hr",sum(liverColdata$Condition=="Exercise24hr")),
                       rep("Liver - 48 hr",sum(liverColdata$Condition=="Exercise48hr")),
                       rep("Adipose - Ctrl",sum(adiposeColdata$Condition=="Control_IPE")),
                       rep("Adipose - 0 hr",sum(adiposeColdata$Condition=="Exercise_IPE")),
                       rep("Adipose - 0.5 hr",sum(adiposeColdata$Condition=="Exercise0.5hr")),
                       rep("Adipose - 1 hr",sum(adiposeColdata$Condition=="Exercise01hr")),
                       rep("Adipose - 4 hr",sum(adiposeColdata$Condition=="Exercise04hr")),
                       rep("Adipose - 7 hr",sum(adiposeColdata$Condition=="Exercise07hr")),
                       rep("Adipose - 24 hr",sum(adiposeColdata$Condition=="Exercise24hr")),
                       rep("Adipose - 48 hr",sum(adiposeColdata$Condition=="Exercise48hr")),
                       rep("Blood - Ctrl",sum(bloodColdata$Condition=="Control_IPE")),
                       rep("Blood - 0 hr",sum(bloodColdata$Condition=="Exercise_IPE")),
                       rep("Blood - 0.5 hr",sum(bloodColdata$Condition=="Exercise0.5hr")),
                       rep("Blood - 1 hr",sum(bloodColdata$Condition=="Exercise01hr")),
                       rep("Blood - 4 hr",sum(bloodColdata$Condition=="Exercise04hr")),
                       rep("Blood - 7 hr",sum(bloodColdata$Condition=="Exercise07hr")),
                       rep("Blood - 24 hr",sum(bloodColdata$Condition=="Exercise24hr")),
                       rep("Blood - 48 hr",sum(bloodColdata$Condition=="Exercise48hr"))))
                     #levels=c("Ctrl","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))

tissueLabels <- factor(c(rep("Gastrocnemius",ncol(gastroData)),
                         rep("Heart",ncol(heartData)),
                         rep("Liver",ncol(liverData)),
                         rep("Adipose",ncol(adiposeData)),
                         rep("Blood",ncol(bloodData))),
                       levels=c("Gastrocnemius","Heart","Liver","Adipose","Blood"))

setwd("/Users/Chigoziri/MoTrPAC")

tiff("tsne_allTissues.tiff", units="in", width=11, height=6, res=300)
tsne(mydata=data4, labels = tissueLabels,
     axistextsize=20,legendtextsize=16,
     dotsize=3,legendtitle="Tissue") +
  xlab("t-SNE1") +
  ylab("t-SNE2") +
  #geom_mark_hull(aes(fill = plotLabels),expand = unit(2, "mm"),concavity = 3) +
  geom_mark_ellipse(aes(fill = tissueLabels, labels = tissueLabels),expand = unit(3, "mm")) +
  labs(fill = "Tissue") +
  #guides(fill = FALSE) +
  theme(legend.title = element_text(size = 20),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()

