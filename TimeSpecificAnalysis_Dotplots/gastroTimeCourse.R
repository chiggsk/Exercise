
rm(list=ls())

library(dplyr)
library(gplots)
library(ggplot2)
library(forcats)
library(tibble)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(biomaRt)
library(sva)
library(varhandle)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(clusterProfiler)
library(DOSE)
library(xlsx)
library(openxlsx)

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

gastroGroups <- gastroPhenotypes$animal.key.anirandgroup
gastroGroupsFactor <- fct_infreq(factor(gastroGroups))

####Split data into experimental groups

ControlIPEIDs <- as.character(intersect(colnames(gastroCounts), ControlIPEPhenotypes$vial_label))
ControlIPEData <- dplyr::select(gastroCounts, ControlIPEIDs)

Control7hrIDs <- as.character(intersect(colnames(gastroCounts), Control7hrPhenotypes$vial_label))
Control7hrData <- dplyr::select(gastroCounts, Control7hrIDs)

ExerciseIPEIDs <- as.character(intersect(colnames(gastroCounts), ExerciseIPEPhenotypes$vial_label))
ExerciseIPEData <- dplyr::select(gastroCounts, ExerciseIPEIDs)

Exercise0.5hrIDs <- as.character(intersect(colnames(gastroCounts), Exercise0.5hrPhenotypes$vial_label))
Exercise0.5hrData <- dplyr::select(gastroCounts, Exercise0.5hrIDs)

Exercise1hrIDs <- as.character(intersect(colnames(gastroCounts), Exercise1hrPhenotypes$vial_label))
Exercise1hrData <- dplyr::select(gastroCounts, Exercise1hrIDs)

Exercise4hrIDs <- as.character(intersect(colnames(gastroCounts), Exercise4hrPhenotypes$vial_label))
Exercise4hrData <- dplyr::select(gastroCounts, Exercise4hrIDs)

Exercise7hrIDs <- as.character(intersect(colnames(gastroCounts), Exercise7hrPhenotypes$vial_label))
Exercise7hrData <- dplyr::select(gastroCounts, Exercise7hrIDs)

Exercise24hrIDs <- as.character(intersect(colnames(gastroCounts), Exercise24hrPhenotypes$vial_label))
Exercise24hrData <- dplyr::select(gastroCounts, Exercise24hrIDs)

Exercise48hrIDs <- as.character(intersect(colnames(gastroCounts), Exercise48hrPhenotypes$vial_label))
Exercise48hrData <- dplyr::select(gastroCounts, Exercise48hrIDs)

###################

#Prepare data for comparisons to control
ControlIPEData_Join <- rownames_to_column(as.data.frame(ControlIPEData), var = "Gene ID")
ExerciseIPEData_Join <- rownames_to_column(as.data.frame(ExerciseIPEData), var = "Gene ID")
Exercise0.5hrData_Join <- rownames_to_column(as.data.frame(Exercise0.5hrData), var = "Gene ID")
Exercise1hrData_Join <- rownames_to_column(as.data.frame(Exercise1hrData), var = "Gene ID")
Exercise4hrData_Join <- rownames_to_column(as.data.frame(Exercise4hrData), var = "Gene ID")
Exercise7hrData_Join <- rownames_to_column(as.data.frame(Exercise7hrData), var = "Gene ID")
Exercise24hrData_Join <- rownames_to_column(as.data.frame(Exercise24hrData), var = "Gene ID")
Exercise48hrData_Join <- rownames_to_column(as.data.frame(Exercise48hrData), var = "Gene ID")

####0hr-Control Comparison####

DataCtrlIPE <- inner_join(ControlIPEData_Join, ExerciseIPEData_Join)
DataCtrlIPE <- column_to_rownames(DataCtrlIPE, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exerciseIPENames <- sprintf("ExerciseIPE%s",seq(1:nrow(ExerciseIPEPhenotypes)))
samples <- c(controlIPENames,exerciseIPENames)
sampleIDs <- colnames(DataCtrlIPE)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("ExerciseIPE",nrow(ExerciseIPEPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrlIPE,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","ExerciseIPE"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

save(res, file = "gastro_IPECtrl.RData")

####0.5hr-Control Comparison####

DataCtrl0.5hr <- inner_join(ControlIPEData_Join, Exercise0.5hrData_Join)
DataCtrl0.5hr <- column_to_rownames(DataCtrl0.5hr, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exercise0.5hrNames <- sprintf("Exercise0.5hr%s",seq(1:nrow(Exercise0.5hrPhenotypes)))
samples <- c(controlIPENames,exercise0.5hrNames)
sampleIDs <- colnames(DataCtrl0.5hr)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("Exercise0.5hr",nrow(Exercise0.5hrPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrl0.5hr,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","Exercise0.5hr"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

save(res, file = "gastro_0.5hrCtrl.RData")

####1hr-Control Comparison####

DataCtrl1hr <- inner_join(ControlIPEData_Join, Exercise1hrData_Join)
DataCtrl1hr <- column_to_rownames(DataCtrl1hr, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exercise1hrNames <- sprintf("Exercise1hr%s",seq(1:nrow(Exercise1hrPhenotypes)))
samples <- c(controlIPENames,exercise1hrNames)
sampleIDs <- colnames(DataCtrl1hr)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("Exercise1hr",nrow(Exercise1hrPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrl1hr,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","Exercise1hr"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

save(res, file = "gastro_1hrCtrl.RData")

####4hr-Control Comparison####

DataCtrl4hr <- inner_join(ControlIPEData_Join, Exercise4hrData_Join)
DataCtrl4hr <- column_to_rownames(DataCtrl4hr, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exercise4hrNames <- sprintf("Exercise4hr%s",seq(1:nrow(Exercise4hrPhenotypes)))
samples <- c(controlIPENames,exercise4hrNames)
sampleIDs <- colnames(DataCtrl4hr)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("Exercise4hr",nrow(Exercise4hrPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrl4hr,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","Exercise4hr"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

save(res, file = "gastro_4hrCtrl.RData")

####7hr-Control Comparison####

DataCtrl7hr <- inner_join(ControlIPEData_Join, Exercise7hrData_Join)
DataCtrl7hr <- column_to_rownames(DataCtrl7hr, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exercise7hrNames <- sprintf("Exercise7hr%s",seq(1:nrow(Exercise7hrPhenotypes)))
samples <- c(controlIPENames,exercise7hrNames)
sampleIDs <- colnames(DataCtrl7hr)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("Exercise7hr",nrow(Exercise7hrPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrl7hr,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","Exercise7hr"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

save(res, file = "gastro_7hrCtrl.RData")

####24hr-Control Comparison####

DataCtrl24hr <- inner_join(ControlIPEData_Join, Exercise24hrData_Join)
DataCtrl24hr <- column_to_rownames(DataCtrl24hr, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exercise24hrNames <- sprintf("Exercise24hr%s",seq(1:nrow(Exercise24hrPhenotypes)))
samples <- c(controlIPENames,exercise24hrNames)
sampleIDs <- colnames(DataCtrl24hr)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("Exercise24hr",nrow(Exercise24hrPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrl24hr,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","Exercise24hr"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05) #Alpha by default is set to 0.1

save(res, file = "gastro_24hrCtrl.RData")

####48hr-Control Comparison####

DataCtrl48hr <- inner_join(ControlIPEData_Join, Exercise48hrData_Join)
DataCtrl48hr <- column_to_rownames(DataCtrl48hr, var = "Gene ID")

##Count matrix input
controlIPENames <- sprintf("ControlIPE%s",seq(1:nrow(ControlIPEPhenotypes)))
exercise48hrNames <- sprintf("Exercise48hr%s",seq(1:nrow(Exercise48hrPhenotypes)))
samples <- c(controlIPENames,exercise48hrNames)
sampleIDs <- colnames(DataCtrl48hr)
condition <- c(rep("ControlIPE",nrow(ControlIPEPhenotypes)),rep("Exercise48hr",nrow(Exercise48hrPhenotypes)))


coldata <- data.frame(samples,sampleIDs,condition)

dds <- DESeqDataSetFromMatrix(countData = DataCtrl48hr,
                              colData = coldata,
                              design = ~ condition)

##Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Factor levels for condition
dds$condition <- factor(dds$condition, levels = c("ControlIPE","Exercise48hr"))
dds$condition <- relevel(dds$condition, "ControlIPE")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

save(res, file = "gastro_48hrCtrl.RData")

####Perform KEGG enrichment for all time points and put data together for dot plots####

setwd("/Users/Chigoziri/MoTrPAC/gastroCtrl")
load("gastro_IPECtrl.RData")

foldchanges_IPE <- res$log2FoldChange
padj_IPE <- res$padj
names(foldchanges_IPE) <- res$entrez
names(padj_IPE) <- res$entrez

#Set padj and logFC cutoffs, and filter genes based on these cutoffs
padj.cutoff <- 0.05
FC.cutoff <- 1

sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Create list with length of all time points
ekeggTimeCourseGeneList <- vector(mode = "list", length = 7)
names(ekeggTimeCourseGeneList) <- c("0hr","0.5hr","1hr","4hr","7hr","24hr","48hr")
sheetTitles <- names(ekeggTimeCourseGeneList)

#Create workbooks to save  enriched pathways
ekeggTimeCourseWB <- createWorkbook()
ekeggEnrichedPathwaysWB <- createWorkbook()
addWorksheet(ekeggEnrichedPathwaysWB, "Group Pathways")

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekeggIPE <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekeggIPEMat <- as.data.frame(ekeggIPE)
ekeggGeneID <- ekeggIPEMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "0hr"
  ekeggTimeCourseList[[name]] <- ekeggIPE
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekeggIPEMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekeggIPEMat$KEGGDescription[j] <- paste(ekeggIPEMat$ID[j], ekeggIPEMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggIPEMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[1])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[1], x = ekeggIPEMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 1)

ekeggIPEMat$medianFC <- medianFC

Time <- as.data.frame(rep("0 hr",nrow(ekeggIPEMat)))
data <- data.frame(ekeggIPEMat$Description,ekeggIPEMat$GeneRatio,ekeggIPEMat$p.adjust,unlist(ekeggIPEMat$medianFC),Time)

#In order to have value for dot plot in case 0 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description=NA,GeneRatio=NA,p.adjust=NA,medianFC=NA,Time="0 hr")
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

####

setwd("/Users/Chigoziri/MoTrPAC/gastro0.5hr")
load("gastro_0.5hrCtrl.RData")

foldchanges_0.5hr <- res$log2FoldChange
padj_0.5hr <- res$padj
names(foldchanges_0.5hr) <- res$entrez
names(padj_0.5hr) <- res$entrez

#Filter genes based on padj and logFC cutoffs
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekegg0.5hr <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekegg0.5hrMat <- as.data.frame(ekegg0.5hr)
ekeggGeneID <- ekegg0.5hrMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "0.5hr"
  ekeggTimeCourseList[[name]] <- ekegg0.5hr
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekegg0.5hrMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekegg0.5hrMat$KEGGDescription[j] <- paste(ekegg0.5hrMat$ID[j], ekegg0.5hrMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekegg0.5hrMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[2])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[2], x = ekegg0.5hrMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 2)

ekegg0.5hrMat$medianFC <- medianFC

Time <- as.data.frame(rep("0.5 hr",nrow(ekegg0.5hr)))
data <- data.frame(ekegg0.5hr$Description,ekegg0.5hr$GeneRatio,ekegg0.5hr$p.adjust,unlist(ekegg0.5hrMat$medianFC),Time)

#In order to have value for dot plot in case 0.5 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description=NA,GeneRatio=NA,p.adjust=NA,medianFC=NA,Time="0.5 hr")
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

####

setwd("/Users/Chigoziri/MoTrPAC/gastro1hr")
load("gastro_1hrCtrl.RData")

foldchanges_1hr <- res$log2FoldChange
padj_1hr <- res$padj
names(foldchanges_1hr) <- res$entrez
names(padj_1hr) <- res$entrez

#Filter genes based on padj and logFC cutoffs
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekegg1hr <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekegg1hrMat <- as.data.frame(ekegg1hr)
ekeggGeneID <- ekegg1hrMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "1hr"
  ekeggTimeCourseList[[name]] <- ekegg1hr
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekegg1hrMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekegg1hrMat$KEGGDescription[j] <- paste(ekegg1hrMat$ID[j], ekegg1hrMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekegg1hrMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[3])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[3], x = ekegg1hrMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 3)

ekegg1hrMat$medianFC <- medianFC

Time <- as.data.frame(rep("1 hr",nrow(ekegg1hr)))
data <- data.frame(ekegg1hr$Description,ekegg1hr$GeneRatio,ekegg1hr$p.adjust,unlist(ekegg1hrMat$medianFC),Time)

#In order to have value for dot plot in case 1 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description=NA,GeneRatio=NA,p.adjust=NA,medianFC=NA,Time="1 hr")
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

####

setwd("/Users/Chigoziri/MoTrPAC/gastro4hr")
load("gastro_4hrCtrl.RData")

foldchanges_4hr <- res$log2FoldChange
padj_4hr <- res$padj
names(foldchanges_4hr) <- res$entrez
names(padj_4hr) <- res$entrez

#Filter genes based on padj and logFC cutoffs
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekegg4hr <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekegg4hrMat <- as.data.frame(ekegg4hr)
ekeggGeneID <- ekegg4hrMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "4hr"
  ekeggTimeCourseList[[name]] <- ekegg4hr
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekegg4hrMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekegg4hrMat$KEGGDescription[j] <- paste(ekegg4hrMat$ID[j], ekegg4hrMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekegg4hrMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[4])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[4], x = ekegg4hrMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 4)

ekegg4hrMat$medianFC <- medianFC

Time <- as.data.frame(rep("4 hr",nrow(ekegg4hr)))
data <- data.frame(ekegg4hr$Description,ekegg4hr$GeneRatio,ekegg4hr$p.adjust,unlist(ekegg4hrMat$medianFC),Time)

#In order to have value for dot plot in case 4 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description=NA,GeneRatio=NA,p.adjust=NA,medianFC=NA,Time="4 hr")
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

####

setwd("/Users/Chigoziri/MoTrPAC/gastro7hr")
load("gastro_7hrCtrl.RData")

foldchanges_7hr <- res$log2FoldChange
padj_7hr <- res$padj
names(foldchanges_7hr) <- res$entrez
names(padj_7hr) <- res$entrez

#Filter genes based on padj and logFC cutoffs
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekegg7hr <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekegg7hrMat <- as.data.frame(ekegg7hr)
ekeggGeneID <- ekegg7hrMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "7hr"
  ekeggTimeCourseList[[name]] <- ekegg7hr
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekegg7hrMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekegg7hrMat$KEGGDescription[j] <- paste(ekegg7hrMat$ID[j], ekegg7hrMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekegg7hrMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[5])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[5], x = ekegg7hrMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 5)

ekegg7hrMat$medianFC <- medianFC

Time <- as.data.frame(rep("7 hr",nrow(ekegg7hr)))
data <- data.frame(ekegg7hr$Description,ekegg7hr$GeneRatio,ekegg7hr$p.adjust,unlist(ekegg7hrMat$medianFC),Time)

#In order to have value for dot plot in case 7 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description=NA,GeneRatio=NA,p.adjust=NA,medianFC=NA,Time="7 hr")
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")
####

setwd("/Users/Chigoziri/MoTrPAC/gastro24hr")
load("gastro_24hrCtrl.RData")

foldchanges_24hr <- res$log2FoldChange
padj_24hr <- res$padj
names(foldchanges_24hr) <- res$entrez
names(padj_24hr) <- res$entrez

#Filter genes based on padj and logFC cutoffs
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekegg24hr <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekegg24hrMat <- as.data.frame(ekegg24hr)
ekeggGeneID <- ekegg24hrMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "24hr"
  ekeggTimeCourseList[[name]] <- ekegg24hr
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekegg24hrMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekegg24hrMat$KEGGDescription[j] <- paste(ekegg24hrMat$ID[j], ekegg24hrMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekegg24hrMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[6])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[6], x = ekegg24hrMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 6)

ekegg24hrMat$medianFC <- medianFC

Time <- as.data.frame(rep("24 hr",nrow(ekegg24hr)))
data <- data.frame(ekegg24hr$Description,ekegg24hr$GeneRatio,ekegg24hr$p.adjust,unlist(ekegg24hrMat$medianFC),Time)

#In order to have value for dot plot in case 24 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description="IL-17 signaling pathway",GeneRatio=NA,p.adjust=NA,medianFC=NA,Time="24 hr")
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

####

setwd("/Users/Chigoziri/MoTrPAC/gastro48hr")
load("gastro_48hrCtrl.RData")

foldchanges_48hr <- res$log2FoldChange
padj_48hr <- res$padj
names(foldchanges_48hr) <- res$entrez
names(padj_48hr) <- res$entrez

#Filter genes based on padj and logFC cutoffs
sig_res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff) %>%
  filter(abs(log2FoldChange) >= FC.cutoff)

#Perform KEGG enrichment and generate dotplot data for each time point
gene <- sig_res$entrez #Input only genes which fell into a particular cluster for that enrichment
ekegg48hr <- enrichKEGG(gene         = gene,
                       organism     = 'rno',
                       pvalueCutoff = 0.05)

ekeggTimeCourseList <- list()

ekegg48hrMat <- as.data.frame(ekegg48hr)
ekeggGeneID <- ekegg48hrMat$geneID
ekeggSplit <- strsplit(ekeggGeneID, "/")

ekeggSymbolList <- vector(mode = "list", length = length(ekeggSplit))
ekeggENTREZList <- vector(mode = "list", length = length(ekeggSplit))
ekeggFCList <- vector(mode = "list", length = length(ekeggSplit))
medianFC <- vector(mode = "list", length = length(ekeggSplit))

for (j in 1:length(ekeggSplit)) {
  name <- "48hr"
  ekeggTimeCourseList[[name]] <- ekegg48hr
  
  ekeggUnlist <- unlist(ekeggSplit[[j]])
  ekeggENTREZList[[j]] <- ekeggUnlist
  ekeggFCList[[j]] <- sig_res$log2FoldChange[sig_res$entrez %in% ekeggUnlist]
  medianFC[[j]] <- median(ekeggFCList[[j]])
  ekeggSymbol = mapIds(org.Rn.eg.db,
                       keys=ekeggUnlist, 
                       column="SYMBOL",
                       keytype="ENTREZID",
                       multiVals="first")
  
  ekeggSymbolList[[j]] <- ekeggSymbol
  geneUnlisted <- unlist(ekeggSymbol)
  ekegg48hrMat$GeneSymbol[j] <- paste(geneUnlisted, collapse = "/")
  ekegg48hrMat$KEGGDescription[j] <- paste(ekegg48hrMat$ID[j], ekegg48hrMat$Description[j])
  
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekegg48hrMat$Description[j])
  ekeggTimeCourseGeneList[[name]] <- c(as.character(ekeggTimeCourseGeneList[[name]]),ekeggSymbol)
}

# Add some sheets to the workbook
addWorksheet(ekeggTimeCourseWB, sheetTitles[7])
# Write the data to the sheets
writeData(ekeggTimeCourseWB, sheet = sheetTitles[7], x = ekegg48hrMat)
writeData(ekeggEnrichedPathwaysWB, sheet="Group Pathways", ekeggTimeCourseGeneList[[name]], colNames = F,rowNames = T,startRow = 2,startCol = 7)

ekegg48hrMat$medianFC <- medianFC

Time <- as.data.frame(rep("48 hr",nrow(ekegg48hr)))
data <- data.frame(ekegg48hr$Description,ekegg48hr$GeneRatio,ekegg48hr$p.adjust,unlist(ekegg48hrMat$medianFC),Time)

#In order to have value for dot plot in case 48 hr shows no enrichment
if (nrow(data)==0) {
  data <- data.frame(Description=c("IL-17 signaling pathway","Chagas disease","Chagas disease"),GeneRatio=NA,p.adjust=NA,medianFC=NA,Time=c("48 hr","24 hr","48 hr"))
}

names(data) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

#Append data to existing dot plot data
ekeggDotPlotData <- rbind(ekeggDotPlotData,data)
names(ekeggDotPlotData) <- c("Description","GeneRatio","p.adjust","medianFC","Time")

ekeggDotPlotData$Time <- factor(ekeggDotPlotData$Time,levels=c("0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))
ekeggDotPlotData$Description <- as.factor(ekeggDotPlotData$Description)

########
setwd("/Users/Chigoziri/MoTrPAC/gastroCtrl")

saveWorkbook(ekeggTimeCourseWB, "gastroTimeCourse_enrichKEGG_FCcutoff1_0.05.xlsx")
saveWorkbook(ekeggEnrichedPathwaysWB,"gastroTimeCourse_enrichedPathways_FCcutoff1.xlsx")

########

#Used to specify enriched pathways specified for supplemental figures
ekeggDotPlotData_filt <- ekeggDotPlotData[ ( !grepl("Epstein-Barr virus infection", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Human cytomegalovirus infection", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Human papillomavirus infection", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Measles", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Salmonella infection", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Yersinia infection", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Hepatitis B", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Hepatitis C", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Transcriptional misregulation in cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Kaposi sarcoma-associated herpesvirus infection" , ekeggDotPlotData$Description)) & 
                                             ( !grepl("Longevity regulating pathway - multiple species", ekeggDotPlotData$Description)) & 
                                             ( !grepl("MicroRNAs in cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Prostate cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Basal cell carcinoma", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Acute myeloid leukemia", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Breast cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Proteoglycans in cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Renal cell carcinoma", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Basal cell carcinoma", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Bladder cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Colorectal cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Basal cell carcinoma", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Cocaine addiction", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Inflammatory bowel disease", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Rheumatoid arthritis", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Malaria", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Legionellosis", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Leishmaniasis", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Thyroid cancer", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Tuberculosis", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Pertussis", ekeggDotPlotData$Description)) & 
                                             ( !grepl("African trypanosomiasis", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Amphetamine addiction", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Chagas disease", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Influenza A", ekeggDotPlotData$Description)) & 
                                             ( !grepl("Amoebiasis", ekeggDotPlotData$Description)), ] 


ekeggDotPlotData_supp <- ekeggDotPlotData[ ( grepl("Epstein-Barr virus infection", ekeggDotPlotData$Description)) | 
                                             ( grepl("Human cytomegalovirus infection", ekeggDotPlotData$Description)) | 
                                             ( grepl("Human papillomavirus infection", ekeggDotPlotData$Description)) | 
                                             ( grepl("Measles", ekeggDotPlotData$Description)) | 
                                             ( grepl("Salmonella infection", ekeggDotPlotData$Description)) | 
                                             ( grepl("Yersinia infection", ekeggDotPlotData$Description)) | 
                                             ( grepl("Hepatitis B", ekeggDotPlotData$Description)) | 
                                             ( grepl("Hepatitis C", ekeggDotPlotData$Description)) | 
                                             ( grepl("Transcriptional misregulation in cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Kaposi sarcoma-associated herpesvirus infection" , ekeggDotPlotData$Description)) | 
                                             ( grepl("Longevity regulating pathway - multiple species", ekeggDotPlotData$Description)) | 
                                             ( grepl("MicroRNAs in cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Prostate cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Basal cell carcinoma", ekeggDotPlotData$Description)) | 
                                             ( grepl("Acute myeloid leukemia", ekeggDotPlotData$Description)) | 
                                             ( grepl("Breast cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Proteoglycans in cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Renal cell carcinoma", ekeggDotPlotData$Description)) | 
                                             ( grepl("Basal cell carcinoma", ekeggDotPlotData$Description)) | 
                                             ( grepl("Bladder cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Colorectal cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Basal cell carcinoma", ekeggDotPlotData$Description)) | 
                                             ( grepl("Cocaine addiction", ekeggDotPlotData$Description)) | 
                                             ( grepl("Inflammatory bowel disease", ekeggDotPlotData$Description)) | 
                                             ( grepl("Rheumatoid arthritis", ekeggDotPlotData$Description)) | 
                                             ( grepl("Malaria", ekeggDotPlotData$Description)) | 
                                             ( grepl("Legionellosis", ekeggDotPlotData$Description)) | 
                                             ( grepl("Leishmaniasis", ekeggDotPlotData$Description)) | 
                                             ( grepl("Thyroid cancer", ekeggDotPlotData$Description)) | 
                                             ( grepl("Tuberculosis", ekeggDotPlotData$Description)) | 
                                             ( grepl("Pertussis", ekeggDotPlotData$Description)) | 
                                             ( grepl("African trypanosomiasis", ekeggDotPlotData$Description)) | 
                                             ( grepl("Amphetamine addiction", ekeggDotPlotData$Description)) | 
                                             ( grepl("Chagas disease", ekeggDotPlotData$Description)) | 
                                             ( grepl("Influenza A", ekeggDotPlotData$Description)) | 
                                             ( grepl("Amoebiasis", ekeggDotPlotData$Description)), ] 

ekeggDotPlotData_supp$GeneRatio <- sapply(ekeggDotPlotData_supp$GeneRatio, function(x) eval(parse(text=x))) #Convert GeneRatio from fraction character to decimal
ekeggDotPlotData_filt$GeneRatio <- sapply(ekeggDotPlotData_filt$GeneRatio, function(x) eval(parse(text=x))) #Convert GeneRatio from fraction character to decimal

setwd("/Users/Chigoziri/MoTrPAC/gastroCtrl")

#Save dot plots as .tiff
tiff("gastroTimeDirection_DotPlot_FCcutoff1_supp_0.05.tiff", units="in", width=11, height=8, res=300)

gastroDotPlot <- ggplot(ekeggDotPlotData_supp, aes(x= Time, y=reorder(Description, dplyr::desc(Description)), size=GeneRatio, color=medianFC, group=Time)) + geom_point(alpha = 0.8) + 
  theme_bw(base_size = 18) + guides(color=guide_colourbar(order=1),size=guide_legend(order=2)) + ylab("Description")
gastroDotPlot = gastroDotPlot+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(ekeggDotPlotData$medianFC),max(ekeggDotPlotData$medianFC)))
gastroDotPlot+scale_size(range = c(2, 8))

dev.off()

####

tiff("gastroTimeDirection_filt_DotPlot_FCcutoff1_0.05.tiff", units="in", width=12, height=8, res=300)

gastroDotPlot <- ggplot(ekeggDotPlotData_filt, aes(x= Time, y=reorder(Description, dplyr::desc(Description)), size=GeneRatio, color=medianFC, group=Time)) + geom_point(alpha = 0.8) + 
  theme_bw(base_size = 18) + guides(color=guide_colourbar(order=1),size=guide_legend(order=2)) + ylab("Description")
gastroDotPlot = gastroDotPlot+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(ekeggDotPlotData_filt$medianFC),max(ekeggDotPlotData_filt$medianFC)))      
gastroDotPlot+scale_size(range = c(2, 8))

dev.off()
