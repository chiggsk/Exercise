rm(list=ls())

library(dplyr)
library(gplots)
library(ggplot2)
library(forcats)
library(tibble)
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)
library(biomaRt)
library(sva)
library(varhandle)
library(AnnotationDbi)
library(org.Rn.eg.db)
library(circlize)

####NOTE: This script requires data generated from LRT analyses for each tissue####

setwd("/Users/Chigoziri/MoTrPAC/gastroCtrl")
load("gastro_LRTData.RData")

#Identify significant genes which meet cutoff
padj.cutoff <- 0.05

sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

# Subset results for faster cluster finding (for classroom demo purposes)
gastro_clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=2000)


#### Obtain mean gene expression data for all pertinent genes at all time points####

# Obtain rlog values for those significant genes
cluster_rld <- rld_mat[gastro_clustering_sig_genes$gene, ]

#Z-score scale transformation
rld_scale <- scale(t(cluster_rld))
rld_scale <- t(rld_scale)

controlIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Control_IPE"]] 
controlIPEMeans <- data.frame(rowMeans(controlIPEVals))
colnames(controlIPEMeans) <- "GControl_IPE"
controlIPEJoin <- rownames_to_column(controlIPEMeans, var="geneNames")

exerciseIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise_IPE"]] 
exerciseIPEMeans <- data.frame(rowMeans(exerciseIPEVals))
colnames(exerciseIPEMeans) <- "GExercise_IPE"
exerciseIPEJoin <- rownames_to_column(exerciseIPEMeans, var="geneNames")

exercise0.5hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise0.5hr"]] 
exercise0.5hrMeans <- data.frame(rowMeans(exercise0.5hrVals))
colnames(exercise0.5hrMeans) <- "GExercise0.5hr"
exercise0.5hrJoin <- rownames_to_column(exercise0.5hrMeans, var="geneNames")

exercise1hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise01hr"]] 
exercise1hrMeans <- data.frame(rowMeans(exercise1hrVals))
colnames(exercise1hrMeans) <- "GExercise1hr"
exercise1hrJoin <- rownames_to_column(exercise1hrMeans, var="geneNames")

exercise4hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise04hr"]] 
exercise4hrMeans <- data.frame(rowMeans(exercise4hrVals))
colnames(exercise4hrMeans) <- "GExercise4hr"
exercise4hrJoin <- rownames_to_column(exercise4hrMeans, var="geneNames")

exercise7hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise07hr"]] 
exercise7hrMeans <- data.frame(rowMeans(exercise7hrVals))
colnames(exercise7hrMeans) <- "GExercise7hr"
exercise7hrJoin <- rownames_to_column(exercise7hrMeans, var="geneNames")

exercise24hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise24hr"]] 
exercise24hrMeans <- data.frame(rowMeans(exercise24hrVals))
colnames(exercise24hrMeans) <- "GExercise24hr"
exercise24hrJoin <- rownames_to_column(exercise24hrMeans, var="geneNames")

exercise48hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise48hr"]] 
exercise48hrMeans <- data.frame(rowMeans(exercise48hrVals))
colnames(exercise48hrMeans) <- "GExercise48hr"
exercise48hrJoin <- rownames_to_column(exercise48hrMeans, var="geneNames")

#Combine data for all time points into a single data frame
dataJoin1 <- inner_join(controlIPEJoin, exerciseIPEJoin)
dataJoin2 <- inner_join(dataJoin1, exercise0.5hrJoin)
dataJoin3 <- inner_join(dataJoin2, exercise1hrJoin)
dataJoin4 <- inner_join(dataJoin3, exercise4hrJoin)
dataJoin5 <- inner_join(dataJoin4, exercise7hrJoin)
dataJoin6 <- inner_join(dataJoin5, exercise24hrJoin)
dataJoin7 <- inner_join(dataJoin6, exercise48hrJoin)
gastro_dataAllTimes <- column_to_rownames(dataJoin7, var = "geneNames") # Remove gene ID column
gastro_dataAllTimes <- as.matrix(gastro_dataAllTimes) 

########

setwd("/Users/Chigoziri/MoTrPAC/heartCtrl")
load("heart_LRTData.RData")

#Identify significant genes which meet cutoff
padj.cutoff <- 0.05

sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

# Subset results for faster cluster finding (for classroom demo purposes)
heart_clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=2000)


#### Obtain mean gene expression data for all pertinent genes at all time points####

# Obtain rlog values for those significant genes
cluster_rld <- rld_mat[heart_clustering_sig_genes$gene, ]

#Z-score scale transformation
rld_scale <- scale(t(cluster_rld))
rld_scale <- t(rld_scale)

controlIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Control_IPE"]] 
controlIPEMeans <- data.frame(rowMeans(controlIPEVals))
colnames(controlIPEMeans) <- "HControl_IPE"
controlIPEJoin <- rownames_to_column(controlIPEMeans, var="geneNames")

exerciseIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise_IPE"]] 
exerciseIPEMeans <- data.frame(rowMeans(exerciseIPEVals))
colnames(exerciseIPEMeans) <- "HExercise_IPE"
exerciseIPEJoin <- rownames_to_column(exerciseIPEMeans, var="geneNames")

exercise0.5hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise0.5hr"]] 
exercise0.5hrMeans <- data.frame(rowMeans(exercise0.5hrVals))
colnames(exercise0.5hrMeans) <- "HExercise0.5hr"
exercise0.5hrJoin <- rownames_to_column(exercise0.5hrMeans, var="geneNames")

exercise1hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise01hr"]] 
exercise1hrMeans <- data.frame(rowMeans(exercise1hrVals))
colnames(exercise1hrMeans) <- "HExercise1hr"
exercise1hrJoin <- rownames_to_column(exercise1hrMeans, var="geneNames")

exercise4hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise04hr"]] 
exercise4hrMeans <- data.frame(rowMeans(exercise4hrVals))
colnames(exercise4hrMeans) <- "HExercise4hr"
exercise4hrJoin <- rownames_to_column(exercise4hrMeans, var="geneNames")

exercise7hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise07hr"]] 
exercise7hrMeans <- data.frame(rowMeans(exercise7hrVals))
colnames(exercise7hrMeans) <- "HExercise7hr"
exercise7hrJoin <- rownames_to_column(exercise7hrMeans, var="geneNames")

exercise24hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise24hr"]] 
exercise24hrMeans <- data.frame(rowMeans(exercise24hrVals))
colnames(exercise24hrMeans) <- "HExercise24hr"
exercise24hrJoin <- rownames_to_column(exercise24hrMeans, var="geneNames")

exercise48hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise48hr"]] 
exercise48hrMeans <- data.frame(rowMeans(exercise48hrVals))
colnames(exercise48hrMeans) <- "HExercise48hr"
exercise48hrJoin <- rownames_to_column(exercise48hrMeans, var="geneNames")

#Combine data for all time points into a single data frame
dataJoin1 <- inner_join(controlIPEJoin, exerciseIPEJoin)
dataJoin2 <- inner_join(dataJoin1, exercise0.5hrJoin)
dataJoin3 <- inner_join(dataJoin2, exercise1hrJoin)
dataJoin4 <- inner_join(dataJoin3, exercise4hrJoin)
dataJoin5 <- inner_join(dataJoin4, exercise7hrJoin)
dataJoin6 <- inner_join(dataJoin5, exercise24hrJoin)
dataJoin7 <- inner_join(dataJoin6, exercise48hrJoin)
heart_dataAllTimes <- column_to_rownames(dataJoin7, var = "geneNames") # Remove gene ID column
heart_dataAllTimes <- as.matrix(heart_dataAllTimes) 

########

setwd("/Users/Chigoziri/MoTrPAC/liverCtrl")
load("liver_LRTData.RData")

#Identify significant genes which meet cutoff
padj.cutoff <- 0.05

sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

# Subset results for faster cluster finding (for classroom demo purposes)
liver_clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=2000)


#### Obtain mean gene expression data for all pertinent genes at all time points####

# Obtain rlog values for those significant genes
cluster_rld <- rld_mat[liver_clustering_sig_genes$gene, ]

#Z-score scale transformation
rld_scale <- scale(t(cluster_rld))
rld_scale <- t(rld_scale)

controlIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Control_IPE"]] 
controlIPEMeans <- data.frame(rowMeans(controlIPEVals))
colnames(controlIPEMeans) <- "LControl_IPE"
controlIPEJoin <- rownames_to_column(controlIPEMeans, var="geneNames")

exerciseIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise_IPE"]] 
exerciseIPEMeans <- data.frame(rowMeans(exerciseIPEVals))
colnames(exerciseIPEMeans) <- "LExercise_IPE"
exerciseIPEJoin <- rownames_to_column(exerciseIPEMeans, var="geneNames")

exercise0.5hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise0.5hr"]] 
exercise0.5hrMeans <- data.frame(rowMeans(exercise0.5hrVals))
colnames(exercise0.5hrMeans) <- "LExercise0.5hr"
exercise0.5hrJoin <- rownames_to_column(exercise0.5hrMeans, var="geneNames")

exercise1hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise01hr"]] 
exercise1hrMeans <- data.frame(rowMeans(exercise1hrVals))
colnames(exercise1hrMeans) <- "LExercise1hr"
exercise1hrJoin <- rownames_to_column(exercise1hrMeans, var="geneNames")

exercise4hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise04hr"]] 
exercise4hrMeans <- data.frame(rowMeans(exercise4hrVals))
colnames(exercise4hrMeans) <- "LExercise4hr"
exercise4hrJoin <- rownames_to_column(exercise4hrMeans, var="geneNames")

exercise7hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise07hr"]] 
exercise7hrMeans <- data.frame(rowMeans(exercise7hrVals))
colnames(exercise7hrMeans) <- "LExercise7hr"
exercise7hrJoin <- rownames_to_column(exercise7hrMeans, var="geneNames")

exercise24hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise24hr"]] 
exercise24hrMeans <- data.frame(rowMeans(exercise24hrVals))
colnames(exercise24hrMeans) <- "LExercise24hr"
exercise24hrJoin <- rownames_to_column(exercise24hrMeans, var="geneNames")

exercise48hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise48hr"]] 
exercise48hrMeans <- data.frame(rowMeans(exercise48hrVals))
colnames(exercise48hrMeans) <- "LExercise48hr"
exercise48hrJoin <- rownames_to_column(exercise48hrMeans, var="geneNames")

#Combine data for all time points into a single data frame
dataJoin1 <- inner_join(controlIPEJoin, exerciseIPEJoin)
dataJoin2 <- inner_join(dataJoin1, exercise0.5hrJoin)
dataJoin3 <- inner_join(dataJoin2, exercise1hrJoin)
dataJoin4 <- inner_join(dataJoin3, exercise4hrJoin)
dataJoin5 <- inner_join(dataJoin4, exercise7hrJoin)
dataJoin6 <- inner_join(dataJoin5, exercise24hrJoin)
dataJoin7 <- inner_join(dataJoin6, exercise48hrJoin)
liver_dataAllTimes <- column_to_rownames(dataJoin7, var = "geneNames") # Remove gene ID column
liver_dataAllTimes <- as.matrix(liver_dataAllTimes) 

########

setwd("/Users/Chigoziri/MoTrPAC/adiposeCtrl")
load("adipose_LRTData.RData")

#Identify significant genes which meet cutoff
padj.cutoff <- 0.05

sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

# Subset results for faster cluster finding (for classroom demo purposes)
adipose_clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=2000)


#### Obtain mean gene expression data for all pertinent genes at all time points####

# Obtain rlog values for those significant genes
cluster_rld <- rld_mat[adipose_clustering_sig_genes$gene, ]

#Z-score scale transformation
rld_scale <- scale(t(cluster_rld))
rld_scale <- t(rld_scale)

controlIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Control_IPE"]] 
controlIPEMeans <- data.frame(rowMeans(controlIPEVals))
colnames(controlIPEMeans) <- "AControl_IPE"
controlIPEJoin <- rownames_to_column(controlIPEMeans, var="geneNames")

exerciseIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise_IPE"]] 
exerciseIPEMeans <- data.frame(rowMeans(exerciseIPEVals))
colnames(exerciseIPEMeans) <- "AExercise_IPE"
exerciseIPEJoin <- rownames_to_column(exerciseIPEMeans, var="geneNames")

exercise0.5hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise0.5hr"]] 
exercise0.5hrMeans <- data.frame(rowMeans(exercise0.5hrVals))
colnames(exercise0.5hrMeans) <- "AExercise0.5hr"
exercise0.5hrJoin <- rownames_to_column(exercise0.5hrMeans, var="geneNames")

exercise1hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise01hr"]] 
exercise1hrMeans <- data.frame(rowMeans(exercise1hrVals))
colnames(exercise1hrMeans) <- "AExercise1hr"
exercise1hrJoin <- rownames_to_column(exercise1hrMeans, var="geneNames")

exercise4hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise04hr"]] 
exercise4hrMeans <- data.frame(rowMeans(exercise4hrVals))
colnames(exercise4hrMeans) <- "AExercise4hr"
exercise4hrJoin <- rownames_to_column(exercise4hrMeans, var="geneNames")

exercise7hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise07hr"]] 
exercise7hrMeans <- data.frame(rowMeans(exercise7hrVals))
colnames(exercise7hrMeans) <- "AExercise7hr"
exercise7hrJoin <- rownames_to_column(exercise7hrMeans, var="geneNames")

exercise24hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise24hr"]] 
exercise24hrMeans <- data.frame(rowMeans(exercise24hrVals))
colnames(exercise24hrMeans) <- "AExercise24hr"
exercise24hrJoin <- rownames_to_column(exercise24hrMeans, var="geneNames")

exercise48hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise48hr"]] 
exercise48hrMeans <- data.frame(rowMeans(exercise48hrVals))
colnames(exercise48hrMeans) <- "AExercise48hr"
exercise48hrJoin <- rownames_to_column(exercise48hrMeans, var="geneNames")

#Combine data for all time points into a single data frame
dataJoin1 <- inner_join(controlIPEJoin, exerciseIPEJoin)
dataJoin2 <- inner_join(dataJoin1, exercise0.5hrJoin)
dataJoin3 <- inner_join(dataJoin2, exercise1hrJoin)
dataJoin4 <- inner_join(dataJoin3, exercise4hrJoin)
dataJoin5 <- inner_join(dataJoin4, exercise7hrJoin)
dataJoin6 <- inner_join(dataJoin5, exercise24hrJoin)
dataJoin7 <- inner_join(dataJoin6, exercise48hrJoin)
adipose_dataAllTimes <- column_to_rownames(dataJoin7, var = "geneNames") # Remove gene ID column
adipose_dataAllTimes <- as.matrix(adipose_dataAllTimes) 

########

setwd("/Users/Chigoziri/MoTrPAC/bloodCtrl")
load("blood_LRTData.RData")

#Identify significant genes which meet cutoff
padj.cutoff <- 0.05

sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

# Subset results for faster cluster finding (for classroom demo purposes)
blood_clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=2000)


#### Obtain mean gene expression data for all pertinent genes at all time points####

# Obtain rlog values for those significant genes
cluster_rld <- rld_mat[blood_clustering_sig_genes$gene, ]

#Z-score scale transformation
rld_scale <- scale(t(cluster_rld))
rld_scale <- t(rld_scale)

controlIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Control_IPE"]] 
controlIPEMeans <- data.frame(rowMeans(controlIPEVals))
colnames(controlIPEMeans) <- "BControl_IPE"
controlIPEJoin <- rownames_to_column(controlIPEMeans, var="geneNames")

exerciseIPEVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise_IPE"]] 
exerciseIPEMeans <- data.frame(rowMeans(exerciseIPEVals))
colnames(exerciseIPEMeans) <- "BExercise_IPE"
exerciseIPEJoin <- rownames_to_column(exerciseIPEMeans, var="geneNames")

exercise0.5hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise0.5hr"]] 
exercise0.5hrMeans <- data.frame(rowMeans(exercise0.5hrVals))
colnames(exercise0.5hrMeans) <- "BExercise0.5hr"
exercise0.5hrJoin <- rownames_to_column(exercise0.5hrMeans, var="geneNames")

exercise1hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise01hr"]] 
exercise1hrMeans <- data.frame(rowMeans(exercise1hrVals))
colnames(exercise1hrMeans) <- "BExercise1hr"
exercise1hrJoin <- rownames_to_column(exercise1hrMeans, var="geneNames")

exercise4hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise04hr"]] 
exercise4hrMeans <- data.frame(rowMeans(exercise4hrVals))
colnames(exercise4hrMeans) <- "BExercise4hr"
exercise4hrJoin <- rownames_to_column(exercise4hrMeans, var="geneNames")

exercise7hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise07hr"]] 
exercise7hrMeans <- data.frame(rowMeans(exercise7hrVals))
colnames(exercise7hrMeans) <- "BExercise7hr"
exercise7hrJoin <- rownames_to_column(exercise7hrMeans, var="geneNames")

exercise24hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise24hr"]] 
exercise24hrMeans <- data.frame(rowMeans(exercise24hrVals))
colnames(exercise24hrMeans) <- "BExercise24hr"
exercise24hrJoin <- rownames_to_column(exercise24hrMeans, var="geneNames")

exercise48hrVals <- rld_scale[,coldata$sampleIDs[coldata$Condition=="Exercise48hr"]] 
exercise48hrMeans <- data.frame(rowMeans(exercise48hrVals))
colnames(exercise48hrMeans) <- "BExercise48hr"
exercise48hrJoin <- rownames_to_column(exercise48hrMeans, var="geneNames")

#Combine data for all time points into a single data frame
dataJoin1 <- inner_join(controlIPEJoin, exerciseIPEJoin)
dataJoin2 <- inner_join(dataJoin1, exercise0.5hrJoin)
dataJoin3 <- inner_join(dataJoin2, exercise1hrJoin)
dataJoin4 <- inner_join(dataJoin3, exercise4hrJoin)
dataJoin5 <- inner_join(dataJoin4, exercise7hrJoin)
dataJoin6 <- inner_join(dataJoin5, exercise24hrJoin)
dataJoin7 <- inner_join(dataJoin6, exercise48hrJoin)
blood_dataAllTimes <- column_to_rownames(dataJoin7, var = "geneNames") # Remove gene ID column
blood_dataAllTimes <- as.matrix(blood_dataAllTimes) 

########

#Select rownames for all significant genes identified above for next step
gastrorownames <- dplyr::select(gastro_clustering_sig_genes, gene)
heartrownames <- dplyr::select(heart_clustering_sig_genes, gene)
liverrownames <- dplyr::select(liver_clustering_sig_genes, gene)
adiposerownames <- dplyr::select(adipose_clustering_sig_genes, gene)
bloodrownames <- dplyr::select(blood_clustering_sig_genes, gene)

#Find genes which are significant in all tissues
names1 <- intersect(gastrorownames$gene, heartrownames$gene)
names2 <- intersect(names1, liverrownames$gene)
names3 <- intersect(names2, adiposerownames$gene)
names4 <- intersect(names3, bloodrownames$gene)

#Convert common significant Ensembl IDs to gene symbols
intersectSymbols <- as.data.frame(mapIds(org.Rn.eg.db,
                                    keys=names4, 
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first"))

gastro_clustering_sig_genes <- column_to_rownames(gastro_clustering_sig_genes,var="gene")
heart_clustering_sig_genes <- column_to_rownames(heart_clustering_sig_genes,var="gene")
liver_clustering_sig_genes <- column_to_rownames(liver_clustering_sig_genes,var="gene")
adipose_clustering_sig_genes <- column_to_rownames(adipose_clustering_sig_genes,var="gene")
blood_clustering_sig_genes <- column_to_rownames(blood_clustering_sig_genes,var="gene")


#Extract expression data for genes which were common in all tissues
gastro_intersect <- gastro_dataAllTimes[names4,]
rownames(gastro_intersect) <- intersectSymbols[,1]
gastro_intersect <- rownames_to_column(as.data.frame(gastro_intersect),var="Gene")

heart_intersect <- heart_dataAllTimes[names4,]
rownames(heart_intersect) <- intersectSymbols[,1]
heart_intersect <- rownames_to_column(as.data.frame(heart_intersect),var="Gene")

liver_intersect <- liver_dataAllTimes[names4,]
rownames(liver_intersect) <- intersectSymbols[,1]
liver_intersect <- rownames_to_column(as.data.frame(liver_intersect),var="Gene")

adipose_intersect <- adipose_dataAllTimes[names4,]
rownames(adipose_intersect) <- intersectSymbols[,1]
adipose_intersect <- rownames_to_column(as.data.frame(adipose_intersect),var="Gene")

blood_intersect <- blood_dataAllTimes[names4,]
rownames(blood_intersect) <- intersectSymbols[,1]
blood_intersect <- rownames_to_column(as.data.frame(blood_intersect),var="Gene")

#Combine data for all tissues into one data frame in order to create heatmap
join1 <- inner_join(gastro_intersect,heart_intersect,by="Gene")
join2 <- inner_join(join1,liver_intersect,by="Gene")
join3 <- inner_join(join2,adipose_intersect,by="Gene")
join4 <- inner_join(join3,blood_intersect,by="Gene")
join4 <- column_to_rownames(join4,var="Gene")

#Assign appropriate labels to each condition
Condition <- rep(c("Control","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"),5)
Tissue <- c(rep("Gastrocnemius",8),rep("Heart",8),rep("Liver",8),rep("Adipose",8),rep("Blood",8))

anno <- as.data.frame(rbind(Tissue,Condition))
anno <- as.data.frame(t(anno))
anno$Condition <- factor(anno$Condition, levels=c("Control","0 hr","0.5 hr","1 hr","4 hr","7 hr","24 hr","48 hr"))
anno$Tissue <- factor(anno$Tissue, levels=c("Gastrocnemius","Heart","Liver","Adipose","Blood"))
rownames(anno) <- colnames(join4)

setwd("/Users/Chigoziri/MoTrPAC")

colours <- list(
  Tissue = c('Gastrocnemius' = '#a6cee3', 'Heart' = '#b2df8a', 'Liver' = '#fb9a99', 'Adipose' = '#fdbf6f', 'Blood' = '#cab2d6'),
  Condition = c('Control' = '#f7fbff', '0 hr' = '#deebf7', '0.5 hr' = '#c6dbef', '1 hr' = '#9ecae1', '4 hr' = '#6baed6', '7 hr' = '#4292c6', '24 hr' = '#2171b5', '48 hr' = '#084594'))

#Set up color scheme and breaks for heatmap
myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)
myBreaks <- seq(-2.3, 2.3, length.out = 100)

colAnn <- HeatmapAnnotation(df=anno,which="col",col= colours,
                            annotation_name_gp = gpar(fontsize = 12, fontface = 'bold'),
                            annotation_height = 0.6,
                            annotation_width = unit(1, 'cm'),
                            gap = unit(1, 'mm'),
                            show_legend = TRUE,
                            annotation_legend_param = list(legend_direction = 'vertical',
                            legend_width = unit(8, 'cm'),
                            legend_height = unit(5.0, 'cm'),
                            title_position = 'topleft',
                            title_gp=gpar(fontsize = 14, fontface = 'bold'),
                            labels_gp=gpar(fontsize = 14)))


hmap <- ComplexHeatmap::Heatmap(as.matrix(join4),
                        col = colorRamp2(myBreaks, myCol),
                        name = 'Mean\nZ-score',
                        # parameters for the colour-bar that represents gradient of expression
                        heatmap_legend_param = list(
                          color_bar = 'continuous',
                          legend_direction = 'vertical',
                          legend_width = unit(8, 'cm'),
                          legend_height = unit(5.0, 'cm'),
                          title_position = 'topcenter',
                          title_gp=gpar(fontsize = 14, fontface = 'bold'),
                          labels_gp=gpar(fontsize = 14)),
                        
                        # row (gene) parameters
                        cluster_rows = TRUE,
                        show_row_dend = TRUE,
                        row_title_side = 'left',
                        row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                        row_title_rot = 90,
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 11),
                        row_names_side = 'right',
                        row_dend_width = unit(25,'mm'),
                        
                        # column (sample) parameters
                        cluster_columns = FALSE,
                        show_column_dend = FALSE,
                        column_title_side = 'bottom',
                        column_title_gp = gpar(fontsize = 16, fontface = 'bold'),
                        column_title_rot = 0,
                        show_column_names = FALSE,
                        column_names_gp = gpar(fontsize = 14, fontface = 'bold'),
                        column_names_max_height = unit(10, 'cm'),
                        column_dend_height = unit(25,'mm'),
                        column_split = anno$Tissue,
                        
                        # specify top and bottom annotations
                        top_annotation = colAnn)

#Save heatmap file as .tiff
tiff("ComplexHeatmap_AllTissues.tiff", units="in", width=18, height=16, res=300)
draw(hmap,
     merge_legends = TRUE,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right',
     row_sub_title_side = 'right')
dev.off()
