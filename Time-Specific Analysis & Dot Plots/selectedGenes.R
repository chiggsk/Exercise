rm(list=ls())

library(dplyr)
library(xlsx)
library(forcats)
library(tibble)

####NOTE: This script requires having previously ran scripts for comparison of each time point to control;
####This script can be adapted to any other tissues you'd like to analyze

####NOTE: Some genes have multiple transcripts which code for the same gene (different Ensembl IDs)
####but map to the same ENTREZID; as it is not possible to distinguish data from multiple
####transcripts once converted to ENTREZ, causative ENTREZIDs must be removed

setwd("/Users/Chigoziri/MoTrPAC")

#Reads list of genes you're interested in; gene symbol is in column 1, Entrez ID is in column 2
selectedGeneList <- read.csv(file="Favorite_Gene_List_Blood.csv", header = TRUE)

####0 hour####

setwd("/Users/Chigoziri/MoTrPAC/bloodCtrl")
load("blood_IPECtrl.RData")

foldchanges_IPE <- res$log2FoldChange
names(foldchanges_IPE) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_IPE <- as.data.frame(foldchanges_IPE[names(foldchanges_IPE) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_IPE)[names(foldchanges_IPE) %in% selectedGeneList[,2]]
selectedGeneFC_IPE$rowname <- geneNames

names(foldchanges_IPE) <- res$symbol
colnames(selectedGeneFC_IPE[1]) <- c("0")

####0.5 hour####

setwd("/Users/Chigoziri/MoTrPAC/blood0.5hr")
load("blood_0.5hrCtrl.RData")

foldchanges_0.5hr <- res$log2FoldChange
names(foldchanges_0.5hr) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_0.5hr <- as.data.frame(foldchanges_0.5hr[names(foldchanges_0.5hr) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_0.5hr)[names(foldchanges_0.5hr) %in% selectedGeneList[,2]]
selectedGeneFC_0.5hr$rowname <- geneNames

names(foldchanges_0.5hr) <- res$symbol
colnames(selectedGeneFC_0.5hr[1]) <- c("0.5")

####1 hour####

setwd("/Users/Chigoziri/MoTrPAC/blood1hr")
load("blood_1hrCtrl.RData")

foldchanges_1hr <- res$log2FoldChange
names(foldchanges_1hr) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_1hr <- as.data.frame(foldchanges_1hr[names(foldchanges_1hr) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_1hr)[names(foldchanges_1hr) %in% selectedGeneList[,2]]
selectedGeneFC_1hr$rowname <- geneNames

names(foldchanges_1hr) <- res$symbol
colnames(selectedGeneFC_1hr[1]) <- c("1")

####4 hour####

setwd("/Users/Chigoziri/MoTrPAC/blood4hr")
load("blood_4hrCtrl.RData")

foldchanges_4hr <- res$log2FoldChange
names(foldchanges_4hr) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_4hr <- as.data.frame(foldchanges_4hr[names(foldchanges_4hr) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_4hr)[names(foldchanges_4hr) %in% selectedGeneList[,2]]
selectedGeneFC_4hr$rowname <- geneNames

names(foldchanges_4hr) <- res$symbol
colnames(selectedGeneFC_4hr[1]) <- c("4")

####7 hour####

setwd("/Users/Chigoziri/MoTrPAC/blood7hr")
load("blood_7hrCtrl.RData")

foldchanges_7hr <- res$log2FoldChange
names(foldchanges_7hr) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_7hr <- as.data.frame(foldchanges_7hr[names(foldchanges_7hr) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_7hr)[names(foldchanges_7hr) %in% selectedGeneList[,2]]
selectedGeneFC_7hr$rowname <- geneNames

names(foldchanges_7hr) <- res$symbol
colnames(selectedGeneFC_7hr[1]) <- c("7")

####24 hour####

setwd("/Users/Chigoziri/MoTrPAC/blood24hr")
load("blood_24hrCtrl.RData")

foldchanges_24hr <- res$log2FoldChange
names(foldchanges_24hr) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_24hr <- as.data.frame(foldchanges_24hr[names(foldchanges_24hr) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_24hr)[names(foldchanges_24hr) %in% selectedGeneList[,2]]
selectedGeneFC_24hr$rowname <- geneNames

names(foldchanges_24hr) <- res$symbol
colnames(selectedGeneFC_24hr[1]) <- c("24")

####48 hour####

setwd("/Users/Chigoziri/MoTrPAC/blood48hr")
load("blood_48hrCtrl.RData")

foldchanges_48hr <- res$log2FoldChange
names(foldchanges_48hr) <- res$entrez

#Select data for all genes in your gene list
selectedGeneFC_48hr <- as.data.frame(foldchanges_48hr[names(foldchanges_48hr) %in% selectedGeneList[,2]])
geneNames <- names(foldchanges_48hr)[names(foldchanges_48hr) %in% selectedGeneList[,2]]
selectedGeneFC_48hr$rowname <- geneNames

names(foldchanges_48hr) <- res$symbol
colnames(selectedGeneFC_48hr[1]) <- c("48")

#Merge data for all time points into a single data frame
selectedGeneFC1 <- inner_join(selectedGeneFC_IPE,selectedGeneFC_0.5hr, by="rowname")
selectedGeneFC2 <- inner_join(selectedGeneFC1,selectedGeneFC_1hr, by="rowname")
selectedGeneFC3 <- inner_join(selectedGeneFC2,selectedGeneFC_4hr, by="rowname")
selectedGeneFC4 <- inner_join(selectedGeneFC3,selectedGeneFC_7hr, by="rowname")
selectedGeneFC5 <- inner_join(selectedGeneFC4,selectedGeneFC_24hr, by="rowname")
selectedGeneFC6 <- inner_join(selectedGeneFC5,selectedGeneFC_48hr, by="rowname")

resGeneIDs <- as.data.frame(res$entrez)
resGeneIDs[2] <- res$symbol
colnames(resGeneIDs) <- c("ENTREZID","GeneSymbol")

#Remove duplicated ENTREZIDs
uniqueResEntrez <- !duplicated(resGeneIDs$ENTREZID)

#Because some ENTREZIDs map to multiple transcripts, must select  gene symbols 
#from a selection of only unique ENTREZIDs
selectedGeneFC6[9] <- resGeneIDs$GeneSymbol[uniqueResEntrez][resGeneIDs$ENTREZID[uniqueResEntrez] %in% selectedGeneFC6$rowname]
timeCourseDF <- column_to_rownames(selectedGeneFC6, var="rowname")
colnames(timeCourseDF) <- c("0","0.5", "1","4","7","24","48","GeneSymbol")

setwd("/Users/Chigoziri/MoTrPAC/bloodCtrl")

write.xlsx2(x = timeCourseDF, 
            file = "blood_selectedGeneList_timeCourse.xlsx", 
            sheetName = "Selected Gene logFC",
            col.names = TRUE,
            row.names = TRUE,
            append = TRUE)
