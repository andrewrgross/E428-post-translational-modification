### Analysis of PTM data -- Andrew R Gross -- 23AUG21
###
### A set of basic preliminary actions to begin analyzing a PTM dataset
### INPUT: PTM quantification tables
### OUTPUT: a differential expression table
####################################################################################################################################################
### 1 - Header
library("DESeq2")

### 1.2 - Functions
####### 1.2.1 - filterIncompleteRows - Remove rows which don't have universal rexpression in iECs
filterIncompleteRows <- function(dataframe){
  iecMinimums <- apply(dataframe[metadata$Shortname[metadata$iEC=='iEC']],1,min)
  dataframe <- dataframe[!is.na(iecMinimums),]
  return(dataframe)
}

####################################################################################################################################################
### 2 - Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Input data/')

ptmAll <- read.csv('PTM-all.csv', fileEncoding="UTF-8-BOM")
ptmLace <- read.csv('PTM-lycine-acetylation.csv', fileEncoding="UTF-8-BOM")
ptmMeth <- read.csv('PTM-methylation.csv', fileEncoding="UTF-8-BOM")
ptmNace <- read.csv('PTM-N-term-acetylation.csv', fileEncoding="UTF-8-BOM")
ptmPhos <- read.csv('PTM-phosphorylation.csv', fileEncoding="UTF-8-BOM")
ptmUnmod <- read.csv('PTM-unmodified.csv', fileEncoding="UTF-8-BOM")

metadata <- read.csv('metadata.csv', fileEncoding="UTF-8-BOM")

####################################################################################################################################################
### 3 - Format

### 3.1 - Create new dataframes separating expression and identity data for each modification row
rowAll  <- ptmAll[1:6]      ; ptmAll <-  ptmAll[7:26]               ; dataset = 'All modifications'
rowLace <- ptmLace[1:6]     ; ptmLace <- ptmLace[7:26]              ; dataset = 'Lysene Acetylation'
rowMeth <- ptmMeth[1:6]     ; ptmMeth <- ptmMeth[7:26]              ; dataset = 'Methylation'
rowNace <- ptmNace[1:6]     ; ptmNace <- ptmNace[7:26]              ; dataset = 'N-terminal Aceytlation'
rowPhos <- ptmPhos[1:6]     ; ptmPhos <- ptmPhos[7:26]              ; dataset = 'Phosphorylation'
rowUnmod <- ptmUnmod[c(1:8)] ; ptmUnmod <-ptmUnmod[9:29]            ; dataset = 'Unmodified'
names(ptmAll) <- names(ptmLace) <- names(ptmMeth) <- names(ptmNace) <- names(ptmPhos) <-  metadata$Shortname

dataContainer <- list(ptmAll, ptmLace, ptmMeth, ptmNace, ptmPhos);  names(dataContainer) <- c('All Modifications', 'L-Acetylation', 'Methylation', 'N-Acetylation', 'Phosphorylation')
rowContainer <- list(rowAll, rowLace, rowMeth, rowNace, rowPhos);   names(rowContainer) <- names(dataContainer)


### 3.2 - Filter out any rows in which iECs aren't expressed in all samples

for(pos in 1:length(dataContainer)) {
  #dataContainer[[pos]] <- filterIncompleteRows(dataContainer[[pos]])
  dataContainer[[pos]] <- round(log(dataContainer[[pos]] + 1)*10,0)
  dataContainer[[pos]][is.na(dataContainer[[pos]])] = 0
  #dataContainer[[pos]] <- round(dataContainer[[pos]],0)
  
  names(dataContainer[[pos]]) <- metadata$Shortname
}

### 3.3 - Format for DESeq
### iEC + HUVEC vs iPSC
(columns1   <- data.frame(row.names = metadata$Shortname, condition = metadata$CellType))   # EC vs not iEC

### iEC vs HUVEC & iPSC
(columns2   <- data.frame(row.names = metadata$Shortname, condition = metadata$iEC))   # iEC vs not iEC

### iEC vs just iPSC
(columns3   <- data.frame(row.names = metadata$Shortname, condition = metadata$iEC)[1:17, , drop = FALSE])   # iEC vs iPSC

####################################################################################################################################################
### 4 - Differential Expression
############################################################################################
### Make our DESeq data sets
"ptmTypeNum = 2
(ptmType <- names(dataContainer)[ptmTypeNum])
ptmData <- dataContainer[[ptmTypeNum]]"

ddsContainer <- list()

for(pos in 1:length(dataContainer)) {
  ptmData <- dataContainer[[pos]]
  dds1 <- DESeqDataSetFromMatrix(countData = as.matrix(ptmData[row.names(columns1)]), colData = columns1, design = ~ condition) ; title1 = 'EC vs iPSC'
  dds2 <- DESeqDataSetFromMatrix(countData = as.matrix(ptmData[row.names(columns2)]), colData = columns2, design = ~ condition) ; title2 = 'iEC vs HUVEC+iPSC'
  dds3 <- DESeqDataSetFromMatrix(countData = as.matrix(ptmData[row.names(columns3)]), colData = columns3, design = ~ condition) ; title3 = 'iEC vs iPSC'
  dds1 <- DESeq(dds1)
  dds2 <- DESeq(dds2)
  dds3 <- DESeq(dds3)
  results1 <- as.data.frame(results(dds1))
  results2 <- as.data.frame(results(dds2))
  results3 <- as.data.frame(results(dds3))
  results1 <- cbind(rowContainer[[pos]][c(1,2,4)], results1, ptmData)
  results2 <- cbind(rowContainer[[pos]][c(1,2,4)], results2, ptmData)
  results3 <- cbind(rowContainer[[pos]][c(1,2,4)], results3, ptmData)
  
  newListElements <- list(results1, results2, results3)
  ptmGroupName <- substring(names(dataContainer[pos]),1,7)
  names(newListElements) <- paste0(ptmGroupName,'-',c('iec+H_versus','iec_versus', 'iec_v_ips'))
  ddsContainer <- c(ddsContainer,newListElements)
}


####################################################################################################################################################
### 5 - Isolate key sets: Upreg in EC v iPSC
############################################################################################
### Reorder dds results by adj. pval and filter out rows where PTMs are DOWN regulated

for(pos in 1:length(ddsContainer)){
  dds <- ddsContainer[[pos]]
  dds <- dds[order(dds$padj,decreasing = FALSE),]
  dds <- dds[dds$log2FoldChange<0,]
  ddsContainer[[pos]] <- dds
}

dds1pos = c(1,4,7,10,13)
dds2pos = dds1pos + 1

iecPlusHuvContainer = ddsContainer[c(1,4,7,10,13)]

for(pos in 1:length(ddsContainer)){
  results <- as.data.frame(results(ddsContainer[[pos]]))
  results <- results1 <- cbind(rowData[1:2], rowData[4:6], results1, ptmData, rowData[3])

  
}



dds1 <- DESeqDataSetFromMatrix(countData = as.matrix(ptmData), colData = columns1, design = ~ condition) ; title1 = 'EC vs iPSC'
dds2 <- DESeqDataSetFromMatrix(countData = as.matrix(ptmData[row.names(columns2)]), colData = columns2, design = ~ condition) ; title2 = 'iEC vs HUVEC'
#dds3 <- DESeqDataSetFromMatrix(countData = as.matrix(ptmData[row.names(columns3)]), colData = columns3, design = ~ condition) ; title = 'iEC vs iPSC'

### Run DESeq
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

### Define results table
results1 <- as.data.frame(results(dds1))
results2 <- as.data.frame(results(dds2))

####################################################################################################################################################
### 5 - Analysis
############################################################################################
### 5.1 - Rejoin the identity info for each row

resultsContainer <- ddsContainer

for(pos in 1:length(rowContainer)) {
  resultsContainer[[pos*3-2]] <- 

}

results1 <- cbind(rowData[1:2], rowData[4:6], results1, ptmData, rowData[3])
results2 <- cbind(rowData[1:2], rowData[4:6], results2, ptmData, rowData[3])

### 5.2 - Remove NAs
results1 <- results1[!is.na(results1$padj),]
results2 <- results2[!is.na(results2$padj),]

### 5.3 Examine Distribution - Remove NAs
summary(results1)
quantile(results1$padj, seq(0,1,0.1))
quantile(results2$padj, seq(0,1,0.1))

### 5.4 - Filter and reorder by padj
results1 <- results1[results1$padj < 0.005,]
results1 <- results1[order(results1$padj, decreasing = FALSE),] ; nrow(results1)
results2 <- results2[results2$padj < 0.005,]
results2 <- results2[order(results2$padj, decreasing = FALSE),] ; nrow(results2)




####################################################################################################################################################
### 6 - Export
############################################################################################

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/DE")

write.csv(results1, paste(dataset, ' - ', title1, '.csv'))
write.csv(results2, paste(dataset, ' - ', title2, '.csv'))

