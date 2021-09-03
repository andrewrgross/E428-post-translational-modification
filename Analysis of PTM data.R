### Analysis of PTM data -- Andrew R Gross -- 23AUG21
###
### A set of basic preliminary actions to begin analyzing a PTM dataset
### INPUT: PTM quantification tables
### OUTPUT: a differential expression table
####################################################################################################################################################
### 1 - Header
library("DESeq2")


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
ptmData <- ptmAll                     ; dataset = 'All modifications'
rowData <- ptmData[1:7]
ptmData <- ptmData[8:27]

ptmData <- ptmLace                    ; dataset = 'Lysene Acetylation'
ptmData <- ptmMeth                    ; dataset = 'Methylation'
ptmData <- ptmNace                    ; dataset = 'N-terminal Aceytlation'
ptmData <- ptmPhos                    ; dataset = 'Phosphorylation'

rowData <- ptmData[1:6]
ptmData <- ptmData[7:26]

ptmData <- ptmUnmod                   ; dataset = 'Unmodified'

rowData <- ptmData[c(1:3,6:10)]
ptmData <- ptmData[11:30]
ptmData[is.na(ptmData)] = 0
ptmData <- round(ptmData,0)

names(ptmData) <- metadata$Shortname

### 3.2 - Format for DESeq

### iEC + HUVEC vs iPSC
(columns1   <- data.frame(row.names = metadata$Shortname, condition = metadata$CellType))   # EC vs iPSC
(columns2   <- data.frame(row.names = metadata$Shortname, condition = metadata$Group)[9:20, , drop = FALSE])   # iEC vs HUVEC
#(columns3   <- data.frame(row.names = metadata$ShortName, condition = metadata$Group)[1:14, , drop = FALSE])
#(columnData.test     <- data.frame(row.names = metadata.data$Shortname, condition = c('rand0', 'rand1', 'rand0', 'rand1', 'rand1', 'rand0', 'rand1', 'rand0', 'rand1', 'rand1', 'rand0', 'rand1'))[1:6, , drop = FALSE])


####################################################################################################################################################
### 4 - Differential Expression
############################################################################################
### Make our DESeq data sets

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

