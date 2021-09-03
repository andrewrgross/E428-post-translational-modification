### Analysis of PTM data -- Andrew R Gross -- 23AUG21
###
### A set of basic preliminary actions to begin analyzing a PTM dataset
### INPUT: PTM quantification tables
### OUTPUT: a differential expression table
####################################################################################################################################################
### 1 - Header
#library("DESeq2")

apply.t.test <- function(dataframe, metadata, group1, group2) {
  statColumns <- data.frame(p.val = numeric(), mean1 = numeric(), mean2 = numeric(), noZero1 = logical(), noZero2 = logical(), noZeroBoth = factor(, levels = c('Null', 'OneFull', 'BothFull')))
  for(rowNum in 1:nrow(dataframe)) {
    row <- dataframe[rowNum,]
    set1 = row[metadata$Shortname[metadata$Group == group1]]
    set2 = row[metadata$Shortname[metadata$Group == group2]]
    testResult <- t.test(set1, set2)
    zCheck1 = min(set1) > 0
    zCheck2 = min(set2) > 0
    zCheck3 = c('Null', 'OneFull', 'BothFull')[zCheck1 + zCheck2 + 1]
    newStatRow  <- data.frame(p.val = testResult$p.value, mean1 = testResult$estimate[[1]], mean2 = testResult$estimate[[2]], noZero1 = zCheck1, noZero2 = zCheck2, noZeroBoth = factor(zCheck3))
    statColumns <- rbind(statColumns, newStatRow)
  }
  return(cbind(dataframe, statColumns))
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
rowAll  <- ptmAll[1:7]      ; ptmAll <-  ptmAll[8:27]               ; dataset = 'All modifications'

rowLace <- ptmLace[1:6]     ; ptmLace <- ptmLace[7:26]              ; dataset = 'Lysene Acetylation'

rowMeth <- ptmMeth[1:6]     ; ptmMeth <- ptmMeth[7:26]              ; dataset = 'Methylation'

rowNace <- ptmNace[1:6]     ; ptmNace <- ptmNace[7:26]              ; dataset = 'N-terminal Aceytlation'

rowPhos <- ptmPhos[1:6]     ; ptmPhos <- ptmPhos[7:26]              ; dataset = 'Phosphorylation'

rowUnmod <- ptmUnmod[c(1:8)] ; ptmUnmod <-ptmUnmod[9:29]            ; dataset = 'Unmodified'

dataContainer <- list(ptmAll, ptmLace, ptmMeth, ptmNace, ptmPhos)
rowContainer <- list(rowAll, rowLace, rowMeth, rowNace, rowPhos)
names(dataContainer) <- c('All', 'Lace', 'Meth', 'Nace', 'Phos')
names(rowContainer) <- names(dataContainer)

#dataContainer <- mapply(log, dataContainer)

for(pos in 1:length(dataContainer)) {
  dataContainer[[pos]] <- round(log(dataContainer[[pos]] + 1),2)
  dataContainer[[pos]][is.na(dataContainer[[pos]])] = 0
  names(dataContainer[[pos]]) <- metadata$Shortname
}

####################################################################################################################################################
### 4 - Perform T-test
############################################################################################
### 4.1 - For each row, perform a t-test using the apply.t.test function and then sort

for(pos in 1:length(dataContainer)) {
  tempRow <- apply.t.test(dataContainer[[pos]], metadata, 'iEC', 'Huvecs')
  dataContainer[[pos]] <- tempRow[order(tempRow$p.val, decreasing = FALSE),]
  print(names(dataContainer)[pos])
  print(dataContainer[[pos]][1,])
}

#hist(as.matrix(ptmData))
### 4.2 - rejoin row data
outputData <- dataContainer
for(pos in 1:length(dataContainer)) {
  rows <- rowContainer[[pos]]
  df <- dataContainer[[pos]]
  tempDF <- cbind(rows[rownames(df),][c(1:4)], df[21:26], df[1:20])
  outputData[[pos]] <- tempDF
}

####################################################################################################################################################
### 5 - Analysis
############################################################################################
### 5.1 - Review the summary of each dF

for(pos in 1:length(outputData)) {
  print(names(outputData)[pos])
  print(summary(outputData[[pos]][4:9]))
}

####################################################################################################################################################
### 6 - Export
############################################################################################

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/DE")

write.csv(outputData[[1]], paste(names(outputData)[[1]], ' - iEC v Huvec.csv'))
write.csv(outputData[[2]], paste(names(outputData)[[2]], ' - iEC v Huvec.csv'))
write.csv(outputData[[3]], paste(names(outputData)[[3]], ' - iEC v Huvec.csv'))
write.csv(outputData[[4]], paste(names(outputData)[[4]], ' - iEC v Huvec.csv'))
write.csv(outputData[[5]], paste(names(outputData)[[5]], ' - iEC v Huvec.csv'))


####################################################################################################################################################
### 6 - Scratchwork
############################################################################################


apply.t.test <- function(dataframe, metadata, group1, group2) {
  statColumns <- data.frame(p.val = numeric(), mean1 = numeric(), mean2 = numeric(), noZero1 = logical(), noZero2 = logical(), noZeroBoth = factor(, levels = c('Null', 'OneFull', 'BothFull')))
  for(rowNum in 1:nrow(dataframe)) {
    row <- dataframe[rowNum,]
    set1 = row[metadata$Shortname[metadata$Group == group1]]
    set2 = row[metadata$Shortname[metadata$Group == group2]]
    testResult <- t.test(set1, set2)
    zCheck1 = min(set1) > 0
    zCheck2 = min(set2) > 0
    zCheck3 = c('Null', 'OneFull', 'BothFull')[zCheck1 + zCheck2 + 1]
    newStatRow  <- data.frame(p.val = testResult$p.value, mean1 = testResult$estimate[[1]], mean2 = testResult$estimate[[2]], noZero1 = zCheck1, noZero2 = zCheck2, noZeroBoth = factor(zCheck3))
    statColumns <- rbind(statColumns, newStatRow)
  }
  return(cbind(dataframe, statColumns))
}

summary(apply.t.test(test, metadata, 'iEC', 'Huvecs'))

