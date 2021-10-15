### E428 - Endothelial protein detection -- Andrew R Gross -- 11OCT21
###
### Reads in a list of endothelial proteins and finds their occurrence in tables of EC proteins
### INPUT: PTM quantification tables; table of EC-associated genes
### OUTPUT: A counts table of the number of PTMs consistently found in each group
### 1 - Header                        #########################################################
######## 1.1 - Libraries  #####################################################################
library("biomaRt")
library(reshape2)

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

######## 1.2 - Functions  #####################################################################
############ 1.2.1 - addGene: adds the gene name based on the Ensembl ID as a new column
addID <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','ensembl_gene_id','external_gene_name'), filters='external_gene_name', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
############ 1.2.3 - add.description: adds the gene description based on the Ensembl ID as a new column
add.description <- function(dataframe) {
  descr <- getBM(attributes=c('external_gene_name','description'), filters='external_gene_name', values=dataframe[1], mart=ensembl)
  descr <- descr[match(dataframe[,1],descr[,1]),]
  descriptions <- c()
  for (rowNumber in 1:length(descr[,1])) {
    newDescr <- descr[rowNumber,][,2]
    newDescr <- strsplit(newDescr, " \\[")[[1]][1]
    descriptions <- c(descriptions, newDescr)
  }
  dataframe[length(dataframe)+1] <- descriptions
  names(dataframe)[ncol(dataframe)] <- "Description"
  return(dataframe)
}
############ 1.2.4 - meltLists: Joins each column into one long dataframe
meltLists <- function(dataframe){
  groups <- names(dataframe)
  listRaw = c()
  df <- data.frame('protein' = c(), 'group' = c())
  for(colNum in 1:ncol(dataframe)){
    col = dataframe[,colNum]
    col = col[!is.na(col)]
    col = col[sapply(col,nchar)>2]
    group <- rep(groups[colNum], length(col))
    
    newRows <- data.frame('protein' = col, 'group' = group)
    df <- rbind(df, newRows)
  }
  return(df)
}

############ 1.2.5 - unmelt: Concatentates all occurences of a gene into one entry with a list of associations
unmelt <- function(dataframe){
  outDf <- dataframe[0,]
  geneList <- unique(dataframe[,1])
  assocList <- rep('',length(geneList))
  for(pos in 1:length(geneList)){
    gene = geneList[pos]
    newAssoc = paste(dataframe$Association[which(dataframe$GeneName == gene)], collapse = '; ')
    assocList[pos] = newAssoc
  }
  return(data.frame(geneList, assocList))
}

####################################################################################################################################################
### 2 - Input
####### 2.1 - Import EC master list
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/')

ecList <- read.csv('Endothelial-Associated Genes/iec-gene-master-list.csv', fileEncoding="UTF-8-BOM")
ecList <- read.csv('Endothelial-Associated Genes/iec-gene-master-list2.csv', fileEncoding="UTF-8-BOM")

ecListGenes <- toupper(gsub('-','',ecList$GeneName))
ecList$GeneName <- ecListGenes
ecListGenes <- sort(unique(ecList$GeneName))
length(ecListGenes)

####### 2.2 - Import PTM detected list
ptmsRaw <- read.csv('Input data/PTM-all.csv', fileEncoding = "UTF-8-BOM")
ptmsAll <- read.csv('Counts/All -PTMs found.csv', fileEncoding = "UTF-8-BOM")[c(1,3,5,7,9,11,13)]
ptmsLace <- read.csv('Counts/Lace -PTMs found.csv', fileEncoding = "UTF-8-BOM")[c(1,3,5,7,9,11,13)]
ptmsMeth <- read.csv('Counts/Meth -PTMs found.csv', fileEncoding = "UTF-8-BOM")[c(1,3,5,7,9,11,13)]
ptmsNace <- read.csv('Counts/Nace -PTMs found.csv', fileEncoding = "UTF-8-BOM")[c(1,3,5,7,9,11,13)]
ptmsPhos <- read.csv('Counts/Phos -PTMs found.csv', fileEncoding = "UTF-8-BOM")[c(1,3,5,7,9,11,13)]

####### 2.3 - Format imported PTMs
ptmsLFull <- meltLists(ptmsLace) ; nrow(ptmsLFull)
ptmsMFull <- meltLists(ptmsMeth) ; nrow(ptmsMFull)
ptmsNFull <- meltLists(ptmsNace) ; nrow(ptmsNFull)
ptmsPFull <- meltLists(ptmsPhos) ; nrow(ptmsPFull)
ptmsAFull <- meltLists(ptmsAll)  ; nrow(ptmsAFull)

ptmsRawGenes <- sort(unique(trimws(ptmsRaw$Gene)))
ptmsRawGenes <- toupper(gsub('-','',ptmsRawGenes))
length(ptmsRawGenes)
#rownames(ecList) <- ecList$GeneName
#ecList <- addID(ecList)

####################################################################################################################################################
### 3 - Cross-referencing
####### 3.1 - First identify any shared entries
shared <- intersect(ptmsRawGenes, ecListGenes)

(sharedL <- intersect(ptmsLFull$protein,ecListGenes))
(sharedM <- intersect(ptmsMFull$protein,ecListGenes))
(sharedN <- intersect(ptmsNFull$protein,ecListGenes))
(sharedP <- intersect(ptmsPFull$protein,ecListGenes))
(sharedA <- intersect(ptmsAFull$protein,ecListGenes))

####### 3.2 - Check overlap when pruning names
"
trimLen = 2
ptmsTrimmed <- substr(ptmsRawGenes,1,trimLen)
ecListTrimmed <- substr(ecListGenes,1,trimLen)

shared <- intersect(ptmsTrimmed, ecListTrimmed)
length(shared)

(abrev1 <- sort(ptmsRawGenes[unique(match(ecListTrimmed,ptmsTrimmed))]))
(abrev2 <- sort(ecListGenes[unique(match(ptmsTrimmed, ecListTrimmed)[!is.na(match(ptmsTrimmed, ecListTrimmed))])]))
"

####### 3.3 - Check overlap when pruning names
matchL <- which(ecList$GeneName %in% ptmsLFull$protein)
(ecListL <- ecList[matchL,])

matchM <- which(ecList$GeneName %in% ptmsMFull$protein)
(ecListM <- ecList[matchM,])

matchN <- which(ecList$GeneName %in% ptmsNFull$protein)
(ecListN <- ecList[matchN,])

matchP <- which(ecList$GeneName %in% ptmsPFull$protein)
(ecListP <- ecList[matchP,])

matchA <- which(ecList$GeneName %in% ptmsAFull$protein)
(ecListA <- ecList[matchA,])

####### 3.4 - Join redundant rows into single row
(ecListL <- unmelt(ecListL))
(ecListM <- unmelt(ecListM))
(ecListN <- unmelt(ecListN))
(ecListP <- unmelt(ecListP))
(ecListA <- unmelt(ecListA))

####### 3.5 - Add protein description
(ecListL <- add.description(ecListL))
(ecListM <- add.description(ecListM))
(ecListN <- add.description(ecListN))
(ecListP <- add.description(ecListP))
(ecListA <- add.description(ecListA))


####################################################################################################################################################
### 4 - Add a column identifying which group it's in, then order by group
####### 4.1 - Add a column
ecListL$group <- ptmsLFull[match(ecListL$geneList,ptmsLFull$protein),]$group
ecListL <- ecListL[order(ecListL$group),]

ecListM$group <- ptmsMFull[match(ecListM$geneList,ptmsMFull$protein),]$group
ecListM <- ecListM[order(ecListM$group),]

ecListN$group <- ptmsNFull[match(ecListN$geneList,ptmsNFull$protein),]$group
ecListN <- ecListN[order(ecListN$group),]

ecListP$group <- ptmsPFull[match(ecListP$geneList,ptmsPFull$protein),]$group
ecListP <- ecListP[order(ecListP$group),]

ecListA$group <- ptmsAFull[match(ecListA$geneList,ptmsAFull$protein),]$group
ecListA <- ecListA[order(ecListA$group),]



####################################################################################################################################################
### 5 - Export files

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Counts/')
write.csv(ecListL, 'Sig-proteins-L-acetylation.csv', row.names = FALSE)
write.csv(ecListM, 'Sig-proteins-Methylation.csv', row.names = FALSE)
write.csv(ecListN, 'Sig-proteins-N-Acetylation.csv', row.names = FALSE)
write.csv(ecListP, 'Sig-proteins-Phosporylation.csv', row.names = FALSE)
write.csv(ecListA, 'Sig-proteins-Any-modification.csv', row.names = FALSE)



