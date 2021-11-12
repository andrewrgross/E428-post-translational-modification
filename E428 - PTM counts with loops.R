### E428 - PTM Counts -- Andrew R Gross -- 22SEP21
###
### A program to count the number of PTMs expressed in an expression table by group
### INPUT: PTM quantification tables
### OUTPUT: A counts table of the number of PTMs consistently found in each group
####################################################################################################################################################
### 1 - Header
####### 1.1 - Libraries
library(ggplot2)
library(gridExtra)
library(cowplot)
#library(grid)
#library(lattice)

####### 1.2 - Functions
########### 1.2.1 - countCheck - Checks whether the items in the specified columns of each row of a dataframe are present or absent above a cutoff
countCheck <- function(dataframe,columns,cutoffNumber){
  rowCheck = c()
  for(rowNum in 1:nrow(dataframe)){
    row = dataframe[rowNum,]
    trues = as.numeric(row[columns]>0)
    if(sum(trues)>=cutoffNumber){
      rowCheck = c(rowCheck,TRUE)
    } else{
      rowCheck = c(rowCheck,FALSE)
    }
  }
  return(rowCheck)
}

########### 1.2.2 - tTestFormatted - Given a DF and column numbers, returns the DF with t-test data columns appended
tTestFormatted <- function(dataframe, columnsX, columnsY, bind = FALSE){
  meansX = c()
  meansY = c()
  FCs    = c()
  pVal    = c()  
  for(rowNum in 1:nrow(dataframe)){
    row = dataframe[rowNum,]
    x = row[columnsX]
    y = row[columnsY]
    raw = t.test(x, y)
    meanX = raw[[5]][[1]]
    meanY = raw[[5]][[2]]
    fc    = meanX/meanY
    pv    = raw[[3]]
    
    meansX = c(meansX, meanX)
    meansY = c(meansY, meanY)
    FCs    = c(FCs, fc)
    pVal   = c(pVal, pv)
    
  }
  dataframe = cbind(dataframe, meansX, meansY, FCs, pVal)
  dataframe = dataframe[order(dataframe$FCs, decreasing = TRUE),]
  if(bind==FALSE){dataframe = dataframe[c('meansX', 'meansY', 'FCs', 'pVal')]}
  return(dataframe)
}

############ 1.2.3 - unmelt: Concatentates all occurences of a gene into one entry with a list of associations
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
####### 2.1 - Import expression data
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Input data/')

ptmAll <- read.csv('PTM-all.csv', fileEncoding="UTF-8-BOM")
ptmLace <- read.csv('PTM-lycine-acetylation.csv', fileEncoding="UTF-8-BOM")
ptmMeth <- read.csv('PTM-methylation.csv', fileEncoding="UTF-8-BOM")
ptmNace <- read.csv('PTM-N-term-acetylation.csv', fileEncoding="UTF-8-BOM")
ptmPhos <- read.csv('PTM-phosphorylation.csv', fileEncoding="UTF-8-BOM")
ptmUnmod <- read.csv('PTM-unmodified.csv', fileEncoding="UTF-8-BOM")

ptmLace <- read.csv('PTM-Lace-norm.csv', fileEncoding="UTF-8-BOM")
ptmMeth <- read.csv('PTM-Meth-norm.csv', fileEncoding="UTF-8-BOM")
ptmNace <- read.csv('PTM-Nace-norm.csv', fileEncoding="UTF-8-BOM")
ptmPhos <- read.csv('PTM-Phos-norm.csv', fileEncoding="UTF-8-BOM")

####### 2.2 - Import metadata
metadata <- read.csv('metadata.csv', fileEncoding="UTF-8-BOM")

####### 2.3 - Import EC gene lists
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/')
ecList <- read.csv('Endothelial-Associated Genes/iec-gene-master-list2.csv', fileEncoding="UTF-8-BOM")

####################################################################################################################################################
### 3 - Format
brk = 4
####### 3.1 - Create new dataframes separating expression and identity data for each modification row
#rowAll  <- ptmAll[1:brk]      ; ptmAll <-  ptmAll[brk+1:26]               ; dataset = 'All modifications'
rowLace <- ptmLace[1:brk]     ; ptmLace <- ptmLace[(brk+1):(brk+20)]              ; dataset = 'Lysene Acetylation'
rowMeth <- ptmMeth[1:brk]     ; ptmMeth <- ptmMeth[(brk+1):(brk+20)]              ; dataset = 'Methylation'
rowNace <- ptmNace[1:brk]     ; ptmNace <- ptmNace[(brk+1):(brk+20)]              ; dataset = 'N-terminal Aceytlation'
rowPhos <- ptmPhos[1:brk]     ; ptmPhos <- ptmPhos[(brk+1):(brk+20)]              ; dataset = 'Phosphorylation'
#rowUnmod <- ptmUnmod[c(1:8)] ; ptmUnmod <-ptmUnmod[9:29]            ; dataset = 'Unmodified'

dataContainer <- list(ptmLace, ptmMeth, ptmNace, ptmPhos)
rowContainer <- list(rowLace, rowMeth, rowNace, rowPhos)
names(dataContainer) <- c('L-Acetylation', 'Methylation', 'N-Acetylation', 'Phosphorylation')
names(rowContainer) <- names(dataContainer)

#dataContainer <- mapply(log, dataContainer)

for(pos in 1:length(dataContainer)) {
  #dataContainer[[pos]] <- round(log(dataContainer[[pos]] + 1),2)
  dataContainer[[pos]][is.na(dataContainer[[pos]])] = 0
  names(dataContainer[[pos]]) <- metadata$Shortname
}

####### 3.2 - Format the Feng EC genes for searching
ecList <- unmelt(ecList)
ecList$geneList <- toupper(ecList$geneList)
ecList <- ecList[order(ecList$geneList),]

####################################################################################################################################################
### 4 - Subset data tables to count based on presence in the specified condition
############################################################################################
####### 4.1 - For each row, check if a set of coditions are met (universal detection, for now)

ptmCountsByCondition = data.frame(allCheck = integer(), iecChec = integer(), huvCheck = integer(), ipsCheck = integer(), iecHuCheck = integer(), iecIpsCheck = integer(), huIpsCheck = integer())

for (ptmTypeNum in 1:length(dataContainer)) {
  ptmType <- names(dataContainer)[ptmTypeNum]
  print(ptmType)
  inputDf <- dataContainer[[ptmTypeNum]]
  
  conditionReport <- data.frame(iecCheck = logical(), ipsCheck = logical(), huvCheck = logical(), allCheck = logical(), iecHuCheck = logical(), iecIpsCheck = logical(), huIpsCheck = logical())
  
  inputDf$iecCheck = countCheck(inputDf,columns = which(metadata$Group=='iEC'), cutoffNumber = 8)
  inputDf$ipsCheck = countCheck(inputDf,columns = which(metadata$Group=='iPSC'), cutoffNumber = 7)
  inputDf$huvCheck = countCheck(inputDf,columns = which(metadata$Group=='Huvecs'), cutoffNumber = 3)
  inputDf$allCheck = as.logical(inputDf$iecCheck * inputDf$ipsCheck * inputDf$huvCheck)
  inputDf$iecHuCheck   = as.logical(inputDf$iecCheck * inputDf$huvCheck)
  inputDf$iecIpsCheck = as.logical(inputDf$iecCheck * inputDf$ipsCheck)
  inputDf$huIpsCheck  = as.logical(inputDf$ipsCheck * inputDf$huvCheck)
  
  newConditionData <- data.frame(t(apply(inputDf[21:27],2,sum)), row.names = ptmType)
  
  ### Subtract out co-occurring rows from the groups themselves to make those exclusive
  newConditionData[c(1,2,3,5,6,7)] = newConditionData[c(1,2,3,5,6,7)]-newConditionData$allCheck
  newConditionData$iecCheck = newConditionData$iecCheck - newConditionData$iecHuCheck - newConditionData$iecIpsCheck
  newConditionData$huvCheck = newConditionData$huvCheck - newConditionData$iecHuCheck - newConditionData$huIpsCheck
  newConditionData$ipsCheck = newConditionData$ipsCheck - newConditionData$huIpsCheck - newConditionData$iecIpsCheck
  
  ptmCountsByCondition <- rbind(ptmCountsByCondition, newConditionData)
}

ptmCountsByCondition

####################################################################################################################################################
### 5 - Generate Bar Graphs
############################################################################################
####### 5.1 - Review the summary of each dF

setwd('C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Counts/')
plotList <- list()
for (ptmTypeNum in 1:length(dataContainer)) {
  (ptmType = names(dataContainer)[ptmTypeNum])
  group = c(rep('iEC',4), rep('HUVEC',4), rep('iPSC',4))
  shared = c('iEC','iEC & HUVEC','iEC & iPSC', 'All', 'HUVEC', 'iEC & HUVEC','HUVEC & iPSC', 'All', 'iPSC', 'iEC & iPSC', 'HUVEC & iPSC', 'All' )
  value = ptmCountsByCondition[ptmType, , drop(FALSE)]
  value = c(value['iecCheck'][[1]], value['iecHuCheck'][[1]], value['iecIpsCheck'][[1]], value['allCheck'][[1]], value['huvCheck'][[1]], value['iecHuCheck'][[1]], value['huIpsCheck'][[1]], value['allCheck'][[1]], value['ipsCheck'][[1]], value['iecIpsCheck'][[1]], value['huIpsCheck'][[1]], value['allCheck'][[1]]  )
  text = rep(NA, 12)
  cutoff = round(max(value)/10)
  text[value>cutoff] <- value[value>cutoff]
  ptmCountsToPlot = data.frame(group, shared, value, text)
  
  valueSummed = c(ptmCountsToPlot$value[1] + ptmCountsToPlot$value[2]+ptmCountsToPlot$value[3]+ptmCountsToPlot$value[4], 
                  ptmCountsToPlot$value[5] + ptmCountsToPlot$value[6]+ptmCountsToPlot$value[7]+ptmCountsToPlot$value[8], 
                  ptmCountsToPlot$value[9] + ptmCountsToPlot$value[10]+ptmCountsToPlot$value[11]+ptmCountsToPlot$value[12])
  ptmCountsToPlot2 = data.frame(group = c(1,2,3), value = valueSummed)
  ptmCountsToPlot$shared <- factor(ptmCountsToPlot$shared, levels = c('iEC', 'HUVEC', 'iPSC', 'iEC & HUVEC', 'iEC & iPSC', 'HUVEC & iPSC', 'All'))
  ptmCountsToPlot$group <- factor(ptmCountsToPlot$group, levels = c('iEC', 'HUVEC', 'iPSC'))
  ptmCountsToPlot2$group <- factor(c('iEC', 'HUVEC', 'iPSC'), levels = c('iEC', 'HUVEC', 'iPSC'))
  
  plot <- ggplot(data = ptmCountsToPlot, aes(x = group, y = value, fill = shared)) +
    geom_bar(position = 'stack', stat = 'identity') + 
    geom_text(aes(x = group, y = value, label = text), color = 'White', size = 5, position = position_stack(vjust = 0.5)) +
    geom_bar(data = ptmCountsToPlot2, aes(x = group, y = value), color = 'Black', fill = NA, position = 'stack', stat = 'identity', size = 1.1) + 
    scale_fill_manual(values=c("#37d51e", "#1e37d5", "#d51e37", "#1e93d5", "#d5bc1e", "#d51e93", "#e3e3e3"), labels = c('iEC only', 'HUVEC only', 'iPSC only', 'iEC & HUVEC', 'iEC & iPSC', 'HUVEC & iPSC', 'Found in all')) +
    #scale_y_continuous(breaks = seq(0,50,5), limits = c(0,45.1), expand = c(0,0)) +
    labs(title = ptmType,
         x = 'Cell Type', 
         y = 'PTMs Detected',
         fill = 'Overlap') +
    theme(plot.title = element_text(size = 14,hjust = 0.5),
          axis.title.x = element_blank(),  #element_text(face="bold", size=12, margin =margin(10,0,0,0)),
          axis.title.y = element_text(face="bold", size=12, margin =margin(0,10,0,0)),
          axis.text.x = element_text(face="italic", size = 12),
          axis.text.y = element_text(size = 12),
          panel.background = element_rect(fill = 'white', color = 'white', size = 1),
          panel.grid = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 2),
          legend.text = element_text(size = 14),
          legend.position = 'none')
  
  plotList[[ptmTypeNum]] = plot
  #png(paste0('counts_bargraph_',ptmType,'.png'), width = 800, height = 600)
  #print(plot)
  #dev.off()
}
#plotList = plotList[2:5]

plotList[[5]] = get_legend(ggplot(data = ptmCountsToPlot, aes(x = group, y = value, fill = shared)) +
                                              geom_bar(position = 'stack', stat = 'identity') + 
                                              scale_fill_manual(values=c("#37d51e", "#1e37d5", "#d51e37", "#1e93d5", "#d5bc1e", "#d51e93", "#e3e3e3"), labels = c('iEC only', 'HUVEC only', 'iPSC only', 'iEC & HUVEC', 'iEC & iPSC', 'HUVEC & iPSC', 'Found in all')) +
                                              labs(fill = 'Overlap') +
                                              theme(legend.title = element_blank(), legend.text = element_text(size = 13),
                                                    legend.position = 'bottom'))

grid.arrange(
  grobs = plotList, ncol = 2, nrow = 3,
  widths = c(1, 1),
  heights = c(1, 1, 0.3),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4),
                        c(5, 5))
)


####################################################################################################################################################
### 6 - Save PTMs to group lists
############################################################################################
####### 6.1 - Define group lists
ptmCountsByCondition
### Define empty lists for proteins and modifications
metaList = list()

####### 6.2 - Parse data based on PTM group to generate lists of which ptms go in what group
for (ptmTypeNum in 1:length(dataContainer)) {
  ### Declare empty lists
  allList = c()
  iecList = c()
  huvList = c()
  ipsList = c()
  iecHuvL = c()
  iecIpsL = c()
  huvipsL = c()
  ### Assign names and define the expression and row data tables for parsing
  ptmType <- names(dataContainer)[ptmTypeNum]
  print(ptmType)
  inputDf <- dataContainer[[ptmTypeNum]]
  rowdata <- rowContainer[[ptmTypeNum]][2:4]
  # Create columns listing what is exclusively present, then use the row numbers to subset DE
  #inputDf$iecCheck = countCheck(inputDf,columns = which(metadata$Group=='iEC'), cutoffNumber = 8)
  #inputDf$ipsCheck = countCheck(inputDf,columns = which(metadata$Group=='iPSC'), cutoffNumber = 7)
  #inputDf$huvCheck = countCheck(inputDf,columns = which(metadata$Group=='Huvecs'), cutoffNumber = 3)
  
  #inputDf$allCheck = as.logical(inputDf$iecCheck * inputDf$ipsCheck * inputDf$huvCheck)

  #inputDf$justIecCheck = (inputDf$iecCheck - inputDf$ipsCheck - inputDf$huvCheck) == 1
  #inputDf$justIpscCheck = (inputDf$ipsCheck - inputDf$iecCheck - inputDf$huvCheck) == 1
  #inputDf$justHuvCheck = (inputDf$huvCheck - inputDf$ipsCheck - inputDf$iecCheck) == 1

  #inputDf$iecHuCheck   = as.logical((inputDf$iecCheck * inputDf$huvCheck) - inputDf$ipsCheck)
  #inputDf$iecIpsCheck = as.logical((inputDf$iecCheck * inputDf$ipsCheck)  - inputDf$huvCheck)
  #inputDf$huIpsCheck  = as.logical((inputDf$ipsCheck * inputDf$huvCheck)  - inputDf$iecCheck)
  
  ### Generate lists of whether a row meets an occurrence criterion 
  iecCheck = countCheck(inputDf,columns = which(metadata$Group=='iEC'), cutoffNumber = 8)
  ipsCheck = countCheck(inputDf,columns = which(metadata$Group=='iPSC'), cutoffNumber = 7)
  huvCheck = countCheck(inputDf,columns = which(metadata$Group=='Huvecs'), cutoffNumber = 3)
  
  ### Generate lists of exclusive occurence
  iecOnly = (iecCheck - ipsCheck - huvCheck) == 1
  ipsOnly = (ipsCheck - huvCheck - iecCheck) == 1
  huvOnly = (huvCheck - iecCheck - ipsCheck) == 1

  ### Generate lists of shared occurence
  allCheck = as.logical(iecCheck * ipsCheck * huvCheck)
  iecHuCheck   = ((iecCheck * huvCheck) - ipsCheck) == 1
  iecIpsCheck  = ((iecCheck * ipsCheck) - huvCheck) == 1
  ipsHuCheck   = ((ipsCheck * huvCheck) - iecCheck) == 1

  ### Generate row data lists
  iecList = rowdata[iecOnly,]
  ipsList = rowdata[ipsOnly,]
  huvList = rowdata[huvOnly,]
  allList = rowdata[allCheck,]
  iecHuvL = rowdata[iecHuCheck,]
  iecIpsL = rowdata[iecIpsCheck,]
  ipsHuvL = rowdata[ipsHuCheck,]
  
  # Add in a column specifying group
  iecList = cbind(iecList, Group_association = 'iEC')
  ipsList = cbind(ipsList, Group_association = 'iPSC')
  huvList = cbind(huvList, Group_association = 'HUVEC')
  allList = cbind(allList, Group_association = 'All')
  iecHuvL = cbind(iecHuvL, Group_association = 'Endothelial')
  iecIpsL = cbind(iecIpsL, Group_association = 'iEC & iPSC')
  ipsHuvL = cbind(ipsHuvL, Group_association = 'iPSC & HUVEC')
    
  # Add in Feng Associations
  fengAssc = ecList$assoc[match(iecList$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  iecList = cbind(iecList, fengAssc)
  fengAssc = ecList$assoc[match(ipsList$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  ipsList = cbind(ipsList, fengAssc)
  fengAssc = ecList$assoc[match(huvList$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  huvList = cbind(huvList, fengAssc)
  fengAssc = ecList$assoc[match(allList$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  allList = cbind(allList, fengAssc)   
  fengAssc = ecList$assoc[match(iecHuvL$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  iecHuvL = cbind(iecHuvL, fengAssc)
  fengAssc = ecList$assoc[match(iecIpsL$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  iecIpsL = cbind(iecIpsL, fengAssc)
  fengAssc = ecList$assoc[match(ipsHuvL$Gene, ecList$geneList)]; fengAssc[is.na(fengAssc)] = ' '
  ipsHuvL = cbind(ipsHuvL, fengAssc)
  
  ### Add FC and p-value to each group using tTestFormatted function
  iecTtest = tTestFormatted(inputDf[iecOnly,], columnsX = which(metadata$Group =='iEC'), columnsY = which(metadata$Group != 'iEC'), bind = TRUE)
  iecList = cbind(iecList, iecTtest)
  
  tTestResults = tTestFormatted(inputDf[huvOnly,], columnsX = which(metadata$Group =='Huvecs'), columnsY = which(metadata$Group != 'Huvecs'), bind = TRUE)
  huvList = cbind(huvList, tTestResults)
  
  tTestResults = tTestFormatted(inputDf[iecHuCheck,], columnsX = which(metadata$CellType =='EC'), columnsY = which(metadata$CellType != 'EC'), bind = TRUE)
  iecHuvL = cbind(iecHuvL, tTestResults)
  
  ### Add empty columns to groups without a T-test so they can be joined with rbind
  tTestResults = data.frame('meansX' = rep(' ',nrow(inputDf)), 'meansY' = rep(' ',nrow(inputDf)), 'FCs' = rep(' ',nrow(inputDf)), 'pVal' = rep(' ',nrow(inputDf)))
  allList = cbind(allList, inputDf[rownames(allList),], tTestResults[1:nrow(allList),])
  ipsList = cbind(ipsList, inputDf[rownames(ipsList),], tTestResults[1:nrow(ipsList),])
  iecIpsL = cbind(iecIpsL, inputDf[rownames(iecIpsL),], tTestResults[1:nrow(iecIpsL),])
  ipsHuvL = cbind(ipsHuvL, inputDf[rownames(ipsHuvL),], tTestResults[1:nrow(ipsHuvL),])
  
  # Assign row data to metaList
  masterList = list('iECs' = iecList, 'HUVECs' = huvList, 'Endotheial (iECs+HUVECs)' = iecHuvL, 'iPSCs' = ipsList, 'iEC+iPSC' = iecIpsL, 'HUVEC+iPSC' = ipsHuvL, 'All Cell Types' = allList)
  metaList[[ptmTypeNum]] = masterList
}

names(metaList) <- names(dataContainer)


####### 6.4A - Format for export. Loop through metaList and generate a single spreadsheet for groups of samples 
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Counts/')

for(listNum in 1:length(metaList)){                ### Loop through the metalist
  ptmType <- names(metaList)[listNum]              ### For each entry define the current PTM type
  masterList <- metaList[[listNum]]                ### Extract the contents back into a master list
  
  outputDf <- masterList[[1]][0,]                  ### Define an empty DF to fill

  ###### 6.4.3 - Extend each DF to the same lenght and bind columns
  for(listNum in 1:length(masterList)){
    df = masterList[[listNum]]
    outputDf <- rbind(outputDf, df)
  }
  write.csv(outputDf, paste(ptmType,'-PTMs found - vert.csv'), row.names = FALSE)
}



####### 6.4 - Format for export. Loop through metaList and generate a single spreadsheet for groups of samples 

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Counts/')
for(listNum in 1:length(metaList)){
  ptmType <- names(metaList)[listNum]
  masterList <- metaList[[listNum]]
  
  ###### 6.4.1 - Find the longest list of PTMs
  maxListLength = 0
  for(list in masterList){
    maxListLength = max(maxListLength, nrow(list))
  }
  
  ###### 6.4.2 - Create an emtpy matrix to write to
  empty = rep('',maxListLength + 1)
  outputDf <- data.frame('All Cell Types_Protein' =empty, 'All Cell Types_Full name'=empty, 'All Cell Types_Modification'=empty, 'All Cell Types_Association'=empty, 'All Cell Types_FC'=empty, 'All Cell Types_pVal'=empty, 'blank'= empty,
                         'iECs_Protein'=empty, 'iECs_Full name'=empty, 'iECs_Modification'=empty, 'iECs_Association'=empty, 'iECs_FC'=empty, 'iECs_pVal'=empty, 'blank'= empty,
                         'HUVECs_Protein'=empty, 'HUVECs_Full name'=empty, 'HUVECs_Modification'=empty, 'HUVECs_Association'=empty, 'HUVECs_FC'=empty, 'HUVECs_pVal'=empty, 'blank'= empty,
                         'iECs+HUVECs_Protein'=empty, 'iECs+HUVECs_Full name'=empty, 'iECs+HUVECs_Modification'=empty, 'iECs+HUVECs_Association'=empty, 'iECs+HUVECs_FC'=empty, 'iECs+HUVECs_pVal'=empty)
  
  ###### 6.4.3 - Extend each DF to the same lenght and bind columns
  for(listNum in 1:length(masterList)){
    df = masterList[[listNum]]
    rowsToAdd <- maxListLength - nrow(df) + 1
    for(addition in 1:rowsToAdd){
      df <- rbind(df,c(' ',' ',' '))
    }
    outputDf[((listNum*7)-6):(listNum*7)] <- df
  }
  write.csv(outputDf, paste(ptmType,'-PTMs found.csv'), row.names = FALSE)
}

####### 6.7 - Export






length(allList)
length(iecList)
length(huvList)
length(ipsList)
length(iecHuvL)
length(iecIpsL)
length(huvipsL)



test <- data.frame('all' = allList, 'iec' = rep(0,103), 'huvec' = rep(0,103))
test$iec = iecList

####################################################################################################################################################
### 5 - Generate Bar Graphs
############################################################################################
####### 5.1 - Review the summary of each dF








fig3a$expandedPercent <- fig3a$Expanded.abnormal/fig3a$expanded.total*100
fig3a$unexpandedPercent <- fig3a$unexpanded.abnormal/fig3a$unexpanded.total*100
fig3a$zScore <- zScoreCalculator(fig3a$Expanded.abnormal, fig3a$unexpanded.abnormal, fig3a$expanded.total, fig3a$unexpanded.abnormal)
zScores <- c()
for(rowNum in 1:nrow(fig3a)) {
  currentRow <- fig3a[rowNum,]
  newZscore = zScoreCalculator(currentRow[,1], currentRow[,3], currentRow[,2], currentRow[,4])
  zScores <- c(zScores, newZscore)
}
fig3a$zScore <- zScores
fig3a$pValue <- 2*(1-pnorm(fig3a$zScore))
plotA <- melt(t(fig3a[5:6]))


(f3aPlot <- ggplot(data = plotA, aes(x = Var2, y = value)) + 
   geom_bar(aes(fill = Var1), stat = 'identity', position=position_dodge()) + 
   scale_fill_manual(values = c('black', 'grey50')) +
   scale_y_continuous(breaks = seq(0,50,5), limits = c(0,45.1), expand = c(0,0)) +
   labs(x = 'Donor Age', 
        y = '% Abnormal') +
   theme(axis.title.x = element_text(face="italic", size=14, margin =margin(10,0,0,0)),
         axis.title.y = element_text(face="italic", size=14, margin =margin(0,10,0,0)),
         axis.text.x = element_text(size = 14),
         axis.text.y = element_text(size = 14),
         panel.background = element_rect(fill = 'white', color = 'white', size = 1),
         panel.grid = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 2),
         legend.position = 'none') ) #+ coord_equal(ratio = 0.08)

sigData <- data.frame(x=c(0.875, 1.875, 2.875, 3.875, 4.875), xend=c(1.125, 2.125, 3.125, 4.125, 5.125),
                      y=c(16, 18, 19, 18, 42), annotation=c('***', '**', '****', ' **** ', ' ** '))

f3aPlot <- f3aPlot + geom_signif(stat="identity", 
                                 data = sigData,
                                 aes(x=x,xend=xend, y=y, yend=y, annotation=annotation),
                                 tip_length = 0,
                                 vjust = 0) 
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

