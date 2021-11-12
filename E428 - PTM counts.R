### E428 - PTM Counts -- Andrew R Gross -- 22SEP21
###
### A program to count the number of PTMs expressed in an expression table by group
### INPUT: PTM quantification tables
### OUTPUT: A counts table of the number of PTMs consistently found in each group
####################################################################################################################################################
### 1 - Header
####### 1.1 - Libraries
library(ggplot2)

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

dataContainer <- list(ptmAll, ptmLace, ptmMeth, ptmNace, ptmPhos)
rowContainer <- list(rowAll, rowLace, rowMeth, rowNace, rowPhos)
names(dataContainer) <- c('All Modifications', 'L-Acetylation', 'Methylation', 'N-Acetylation', 'Phosphorylation')
names(rowContainer) <- names(dataContainer)

#dataContainer <- mapply(log, dataContainer)

for(pos in 1:length(dataContainer)) {
  dataContainer[[pos]] <- round(log(dataContainer[[pos]] + 1),2)
  dataContainer[[pos]][is.na(dataContainer[[pos]])] = 0
  names(dataContainer[[pos]]) <- metadata$Shortname
}



test <- countCheck(inputDf, columns= 1:8, cutoffNumber = 7)
####################################################################################################################################################
### 4 - Subset data tables to count based on presence in the specified condition
############################################################################################
### 4.1 - For each row, check if a set of coditions are met (universal detection, for now)

ptmCountsByCondition = data.frame(allCheck = integer(), iecChec = integer(), huvecCheck = integer(), ipscCheck = integer(), iecHuCheck = integer(), iecIpscCheck = integer(), huIpscCheck = integer())

for (ptmTypeNum in 1:length(dataContainer)) {
  ptmType <- names(dataContainer)[ptmTypeNum]
  print(ptmType)
  inputDf <- dataContainer[[ptmTypeNum]]
  
  conditionReport <- data.frame(allCheck = logical(), iecChec = logical(), huvecCheck = logical(), ipscCheck = logical(), iecHuCheck = logical(), iecIpscCheck = logical(), huIpscCheck = logical())
  
  inputDf$iecCheck = countCheck(inputDf,columns = which(metadata$Group=='iEC'), cutoffNumber = 8)
  inputDf$ipsCheck = countCheck(inputDf,columns = which(metadata$Group=='iPSC'), cutoffNumber = 7)
  inputDf$huvCheck = countCheck(inputDf,columns = which(metadata$Group=='Huvecs'), cutoffNumber = 3)
  inputDf$allCheck = as.logical(inputDf$iecCheck * inputDf$ipsCheck * inputDf$huvCheck)
  inputDf$iecHuCheck   = as.logical(inputDf$iecCheck * inputDf$huvCheck)
  inputDf$iecIpscCheck = as.logical(inputDf$iecCheck * inputDf$ipsCheck)
  inputDf$huIpscCheck  = as.logical(inputDf$ipsCheck * inputDf$huvCheck)
  
  ### For each row, check if the PTM is listed in all of each cell type
  for(rowNum in 1:nrow(inputDf)) {
    row = inputDf[rowNum,]                                       # Loop through row numbers
    ### Subset each cell type
    ipsc = row[metadata[metadata$Group=='iPSC',]$Shortname]      
    iec = row[metadata[metadata$Group=='iEC',]$Shortname]
    huvec = row[metadata[metadata$Group=='Huvecs',]$Shortname]
    ### Check if all samples in a group have a non-zero value
    
    allCheck = min(row)       != 0
    iecCheck = countCheck(iec)
    huvecCheck = min(huvec)   != 0
    ipscCheck = min(ipsc)     != 0
    ### Check if the row is universally present in multiple groups
    iecHuCheck = as.logical(iecCheck * huvecCheck)
    iecIpscCheck = as.logical(iecCheck * ipscCheck)
    huIpscCheck = as.logical(huvecCheck * ipscCheck)
    ### Assign the result to the condition report file
    conditionReport = rbind(conditionReport,data.frame(allCheck, iecCheck, huvecCheck, ipscCheck, iecHuCheck, iecIpscCheck, huIpscCheck))
  }
  
  newConditionData <- data.frame(t(apply(conditionReport,2,sum)), row.names = ptmType)
  ### Subtract out co-occurring rows from the groups themselves to make those exclusive
  newConditionData[2:7] = newConditionData[2:7]-newConditionData$allCheck
  newConditionData$iecCheck = newConditionData$iecCheck - newConditionData$iecHuCheck - newConditionData$iecIpscCheck
  newConditionData$huvecCheck = newConditionData$huvecCheck - newConditionData$iecHuCheck - newConditionData$huIpscCheck
  newConditionData$ipscCheck = newConditionData$ipscCheck - newConditionData$huIpscCheck - newConditionData$iecIpscCheck
  
  
  ptmCountsByCondition <- rbind(ptmCountsByCondition, newConditionData)
}

####################################################################################################################################################
### 5 - Generate Bar Graphs
############################################################################################
### 5.1 - Review the summary of each dF

(ptmType = names(dataContainer)[5])

group = c(rep('iEC',4), rep('HUVEC',4), rep('iPSC',4))
shared = c('iEC','iEC & HUVEC','iEC & iPSC', 'All', 'HUVEC', 'iEC & HUVEC','HUVEC & iPSC', 'All', 'iPSC', 'iEC & iPSC', 'HUVEC & iPSC', 'All' )
value = ptmCountsByCondition[ptmType, , drop(FALSE)]
value = c(value['iecCheck'][[1]], value['iecHuCheck'][[1]], value['iecIpscCheck'][[1]], value['allCheck'][[1]], value['huvecCheck'][[1]], value['iecHuCheck'][[1]], value['huIpscCheck'][[1]], value['allCheck'][[1]], value['ipscCheck'][[1]], value['iecIpscCheck'][[1]], value['huIpscCheck'][[1]], value['allCheck'][[1]]  )
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

ggplot(data = ptmCountsToPlot, aes(x = group, y = value, fill = shared)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  geom_text(aes(x = group, y = value, label = text), color = 'White', size = 8, position = position_stack(vjust = 0.5)) +
  geom_bar(data = ptmCountsToPlot2, aes(x = group, y = value), color = 'Black', fill = NA, position = 'stack', stat = 'identity', size = 1.1) + 
  scale_fill_manual(values=c("#37d51e", "#1e37d5", "#d51e37", "#1e93d5", "#d5bc1e", "#d51e93", "#e3e3e3"), labels = c('iEC only', 'HUVEC only', 'iPSC only', 'iEC & HUVEC', 'iEC & iPSC', 'HUVEC & iPSC', 'Found in all')) +
  #scale_y_continuous(breaks = seq(0,50,5), limits = c(0,45.1), expand = c(0,0)) +
  labs(title = ptmType,
       x = 'Cell Type', 
       y = 'PTMs Detected',
       fill = 'Overlap') +
  theme(plot.title = element_text(size = 24,hjust = 0.5),
        axis.title.x = element_text(face="bold", size=14, margin =margin(10,0,0,0)),
        axis.title.y = element_text(face="bold", size=14, margin =margin(0,10,0,0)),
        axis.text.x = element_text(face="italic", size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = 'white', color = 'white', size = 1),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 2))



####################################################################################################################################################
### 6 - Save PTMs to group lists
############################################################################################
### 6.1 - Define group lists
ptmCountsByCondition

allList = c() ; allListM = c()
iecList = c() ; iecListM = c()
huvList = c() ; huvListM = c()
ipsList = c() ; ipsListM = c() 
iechuvL = c() ; iechuvLM = c()
iecipsL = c() ; iecipsLM = c()
huvipsL = c() ; huvipsLM = c()

### 6.2 - Call data based on PTM group

ptmTypeNum = 5
(ptmType <- names(dataContainer)[ptmTypeNum])
inputDf <- dataContainer[[ptmTypeNum]]
proteins = rowContainer[[ptmTypeNum]]$Gene
mods =  rowContainer[[ptmTypeNum]]$Modified.Sequence

conditionReport <- data.frame(allCheck = logical(), iecChec = logical(), huvecCheck = logical(), ipscCheck = logical(), iecHuCheck = logical(), iecIpscCheck = logical(), huIpscCheck = logical())
### 6.3 - For each row, check if the PTM is listed in all of each cell type
for(rowNum in 1:nrow(inputDf)) {
  row = inputDf[rowNum,]                                       # Loop through row numbers
  ### Subset each cell type
  ipsc = row[metadata[metadata$Group=='iPSC',]$Shortname]      
  iec = row[metadata[metadata$Group=='iEC',]$Shortname]
  huvec = row[metadata[metadata$Group=='Huvecs',]$Shortname]
  ### Check if all samples in a group have a non-zero value
  allCheck = min(row)       != 0
  iecCheck = min(iec)       != 0
  huvecCheck = min(huvec)   != 0
  ipscCheck = min(ipsc)     != 0
  ### Check if the row is universally present in multiple groups
  iecHuCheck = as.logical(iecCheck * huvecCheck)
  iecIpscCheck = as.logical(iecCheck * ipscCheck)
  huIpscCheck = as.logical(huvecCheck * ipscCheck)
  conditionReport = rbind(conditionReport,data.frame(allCheck, iecCheck, huvecCheck, ipscCheck, iecHuCheck, iecIpscCheck, huIpscCheck))
  ### Add protein to appropriate list and skip others
  if(allCheck){                                             # If all check is true, add to allList
    allList = c(allList, proteins[rowNum]) ; allListM = c(allListM, mods[rowNum]) 
  } else if(iecHuCheck) {                             # If iecHuCheck is true add to iechuvL
    iechuvL = c(iechuvL, proteins[rowNum]) ; iechuvLM = c(iechuvLM, mods[rowNum])
  } else if(iecIpscCheck) {                          # If iecIpscCheck is true add to iecipsL 
    iecipsL = c(iecipsL, proteins[rowNum]) ; iecipsLM = c(iecipsLM, mods[rowNum])
  } else if(huIpscCheck) {                           # If huIpscCheck is tru add to huvipsL
    huvipsL = c(huvipsL, proteins[rowNum]) ; huvipsLM = c(huvipsLM, mods[rowNum])
  } else if(iecCheck) {                        # If iecCheck is true add to iecList
    #print(paste('iEC:', proteins[rowNum]))
    iecList = c(iecList, proteins[rowNum]) ; iecListM = c(iecListM, mods[rowNum])
  } else if(huvecCheck) {                      # If huvecCheck is true add to huvList
    #print(paste('HUVEC:', proteins[rowNum]))
    huvList = c(huvList, proteins[rowNum]) ; huvListM = c(huvListM, mods[rowNum])
  } else if(ipscCheck) {                       # If ipscCheck is true add ipsList
    #print(paste('iPSC:', proteins[rowNum]))
    ipsList = c(ipsList, proteins[rowNum]) ; ipsListM = c(ipsListM, mods[rowNum])
  }
}
masterList = list('All Cell Types' = allList,'iECs' = iecList, 'HUVECs' = huvList, 'iPSCs' = ipsList, 'Endotheial (iECs+HUVECs)' = iechuvL, 'Induced (iECs+iPSCs)' = iecipsL, 'HUVECs+iPSCs' = huvipsL)
modsList = list('All Cell Types' = allListM,'iECs' = iecListM, 'HUVECs' = huvListM, 'iPSCs' = ipsListM, 'Endotheial (iECs+HUVECs)' = iechuvLM, 'Induced (iECs+iPSCs)' = iecipsLM, 'HUVECs+iPSCs' = huvipsLM)

### 6.4 - Find the longest list
maxListLength = 0
for(list in masterList){
  print(length(list))
  maxListLength = max(maxListLength, length(list))
}

### 6.5 - Create an emtpy matrix to write to
outputMatrix <- matrix('', ncol = 14, nrow = maxListLength)
colnames(outputMatrix) <- c('All Cell Types - Protein', 'All - Mod', 'iECs - Protein', 'iECs - Mod', 'HUVECs - Protein', 'HUVECs - Mod', 'iPSCs - Protein', 'iPSC - Mod', 'iECs+HUVECs - Protein', 'iEC+HUVECs - Mod', 'iECs+iPSCs - Protein', 'iEC+iPSCs - Mod', 'HUVECs+iPSCs - Protein', 'HUVECs+iPSCs - Mod')

### 6.6 - Assign each list to a column
for(listNum in 1:length(masterList)){
  selectedList = masterList[[listNum]]
  selectedListM = modsList[[listNum]]
  #outputMatrix[,listNum][1:length(selectedList)] <- selectedList
  assignmentPos = listNum*2-1
  outputMatrix[0:length(selectedList),assignmentPos] <- selectedList
  outputMatrix[0:length(selectedListM),assignmentPos+1] <- selectedListM
}


setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Counts/')
write.csv(outputMatrix, paste(ptmType,'-PTMs found.csv'), row.names = FALSE)





length(allList)
length(iecList)
length(huvList)
length(ipsList)
length(iechuvL)
length(iecipsL)
length(huvipsL)



test <- data.frame('all' = allList, 'iec' = rep(0,103), 'huvec' = rep(0,103))
test$iec = iecList

####################################################################################################################################################
### 5 - Generate Bar Graphs
############################################################################################
### 5.1 - Review the summary of each dF








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

