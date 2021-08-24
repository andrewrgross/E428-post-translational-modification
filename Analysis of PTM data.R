### Analysis of PTM data -- Andrew R Gross -- 23AUG21
###
### A set of basic preliminary actions to begin analyzing a PTM dataset
### INPUT: PTM quantification tables
### OUTPUT: a differential expression table
####################################################################################################################################################
### Header


####################################################################################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E428 - PTM of iECs/Input data/')

ptmAll <- read.csv('PTM-all.csv')
ptmLace <- read.csv('PTM-lycine-acetylation.csv')
ptmMeth <- read.csv('PTM-methylation.csv')
ptmNace <- read.csv('PTM-N-term-acetylation.csv')
ptmPhos <- read.csv('PTM-phosphorylation.csv')
ptmUnmod <- read.csv('PTM-unmodified.csv')

