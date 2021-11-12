# -*- coding: utf-8 -*-
"""
E428-PTM-Normalization.py

This program will divide PTM expression values based on the detected value of total protein.
INPUT: A PTM expression table; a total protein expression table;
OUTPUT: A normalized PTM expresison table

Created on Fri Nov  5 11:13:43 2021

@author: GrossAR
"""

##############################################################################
### 1. Header
####### 1.1 Libraries

import os
import numpy as np
import pandas as pd

####### 1.2 Functions
def countNaNs(df) :
    nNonZero = []
    for rowNum in range(0,len(df)):
        row = df.iloc[rowNum]
        nNonZero.append(np.count_nonzero(row))
    return(nNonZero)


####### 1.3 Set display parameters
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
np.set_printoptions(edgeitems=3, infstr='inf',linewidth=200, nanstr='nan', precision=8,suppress=False, threshold=1000, formatter=None)


##############################################################################
### 2. Import data

dfPTM = pd.read_excel('C://Users//grossar//Box//Sareen Lab Shared//Data//Proteomics Datasets//2021-07-15_Sarah_Parker_PTMs-of-iECs//Sareen_FraggerDIANNlibrary_PTMS_MS1CorrFiltered0.8.xlsx')
dfTotalProt = pd.read_excel('C://Users/grossar/Box/Sareen Lab Shared//Data//Proteomics Datasets//DIANN_iPSC_EC_TotalProteome_ProteotypicIdentifications.xlsx')
metadata = pd.read_csv('C://Users//grossar//Box//Sareen Lab Shared//Data//Andrew//E428 - PTM of iECs//Input data//metadata.csv')

##############################################################################
### 3. Format

####### 3.1 - Format columns
########## 3.1.1 - Assign IDs as row indecies
dfTotalProt.index = dfTotalProt['Protein.Group']
dfPTM.index = dfPTM['Protein.Group']

########## 3.1.2 - Remove extraneous columns and reorder
dfTotalProt = dfTotalProt[metadata['Name2']]
dfPTM       = dfPTM[metadata['Name2']]

####### 3.2 - Check whether a row has enough values for each group
dfPTM = dfPTM.fillna(0)

colIEC = metadata['Name2'][metadata['Group'] == 'iEC']
colIPSC = metadata['Name2'][metadata['Group'] == 'iPSC']
colHUVEC = metadata['Name2'][metadata['Group'] == 'Huvecs']

nonzeroesIEC = countNaNs(dfPTM[colIEC])
nonzeroesIPSC = countNaNs(dfPTM[colIPSC])
nonzeroesHUVEC = countNaNs(dfPTM[colHUVEC])

nonzeroesIEC = np.array(nonzeroesIEC)>=8
nonzeroesIPSC = np.array(nonzeroesIPSC)>=7
nonzeroesHUVEC = np.array(nonzeroesHUVEC)>=3
rowsToKeep = nonzeroesIEC + nonzeroesIPSC + nonzeroesHUVEC

####### 3.3 - Filter by nonzeroes
########## 3.3.1 - Subset PTMs
dfPTM = dfPTM.iloc[rowsToKeep]

########## 3.3.2 - Identify shared proteins
shared = set.intersection(set(dfPTM.index), set(dfTotalProt.index))
missing = set.difference(set(dfPTM.index), shared)

########## 3.3.3 - Drop missing proteins from PTM list and reorder
dfPTM = dfPTM.drop(missing)
reorder = np.argsort(dfPTM.index)
dfPTM = dfPTM.iloc[reorder]

########## 3.3.4 - Find index of each matching row for dfPTM
newTotProtOrder = []

for ID in dfPTM.index:
    try:
        newPos = list(dfTotalProt.index).index(ID)
        newTotProtOrder.append(newPos)
    except:
        pass

dfTotalProt = dfTotalProt.iloc[newTotProtOrder]


##############################################################################
### 4. Normalize

dfPTMnorm = dfPTM/dfTotalProt
dfPTMnorm = dfPTMnorm.fillna(0)

##############################################################################
### 5. Output
dfPTM = pd.read_excel(
dfPTMnorm.to_excel('C://Users//grossar//Box//Sareen Lab Shared//Data//Andrew//E428 - PTM of iECs//Input data//PTM-expression-normalized-to-protein-expression.xls')
