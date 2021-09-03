### Plotting DE in boxplots -- Andrew R Gross -- 2018-08-14
### This script will generate boxplots from differential expression data. 
### This script operates downstream of the standard DESeq2 Pipeline
### INPUT: Expression data; differential expression data is typical
### OUTPUT: Boxplot figures

############################################################################################
### Header
library(ggplot2)

############################################################################################
### Functions

plotPCA <- function(data, xAxisPC = 1, yAxisPC = 2, subtitle = '', textSize = 16, xPadding = 0, yPadding = 0, Shortname = data$Shortname, Group = data$Group) {
  return(
    ggplot(data = data, aes(x = data[,xAxisPC + 1], y = data[,yAxisPC + 1], label = Shortname)) +
      geom_point(size = 5, aes(fill = Group), color = 'black', pch = 21, stroke = 3) +
      #geom_text(hjust=-0.2,vjust=0.5) + 
      xlim(min(data[xAxisPC + 1]) - xPadding, max(data[xAxisPC + 1]) + xPadding) +
      ylim(min(data[yAxisPC + 1]) - yPadding, max(data[yAxisPC + 1]) + yPadding) +
      scale_fill_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
      labs(       subtitle = subtitle,
                  x = paste('PC', xAxisPC, ' - ', pca.data.var.per[xAxisPC], '%', sep = ''), 
                  y = paste('PC', yAxisPC, ' - ', pca.data.var.per[yAxisPC], '%', sep = '')) +
      theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
            axis.title.x = element_text(face="bold", size= textSize + 2 ,margin =margin(10,0,10,0)),
            axis.title.y = element_text(face="bold", size= textSize + 2 ,margin =margin(0,10,0,10)),
            panel.background = element_rect(fill = 'white', color = 'black', size = 3),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = textSize)) 
  )
}
############################################################################################
### Input

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/iECs/')
list.files()
ang.data <- read.csv('RDSS-12511--04--28--2021_FPKM.csv', row.names = 1, fileEncoding="UTF-8-BOM")
(ang.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/metadata.csv', fileEncoding="UTF-8-BOM"))

############################################################################################
### Format

### Rename columns
data.frame(names(ang.data),ang.metadata)
names(ang.data) <- ang.metadata$Shortname
summary(ang.data)

### Select data to plot
title = 'E417-all data'
title = 'E417- wo iPSCs'; ang.metadata <- ang.metadata[1:9,]      # Everything but iPSCs

ang.data.f <- ang.data[ang.metadata$Shortname]

### Replace NAs with 0
ang.data.f[is.na(ang.data.f)] <- 1
ang.data.max <- apply(ang.data.f, 1, max)
(cutoff <- quantile(ang.data.max, 0.5))
rows.to.keep <- ang.data.max > cutoff
summary(rows.to.keep)
### Reorder columns
results.ang <- ang.data.f[rows.to.keep,]

### Convert to matrix
results.ang <- as.matrix(results.ang)
summary(results.ang)
############################################################################################
### Calculate Principle Components

### Calculate the actual components
pca.ang <- prcomp(t(results.ang), scale = TRUE)
#pca.ang <- prcomp(t(results.ang), scale = FALSE)

### Calculate the percent variation accounted for by each component
pca.data.var <- pca.ang$sdev^2
pca.data.var.per <- round(pca.data.var/sum(pca.data.var)*100, 1)



############################################################################################
### Plot Data
barplot(pca.data.var.per, main = 'Scree Plot', xlab = 'Principle Component', ylab = 'Percent Variation')

### Define data frame for ggplot
pca.data.to.plot <- data.frame(Sample = rownames(pca.ang$x), 
                               PC1 = pca.ang$x[,1],
                               PC2 = pca.ang$x[,2],
                               PC3 = pca.ang$x[,3],
                               PC4 = pca.ang$x[,4],
                               PC5 = pca.ang$x[,5])

(pca.data.to.plot <- cbind(pca.data.to.plot, ang.metadata[-1]))
#pca.data.to.plot$Group <- pca.data.to.plot$CellLine

### Basic plot of PC1 v PC2
(pca.plot <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC2, label = Condition, color = Condition)) +
    geom_point() + 
    geom_text() + 
    xlab(paste('PC1 - ', pca.data.var.per[1], '%', sep = '')) +
    ylab(paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
    #xlim(c(-50,70)) +
    theme_bw() +
    ggtitle('PCA of E417'))

############################################################################################
### Generate formatted pca plots using function

pca.1.v.2 <- plotPCA(data = pca.data.to.plot, Group = pca.data.to.plot$Condition, xPadding = 20, yPadding = 20)
pca.1.v.2

pca.1.v.2 <- plotPCA(data = pca.data.to.plot, Group = pca.data.to.plot$Condition, xPadding = 10, yPadding = 10)



############################################################################################
### Write to folder
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/PCA/')

tiff(filename= paste0(title, ' - PCA 1v2.tiff'), width = 2000, height = 1600, units = "px", pointsize = 12, res = 250)
pca.1.v.2
dev.off()

tiff(filename= paste0(title, ' - PCA 1v3.tiff'), width = 2000, height = 1600, units = "px", pointsize = 12, res = 250)
pca.1.v.3
dev.off()

