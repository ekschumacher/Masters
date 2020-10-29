###############################################
############### Libraries #####################
###############################################

library(adegenet)
library(poppr)
library(ggplot2)
library(poppr)
library(hierfstat)
library(car)
library(PMCMR)
library(Demerelate)
library(diveRsity)

###############################################
############### Loading Files #################
###############################################
setwd("G:\\My Drive\\Masters")

master_gendiv <- read.csv("master_gendiv.csv")
##set as factor
master_gendiv$Accession <- as.factor(master_gendiv$Accession)

##create new rows 
ne_matrix <- matrix(ncol = 2, nrow = length(master_gendiv[,1]))

##create new accession names
ne_matrix <- cbind(paste0(master_gendiv$Species,master_gendiv$Accession), as.numeric(master_gendiv$Ne))
ne_matrix[,2] <- as.numeric(ne_matrix[,2])

##accession list
accession_list <- unique(ne_matrix[,1])

##make ne matrix results
gendiv_df <- matrix(nrow = length(accession_list), ncol = 5)

##write loop to calculate means 
for(a in 1:length(accession_list)){
  
  gendiv_df[a,1] <- mean(as.numeric(ne_matrix[ne_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
  
  
}


ggplot(data = master_gendiv, aes(x=Species, y=Ho)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) 

