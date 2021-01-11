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
library(plyr)
library(dplyr)
library(tidyr)

###############################################
############## Working Directory ##############
###############################################
setwd("G:\\My Drive\\Masters")

##load file 
master_gendiv <- read.csv("master_gendiv.csv")

##set as factor
master_gendiv$Accession <- as.factor(master_gendiv$Accession)

##create new accession names
ne_matrix <- cbind(paste0(master_gendiv$Species,"_",master_gendiv$Accession), as.numeric(master_gendiv$Ne))

##accession list
accession_list <- unique(ne_matrix[,1])

##make ne matrix results
ne_df <- matrix(nrow = length(accession_list), ncol = 5)

##se column matrix
se_matrix <- matrix(nrow = length(ne_matrix[,1]), ncol = 2)
colnames(se_matrix) <- c("Min","Max")

##write loop to calculate means 
for(a in 1:length(accession_list)){
  
  ne_df[a,3] <- mean(as.numeric(ne_matrix[ne_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
  ne_df[a,4] <- (ne_df[a,3]) + std.error(as.numeric(ne_matrix[ne_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
  ne_df[a,5] <- (ne_df[a,3]) - std.error(as.numeric(ne_matrix[ne_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
}  

##update doc
ne_df[,1] <- gsub("\\_.*","",unique(ne_matrix[,1]))
ne_df[,2] <- gsub("^.*\\_","", unique(ne_matrix[,1]))
colnames(ne_df) <- c("Species", "Accession", "Mean", "SE_1","SE_2")
ne_df <- data.frame(ne_df)
ne_df[,2] <- as.factor(ne_df[,2])

##replicate 

# for(r in 1:length(ne_matrix[,1])){

#  se_matrix[r,1] <- replicate(length(ne_matrix[ne_matrix[,1] == paste0(accession_list[[1]]),][,2]), ne_df[a,4])

# }


###############################################
############# Plot Gendiv #####################
###############################################
##na data frame
na_df <- cbind(master_gendiv[,c(1:2,5)])

##barplot
pdf("barplot_accession.pdf")

par(mfrow = c(2,2))

ggplot(master_gendiv, aes(x = Species, y = Na, fill = Accession)) + ylim(c(0,27)) +
  geom_bar(stat = "identity", position = "dodge") + xlab("Species") + 
  ylab("Actual Number of Alleles") + ggtitle("Actual Number of Alleles Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  stat_summary(fun.data = master_gendiv, geom = "errorbars")
    
    
ggplot(data = master_gendiv, aes(x=Species, y=Ne)) + 
  geom_bar(aes(fill = Accession), stat = "identity", position = "dodge") + xlab("Species") + 
  ylab("Effective Number of Alleles") + ggtitle("Effective Number of Alleles by Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) + 
  ylim(c(0,20))

ggplot(data = master_gendiv, aes(x=Species, y=He)) + 
  geom_bar(aes(fill = Accession), stat = "identity", position = "dodge") + xlab("Species") + 
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77"))

ggplot(data = master_gendiv, aes(x=Species, y=Ho)) + 
  geom_bar(aes(fill = Accession), stat = "identity", position = "dodge") + xlab("Species") + 
  ylab("Observed Heterozygosity") + ggtitle("Observed Heterozygosity") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77"))

dev.off()
######
pdf("barplots_species.pdf")

par(mfrow = c(2,2))

barplot(ne_df[,3], beside=T, col = color_names)

dev.off()


