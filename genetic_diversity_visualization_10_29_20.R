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
ne_matrix <- cbind(paste0(master_gendiv$Species,"_",master_gendiv$Accession), as.numeric(master_gendiv$Ne))

##na matrix
na_matrix <- cbind(paste0(master_gendiv$Species,master_gendiv$Accession), master_gendiv$Na)

##he matrix
he_matrix <- cbind(paste0(master_gendiv$Species,master_gendiv$Accession), master_gendiv$He)

##
ho_matrix <- cbind(paste0(master_gendiv$Species,master_gendiv$Accession), master_gendiv$Ho)

#F matrix
f_matrix <- cbind(paste0(master_gendiv$Species,master_gendiv$Accession), master_gendiv$F)

##accession list
accession_list <- unique(ne_matrix[,1])

##make ne matrix results
ne_df <- matrix(nrow = length(accession_list), ncol = 3)

##colnames(gendiv_df) <- c("Ne","Na","He","Ho","F")

##write loop to calculate means 
for(a in 1:length(accession_list)){
  
  ne_df[a,3] <- mean(as.numeric(ne_matrix[ne_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
  
}  
#  gendiv_df[a,2] <- mean(as.numeric(na_matrix[na_matrix[,1] == paste0(accession_list[[a]]),][,2]))

ne_df <- data.frame(ne_df)
ne_df[,1] <- accession_list

##now separate
ne_df[,1] <- gsub("\\_.*","",ne_df[,1])
ne_df[,2] <- gsub("^.*\\_","", ne_df[,1])


  
 # gendiv_df[a,3] <- mean(as.numeric(he_matrix[he_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
 # gendiv_df[a,4] <- mean(as.numeric(ho_matrix[ho_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
#  gendiv_df[a,5] <- mean(as.numeric(f_matrix[f_matrix[,1] == paste0(accession_list[[a]]),][,2]))
  
#}

##reformat loop

ne_matrix <- cbind(paste0(master_gendiv$Species,master_gendiv$Accession), as.numeric(master_gendiv$Ne))


##write to data frame
##write loop to calculate means 
for(a in 1:length(accession_list)){
  
  gendiv_df[a,1] <- mean(as.numeric(ne_matrix[ne_matrix[,1] == paste0(accession_list[[a]]),][,2]))

}

##create df 
ne_df <- data.frame(ne_df)
colnames(ne_df) <- c("Species", "Accession","NE")

ne_df[,2] <- as.factor(ne_df[,2])

###Plot 
ggplot(data = ne_df, aes(x=Species, y=Ho)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Inbreeding Coefficient") + ggtitle("Inbreeding by Species") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) 

