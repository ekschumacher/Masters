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
############## Working Directory ##############
###############################################
setwd("G:\\Shared drives\\Emily_Schumacher\\masters")

##load file 
master_gendiv <- read.csv("master_gendiv.csv")

gendiv_species <- read.csv("gendiv_byspecies.csv")

gendiv_species <- gendiv_species[-21,-9]

##Stacked barplot
#pdf("stack_bp.pdf")
#ggplot(gendiv_species, aes(x = Ne)) + 
 # geom_histogram(binwidth = 1, aes(fill = Species)) + 
  #xlab("Effective Number of Alleles") + xlim(0,6) +
 # scale_color_brewer(palette="Blues") +
 # scale_fill_brewer(palette="Blues")
#dev.off()

###
na_lm <- lm(master_gendiv$Na ~ master_gendiv$Species*master_gendiv$Accession)

####
setwd("G:\\Shared drives\\Emily_Schumacher\\masters\\masters_gen")

##genind list 
sos_genind_list <- list.files(pattern = ".gen$")

##create a list to store genind files
sos_genind <- list()

##poppr list
sos_poppr <- list()

##all rich list
sos_allrich_list <- list()

##make a list for hexp
sos_hexp_list <- list()

##now read in a genepop files
for(i in 1:length(sos_genind_list)){
  
  sos_genind[[i]] <- read.genepop(sos_genind_list[[i]], ncode = 3)
  
  sos_poppr[[i]] <- poppr(sos_genind[[i]])
  
  sos_allrich_list[[i]] <- colMeans(allelic.richness(sos_genind[[i]])$Ar)
  
  sos_hexp_list[[i]] <- sos_poppr[[i]]$Hexp
  
}

##calculate allelic richness 

##species name list 
species_names <- unique(master_gendiv$Species)

#trim matrix
master_gendiv <- master_gendiv[,-c(7,10)]

##


##matrix 
pvalue_matrix_accessions <- matrix(nrow = length(species_names), ncol = length(master_gendiv[,5:9]))
rownames(pvalue_matrix_accessions) <- species_names
colnames(pvalue_matrix_accessions) <- names(master_gendiv[,5:9])

##tests 
##set up p value test
for(p in 1:length(species_names)){
  
  pvalue_matrix_accessions[p,1] <- kruskal.test(master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,5]~master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,2])[[3]]
  
  pvalue_matrix_accessions[p,2] <- kruskal.test(master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,6]~master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,2])[[3]]
  
  pvalue_matrix_accessions[p,3] <- kruskal.test(master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,7]~master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,2])[[3]]
  
  pvalue_matrix_accessions[p,4] <- kruskal.test(master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,8]~master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,2])[[3]]
  
  pvalue_matrix_accessions[p,5] <- kruskal.test(master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,9]~master_gendiv[master_gendiv$Species == paste0(species_names[[p]]),][,2])[[3]]
  
  
}
write.csv(pvalue_matrix_accessions,"pvalue_matrix_accessions.csv")

##now create a matrix for by species 
pvalue_matrix_species <- matrix(nrow = 1, ncol = length(master_gendiv[,5:9]))
rownames(pvalue_matrix_species) <- "pvalues"
colnames(pvalue_matrix_species) <- c(names(master_gendiv[,5:9]))

pvalue_matrix_species[,1] <- kruskal.test(master_gendiv[,5]~master_gendiv[,1])[[3]]
pvalue_matrix_species[,2] <- kruskal.test(master_gendiv[,6]~master_gendiv[,1])[[3]]
pvalue_matrix_species[,3] <- kruskal.test(master_gendiv[,7]~master_gendiv[,1])[[3]]
pvalue_matrix_species[,4] <- kruskal.test(master_gendiv[,8]~master_gendiv[,1])[[3]]
pvalue_matrix_species[,5] <- kruskal.test(master_gendiv[,9]~master_gendiv[,1])[[3]]
##do multiple comparison tests 
posthoc.kruskal.nemenyi.test(master_gendiv[,5]~as.factor(master_gendiv[,1]), data = master_gendiv)
posthoc.kruskal.nemenyi.test(master_gendiv[,6]~as.factor(master_gendiv[,1]), data = master_gendiv)
posthoc.kruskal.nemenyi.test(master_gendiv[,7]~as.factor(master_gendiv[,1]), data = master_gendiv)
posthoc.kruskal.nemenyi.test(master_gendiv[,8]~as.factor(master_gendiv[,1]), data = master_gendiv)
posthoc.kruskal.nemenyi.test(master_gendiv[,9]~as.factor(master_gendiv[,1]), data = master_gendiv)

write.csv(pvalue_matrix_species, "pvalue_byspecies.csv")


##calculate NA by species 
na_list <- list()
ne_list <- list()


for(s in 1:length(species_names)){
  
  na_list[[s]] <- mean(master_gendiv[master_gendiv$Species == paste0(species_names[[s]]),][,5])
  
  ne_list[[s]] <- mean(master_gendiv[master_gendiv$Species == paste0(species_names[[s]]),][,6])
  
}



###############################################
############### Loading Files #################
###############################################
setwd("G:\\My Drive\\Masters")

master_gendiv <- read.csv("master_gendiv.csv")
##set as factor
master_gendiv$Accession <- as.factor(master_gendiv$Accession)


#na_comparisons <- list(c("ACMI", "ASCA"), c("ERUM", "ASCA"),c("GRSQ", "ASCA"))
#ne_comparisons <- list(c("ASCA", "ACMI"), c("ASCA", "CHVI"), c("ERUM", "ASCA"),
                     #  c("GRSQ", "ASCA"))
#ho_comparisons <- list(c("ACMI", "ASCA"), c("ACMI", "CHVI"), c("ACMI", "GRSQ"),
                   #    c("ERUM", "ASCA"), c("CHVI","ERUM"), c("ERNA", "ERUM"),
                     #  c("ERUM","GRSQ"))
#he_comparisons <- list(c("ACMI", "ASCA"), c("ACMI", "CHVI"), c("ERUM", "ASCA"),
                      # c("GRSQ", "ASCA"))

#f_comparisons <- list(c("ACMI", "CHVI"), c("ACMI","GRSQ"), c("ASCA", "GRSQ"),
                     # c("CHVI", "ERUM"), c("CHVI", "ERUM"), c("ERUM", "GRSQ"))

####By Species comparisons
##NA
pdf("na_species.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=Na)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Actual Number of Alleles") + ggtitle("Actual Number of Alleles Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77"))


##Ne
pdf("ne_species.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=Ne)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Effective Number of Alleles") + ggtitle("Effective Number of Alleles Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  stat_compare_means(comparisons = ne_comparisons)
dev.off()

##Ho
pdf("ho_species.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=Ho)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Observed Heterozygosity") + ggtitle("Observed Heterozygosity") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  stat_compare_means(comparisons = ho_comparisons)
dev.off()

##He
pdf("he_species.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=Ho)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  stat_compare_means(comparisons = he_comparisons)
dev.off()

##F
pdf("f_species.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=Ho)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Inbreeding Coefficient") + ggtitle("Inbreeding by Species") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  ylim(-1,1)
dev.off()

#########By accession
##NA
pdf("na_accession.pdf", width = 8, height = 6)

ggplot(data = master_gendiv, aes(x=Species, y=Na)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Actual Number of Alleles") + ggtitle("Actual Number of Alleles Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) 

dev.off()


###NE

pdf("ne_accession.pdf", width = 8, height = 6)

ggplot(data = master_gendiv, aes(x=Species, y=Ne, fill = Accession)) + 
  geom_boxplot() + xlab("Species") + 
  ylab("Effective Number of Alleles") + ggtitle("Effective Number of Alleles Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) + 
  stat_compare_means(aes(group = Accession), label = "p.signif") + 
  ylim(c(0,20)) 

dev.off()

####HO 
pdf("ho_accession.pdf", width = 8, height = 6)

ggplot(data = master_gendiv, aes(x=Species, y=Ho)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Observed Heterozygosity") + ggtitle("Observed Heterozygosity Between Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) + 
  stat_compare_means(aes(group = Accession), label = "p.signif")

dev.off()


####He
pdf("he_accession.pdf", width = 8, height = 6)

ggplot(data = master_gendiv, aes(x=Species, y=He)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity Between Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) + 
  stat_compare_means(aes(group = Accession), label = "p.signif")

dev.off()

##F statistic
pdf("f_accession.pdf", width = 8, height = 6)

ggplot(data = master_gendiv, aes(x=Species, y=F)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Inbreeding Coefficient") + ggtitle("Inbreeding Coefficient Between Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) + 
  stat_compare_means(aes(group = Accession), label = "p.signif")

dev.off()

##HEXP
pdf("boxplot_hexp_allspecies.pdf", width = 8, height = 6)

ggplot(data = master_gendiv, aes(x=Species, y=He)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + ylim(0,1) +
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity Among Species") 

dev.off()

##plot within species 
pdf("hexp_byaccession.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=He)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  ylim(c(0,1))
dev.off()  

##plot within species 
pdf("hexp_byspecies.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=He)) + 
  geom_boxplot() + xlab("Species") + 
  ylab("Expected Heterozygosity") + ggtitle("Expected Heterozygosity Among Species") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  ylim(c(0,1)) 
dev.off()  


###Inb
pdf("inb_byaccession.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=F)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Inbreeding Coefficient (F)") + ggtitle("Inbreeding Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  ylim(c(-1,1)) 
dev.off()  

#by species
pdf("inb_byspecies.pdf", width = 8, height = 6)
ggplot(data = master_gendiv, aes(x=Species, y=F)) + 
  geom_boxplot() + xlab("Species") + 
  ylab("Inbreeding Coefficient (F)") + ggtitle("Inbreeding Among Species") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  ylim(c(-1,1)) 
dev.off()  

##NA graph 
ggplot(data = master_gendiv, aes(x=Species, y=Na)) + 
  geom_boxplot(aes(fill=Accession)) + xlab("Species") + 
  ylab("Na") + ggtitle("Na Among Accessions") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77"))


##load in diversity simulations
diversity_simulations_path <- "G:\\My Drive\\Masters\\Diversity_Simulations"
setwd(diversity_simulations_path)

##
eff_all <- read.csv("eff_num_all_gendiv.csv")

plot(eff_all[,2]~eff_all[,1])
plot(eff_all[,4]~eff_all[,1], add = TRUE, pch = 1)


#################################################
############# Allelic Richness ##################
#################################################
##load in gen documents
setwd("G:\\My Drive\\Masters\\masters_gen")

##load in genind files 
sos_list <- list.files(pattern = ".gen$")

##genind list 
sos_genind_list <- list()

species_list <- c("ACMI", "ARTR", "ASCA", "CHVI", "ERNA",
                  "ERUM", "GRSQ")

##write loop to load in 
for(g in 1:length(sos_list)){
  
  sos_genind_list[[g]] <- read.genepop(sos_list[[g]], ncode = 3)
  
}

##allelic richness 
##write out new documents with allelic richness 
for(a in 1:length(sos_list)){
  
  temp_ar <- allelic.richness(sos_genind_list[[a]])$Ar
  
  write.csv(temp_ar, paste0(species_list[[a]], "_ar.csv"))
  
}


##calculate allelic richness 
for(a in 1:length(allrich_list)){
  
  temp_ar <- colSums(allelic.richness(sos_genind_list[[a]])$Ar)/length(sos_genind_list[[a]]@loc.n.all)
 
  write.csv(temp_ar, paste0(species_list[[a]], "_allrich_sum.csv"))
  
}





#by species
pdf("inb_byspecies.pdf", width = 8, height = 6)
ggplot(data = allrich_df, aes(x=Species, y=AllRich)) + 
  xlab("Species") + geom_point(aes(Accession))
  ylab("Allelic Richness") + ggtitle("Allelic Richness By Species") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77")) +
  ylim(c(-1,1)) 
dev.off() 

##load arp file 
setwd("G:\\My Drive\\Masters\\masters_gen")

##ERNA load in 
erna_arp <- arp2gen("erna.arp")

##write out all rich
erna_genind <- read.genepop("erna.gen", ncode = 3)

##calculate allrich
erna_ar <- allelic.richness(erna_genind)$Ar

##write ar
write.csv(erna_ar, "erna_ar.csv")

##UPLOAD 
setwd("G:\\My Drive\\Masters")

##upload allelic richness
allrich_species <- read.csv("all_species_ar.csv")

##
pdf("allrich_by_accession.pdf", width = 8, height = 6)

ggplot(data = allrich_species, aes(x=Species, y=AllRich)) + 
 xlab("Species") + geom_boxplot(aes(fill = Accession)) +
  ylab("Allelic Richness") + ggtitle("Allelic Richness By Accession") +
  scale_fill_manual(values = c("gray33", "gray48", "gray77"))

dev.off() 
