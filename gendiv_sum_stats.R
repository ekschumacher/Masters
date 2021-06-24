#Code below was the first pass at genetic analyses with the increased 
#sampling from butternut's northern population range. 
#These are the next steps in calculating genetic diversity 
#First, we tested linkage disequilibrium, null alleles, 
#and Hardy Weinberg equilibrium.
#Next, a data frame including expected heterozygosity, allelic richness,
#number of alleles, mean longtiude and latitude by population, and 
#individual numbers. This table is included in full in the supplemental text 
#of this manuscript.

######################################
########### Load Libraries ###########
######################################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)

##########################################
############# Set directories ############
##########################################
##set directory to all butternut files 
SOS_drive <- "G:\\My Drive\\Masters\\masters_gen"

##########################################
############ Load in files ###############
##########################################
setwd(SOS_drive)

##load genind files for all species 
species_gen_list <- list.files(path = SOS_drive, pattern = ".gen$")

##create a list to store all the genind files 
species_genind_list <- list()

##write a loop to load these files in 
for(s in 1:length(species_gen_list)){
  
  species_genind_list[[s]] <- read.genepop(species_gen_list[[s]], ncode = 3)
  
}

##list out species names 
species_names <- c("ACMI","ARTR","ASCA","CHVI","ERNA","ERUM","GRSQ")


###########################################
####### Linkage Disequilibrium  ###########
###########################################
##run summary stats over the species 
#linkage disequilibrium
species_ld_list <- list()

##
for(s in 1:length(species_genind_list)){
  
  species_ld_list[[s]] <- pair.ia(species_genind_list[[s]], sample = 1000)
  
}

###combine everything into one table

##test for LD
ld_comp_df <- round(rbind(species_ld_list[[1]],species_ld_list[[2]],species_ld_list[[3]],species_ld_list[[4]],
                          species_ld_list[[5]],species_ld_list[[6]],species_ld_list[[7]],),2)

##write out data files 
write.csv(ld_comp_df, "G:\\My Drive\\Masters\\ld_loci.csv")

######################################
############# basic stats ############
######################################
##reorg data file 
#bn_sumstats <- summary(butternutgen_reorg)

##create poppr file 
#BN_poppr <- poppr(butternutgen_reorg)

##expected heterozygosity 
##BN_hexp <- BN_poppr[1:24, 10]
####allele numbers by pop 
##BN_nall <- bn_sumstats$pop.n.all
####individual numbers
##BN_ind <- BN_poppr[1:24, 2:3]
####allelic richness code 
##BN_alleles <-bn_sumstats$pop.n.all/length(butternutgen_reorg@loc.n.all)
##BN_all_rich <- colMeans(allelic.richness(butternutgen_reorg)$Ar)	
##
####create data frame 
##butternut_stat_df <- signif(cbind(butternut_mean_lon_lat[,2:3], BN_ind, BN_nall, BN_all_rich, BN_hexp),3)
##
####name columns and rows 
##rownames(butternut_stat_df) <- butternut_24pop_names
##colnames(butternut_stat_df) <- c("Mean Longitude", "Mean Latitude", "Number of Individuals", "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")
##
####write out csv 
##write.csv(butternut_stat_df, "genetic_analyses_results\\butternut_stat_df.csv")