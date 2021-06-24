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
sos_drive <- "G:\\Shared drives\\Emily_Schumacher\\masters\\masters_gen"
gendiv_output_drive <- "G:\\Shared drives\\Emily_Schumacher\\masters\\gendiv_ouput"

##
erna_new <- list.files(pattern = ".arp$")
arp2gen(erna_new)

##########################################
############ Load in files ###############
##########################################
setwd(sos_drive)

##load genind 
sos_genind_list <- list.files(pattern = ".gen$")

##create a list to store genind files
sos_genind <- list()

for(i in 1:length(sos_genind_list)){
  
  sos_genind[[i]] <- read.genepop(sos_genind_list[[i]], ncode = 3)
}

##create list of species names 
species_names <- c("ACMI", "ARTR", "ASCA", "CHVI", "ERNA", "ERUM", "GRSQ")

############################################################################
####### Run Genetic Diversity Checks like LD, HWE, Null Alleles  ###########
############################################################################
##create list to store null alleles
sos_null_all_list <- list()

##list to store all the data frames generated from this analysis 
sos_null_all_df <- c(list(), list(), list(), list(), list(), list(), list())

##create list to store hwe 
sos_hwe_df <- c(list(), list(), list(), list(), list(), list(), list())

##create list to store linkage disequilbrium
sos_ld_list <- c(list(), list(), list(), list(), list(), list(), list())

##run summary stats in one loop
for(s in 1:length(sos_genind)){

  ##calculate null alleles 
  sos_null_all_list[[s]] <- null.all(sos_genind[[s]])
  
  ##create null allele table
  sos_null_all_df[[s]] <- matrix(nrow = length(rownames(sos_null_all_list[[s]]$null.allele.freq$summary1)),
                          ncol = length(colnames(sos_null_all_list[[s]]$null.allele.freq$summary1)))
  
  ##save null allele frequency by loci for all species
  sos_null_all_df[[s]] <- sos_null_all_list[[s]]$null.allele.freq$summary1
  
  ##write out data frames 
  write.csv(sos_null_all_df[[s]], paste0(gendiv_output_drive, "\\", species_names[[s]], "_null_all_df.csv"))
  
  ##separate by population
  sos_hwe_df[[s]] <- seppop(sos_genind[[s]]) %>% lapply(hw.test, B = 1000)
  
  ##write out data frames
  write.csv(sos_hwe_df[[s]], paste0(gendiv_output_drive, "\\", species_names[[s]], "_hwe_df.csv"))
  
  ##test for LD
  sos_ld_list[[s]] <- pair.ia(sos_genind[[s]], sample = 1000)
  
  ##write out data frames
  write.csv(sos_ld_list[[s]], paste0(gendiv_output_drive, "\\", species_names[[s]], "_ld_df.csv"))
  
  
}