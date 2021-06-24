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

############################################################################
####### Run Genetic Diversity Checks like LD, HWE, Null Alleles  ###########
############################################################################
##create list to store null alleles
sos_null_all_list <- list()

##list to store all the data frames generated from this analysis 
sos_null_all_df <- c(list(), list(), list(), list(), list(), list(), list())

##create list to store hwe 
sos_hwe_list <- list()

##create list to store linkage disequilbrium
sos_ld_list <- list()

##run summary stats in one loop
for(s in 1:length(sos_genind)){

  ##calculate null alleles 
  sos_null_all_list[[s]] <- null.all(sos_genind[[s]])
  
  ##create null allele table
  sos_null_all_df[[s]] <- matrix(nrow = length(rownames(sos_null_all_list[[s]]$null.allele.freq$summary1)),
                         ncol = length(colnames(sos_null_all_list[[s]]$null.allele.freq$summary1)))
  
  ##save null allele frequency by loci for all species
  sos_null_all_df[[s]] <- sos_null_all_list[[s]]$null.allele.freq$summary1
  
  #bn_hwe <- hw.test(butternutgen_reorg, B = 1000)
  
}

##create null allele table
#bn_null_all_df <- matrix(nrow = length(rownames(bn_null_all$null.allele.freq$summary1)),
 #                        ncol = length(colnames(bn_null_all$null.allele.freq$summary1)))

#bn_null_all_df <- bn_null_all$null.allele.freq$summary1

##bn HWE test
#bn_hwe <- hw.test(butternutgen_reorg, B = 1000)
#bn_hwe_pop <- seppop(butternutgen_reorg) %>% lapply(hw.test, B = 0)

##create table by populations
#bn_hwe_bypop <- sapply(bn_hwe_pop, "[", i = TRUE, j = 3)

##name columns
#colnames(bn_hwe_bypop) <- butternut_24pop_names

##test for LD
#ld_comp <- pair.ia(butternutgen_reorg, sample = 1000)
#ld_comp_df <- data.frame(round(ld_comp,digits = 2))

##write out data files 
#write.csv(bn_hwe, "genetic_analyses_results\\bn_hwe_overall.csv")
#write.csv(bn_hwe_bypops, "genetic_analyses_results\\bn_hwe_bypop.csv")
#write.csv(ld_comp_df, "genetic_analyses_results\\ld_loci.csv")

######################################
############# basic stats ############
######################################
##reorg data file 
bn_sumstats <- summary(butternutgen_reorg)

##create poppr file 
BN_poppr <- poppr(butternutgen_reorg)

##expected heterozygosity 
BN_hexp <- BN_poppr[1:24, 10]
##allele numbers by pop 
BN_nall <- bn_sumstats$pop.n.all
##individual numbers
BN_ind <- BN_poppr[1:24, 2:3]
##allelic richness code 
BN_alleles <-bn_sumstats$pop.n.all/length(butternutgen_reorg@loc.n.all)
BN_all_rich <- colMeans(allelic.richness(butternutgen_reorg)$Ar)	

##create data frame 
butternut_stat_df <- signif(cbind(butternut_mean_lon_lat[,2:3], BN_ind, BN_nall, BN_all_rich, BN_hexp),3)

##name columns and rows 
rownames(butternut_stat_df) <- butternut_24pop_names
colnames(butternut_stat_df) <- c("Mean Longitude", "Mean Latitude", "Number of Individuals", "MLG","Number of Alleles", "Allelic Richness", "Expected Heterozygosity")

##write out csv 
write.csv(butternut_stat_df, "genetic_analyses_results\\butternut_stat_df.csv")