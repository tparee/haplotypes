require(limSolve)
require(dplyr)
library(data.table)
library(readr)
source("/Users/tomparee/Desktop/happlotype_TUC/Functions_SearchFounderHaplo.R")

#########################################
######### Impor haplotypes #############


genomatrix = read_tsv("~/Desktop/happlotype_TUC/Genotype_Matrix.tsv")
genomatrix = as.data.frame(genomatrix)
genomatrix = subset(genomatrix, CHROM == "I")
snps = genomatrix[,1:3]
colnames(snps) = tolower(colnames(snps))
genomatrix = genomatrix[,5:ncol(genomatrix)]
colnames(genomatrix)

foundernames = c("AB1","CB4507","CB4852","CB4855","CB4856","CB4858","JU319","JU345","JU400","MY1","MY16","N2","PB306","PX174","PX179","RC301")

# Change in -1/0/1 in 0/0.5/1 (I think it should wroks bot way but not sure)
allnames = colnames(genomatrix)
genomatrix = do.call(cbind,lapply(1:ncol(genomatrix), function(i){as.numeric(genomatrix[,i])}))
genomatrix[genomatrix==0] = 0.5
genomatrix[genomatrix==-1] = 0
colnames(genomatrix)=allnames

founders = genomatrix[,colnames(genomatrix) %in% foundernames]
rils = genomatrix[,!(colnames(genomatrix) %in% foundernames)]

# Let's do a subset of 15 rils for example to reduce the computing time
#rils = rils[,sample(1:ncol(rils), 15)]


#######################################
############## RUN ####################
#######################################



#####################################
#### Call haplotype for the lines ###
haplotype = do.call(rbind,lapply(1:ncol(rils), function(thisril){
  print(thisril)
  ril = rils[,thisril]
  ril[ril!=1 & ril!=0]=NA
  output = haplosearch(ril, founders, snps)
  output$rilname = thisril
  output
}))

haplotype=haplotype[!is.na(haplotype[,3]),]

# pos1,pos1 = interval
# rilname = the number of the ril in the order of rils columns
# for each founder a number indicate if there is a match
# 0 no match with the corresponding founder
# 1 = match only with this founder
# 1/n can match with n different founder, which are identical for the corresponding interval
# For these interval with several possible founder:
# => Assume/choose the founder the most frequent among the rils
haplotype = KeepMostFrequentHaplo(haplotype, rils)


#Save
save(haplotype, file = "/Users/tomparee/Desktop/happlotype_TUC/rils_haplotype.Rdata")


