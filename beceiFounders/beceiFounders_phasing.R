library(readr)

rils <- read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/041724_becei_chr_II_Cross_AB_Founder_RIL_GT_table.csv")
#####################################
#### format the input data ##########
founders = rils[,c('FM', 'FA', 'FB')] # founder genotype infered by poolseq
info = rils[,1:4] # df with ID, CHROM, POS, cM

rils = rils[,grepl("_QG", colnames(rils))] # rils
rilnames =  colnames(rils)

# format rils in a numeric genotype matrix with 0 or 1 (or NA value)
rils[rils=="1/1"] = "1"
rils[rils=="0/0"] = "0"
rils[rils=="0/1"]=NA
rils = matrix(as.numeric(unlist(c(rils))), ncol=ncol(rils))
colnames(rils) = rilnames


# Filter out fixed snps or snps with more than 5% NA values
fixed = apply(rils, 1, function(x){x = mean(x, na.rm=T); x %in% c(0,1)})
nafreq = apply(rils,1, function(x){sum(is.na(x))})/ncol(rils)

rils = rils[!fixed & nafreq<0.05,]
info = info[!fixed & nafreq<0.05,]
founders = founders[!fixed & nafreq<0.05,]

# Filter out	SNP that are in high LD  (r > 0.9 or < -0.9) with less than than 3 other SNP within a 500 SNP window.  
whichlowLDSNPs = Search.LowLDSNP(rils, winsizeLD = 500, LDth=0.9, min.nHighlink=3)
rils = rils[-whichlowLDSNPs,]
info = info[-whichlowLDSNPs,]
founders = founders[-whichlowLDSNPs,]


# Split rils by cross A and cross B
rilsA = rils[,grepl("A_QG", colnames(rils))]
rilsB = rils[,grepl("B_QG", colnames(rils))]


#founder A, B and M in format snp x two haplotype (phase don't matter here)
founderM = founders$FM
founderA = founders$FA
founderB = founders$FB
founderM = do.call(rbind, lapply(strsplit(founderM, "/"), function(x){as.numeric(x)}))
founderA = do.call(rbind, lapply(strsplit(founderA, "/"), function(x){as.numeric(x)}))
founderB = do.call(rbind, lapply(strsplit(founderB, "/"), function(x){as.numeric(x)}))


# First, InferFoundersHaploBlocks look at the different possible haplotypes in each genomic window
# the most frequent haplotypes are assumed to be acestral
# each window need to be big enough to diferenciate founders  (default is 100SNP)
# And not too large in genetid distance => if too much recombination, the ancestral haplotype are too broken to be recognize
# If it is the case, the window is subdivided
# Then the major haplotypes are attributed to each founders in the way that minimize the change from the initially inferred genotypes (given in founder1 and founder2)
# The function return the four founder hap order for each window, in a list
founderhaplotypes = InferFoundersHaploBlocks(founderA, founderB, founderM, rilsA, rilsB, info,nsnpWin = 100, maxSizeCM=3)
#This function phase the haplotype block inferred above in a way that minimize the number of breakpoints
phasedfounderhaplotypes = phaseHaplotypes(founderhaplotypes, info, founderA, founderB, founderM)

#=> More details of functyion are provided in "functions_beceiFounders_phasing.R"



# Here is a bit of code that try to remove the snps which are divergent
# witth the poolseq within some window win
# and return some statistics
#foundersAreMajorHaplo: T/F, founders haplotypes are the major haplotypes within the window
#nbroken: number of broken haplotypes
#nfounderhaplo: the number of different founder haplo
#ndivergent: the number of divergent snps with poolseq genotypes
#excludedDivergent: T/F, wether we excluedd the divergent snps

win =22000:23024



lapply(c(F,T), function(excludeDivergingWithPoolseq){
  gisze=diff(info$cM[phasedfounderhaplotypes[c(min(win),max(win)),1]])

  p = phasedfounderhaplotypes[win,]

  fgt = cbind(A=(founderA[p$whichsnp,1]+founderA[p$whichsnp,2])/2,
              B=(founderB[p$whichsnp,1]+founderB[p$whichsnp,2])/2,
              M= (founderM[p$whichsnp,1]+founderM[p$whichsnp,2])/2)

  pgt = cbind(A=(p[,2]+p[,3])/2,
              B=(p[,4]+p[,5])/2,
              M= (p[,6]+p[,7])/2)

  genodist=apply(abs(pgt-fgt),1,sum)

  gdistth=0
  if(excludeDivergingWithPoolseq){
    p=p[!(genodist>gdistth),]
  }

  gtrils=cbind(rilsA, rilsB)[p$whichsnp,]
  gx = cbind(p[,-1], gtrils)

  groups = cutree(hclust(dist(t(gx))), h = 0)
  npergroups = table(groups)
  foundergroups = unique(groups[1:(ncol(p)-1)])
  nfounderhaplo=length(foundergroups)

  #Here we verify that the founders are the major haplotypes
  foundersAreMajorHaplo=sum((names(sort(npergroups, decreasing = T)) %in% foundergroups)[1:nfounderhaplo])==nfounderhaplo
  nbroken=sum(npergroups[!(names(npergroups) %in% foundergroups)])

  data.frame(foundersAreMajorHaplo=foundersAreMajorHaplo,
             nbroken=nbroken,
             nfounderhaplo=nfounderhaplo,
             ndivergent = sum(genodist>gdistth),
             excludedDivergent=excludeDivergingWithPoolseq,
             gisze=gisze)
})




###For chromosome I: snps 30752:30762 were detected as problematic
###phasedfounderhaplotypes = phasedfounderhaplotypes[-(30752:30762),]


# SAVE, and add genotype distance an info do the phased haplotypes
fgtpool=cbind((founderA[phasedfounderhaplotypes$whichsnp,1]+founderA[phasedfounderhaplotypes$whichsnp,2])/2,
      (founderB[phasedfounderhaplotypes$whichsnp,1]+founderB[phasedfounderhaplotypes$whichsnp,2])/2,
      (founderM[phasedfounderhaplotypes$whichsnp,1]+founderM[phasedfounderhaplotypes$whichsnp,2])/2)

fgtphase = cbind((phasedfounderhaplotypes[,2]+phasedfounderhaplotypes[,3])/2,
                 (phasedfounderhaplotypes[,4]+phasedfounderhaplotypes[,5])/2,
                 (phasedfounderhaplotypes[,6]+phasedfounderhaplotypes[,7])/2)

genotypedistance = apply(abs(fgtpool-fgtphase), 1,sum)
sum(genotypedistance>0)/length(genotypedistance)

# phasedfounderhaplotypes=cbind(info[phasedfounderhaplotypes$whichsnp,], phasedfounderhaplotypes[,-1],genotypedistance)
# 
# save(phasedfounderhaplotypes, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240502_beceiFounderPhasedHaplotypes_chrIII.Rdata")



#(founderA[phasedfounderhaplotypes$whichsnp,1]+founderA[phasedfounderhaplotypes$whichsnp,2])/2






















