
#This function return a list of "which" window
# e.g, list(1:50, 51:100, ..., 1000:1050, 1051:1060)
# size is the total size (1060 in our example)
# winsize is the windowsize (50)
# minsize is the minimal window size which is applied to the last window
# => if the last is smaller than minsize, the two last window are fused:  1000:1050, 1051:1060) =>  1000:1060)
get.win = function(size, winsize, minsize){
  
  if(winsize < size){
    
    nblocks = size/winsize
    fullblocks = floor(nblocks)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:(fullblocks*winsize), rep(1:fullblocks, each = winsize))
    if(nblocks>fullblocks){ 
      dblock = (nblocks-fullblocks)*winsize
      SPLIT[[length(SPLIT)+1]] = 1:dblock + max(fullblocks*winsize)
    }
    
    if(length(SPLIT[[length(SPLIT)]]) < minsize & length(SPLIT) > 1){
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)-1]],SPLIT[[length(SPLIT)]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}


# Count the number of broken founder haplotype in rils 
# = haplotype which not match with one of the founder gt
# return a data.farme with some info (number of braks, genetic size, the start and the end of the window)
# Right now the function is ugly and inefficient
# I will modify it
# (it works though)
countBreaks = function(windows, rils, foundergt,info){
  do.call(rbind, lapply(1:(length(windows)-1), function(i){
    
    win = windows[i]:windows[i+1]
    if(i %% 1000 ==0 ) print(i)
    fx = foundergt[win,]
    gx = rils[win,]
    
    nbreak = sum(apply(gx, 2, function(x){ ifelse(sum(apply(fx!=x,2,mean,na.rm=T) == 0) == 0, 1,0) }), na.rm=T)
    data.frame(nbreak = nbreak, gsize = info$cM[win[length(win)]] - info$cM[win[1]], start = win[1], end = win[length(win)])
  }))
}




#rils = cbind(rilsA,rilsB)
Search.LowLDSNP = function(rils, winsizeLD = 500, LDth=0.9, min.nHighlink=3){
  
  steps = seq(1, nrow(rils), winsizeLD/2)
  steps = steps[-1]
  
  whichlowLDSNPs = unlist(lapply(steps, function(step){
    
    print(step)
    win = (step-winsizeLD):(step+winsizeLD)
    win = win[win >= 1 & win <= nrow(rils)]
    rx = cor(t(rils[win,]), use = "pairwise.complete.obs")
    
    #For each SNP count the number of other SNP it is in LD with
    nHighlink = apply(rx^2, 2, function(x){sum(x>LDth, na.rm=T)}) -1 #-1 because to account for the diagonal
    
    win[which(nHighlink<min.nHighlink)]
  }))
  
  whichlowLDSNPs=sort(unique(whichlowLDSNPs))
  whichlowLDSNPs
}


# Change the swap value (T|F), which tell us about the phase between window i and window j
# to phases value (1 | -1) which tell us about the phase of all blocks relative to each other
# example:  F,F,F,T,F,T,F => 1, 1, 1,-1,-1, 1, 1
# this function is used within phaseHaplotypes function
SwapToPhase = function(swapvector){
  #swapvector=c(rep(F,10), T, rep(F,10), T,rep(F,10))
  phasevector = rep(1, length(swapvector))
  
  for(swaphere in which(swapvector)){
    phasevector[swaphere:length(phasevector)] = phasevector[swaphere:length(phasevector)]*-1
  }
  
  phasevector
}


# return the percent of match between haplotypes in hapmatrix1 and hapmatrix2
# (hapmatrix1 and hapmatrix2 can be the same)
getHaplotypeSimilarirty = function(hapmatrix1, hapmatrix2){
  
  hapsim = setNames(as.data.frame(t(apply(expand.grid(1:ncol(hapmatrix1), 1:ncol(hapmatrix2)), 1, function(pair){
    h1 = hapmatrix1[,pair[1]]
    h2 = hapmatrix2[,pair[2]]
    percentMatch = sum(h1==h2, na.rm=T)/sum(!is.na(h1) & !is.na(h2))
    c(pair,percentMatch)
  }))), c("h1", "h2", "percentMatch"))
  
  hapsim=hapsim[order(hapsim$h1, hapsim$h2),]
  
  hapsim
  
}



checkPhaseOverlapWindown = function(fgt1,fgt2){
  #fgt1 and fgt2, are the two haplotype matrix [snp,  c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2')]
  # corresponding to the same set of snp
  # for each founder, .g1 and .g2 can be swaped (trans) or not (cis) between fgt1 or fgt2
  # ex: 
  #     fgt1:                          fgt2:
  # FA.g1   FA.g2  ...            FA.g1   FA.g2  
  #   0       1                     1       0
  #   0       1                     1       0
  #   0       1                     1       0
  # => the two matrix fgt1 and fgt2 are compatible in trans
  
 
  
  #get the pariwise similarity data.frame between fgt1 and fgt2
  hapsim = getHaplotypeSimilarirty(hapmatrix1=fgt1, hapmatrix2=fgt2)
  
  # iterate in c(1,3,5) = the "genome1" of founder A, B and M
  # for n in c(1,3,5): c(n,n+1) return c(1,2), c(3,4), c(5:6) => the two haplotype of founder A, B and M
  phase = do.call(cbind,lapply(c(1,3,5), function(n){
    
    
    ## match of the genome 1 of the founder under focus in fgt1 with the genome 1 and genome 2 in fgt2
    matchg1=hapsim[hapsim$h1==n & hapsim$h2 %in% c(n,n+1),] 
    
    ## match of the genome 2 of the founder under focus in fgt1 with the genome 1 and genome 2 in fgt2
    matchg2=hapsim[hapsim$h1==n+1 & hapsim$h2 %in% c(n,n+1),]
    
    # Genome 1 of fgt1 corresponds to genome 1 of fgt2
    # &  genome 2 of fgt1 corresponds to genome 2 of fgt2
    cis = matchg1$percentMatch[1] == 1 & matchg2$percentMatch[2] == 1
    
    # Genome 1 of fgt1 corresponds to genome 2 of fgt2 
    # &  genome 2 of fgt1 corresponds to genome 1 of fgt2
    trans = matchg1$percentMatch[2] == 1 & matchg2$percentMatch[1]  == 1
    
    c(cis=cis, trans=trans)
    
  }))
  
  # the output of the function is a T/F matrix with rows correspond to cis and trans and columns to founders.
  # example:
  ##       FA    FB   FM
  #cis    TRUE  TRUE TRUE
  #trans FALSE FALSE TRUE
  colnames(phase) = data.table::tstrsplit(colnames(fgt1)[c(1,3,5)], ".g")[[1]]
  
  phase
  
}





InferFoundersHaploBlocksWIN = function(win,founderA, founderB, founderM, rilsA, rilsB, info){
  #fgt = founder genotypes infered by poolseq within the target window (1=hom alt; 0.5=het; 0 = hom ref)
  fgt = list(A=(founderA[win,1]+founderA[win,2])/2, 
             B=(founderB[win,1]+founderB[win,2])/2,
             M= (founderM[win,1]+founderM[win,2])/2)
  
  # List crontaining rils genotype within the window for cross A and cross B
  rilgeno = list(A=rilsA[win,], B=rilsB[win,])
  
  # Infer founder haplotype in each cross (=major haplotypes in rils)
  rilshap = lapply(rilgeno, function(gx) {
    
    # Remove RILs with a high proportion of missing values
    gx = gx[, apply(gx, 2, function(x) { sum(is.na(x)) }) / nrow(gx) < 0.3]
    
    # Identify major haplotypes in rils
    grouprils = cutree(hclust(dist(t(gx))), h = 0) # Cluster RILs based on their haplotypes
    npergroup = table(grouprils) # Count the number of RILs in each cluster
    npergroup = npergroup[npergroup > 10] # Keep only clusters with more than 10 RILs (I don't have a good rational for this number, could be more)
    npergroup = sort(npergroup, decreasing = TRUE)  # Sort the clusters by size in decreasing order
    if (length(npergroup) > 4) npergroup = npergroup[1:4] # Limit to the top 4 major haplotypes, if more than 4 are present (we expect four founder haplo at max)
    
    # Gat each major haplotype 
    # (Use info from all RILs from a cluster rather than one, though they are identical, because they are some NA)
    rilshapx = do.call(cbind, lapply(as.numeric(names(npergroup)), function(x) {
      apply(gx[, grouprils == x], 1, mean, na.rm = TRUE)
    }))
    
    rilshapx
  })
  
 # rilshap is a list of two genotype matrix [snp,haplo] corresponding to cross A or B
  
  
  
  # First, let's find the common RILs haplotypes between cross A and cross B (should correspond to the founder M)
  # Here we generate a df with hA the ith column number in cross A, hB the ith column number in cross B, 
  # and percentMatch the % of match between the two (0=>1)
  sharedRilsHaplo = setNames(as.data.frame(t(apply(expand.grid(1:ncol(rilshap[["A"]]), 1:ncol(rilshap[["B"]])), 1, function(pair){
    hA = rilshap[["A"]][,pair[1]]
    hB = rilshap[["B"]][,pair[2]]
    percentMatch = sum(hA==hB, na.rm=T)/sum(!is.na(hA) & !is.na(hB))
    c(pair,percentMatch)
  }))), c("hA", "hB", "percentMatch"))
  
  # Only keep pairs of haplotypes with percentMatch == 1, i.e, haplotype that are the same between cross A and B
  sharedRilsHaplo = sharedRilsHaplo[sharedRilsHaplo$percentMatch==1,]
  
  
  # Classify the RILs haplotype in shared and unshared: 
  # i.e., classify the haplotypes with percentMatch==1 as shared, if not in onlyA or onlyB
  rilshap = list(onlyA=matrix(rilshap[["A"]][,-sharedRilsHaplo$hA], nrow=nrow(rilshap[["A"]])), #Rils haplotypes only found in cross A
                 onlyB=matrix(rilshap[["B"]][,-sharedRilsHaplo$hB], nrow=nrow(rilshap[["B"]])), #Rils haplotypes only found in cross B
                 shared = matrix(rilshap[["B"]][,sharedRilsHaplo$hB], nrow=nrow(rilshap[["B"]]))) # Shared Rils haplotypes
  
  # Transform the list in a matrix with column names indicating shared, onlyA or onlyB
  rilshap = do.call(cbind, lapply(1:length(rilshap), function(i){
    x=rilshap[[i]]
    if(ncol(x)>0) colnames(x) = paste0(names(rilshap)[i], 1:ncol(x)) 
    x
  }))
  
  
  
  # For each founder, calulcate the genotype distance (dist) 
  # Between the genotype infered from poolseq (stored within fgt) and the possible pairs of haplotypes
  # We focus on the haplotype that are relevant for each founder (i.e, founder A, we only look at haplo identified in cross A or shared)
  # A result for each founder in returned into a list with genotype ranked  by genetic distances
  gdistance = lapply(c("A", "B", "M"), function(thisfounder){
    
    # fx = the genotype inferred by poolseq for the founder under consideration
    fx = fgt[[thisfounder]]
    
    # In rilsHhap, restrict the haplotypes that corresponds to this founder
    # If founder M = shared haplo; if A/B, only cross A/B or shared
    targetHaploRils = c("onlyA", "onlyB", "shared")[which(c("A", "B", "M")==thisfounder)]
    targetHaploRils = which(grepl(paste0(targetHaploRils,"|shared"), colnames(rilshap)))
    
    # Get all the possible unique pairwise combinations of the targetted haplotypes
    # i.e. two haplotype within a founder
    comb = expand.grid(targetHaploRils, targetHaploRils)
    comb = comb[comb$Var2 >= comb$Var1,]
    colnames(comb)=c("g1","g2")
    
    # For each pair, caluclate the genotype distance for each n snp with the poolseq infered ones
    # snpdist is a matrix where column are snp, rows corresponds to the possible combinations in combs
    # and value difference between diploid genotype
    snpdist=do.call(rbind,(lapply(1:nrow(comb), function(i){
      g1 = comb[i,1] 
      g2 = comb[i,2] 
      cgt = (rilshap[,g1]+rilshap[,g2])/2 #
      
      #dist = sum(abs(cgt-fx), na.rm=T)
      matrix(abs(cgt-fx), nrow=1)
    })))
    
    # get the genotype distance summed across all sites within the window for each combinations
    windist = apply(snpdist, 1, sum, na.rm=T)
    
    #cbind comb (=g1,g2) with windist (the genotype distance summed across all snp for a given comb), and snpdist (for each snp)
    comb = cbind(comb,windist = windist , snpdist)
    colnames(comb)[4:ncol(comb)] = paste0("snpdist", colnames(comb)[4:ncol(comb)])
    #Order by increasing distance
    # then by ifelse(comb[,"g2"]-comb[,"g1"] ==0,0,1) => This just ensure that if different combinations have the same distance,
    # The combinations with different rils haplotype will be classed higher => maximize the number of rils haplotype represented
    comb = comb[order(comb[,"windist"], ifelse(comb[,"g2"]-comb[,"g1"] ==0,0,1), decreasing = c(F,T)),]
    
    # matrix with colum = g1,g2: the rils haplotype combination (value = ith column in rilsname)
    # windist the genotype distance summed across the windo
    #snpsdist1,2,.. = snp distance for each individual snp
    comb
    
  })
  
  # Now we want to choose the haplotype attribution that minimize the total genetic distance 
  # i.e. taking the first combination (already ranked by distance) for each founder
  
  fhap = unlist(lapply(gdistance, function(x){x[1,c(1,2)]}))
  totdist = sum(unlist(lapply(gdistance, function(x){x[1,3]})))
  
  #do.call(rbind,lapply(gdistance, function(x){x[1,4:ncol(x)]}))
  
  names(fhap) = c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2')
  
  # Verify that all the major haplotypes identified within RILs are attributed to the founders
  representedRilsHap = sort(unique(fhap))
  allRepresented = sum(1:ncol(rilshap) %in% representedRilsHap)==ncol(rilshap)
  
  # If sole of the RILs major haplotypes are not atrributed to one founder,
  # try to see if the missing haplotype only differ of a single snp with other haplotypes
  # if yes, this single snp may be a genotyping error => filter out
  if(!allRepresented){
    
    #The non represented haplotype(s)
    missingHaplo = which(is.na(match(1:ncol(rilshap), representedRilsHap)))
    
    # Generate a table with the number of different snp between the missing haplotypes and other (non-missing) haplotypes
    hapDiff = setNames(as.data.frame(t(apply(expand.grid(missingHaplo, (1:ncol(rilshap))[-missingHaplo]), 1, function(pair){
      h1 = rilshap[,pair[1]]
      h2 = rilshap[,pair[2]]
      ndiff = sum(h1!=h2, na.rm=T)
      c(pair,ndiff)
    }))), c("missHap", "otherHap", "ndifferent"))
    
    # If only one different SNP between the missing hap and another hap,
    # consider that the two haplotype are the same and that the divergent snp is problematic and exclude it
    ## Note: Here we try to infer a robust base for the founder haplotypes with only SNP that we are very confident with
    ## there will be time later to add low quality SNP and other suspicious SNP, so SNP exluded here are not lost forever.
    hapDiff = hapDiff[hapDiff$ndifferent==1,]
    
    # if there are some haplotypes with a divergence of a signle snp with the missing snp,
    # identify the divergent snp ("badsnp")
    if(nrow(hapDiff)>0){
      badsnp = unlist(lapply(1:nrow(hapDiff), function(i){
        x=hapDiff[i,]
        h1 = rilshap[,unlist(x[1])]
        h2 = rilshap[,unlist(x[2])]
        which(h1!=h2)
      }))
      
      # Just a line of code that would indicate a bug
      if(nrow(hapDiff) != length(badsnp)){stop("DEBUG ME")}
      
      # Update after filtering out the suspicious snp (badsnp): 
      allRepresented = T
      win = win[-badsnp]
      rilshap = rilshap[-badsnp,]
      rilgeno = lapply(rilgeno, function(x){x[-badsnp,]})
      totdist=sum(unlist(lapply(gdistance, function(x){x[1,-(3+badsnp)]})))
    }
    
  }
  
  # output as a list:
  ### $win: the window (note: this output window can be slightly different than the input window because algorithm can exclude some suspicious snps)
  ### $foundergt: matrix [snp,founder haplotype] with each of the six column contains a haplotype c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2'),
  ###### => FA,FB,FM = the founders; g1/g2 stands for genome 1 and genome 2 (here which one is genome 1 or genome 2 is random)
  ### $heterozygous: a [snp, founder] matrix where value T | F indicates if the snp is heterozygous among each founders
  ### $rilshap: a vector where value is a number corresponding to unique haplotype, useful to see if a haplotype is reapeated amoing founders
  ### $nbreaks: the number of rils where the founder haplotype blocks are broken
  ### $totdist: the total genotype distanc between those haplotypes and the poolseq genotypes, summed among all founders among all snps
  ### $allRepresented= T or F, wether all the major haplotypes found in RILs were attributed to founders
  ### $gsize = the genetic size of the window
  whichrilshap = fhap
  het = ifelse(fhap[c(1,3,5)]-fhap[c(2,4,6)] == 0,F,T)
  fhap = rilshap[,fhap]
  colnames(fhap) = c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2')
  rilgeno=do.call(cbind, lapply(rilgeno, function(x){x}))
  nbreaks=countBreaks(windows=c(1,length(win)), rilgeno, fhap, info[win,])$nbreak
 
  return(list(win=win,foundergt=fhap, heterozygous = het, rilshap=whichrilshap, nbreaks = nbreaks, totdist = totdist, allRepresented=allRepresented, gsize=diff(range(info$cM[win]))))
}


InferFoundersHaploBlocks = function(founderA, founderB, founderM, rilsA, rilsB, info,nsnpWin = 100, maxSizeCM=3){
  
  
  #Obtain a list of window where value corresponds to the ordered SNP position
  # list(1:100,101:200, ...)
  windows = get.win(nrow(info), nsnpWin, nsnpWin)
  
  # If one of the window is larger than maxSizeCM, split it in two
  windows = do.call(c, lapply(1:length(windows), function(i){
    win = windows[[i]]
    gsize = diff(range(info$cM[win]))
    
    if (gsize < maxSizeCM) {
      # If within limits, return the first marker index of the window
      return(list(win))
    } else {
      # If not, return splitted window
      wheretosplit = which.min(abs(info$cM[win] - median(info$cM[win]) ))
      return(list(win[1:wheretosplit], win[(wheretosplit+1):length(win)]))
    }
    
  }))
  
  # Add intermediate windows which overlap with the adjacent windows
  # list(1:100,101:200, ...) => list(1:100,51:150, 101:200, ...)
  windows = c(do.call(c, lapply(1:(length(windows)-1), function(i){
    wini = windows[[i]]
    winj = windows[[i+1]]
    midwin = c(wini[ceiling(0.5*length(wini)+1):length(wini)], winj[1:floor(0.5*length(winj)-1)])
    
    list(wini,midwin)
  })), list(windows[[length(windows)]]))
  
  
  
  # For each window infer haplotype and attribute two to each founders
  founderhaplotypes = lapply(1:length(windows), function(step){
    #print(step)
    if(step %in% 10 == 0) print(step)
    
    # The window under consideration (this step to the next one -1)
    win =windows[[step]]
    
    # Infer founder haplotype blocks in the window + some stat about inference
    # output as a list;
    ### $win: the window (note: this output window can be slightly different than the input window because algorithm can exclude some suspicious snps)
    ### $foundergt: matrix [snp,founder haplotype] with each of the six column contains a haplotype c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2'),
    ###### => FA,FB,FM = the founders; g1/g2 stands for genome 1 and genome 2 (here which one is genome 1 or genome 2 is random)
    ### $heterozygous: a [snp, founder] matrix where value T | F indicates if the snp is heterozygous among each founders
    ### $rilshap: a vector where value is a number corresponding to unique haplotype, useful to see if a haplotype is reapeated amoing founders
    ### $nbreaks: the number of rils where the founder haplotype blocks are broken
    ### $totdist: the total genotype distanc between those haplotypes and the poolseq genotypes, summed among all founders among all snps
    ### $allRepresented= T or F, wether all the major haplotypes found in RILs were attributed to founders
    ### $gsize = the genetic size of the window
    ###### Note: most of these stat are not used by the following steps of the methods but have been useful to
    ###### develop the algorithm or to debug it
    ###### see the function "InferFoundersHaploBlocksWIN" for details 
    InferFoundersHaploBlocksWIN(win=win,
                                founderA=founderA, founderB=founderB, founderM=founderM,
                                rilsA=rilsA, rilsB=rilsB, 
                                info)
  })
  
  
  # Verify that all the overlaping window are compatible
  # Meaning that, for a given founder, at shared SNP genome 1 of win i should corresonds to genome 1 or 2 of win j
  # If it is not the case, it means that the two different window result in inferred haplotypes that are cannot be the same
  # for each  window i
  compatible = as.data.frame(do.call(rbind, lapply(1:(length(founderhaplotypes)-1), function(i){ 
    #print(i)
    j=i+1 # the next window j 
    fi = founderhaplotypes[[i]]$foundergt # the founder haplotypes for window i
    fj = founderhaplotypes[[j]]$foundergt  # the founder haplotypes for window j
    wini = founderhaplotypes[[i]]$win  # window i
    winj = founderhaplotypes[[j]]$win # window j
    
    winoverlap = wini[wini %in% winj] #common snps between the two windows
    
    #founder haplotypes at the common snps
    ovfi = fi[wini %in% winoverlap,] 
    ovfj = fj[winj %in% winoverlap,]
    
    #  checkPhaseOverlapWindown check that for each founder, 
    #  the phase between the two genotype matrix (cis or trans)
    # ex: 
    #     fgt1:                          fgt2:
    # FA.g1   FA.g2  ...            FA.g1   FA.g2  
    #   0       1                     1       0
    #   0       1                     1       0
    #   0       1                     1       0
    # => the two matrix fgt1 and fgt2 are compatible in trans
    # the output is a T/F matrix with rows correspond to cis and trans and columns to founders.
    # example:
    ##       FA    FB   FM
    #cis    TRUE  TRUE TRUE
    #trans FALSE FALSE TRUE
    iscompatible = checkPhaseOverlapWindown(fgt1=ovfi, fgt2=ovfj)
    
    # compatible if for a founder we have at least a T is cis or trans
    iscompatible = unlist(apply(iscompatible, 2, function(x){ x[1]|x[2]}))
    iscompatible
    
  })))
  
  #identify the incompatible windows
  # ex: incompatibility = 10 => means that window 10 and 11 and incompatible 
  incompatibility = which(apply(compatible,1,function(x){sum(x)<length(x)}) > 0)
  
  
  if( length(incompatibility)>0 ){ # if there are some incompatible window
    
    #For each incompatible window, we are going to to try to redo the haplotype inference
    # by fusing the window ith and jth and see if it helps
    
    # if incompatible window are adjacent (=> diff(incompatibility)<=2)
    # group them together in a list
    # incompatiblility = 10 11 12 40  => incompatiblility = list(10, c(11,12), 40)
    splithere = c(0,which(diff(incompatibility)>2), length(incompatibility))
    incompatibility=lapply(1:(length(splithere)-1), function(i){
      incompatibility[(splithere[i]+1):(splithere[i+1])]
    })
    
    # for each incompatibility groups
    compatibleHaplo = lapply(incompatibility, function(i){
    
      whichwindows = min(i):(max(i)+1) # target the window within the incompatibility groups and the next window as well
      #vector of the target snps
      winall = unlist(lapply(founderhaplotypes[whichwindows], function(x){x$win}))
      winall = sort(unique(winall))
      
      # Redo the haplotype inefrence with the extended window
      newfhap = InferFoundersHaploBlocksWIN(win=winall,
                                            founderA=founderA, founderB=founderB, founderM=founderM,
                                            rilsA=rilsA, rilsB=rilsB, 
                                            info)
      
      
      # Previously, block i was compatible with block h & block j with k
      # So we want to check that we did not change that by with the reestimation of the founder haplotypes
      
      h=min(whichwindows)-1
      k=max(whichwindows)+1
      winh=founderhaplotypes[[h]]$win
      wink=founderhaplotypes[[k]]$win
      
      
      iscompatible = cbind(#Check left compatibility
        checkPhaseOverlapWindown(fgt1=founderhaplotypes[[h]]$foundergt[winh %in% newfhap$win,], 
                                 fgt2=newfhap$foundergt[newfhap$win %in%  winh,]),
        #Check right compatibility
        checkPhaseOverlapWindown(fgt1=founderhaplotypes[[k]]$foundergt[wink %in% newfhap$win,], 
                                 fgt2=newfhap$foundergt[newfhap$win %in%  wink,]))
      
      iscompatible = c(unlist(apply(iscompatible, 2, function(x){ x[1]|x[2]})))
      iscompatible = sum(iscompatible)==length(iscompatible)
      
      # For now, I did not encountered a case where extending the window did not resolve the problem
      # If this happens, the two best solution would be to try to extend a bit more the window
      # and if does not work, just filter out the problematic window
      if(!iscompatible) stop("Code needs to be improved here")
      
      return(newfhap)
      
    })
    
    
    # Just replace the problematic windows with the extended windows
    for(i in 1:length(incompatibility)){
      founderhaplotypes[[ sapply(incompatibility,min)[i] ]] <- compatibleHaplo[[i]]
    }
    
    founderhaplotypes = founderhaplotypes[-unlist(lapply(incompatibility, function(i){i=c(min(i):(max(i)+1))[-1]}))]
    
  }
  
  
  return(founderhaplotypes)
}




phaseHaplotypes = function(founderhaplotypes, info){
  # each founderhaplotypes[[i]]$foundegt, matrix with haplotypes inferred for c('FA.g1', 'FA.g2', 'FB.g1', 'FB.g2', 'FM.g1', 'FM.g2')
  # where FA,FB,FM are the three differents founders and g1 and g2, their diploid genomes
  # g1 and g2 where attributed randomly
  # When considering two windows, g1 and g2 of window i can be in the same phase as g1 and g2 in window j (cis) or swaped (trans)
  # We want to swap g1 and g2 within founders when it is needed so all the window are in the same phase
  # This phasing occurs in to main step
  # 1) The first step is looking the phase between overlapping window by matching genotypes of shared snps
  # This step is only informative when the window i is heterozygous, so this step can only gives us blocks of
  # phased windows that are inetrupted by regions of homozygosity (the blocks not being phased with each other)
  # 2) In the second step, for each founder we consider that the correct phase correspond to the most common
  # combination that still exist in the rils (= higher linkage, = phase that minimize the number of breaks in RILs)
  # Most of the case a phase is really obvious, but sometimes, for one founder two heterozygous chunk of the chromosome
  # are separated by a big homozygous part. In this case, it is difficult to phase the two heterozygous chunk because
  # they are very distant and linkage was highly broken during the RILs. We still chose the phase that is the more frequent
  # in the RILs as the most likely founder phases (minimizing the breaks), but with less confidence
  
  
  # 1) We check the phases between overlapping blocks by checking correspondance of the inferred haplotypes
  #  =>  return a SWAP matrix with value T/F, row = windows in founderhaplotypes and cols = the three founders
  # Each value T in SWAP[i,thisfounder] means that the value window j is swapped relative to window i
  # When one of the two window is homozygous, both phase are correct between window i and j
  # However, only one phase is correct we previous heterozygous windows
  # example:
  # three window 1,2,3,4,... are het,het,hom,het. 
  # We iterate from 1 to nwindow - 1  to phase the window
  # We see that 2 is in cis to 1;
  # for 3 it does not matters (both cis and trans are correct)
  # for 4, the phase is important relative to 1,2 which are heterozygous
  # However, at this step, we compare overlaping windows 3 and 4. As 3 is homozygous, we cannot infer the phase
  # When we encounter these cases where win i is homozygous and win j is heterozygous,
  # we indicate SWAP = NA
  # (What stand between two NA value will be a block of phased windows) 
  
  SWAP = as.data.frame(do.call(rbind, lapply(1:(length(founderhaplotypes)-1), function(i){
    #print(i)
    j=i+1
    #founder haplotypes at windows i and j
    fi = founderhaplotypes[[i]]$foundergt
    fj = founderhaplotypes[[j]]$foundergt
    
    # windows i and j
    wini = founderhaplotypes[[i]]$win
    winj = founderhaplotypes[[j]]$win
    
    winoverlap = wini[wini %in% winj] # overlaping snp between the adjacent windows
    
    # The founder haplotypes at the two overlaping snps
    ovfi = fi[wini %in% winoverlap,]
    ovfj = fj[winj %in% winoverlap,]
    
    phasepair = checkPhaseOverlapWindown(fgt1=ovfi, fgt2=ovfj)
    
    #  checkPhaseOverlapWindown check the phase for each founder based on the match of genotypes
    #  the phase between the two genotype matrix (cis or trans)
    # ex: 
    #     fgt1:                          fgt2:
    # FA.g1   FA.g2  ...            FA.g1   FA.g2  
    #   0       1                     1       0
    #   0       1                     1       0
    #   0       1                     1       0
    
    #phasepair look like this:
    #        FA    FB   FM
    # cis    TRUE  TRUE TRUE
    # trans FALSE FALSE TRUE
    # cis=T: the two windows are already in the same phase
    # trans=T: the two windows are in opposite phase
    # cis & trans = T; happen when fi & fj are homozygous within 
    # cis & trans = F, means that code is not robust to some case and need to be revised
    
    # Here, we transform the cis/trans matrix in a "swap vector" where T = swap (trans); F = no swap (cis)
    sapply(1:ncol(phasepair), function(n){
      swapx = c(F,T)[which(phasepair[,n])]
      
      # If the two swap are possible, check if fj is homozygous
      # if fj is hom, the phase does not matter so swap = T or F would be good, we choose F
      # In other case (ex: fi/fj are hom within winoverlap but fj is still heterozygous winthin winh ):
      # => we cannot infer the phase like this
      # => And will need to use other info later
      if(length(swapx)==2 & founderhaplotypes[[j]]$het[n]){swapx=NA}
      if(length(swapx)==2 & !founderhaplotypes[[j]]$het[n]){swapx=F}
      if(length(swapx)==0) stop("This should not happen. Need some revision")
      
      swapx
    })
    
  })))
  
  SWAP = rbind(rep(F,ncol(SWAP)), SWAP) # Just add a first row with c(F,F,F) so row i corresponds to window i
  
  
  
  # Now, for a given founder, we want to fuse all the window within a phased block
  # and swaping g1 and g2 according to the phases contained in SWAP
  # and accounting the overlaping SNPs (= deleting them)
  # phasedfounderhaplotypes is a list of 3 founders
  # each founders is a data.frame, with columns = block, whichsnp, FX.g1, FX.g2
  phasedfounderhaplotypes = lapply(1:ncol(SWAP), function(thisfounder){
    
    #swapx = the swap vector for thisfounder
    swapx = SWAP[,thisfounder] 
    # breaks = where to split the different phased blocks
    # = where we have NA values in the swap vector
    breaks = c(1,which(is.na(swapx)), length(swapx)+1)
    
    # For each block (i.e, between two breaks):
    PHASED = do.call(rbind,lapply(1:(length(breaks)-1), function(b){
      #print(b)
      wins = breaks[b]:(breaks[b+1]-1) # The windows under considerations
      subswapx = swapx[wins] # The SWAP value for these windows for this founder
      
      # Now we want to change the SWAP value (T|F), which tell us about the phase between window i and window j
      # To phases value (1 | -1) which tell us about the phase of all window relative to each other
      # i.e., subswapx = F,F,F,T,F,T,F => phasex= 1, 1, 1,-1,-1, 1, 1
      phasex = rep(1, length(wins)) 
      for(swaphere in which(subswapx)){
        phasex[swaphere:length(phasex)] = phasex[swaphere:length(phasex)]*-1
      }
      
  
      # Now, we just swap g1 and g2 of thisfounder whenever phasex = -1
      phasedhaplo = lapply(1:length(wins), function(whichwin){
        x = founderhaplotypes[[ wins[whichwin] ]] 
        whichgt = (thisfounder*2)-c(1,0)
        cnames = colnames(x$foundergt[,whichgt])
        
        if(phasex[whichwin] == -1){
          whichgt=whichgt[c(2,1)]
          fgt = x$foundergt[,whichgt]
          colnames(fgt)=cnames
        }else{
          fgt = x$foundergt[,whichgt]
        }
        
        list(win=x$win,foundergt=fgt, gsize=x$gsize)
      })
      
      #Now we want to fuse haplotypes knowing there is an overlap (we want to delete)
      # Will will start deleting overlaping from the second window of this block. We store the first window here.
      firsthap = cbind(block=b,whichsnp=phasedhaplo[[1]]$win, phasedhaplo[[1]]$foundergt) 
      
      # wid is the a vector of iterations between window 2:n
      # if n window <= 2, need to be a bit modified
      if(length(phasedhaplo)>2) wid = 2:length(phasedhaplo)
      if(length(phasedhaplo)==2) wid = 2
      if(length(phasedhaplo)<2) wid = NULL
      
      
      if(!is.null(wid)){
        # for each window supress the snp overlaping with the previous window
        # Then bind all the window together
        # This gives phased haplotype (for this block)
        phasedhaplo=as.data.frame(do.call(rbind, lapply(wid, function(whichwin){
          print(whichwin)
          
          # Phased haplo i and j
          pi= phasedhaplo[[ whichwin-1 ]]$foundergt
          pj = phasedhaplo[[ whichwin ]]$foundergt
          
          # window i and j
          wini = phasedhaplo[[ whichwin-1 ]]$win
          winj = phasedhaplo[[ whichwin ]]$win
          
          winoverlap = wini[wini %in% winj] #overlaping snps between i and j
          
          #Verification that the window are in good phase to be sure (their overlaping genotype should be the same)
          if(sum(abs(c(pi[wini %in% winoverlap,] - pj[winj %in% winoverlap,]))) > 0){
            stop("Windows not correctly phased. Should not happens. There is a bug.")}
          
          #return the founder genotype and window for j, exlucing the overlaping snps with i
          cbind(block=b,whichsnp=winj[!(winj %in% winoverlap)], pj[!(winj %in% winoverlap),])
        })))
        
        
        phasedhaplo=rbind(firsthap,phasedhaplo) # bind with the first window we previously stored
        
      }else{
        
        phasedhaplo=firsthap
        
      }
      
      
      
    }))
    
  })
  
  names(phasedfounderhaplotypes) = c("FA", "FB", "FM")
  
 
  # phasedfounderhaplotypes[[thisfounder]] looks like this:
  # block, whichsnp, FX.g1, FX.g2
  #   1       1        0     0
  #   1       2        0     1
  #   1       3        0     1
  #   ...     
  #   2     1230       1     0
  #   2     1231       0     1
  #   ...
  # Haplotypes being phased within a given block, but not between different blocks
  # Here, we are going to phase the different block by looking at the different adjacent block combinations
  # that are existing in the rils. We assume that the most ferquent combination is the founder phase
  phasedfounderhaplotypes=lapply(1:3, function(thisfounder){
    
    #phap = the phased haplotype data.frame for thisfounder
    phap = as.data.frame(phasedfounderhaplotypes[[thisfounder]]) 
    
    # For each block, we wan the phase 
    blockphase=do.call(rbind, lapply(1:(max(phap$block)-1), function(i){
      #print(i)
      j=i+1
      #phased haplo at block i and j
      pi=phap[phap$block==i,]
      pj=phap[phap$block==j,]
      
      # genetic size of block i + j
      glength = diff(range(c(info$cM[pi$whichsnp],info$cM[pj$whichsnp])))
      
      # If genetic size too long (>10 cM), we are going to trim a bit the blocks
      # We trim the block i and j to reduce their size, while ensuring we keep enough 
      # heterozygous snp in each block to be able to infer the phase
      # (if we trim out all the heterozygous snps, we cannot tell anything about the phase)
      if(glength>10){
        #heterozygous snps in block i and block j
        hetsnpsi = which(pi[,3]-pi[,4] != 0) 
        hetsnpsj = which(pj[,3]-pj[,4] != 0)
        
        # trim pi and pj at the "exterior sides' , keeping 10 heterozygous snps for each
        pi=pi[hetsnpsi[length(hetsnpsi)-10]:nrow(pi),]
        pj=pj[1:hetsnpsj[10],]
      }
      
      
      # Now we count the number of combinations between blocki g1/g2 and block j g1/g2 in rils
      # We only keep the relevant rils panel for the founder
      if(names(phasedfounderhaplotypes)[thisfounder]=="FA") rx = rilsA #rils A if founder A
      if(names(phasedfounderhaplotypes)[thisfounder]=="FB") rx = rilsB # rils B if founder B
      if(names(phasedfounderhaplotypes)[thisfounder]=="FM") rx = cbind(rilsA,rilsB) # rils A and B if founder M
      
      #gi and gj and the founder haplotypes and rils binded for block i and j
      gi = cbind(pi[,3:4], rx[pi$whichsnp,])
      gj = cbind(pj[,3:4], rx[pj$whichsnp,])
      
      # Just exclude lines with high proportion of NA value
      gi=gi[,apply(gi,2,function(x){sum(is.na(x))})/ncol(gi) < 0.3]
      gj=gj[,apply(gj,2,function(x){sum(is.na(x))})/ncol(gj) < 0.3]
      
      #Just keep the common rils (as some rils with high NA value can be filtered in only one block)
      commonrils = colnames(gi)[colnames(gi) %in% colnames(gj)]
      gi = gi[,commonrils]
      gj = gj[,commonrils]
      
      # Group founder and rils by haplotypes
      groupsi = cutree(hclust(dist(t(gi))), h = 0) 
      groupsj = cutree(hclust(dist(t(gj))), h = 0) 
      # here, FM.g1 and FM.g2 are always groups 1 and 2
      
      # Just something needed when the first block is homozgous
      if(i==1 & groupsi[1]==groupsi[2]){return(data.frame(block=j, swap=F))}
      
    
      #fcomn = a 2x2 matrix with the number of rils with the different combination of g1 and g2 at the two block
      # 4 possible combinations: cis: g1_i /g1_j; g2_i /g2_j; trans: g1_i /g2_j; g2_i /g1_j
      fcomb = table(groupsi[groupsi %in% groupsi[1:2] & groupsj %in% groupsj[1:2]], 
                    groupsj[groupsi %in% groupsi[1:2] & groupsj %in% groupsj[1:2]])

      # Count the number of cis and trans
      ncis=sum(fcomb[cbind(c(1,2), c(1,2))])
      ntrans=sum(fcomb[cbind(c(1,2), c(2,1))])
      
      # Some threshold to warn that there is a case ncis and ntrans are close
      if(abs(ncis-ntrans)<5){ warning("CHECK HERE: ODD CASE?")}
      
      # swap = F if ncis > ntrans, and vice versa
      swap=c(F,T)[which.max(c(ncis,ntrans))] 
      
      return(data.frame(block=j, swap=swap))
      
    }))
    
    
    # Change the swap value (T|F), which tell us about the phase between window i and window j
    # to phases value (1 | -1) which tell us about the phase of all blocks relative to each other
    # example:  F,F,F,T,F,T,F => 1, 1, 1,-1,-1, 1, 1
    blockphase$phase=SwapToPhase(swapvector=blockphase$swap)
    
    # Bind all blocks together, accounting for the phases and deleting redundantn snps
    phap=do.call(rbind, lapply(split(phap, phap$block), function(x){
      cnames=colnames(x)
      if(x$block[1]==1){phasex = 1}else{phasex=blockphase$phase[blockphase$block==x$block[1]]}
      if(phasex==-1){x[,3:4]=x[,4:3]}
      x=x[,-which(colnames(x)=="block")]
      return(x)
    }))
    
    phap=phap[!duplicated(phap$whichsnp),]
    phap=phap[order(phap$whichsnp),]
    
    
    as.data.frame(phap)
    
  })
  
  
  
  phasedfounderhaplotypes=cbind(whichsnp=phasedfounderhaplotypes[[1]]$whichsnp, 
                                do.call(cbind, lapply(phasedfounderhaplotypes, function(x){x[,-1]})))
  
  
  
  return(phasedfounderhaplotypes)
}



# findDivergentSnps = function(hapmatrix1, hapmatrix2){
#   
#   phasex=checkPhaseOverlapWindown(fgt1=hapmatrix1,fgt2=hapmatrix2)
#   iscompatible = unlist(apply(phasex, 2, function(x){ x[1]|x[2]}))
#   
#   if(!sum(iscompatible)==length(iscompatible)){
#     
#     whichincomp = which(iscompatible==F)
#     
#     divergentsnps = unlist(lapply(whichincomp, function(wx){
#       f1x = hapmatrix1[,(2*wx)-c(1,0)]
#       f2x = hapmatrix2[,(2*wx)-c(1,0)]
#       
#       hsim=getHaplotypeSimilarirty(hapmatrix1=f1x, hapmatrix2=f2x)
#       swap = sum(hsim$percentMatch[c(1,4)]) < sum(hsim$percentMatch[c(2,3)])
#       if(swap){f2x = f2x[,c(2,1)]}
#       divsnp = which(apply(f1x-f2x, 1, function(g){sum(g)!=0}))
#       #divsnp = overlapwin[divsnp]
#       divsnp
#     }))
#     
#     divergentsnps=sort(unique(divergentsnps))
#     return(divergentsnps)
#   }else{
#     return(NULL)
#   }
#   
# }
# 

# return the broken founder haplotype
# Not used 
#whichRilsBroken = function(rilsgt, foundergt){
#  which(apply(rilsgt, 2, function(x){ ifelse(sum(apply(foundergt!=x,2,mean,na.rm=T) == 0) == 0, T,F) }))
#}


# highdistwin=which(unlist(lapply(founderhaplotypes, function(x){ x$totdist/length(x$win) }))>0)
# splithere=c(0,which(diff(highdistwin)>1))
# highdistwin=lapply(1:(length(splithere)-1), function(i){
#   highdistwin[(splithere[i]+1):(splithere[i+1])]
# })
# 
# 
# i=75
# 
# whichwin = highdistwin[[i]]
# extwhichwin=c(whichwin[1]-1, whichwin, whichwin[length(whichwin)]+1)
# #extwhichwin = whichwin
# 
# highdisthap = founderhaplotypes[ extwhichwin ]
# extwin=sort(unique(unlist(lapply(highdisthap, function(x){x$win}))))
# #extwin = extwin[-which(extwin==33850)]
# sum(unique(unlist(lapply(highdisthap, function(x){x$totdist}))))
# sum(unique(unlist(lapply(highdisthap, function(x){x$nbreaks}))))
# 
# exthap = InferFoundersHaploBlocksWIN(win=extwin,
#                                      founderA=founderA, founderB=founderB, founderM=founderM,
#                                      rilsA=rilsA, rilsB=rilsB, 
#                                      info)
# 
# 
# diff(DivergentsSnps)
# 
# DivergentsSnps = sort(unique(unlist(lapply(highdisthap, function(x){
#   
#   overlapwin = exthap$win[exthap$win %in% x$win]
#   f1 = exthap$foundergt[exthap$win %in% overlapwin,]
#   f2 = x$foundergt[x$win %in% overlapwin,]
#   
#   divsnps = findDivergentSnps(hapmatrix1=f1, hapmatrix2=f2)
#   overlapwin[divsnps]
#   
# }))))
# 
# 
# if(length(DivergentsSnps)>0){
#   extwin = extwin[-which(extwin %in% DivergentsSnps)]
#   exthap = InferFoundersHaploBlocksWIN(win=extwin,
#                                        founderA=founderA, founderB=founderB, founderM=founderM,
#                                        rilsA=rilsA, rilsB=rilsB, 
#                                        info)
#   DivergentsSnps = sort(unlist(lapply(highdisthap, function(x){
#     
#     overlapwin = exthap$win[exthap$win %in% x$win]
#     f1 = exthap$foundergt[exthap$win %in% overlapwin,]
#     f2 = x$foundergt[x$win %in% overlapwin,]
#     
#     divsnps = findDivergentSnps(hapmatrix1=f1, hapmatrix2=f2)
#     overlapwin[divsnps]
#     
#   })))
#   
# }
# 
# 
# 
# 

# win =30000:32000
# 
# #win=35000:36548
# 
# lapply(c(F,T), function(excludeDivergingWithPoolseq){
#   gisze=diff(info$cM[phasedfounderhaplotypes[c(min(win),max(win)),1]])
#   
#   p = phasedfounderhaplotypes[win,]
#   
#   fgt = cbind(A=(founderA[p$whichsnp,1]+founderA[p$whichsnp,2])/2, 
#               B=(founderB[p$whichsnp,1]+founderB[p$whichsnp,2])/2,
#               M= (founderM[p$whichsnp,1]+founderM[p$whichsnp,2])/2)
#   
#   pgt = cbind(A=(p[,2]+p[,3])/2, 
#               B=(p[,4]+p[,5])/2, 
#               M= (p[,6]+p[,7])/2)
#   
#   genodist=apply(abs(pgt-fgt),1,sum)
#   
#   
#   if(excludeDivergingWithPoolseq){
#     gdistth=0
#     p=p[!(genodist>gdistth),]
#   }
#   
#   gtrils=cbind(rilsA, rilsB)[p$whichsnp,]
#   gx = cbind(p[,-1], gtrils)
#   
#   groups = cutree(hclust(dist(t(gx))), h = 0) 
#   npergroups = table(groups)
#   foundergroups = unique(groups[1:(ncol(p)-1)])
#   nfounderhaplo=length(foundergroups)
#   
#   #Here we verify that the founders are the major haplotypes
#   foundersAreMajorHaplo=sum((names(sort(npergroups, decreasing = T)) %in% foundergroups)[1:nfounderhaplo])==nfounderhaplo
#   nbroken=sum(npergroups[!(names(npergroups) %in% foundergroups)])
#   
#   data.frame(foundersAreMajorHaplo=foundersAreMajorHaplo,
#              nbroken=nbroken,
#              nfounderhaplo=nfounderhaplo,
#              ndivergent = sum(genodist>gdistth),
#              excludedDivergent=excludeDivergingWithPoolseq,
#              gisze=gisze)
# })
# 
# 
# phasedfounderhaplotypes = phasedfounderhaplotypes[-(30752:30762),]
# 
# #Problematic snp for chromosome I: 30752:30762
# 
# 
# phasedfounderhaplotypes=cbind(info[phasedfounderhaplotypes$whichsnp,], phasedfounderhaplotypes[,-1])

# save(founderhaplotypes, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/240501_beceiFounderPhasedHaplotypes_chrI.Rdata")
# 
# 
# 
# hapx = apply(gtrils[, groups[-(1:(ncol(p)-1))] == 12], 1, mean, na.rm = TRUE)
# 
# getHaplotypeSimilarirty(hapmatrix1=p[,2:7], hapmatrix2=matrix(hapx, ncol=1))
# hapsim=getHaplotypeSimilarirty(hapmatrix1=p[,2:7], hapmatrix2=matrix(hapx, ncol=1))
# 
# View(cbind(p,hapx,genodist))
# which(p[,1+6] - hapx != 0)
# which(p[,1+1] - p[,1+5] != 0)
# 
# p$whichsnp[which(p[,1+6] - hapx != 0)]
# 
# 
# which(unlist(lapply(founderhaplotypes, function(x){30752 %in% x$win})))
# 
# do.call(cbind,lapply(1:3, function(n){
#   
#   gfix = groupspi[(2*n)-c(1,0)]
#   gfjx = groupspj[(2*n)-c(1,0)]
#   
#   
#   if(n == 1){rx = rilsA[winall,]}
#   if(n == 2){rx = rilsB[winall,]}
#   if(n == 3){rx = cbind(rilsA,rilsB)[winall,]}
#   
#   grix = groupsi[colnames(rx)]
#   grjx = groupsj[colnames(rx)]
#   
#   targetcomb = cbind(gfix, gfjx)
#   
#   do.call(cbind,lapply(1:2, function(g){
#     gtrils = rx[,grix==targetcomb[g,1] & grjx==targetcomb[g,2]]
#     
#     grpx = cutree(hclust(dist(t(gtrils))), h = 0) 
#     ngrpx = table(c(grpx))
#     ngrpx = ngrpx[ngrpx>1]
#     ngrpx = as.numeric(names(ngrpx))
#     
#     hapx = do.call(cbind, lapply(ngrpx, function(thisgroup){
#       apply(gtrils[, grpx == thisgroup], 1, mean, na.rm = TRUE)
#     }))
#     
#     which(hapx[,1]-phasedfounderhaplotypes[min(wini):max(winj),1+(2*n)-abs(2-g)] != 0)
#     
#     
#     
#     founder
#     
#     
#     
#     hapx
#   }))
#   
# }))
# 
# 
# 
# fcomb = table(groupsi[groupsi %in% groupsi[1:6] & groupsj %in% groupsj[1:6]], 
#               groupsj[groupsi %in% groupsi[1:6] & groupsj %in% groupsj[1:6]])
# 
# 
# table(groupsi, 
#       groupsj)
# 
# 
# 
# View(cbind(cbind(rilsA,rilsB)[extwin,names(which(groupsi == 3 & groupsj==3)[-1])], exthap$foundergt))
# 
# View(cbind(cbind(rilsA,rilsB)[extwin[-which(extwin %in% DivergentsSnps)],names(which(groupsi == 3 & groupsj==3)[-1])], exthap$foundergt[-which(extwin %in% DivergentsSnps),]))
# 
# 
# 
# 
# # i=650
# # j=700
# # fi = founderhaplotypes[[i]]$foundergt
# # wini = founderhaplotypes[[i]]$win
# # fj = founderhaplotypes[[j]]$foundergt
# # winj = founderhaplotypes[[j]]$win
# # 
# # winall = sort(unique(unlist(lapply(founderhaplotypes[i:j], function(x){x$win}))))
# # 
# # gi = cbind(fi, cbind(rilsA, rilsB)[wini,])
# # gj = cbind(fj, cbind(rilsA, rilsB)[winj,])
# # 
# # 
# # groupsi = cutree(hclust(dist(t(gi))), h = 0) 
# # table(groupsi)
# # 
# # groupsj = cutree(hclust(dist(t(gj))), h = 0) 
# # table(groupsj)
# # 
# # 
# # groupsfi = groupsi[1:ncol(fi)]
# # groupsfj = groupsj[1:ncol(fj)]
# # 
# # 
# # do.call(cbind,lapply(1:3, function(n){
# #   
# #   gfix = groupsfi[(2*n)-c(1,0)]
# #   gfjx = groupsfj[(2*n)-c(1,0)]
# #   
# #   
# #   if(n == 1){rx = rilsA[winall,]}
# #   if(n == 2){rx = rilsB[winall,]}
# #   if(n == 3){rx = cbind(rilsA,rilsB)[winall,]}
# #   
# #   grix = groupsi[colnames(rx)]
# #   grjx = groupsj[colnames(rx)]
# #   
# #   if( length(c(unique(gfix), unique(gfjx))) == 4){
# #     fcomb = table(groupsri[grix %in% gfix & grjx %in% gfjx],
# #                   groupsrj[grix %in% gfix & grjx %in% gfjx])
# #     
# #     swap=c(F,T)[which.max(c(sum(fcomb[cbind(c(1,2), c(1,2))]),  #cis numbers
# #                             sum(fcomb[cbind(c(1,2), c(2,1))])))] # trans numbers
# #     
# #     if(swap){gfjx=gfjx[c(2,1)]}
# #   }
# #   
# #   
# #   targetcomb = cbind(gfix, gfjx)
# #   
# #   do.call(cbind,lapply(1:2, function(g){
# #     gtrils = rx[,grix==targetcomb[g,1] & grjx==targetcomb[g,2]]
# #     
# #     grpx = cutree(hclust(dist(t(gtrils))), h = 0) 
# #     ngrpx = table(grpx)
# #     ngrpx = ngrpx[ngrpx>1]
# #     
# #     if(length(ngrpx)==1){
# #       hapx = apply(gtrils[, grpx == as.numeric(names(ngrpx))], 1, mean, na.rm = TRUE)}else{
# #         hapx=rep(NA,nrow(gtrils))
# #       }
# #     
# #     hapx
# #   }))
# #   
# # }))
# # 
# # 
# # 
# # fcomb = table(groupsi[groupsi %in% groupsi[1:6] & groupsj %in% groupsj[1:6]], 
# #               groupsj[groupsi %in% groupsi[1:6] & groupsj %in% groupsj[1:6]])
# # 
# # 
# # table(groupsi, 
# #       groupsj)
# # 
# # 
# # 
# # View(cbind(cbind(rilsA,rilsB)[extwin,names(which(groupsi == 3 & groupsj==3)[-1])], exthap$foundergt))
# # 
# # View(cbind(cbind(rilsA,rilsB)[extwin[-which(extwin %in% DivergentsSnps)],names(which(groupsi == 3 & groupsj==3)[-1])], exthap$foundergt[-which(extwin %in% DivergentsSnps),]))
# # 
# # 
# # 
# # 
