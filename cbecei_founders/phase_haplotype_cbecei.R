library(readr)

#########################################################################
################ FUNCTIONS ##############################################



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


InferFoundersHaploBlocks = function(founder1, founder2, rils, info,nsnpWin = 50, maxSizeCM=3){
  
  
  steps = seq(1,nrow(rils),nsnpWin)
  steps[length(steps)]=nrow(rils)+1
  
  # Ensure that genetic distance is not too large
  # If too high => lot of recombination => Though to infer founder haplotype because too much broken
  steps = c(unlist(lapply(1:(length(steps)-1), function(step){
    win = steps[step]:(steps[step+1]-1)
    gsize = diff(range(info$cM[win]))
    if(gsize<maxSizeCM){
      return(win[1])
    }else{
      return(win[c(1, which.max(diff(info$cM[win]))+1 )])
    }
  })), steps[length(steps)])
  
  haplotype = lapply(1:(length(steps)-1), function(step){
    #print(step)
    #if(step %in% 10 == 0) print(step)
    
    win =steps[step]:(steps[step+1]-1)
    fgt1 = (founder1[win,1]+founder1[win,2])/2
    fgt2 = (founder2[win,1]+founder2[win,2])/2
    
    gx = rils[win,]
    gx = gx[,apply(gx, 2, function(x){sum(is.na(x))})/nrow(gx) < 0.3]
    grouprils = cutree(hclust(dist(t(gx))),h=0)
    npergroup = table(grouprils)
    npergroup = npergroup[npergroup>10]
    npergroup = sort(npergroup, decreasing = T)
    if(length(npergroup)>4) npergroup=npergroup[1:4]
    
    #rilshap = gx[,sapply(as.numeric(names(npergroup)), function(x){which(grouprils==x)[1]})]
    
    rilshap=do.call(cbind, lapply(as.numeric(names(npergroup)), function(x){ apply(gx[,grouprils==x],1, mean, na.rm=T) }))
    
    
    comb = expand.grid(1:ncol(rilshap), 1:ncol(rilshap))
    comb = comb[comb$Var2 >= comb$Var1,]
    colnames(comb)=c("g1","g2")
    
    comb = cbind(comb, do.call(rbind, lapply(1:nrow(comb), function(i){
      g1 = comb[i,1] 
      g2 = comb[i,2] 
      cgt = (rilshap[,g1]+rilshap[,g2])/2
      data.frame(dist1 = sum(abs(cgt-fgt1), na.rm=T),
                 dist2 = sum(abs(cgt-fgt2), na.rm=T))
    })))
    
    rank1 = order(comb$dist1)
    rank2 = order(comb$dist2)
    
    fhap = comb[ c(rank1[1], rank2[1]),]
    representedRilsHap = sort(unique(unlist(c(fhap[,1:2]))))
    if(sum(1:ncol(rilshap) %in% representedRilsHap)!=ncol(rilshap)){
      
      combcomb = expand.grid(1:nrow(comb), 1:nrow(comb))
      
      combcomb=cbind(combcomb,do.call(rbind, lapply(1:nrow(combcomb), function(i){
        
        i1 = combcomb[i,1] 
        i2 = combcomb[i,2] 
        c1 = comb[i1,]
        c2 = comb[i2,]
        
        totdist = c1[,"dist1"]+c2[,"dist2"]
        nrepresented = sum(1:ncol(rilshap) %in% sort(unique(unlist(c(c1[,1:2],c2[,1:2])))))
        data.frame(totdist = totdist, nrepresented=nrepresented)
      })))
      
      combcomb = combcomb[combcomb$nrepresented==ncol(rilshap), ]
      combcomb = combcomb[which.min(combcomb$totdist),]
      
      fhap = comb[ c(combcomb[,1], combcomb[,2]),]
      
    }
    
    totdist = fhap[1,3] + fhap[2,4]
    fhap1 = rilshap[,unlist((c(fhap[1,1:2])))]
    fhap2 = rilshap[,unlist((c(fhap[2,1:2])))]
    fgt = cbind(fhap1,fhap2)
    colnames(fgt) = c("founder1.g1","founder1.g2","founder2.g1","founder2.g2")
    
    nbreaks1=countBreaks(windows=c(1,length(win)), gx, fgt, info[win,])$nbreak
    #nbreaks2=sum(!(grouprils %in% (as.numeric(names(npergroup)))))
    #if(nbreaks1 > nbreaks2){print("PROBLEM");print(step)}
    
    return(list(win=win,foundergt=fgt, nbreaks = nbreaks1, totdist = totdist))
    
  })
  
  #bb = do.call(rbind,lapply(test, function(x){
  #  
  #  data.frame(gdis=diff(range(info$cM[x$win])),
  #             nbreaks = x$nbreaks)
  #}))
  
  return(haplotype)
}


#founderhaplotypes = founderhaplotypessave
phaseHaplotypes = function(founderhaplotypes){
  
  i = 1
  #while(i < 379){
  while(i < length(founderhaplotypes)-1 ){
    #print(i)
    
    fi = founderhaplotypes[[i]]$foundergt
    fj = founderhaplotypes[[i+1]]$foundergt
    wini = founderhaplotypes[[i]]$win
    winj = founderhaplotypes[[i+1]]$win
    midwin = c(wini[floor(length(wini)/2):length(wini)], winj[1:floor(length(winj)/2)])
    winall = c(wini, winj)
    #winall = c(winall, max(winall)+1) # Just some add on to compensate a dumb thing that I am to lazy to correct nicely right now
    
    potentialShift = setNames(expand.grid(c(1,-1), c(1,-1)), c("shiftF1", "shiftF2"))
    
    nbreaks = unlist(lapply(1:nrow(potentialShift), function(ii){
      
      if(potentialShift[ii,1]==1){shiftF1=c(1,2)}else{shiftF1=c(2,1)}
      if(potentialShift[ii,2]==1){shiftF2=c(3,4)}else{shiftF2=c(4,3)}
      shift = c(shiftF1,shiftF2)
      
      fused = rbind(fi[wini %in% midwin,], fj[winj %in% midwin,shift])
      countBreaks(windows=c(1,nrow(fused)), rils=rils[midwin,], foundergt=fused, info=info[midwin,])$nbreak
    }))
    
    shift = as.data.frame(potentialShift[which.min(nbreaks),])
    shift = c(apply(shift, 2, function(x){if(x==1){c(1,2)}else{c(2,1)}})) + c(0,0,2,2)
    
    haploall_phased = rbind(fi, fj[,shift])
    nbreaksPhased = nbreaks[which.min(nbreaks)]
    
    haploall_infered = InferFoundersHaploBlocks(founder1=founder1[winall,],
                                                founder2=founder2[winall,],
                                                rils=rils[winall,],
                                                info=info[winall,], 
                                                nsnpWin = length(winall)-1, 
                                                maxSizeCM = 1e6)[[1]] # Just an absuerdely big number so no window subdivision 
    
    nbreaksInfered = haploall_infered$nbreaks
    
    haploall_infered = haploall_infered$foundergt
    
    #hapdist = sum(abs(haploall_infered - haploall_phased),na.rm=T)
    
    countBreaks(windows=c(1,nrow(haploall_infered)), rils=rils[winall,], foundergt=haploall_infered, info=info[winall,])$nbreak
    countBreaks(windows=c(1,nrow(haploall_infered)), rils=rils[winall,], foundergt=haploall_phased, info=info[winall,])$nbreak
    
    choosePhased =  nbreaksInfered >= nbreaksPhased
    
    if(choosePhased){
      founderhaplotypes[[i+1]]$foundergt = fj[,shift]
      i=i+1
    }else{
      founderhaplotypes[[i+1]]$foundergt = haploall_infered
      founderhaplotypes[[i+1]]$win =  c(wini, winj)
      founderhaplotypes[[i+1]]$nbreaks =  nbreaksInfered
      founderhaplotypes = founderhaplotypes[-i]
      i=i-1
    }
    
    
  }
  
  #founderhaplotypes[[380]]
  
  phased = do.call(rbind, lapply(founderhaplotypes, function(x){x$foundergt}))
  
  return(phased)
  
  
}

#########################################################################
#########################################################################
founders <- read_csv("~/Downloads/Cross_ABC_I_Pool_GT_table.csv")
rils <- read_csv("~/Downloads/becei_chr_I_Cross_A_filtered_genotype_table.csv")
founders = founders[match(rils$ID, founders$ID),]

# some formating
info = rils[,1:4]
rils = rils[,7:ncol(rils)]
rils[rils=="1/1"] = "1"
rils[rils=="0/0"] = "0"
rils[rils=="0/1"]=NA
rils = matrix(as.numeric(unlist(c(rils))), ncol=ncol(rils))

# remove fixed snps
fixed = apply(rils, 1, function(x){x = mean(x, na.rm=T); x %in% c(0,1)})
rilsall=rils
infoall = info
rils = rils[!fixed,]
info = info[!fixed,]
founders = founders[!fixed,]

#founder 1 & 2 in format snp x two haplotype
founder1 = founders$FM
founder2 = founders$FA
founder1 = do.call(rbind, lapply(strsplit(founder1, "/"), function(x){as.numeric(x)}))
founder2 = do.call(rbind, lapply(strsplit(founder2, "/"), function(x){as.numeric(x)}))


# First, InferFoundersHaploBlocks look at the different possible haplotypes in each genomic window
# the most frequent haplotypes are assumed to be acestral
# each window need to be big enough to diferenciate founders  (default is 50SNP)
# And not too large in genetid distance => if too much recombination, the ancestral haplotype are too broken to be recognize
# If it is the case, the window is subdivided
# Then the major haplotypes are attributed to each founders in the way that minimize the change from the initially inferred genotypes (given in founder1 and founder2)
# The function return the four founder hap in  c("founder1.g1","founder1.g2","founder2.g1","founder2.g2") order for each window, in a list
founderhaplotypes = InferFoundersHaploBlocks(founder1=founder1, founder2=founder2, rils=rils, info=info)

#This function phase the haplotype block inferred above in a way that minimize the number of breakpoints
founderhaplotypes = phaseHaplotypes(founderhaplotypes)

#As a "control" we can count the number of breaks that happened in different intervals
#(= the haplotype that are not identical to one founder for a given window)
# window by snps:
windows = seq(1,nrow(founderhaplotypes2), 100)
windows[length(windows)]=nrow(founderhaplotypes2)
breakBins = countBreaks(windows, rils, founderhaplotypes2, info)

# we are very close from the number of recombination event one may expect from 5 generation of oucrossing with 1 CO per meiosis
# expected: 148 lines * 5 generation * 0.5 CO per chrom = 370
# observed: 374
sum(breakBins$nbreak)

ggplot()+theme_minimal()+
  geom_point(data=breakBins, aes(gsize, nbreak))+
  geom_point(data=subset(breakBins, nbreak>20 & gsize < 2), aes(gsize, nbreak), color='red')
#=> there is a few suspicious windows (a lot of breakpoint for small genetic distances)
#   I looked at the two worst (in red in the gglot above)
#   There is no obvious pattern indicating that there are not legit
#   => Important haplotype diversity within the rils at this points
#   => breakpoints that are relatively dispered (a phase error would more likely focus the breakpoints at a unique point)
#   If they are indeed legit, it might indicate that the genetic map underestimate genetic distance there
#   or that an higher heterozygosity is maintained there for some reason (overdominant selection), increasing effective recombination
#   or simply selection on recombinant haplotypes
#   Many thing to explore, but interesting. And not so unusual to observe these strange pattern.

#window by genetic distances:
windows = sapply(seq(0,max(info$cM), 5), function(x){ which.min(abs(info$cM-x))})
breakBins = countBreaks(windows, rils, founderhaplotypes2, info)


ggplot()+theme_minimal()+theme(legend.position = "none")+
  geom_point(data=breakBins, aes(x=1:nrow(breakBins), y=nbreak, color = ifelse(nbreak>50, "suspicious?", "OK")))+
  scale_color_manual(values = c("black", "red"))



# save

founderCrossA = matrix(NA, ncol=4, nrow=nrow(rilsall))
colnames(founderCrossA) = c("founderM.g1","founderM.g2", "founderA.g1","founderA.g2")
founderCrossA[fixed] = matrix(apply(rilsall[fixed,], 1, mean, na.rm=T), ncol=1)[1:sum(fixed), rep(1,4)]
founderCrossA[!fixed] = founderhaplotypes
founderCrossA = cbind(infoall,founderCrossA)
save(founderCrossA, file="/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/cbecei_haplo/founderHaplotypesCrossA.Rdata")


#20301
#20301

#21601
#21701

# countBreaks(windows, rils, foundergt, info)
# founderhaplotypes2[21601:21701,]
# breakBins = countBreaks(1:length(21601:21701), rils[21601:21701,], founderhaplotypes2[21601:21701,], info[21601:21701,])
# 
# countBreaks(1:length(20301:20401), rils[20301:20401,], founderhaplotypes2[20301:20401,], info[20301:20401,])
# 
# test = founderhaplotypes2[18890:18910,]
# 
# View(rils[21601:21701,])
# founderhaplotypes2[21601:21701,]
# 
# 
# cutree(hclust(dist(t(rils[20301:20401,]))),h=0)
# 
# grouprils = cutree(hclust(dist(t(rils[20340:20380,]))),h=0)
# npergroup = table(grouprils)
# npergroup = npergroup[npergroup>10]
# npergroup = sort(npergroup, decreasing = T)
# if(length(npergroup)>4) npergroup=npergroup[1:4]

# inferPhaseFromRilsLD = function(founder1, founder2, rils, info, LORth = 3){
#   #founder1 = founders$FM
#   #founder2 = founders$FA
#   
#   heterozygous = founder1 == "0/1"
#   freq = apply(rils, 1, function(x){x = mean(x, na.rm=T)})
#   nmiss = apply(rils, 1, function(x){x = sum(is.na(x))})
#   
#   KEEP = heterozygous & freq>0.075 & nmiss < 30
#   
#   infoall = info
#   founder1all = founder1
#  
#   info = info[KEEP,]
#   rils = rils[KEEP,]
#   founder1 = founder1[KEEP]
#   founder2 = founder2[KEEP]
#   founder1 = do.call(rbind, lapply(strsplit(founder1, "/"), function(x){as.numeric(x)}))
#   founder2 = do.call(rbind, lapply(strsplit(founder2, "/"), function(x){as.numeric(x)}))
#   founder1all = do.call(rbind, lapply(strsplit(founder1all, "/"), function(x){as.numeric(x)}))
#   
#   
#   
#   SHIFT = unlist(lapply(1:(nrow(rils)-1), function(wsnp){
#   #SHIFT = unlist(lapply(1:2000, function(wsnp){
#     if(wsnp %% 1000 ==0) print(wsnp)
#     out = try(phaseshift(i=wsnp,j=wsnp+1, rils=rils, founder2=founder2, info=info), silent = T)
#     if(class(out)=='try-error') out = NA
#     out
#   }))
#   
#  
#   founder1 = cbind(founder1[,c(1,2)], c(SHIFT,NA))
#   
#   
#   nsnp = nrow(founder1)
#   #nsnp = 2000
#   i=0
#   haplotype = list()
#   
#   while(i<(nsnp-1)){
#     
#     start=as.integer(i+1)
#     hap = matrix(founder1[start,], nrow=1)
#     
#     i = start
#     end = F
#     
#     while(!end){
#       shiftx = founder1[i,3]
#       
#       if(i<nsnp){x2 = founder1[i+1,]}else{end = T}
#       
#       if(!end & is.na(shiftx)) end = T
#       #if(!end & (shiftx > -LORth & shiftx < LORth)) end = T
#       if(!end & shiftx>0) hap = rbind(hap, c(hap[nrow(hap),c(1,2)], x2[3]))
#       if(!end & shiftx < 0) hap = rbind(hap, c(hap[nrow(hap),c(2,1)], x2[3]))
#       if(!end) i=i+1
#     }
#     
#     hap = cbind(start:i, hap)
#     colnames(hap) = c("wsnp", "g1", "g2", "shift")
#     haplotype = c(haplotype, list(hap))
#   }
#   #haplotype2=haplotype
#   #haplotype=haplotype2
#   i = 1
#   while(i<length(haplotype)){
#     h1 = haplotype[[i]]
#     h2 = haplotype[[i+1]]
#     hp = hapPhase(hap1=h1, hap2=h2, LORth = 2)
#     if(!is.na(hp)){
#       
#       if(hp == -1){h2 = matrix(h2[,c(1,3,2,4)], ncol=4);colnames(h2)=colnames(h1)}
#       
#       haplotype[[i]] = rbind(h1,h2)
#       haplotype = haplotype[-(i+1)]
#       
#       i=ifelse(i-1 >= 1, i-1, 1)
#       
#     }else{
#       i=i+1
#     }
#   }
#   
#   print(length(haplotype))
#   
#   if(length(haplotype)>1){
#     # At this stages, we trash "single nucleotide haplotype"
#     # Not informative about phase and/or genotyping error/missing genotypes
#     haplotype = haplotype[which(unlist(lapply(haplotype, function(hap){nrow(hap)>1})))]
#     hapcomb =  expand.grid(1:length(haplotype), 1:length(haplotype))
#     hapcomb = hapcomb[hapcomb$Var2 > hapcomb$Var1,]
#     hapcomb = hapcomb[hapcomb$Var1 - hapcomb$Var2 <= 5,]
#     
#     hapcomb = do.call(rbind, lapply(1:nrow(hapcomb), function(i){
#       h1 = haplotype[[ hapcomb[i,1] ]]
#       h2 = haplotype[[ hapcomb[i,2] ]]
#       phase = hapPhase(hap1=h1, hap2=h2)
#       gap = h2[1,"wsnp"] - h1[nrow(h1),"wsnp"]
#       gdis = info$cM[h2[1,"wsnp"]] - info$cM[h1[nrow(h1),"wsnp"]]
#       data.frame(hap1= hapcomb[i,1], hap2 = hapcomb[i,2], phase = phase, snpdis = gap, gdis = gdis)
#     }))
#     
#     hapcomb = subset(hapcomb, gdis < 4 & !is.na(phase))
#     hapcomb = hapcomb[order(hapcomb$gdis), ]
#     
#     happhase = rep(NA, length(haplotype))
#     
#     if(nrow(hapcomb)>0){
#       for(i in 1:nrow(hapcomb)){
#         hcx = hapcomb[i,]
#         h1 = hcx[,1]
#         h2= hcx[,2]
#         if( sum(is.na(happhase[c(h1,h2)])) == 2 ) happhase[c(h1,h2)] = c(1, hcx$phase)
#         if( sum(is.na(happhase[c(h1,h2)])) == 1 ){
#           isna = is.na(happhase[c(h1,h2)])
#           happhase[c(h1,h2)[isna]] = happhase[c(h1,h2)[!isna]]*hcx$phase
#         }
#       }
#     }
#     
#     
#     happhase[is.na(happhase)] = 1
#     
#   }else{happhase=1}
#   
#   
#   
#   haplotype = do.call(rbind, lapply(1:length(haplotype), function(i){
#     
#     phase = happhase[i]
#     hap = haplotype[[i]][,1:3]
#     colnames(hap) = NULL
#     
#     if(phase == -1){hap = hap[,c(1,3,2)]}
#     
#     return(hap)
#     
#   }))
#   
#   colnames(haplotype) = c("wsnp", "g1", "g2")
#   
#   snpmatch = which(KEEP)[haplotype[,1]]
#   
#   founder1all[snpmatch,] = haplotype[,2:3]
#   founder1all = as.data.frame(founder1all)
#   founder1all$tag = "unphased"
#   founder1all$tag[snpmatch] = "phased"
#   founder1all$tag[founder1all[,1] == founder1all[,2]] = "phased" #consider phased if homozygous
#   colnames(founder1all) = c("g1", "g2", "tag")
#   founder1all = cbind(infoall, founder1all)
#   return(founder1all)
#   
# }
# 
# 
# FMphased = inferPhaseFromRilsLD(founder1=founders$FM, founder2=founders$FA, rils=rils, info=info, LORth = 3)
# FAphased = inferPhaseFromRilsLD(founder1 = founders$FA, founder2 = founders$FM, rils=rils, info=info, LORth = 3)
# 
# foundergt = cbind(FMphased[,c("g1","g2")], FAphased[,c("g1","g2")])
# foundergt = foundergt[FAphased$tag == "phased" & FMphased$tag == "phased",]
# rils = rils[FAphased$tag == "phased" & FMphased$tag == "phased",]
# info = info[FAphased$tag == "phased" & FMphased$tag == "phased",]
# 
# 
# 
# 
# win = 11051:11101
# #Breaks between adjacent snp
# breakAdjSNP = countBreaks(windows=win, rils, foundergt, info)
# breakAdjSNP  = breakAdjSNP[breakAdjSNP$nbreak>3,]
# 
# BAD = NULL
# for(b in 1:nrow(breakAdjSNP)){
#   
#   i = breakAdjSNP[b,"start"]
#   
#   fi = foundergt[c(i,i+1),]
#   multiphase = sum(apply(fi[,1:2], 1, sum)==1)==2 | sum(apply(fi[,3:4], 1, sum)==1)==2
#   if(multiphase) print(i)
#   
#   ix0 = c((i-5):(i+5)) # around i
#   ix1 = c((i-5):(i-1),(i+1):(i+5)) # around i without i
#   ix2 =  c((i-5):(i),(i+2):(i+5)) # around i without i+1
#   ix3 =   c((i-5):(i-1),(i+2):(i+5)) # around i without i & i+1
#   
#   nbreaks_base = countBreaks(windows=c(1,length(ix0)), rils[ix0,], foundergt[ix0,], info[ix0,])$nbreak
#   nbreaksWithout_ij = countBreaks(windows=c(1,length(ix3)), rils[ix3,], foundergt[ix3,], info[ix3,])$nbreak
#   
#   bad = NULL
#   if(nbreaks_base > nbreaksWithout_ij){
#     nbreaksWithout_i = countBreaks(windows=c(1,length(ix1)), rils[ix1,], foundergt[ix1,], info[ix1,])$nbreak
#     nbreaksWithout_j = countBreaks(windows=c(1,length(ix2)), rils[ix2,], foundergt[ix2,], info[ix2,])$nbreak
#     bad = i + c(0,1)[nbreaks_base > c(nbreaksWithout_i,nbreaksWithout_j)]
#     if(length(bad)==0) bad = c(i,i+1)
#   }
#   
#   BAD = c(BAD,bad)
#   
# }
# 
# BAD = unique(BAD)
# 
# win2 = win[-which(win %in% BAD)]
# 
# win2 =5452:5501
# #win2 = c(2145:2150, win2)
# countBreaks(windows=c(1,length(win2)), rils[win2,], foundergt[win2,], info[win2,])
# 
# test = cor(t(foundergt[win2,]),use="pairwise.complete.obs")
# test2 = cor(t(rils[win2,]),use="pairwise.complete.obs")
# r_ratio = cor(t(rils[win2,]),use="pairwise.complete.obs")/cor(t(foundergt[win2,]),use="pairwise.complete.obs")
# colnames(r_ratio) = rownames(r_ratio) = win2
# colnames(test) = rownames(test) = win2
# colnames(test2) = rownames(test2) = win2
# 
# r_shift = sign(r_ratio+1)
# nshift = apply(r_shift, 1, function(x){sum(x==-1,na.rm=T)})
# 
# while(max(nshift,na.rm=T)>0){
#   i = which.max(nshift)
#   r_shift = r_shift[-i,-i]
#   nshift = apply(r_shift, 1, function(x){sum(x==-1,na.rm=T)})
# }
# 
# win3 = as.numeric(rownames(r_shift))
# 
# countBreaks(windows=c(1,length(win3)), rils[win3,], foundergt[win3,], info[win3,])$nbreak
# 
# countBreaks(windows=1:length(win3), rils[win3,], foundergt[win3,], info[win3,])
# 
# bad = info[win2[!(win2 %in% win3)],]
# 
# for(b in 1:nrow(breakAdjSNP)){
#   
#   r_shift2 = r_shift[-i,-i]
#   nshift2 = apply(r_shift2, 1, function(x){sum(x==-1)})
#   
#   ix0 = c((i-5):(i+5)) # around i
#   ix1 = c((i-5):(i-1),(i+1):(i+5)) # around i without i
#   ix2 =  c((i-5):(i),(i+2):(i+5)) # around i without i+1
#   ix3 =   c((i-5):(i-1),(i+2):(i+5)) # around i without i & i+1
#   
#   nbreaks_base = countBreaks(windows=c(1,length(ix0)), rils[ix0,], foundergt[ix0,], info[ix0,])$nbreak
#   nbreaksWithout_ij = countBreaks(windows=c(1,length(ix3)), rils[ix3,], foundergt[ix3,], info[ix3,])$nbreak
#   
#   bad = NULL
#   if(nbreaks_base > nbreaksWithout_ij){
#     nbreaksWithout_i = countBreaks(windows=c(1,length(ix1)), rils[ix1,], foundergt[ix1,], info[ix1,])$nbreak
#     nbreaksWithout_j = countBreaks(windows=c(1,length(ix2)), rils[ix2,], foundergt[ix2,], info[ix2,])$nbreak
#     bad = i + c(0,1)[nbreaks_base > c(nbreaksWithout_i,nbreaksWithout_j)]
#     if(length(bad)==0) bad = c(i,i+1)
#   }
#   
#   BAD = c(BAD,bad)
#   
# }
# 
# 
# 
# win[-which(ix0==bad)]
# 
# 
# i=6894
# ix = i:(i+1)
# #View(rils[ix,])
# foundergt[ix,]
# 
# info[ix,]
# 
# 
# View(rils[(i-2):(i+2),])
# foundergt[(i-2):(i+2),]
# 
# cor(t(rils[ix,]), use = "pairwise.complete.obs")
# 
# 
# 
# #i=6750
# 
# ix0 = c((i-5):(i+5))
# ix1 = c((i-5):(i-1),(i+1):(i+5))
# ix2 =  c((i-5):(i),(i+2):(i+5))
# 
# ix3 =   c((i-5):(i-1),(i+2):(i+5))
# 
# v = rils[ix0,]
# rownames(v)=ix0
# foundergt[ix0,]
# View(v)
# cor(t(v),use="pairwise.complete.obs")/cor(t(foundergt[ix0,]),use="pairwise.complete.obs")
# 
# nbreaks_base = countBreaks(windows=c(1,length(ix0)), rils[ix0,], foundergt[ix0,], info[ix0,])$nbreak
# nbreaks_WithoutSNP1 = countBreaks(windows=c(1,length(ix1)), rils[ix1,], foundergt[ix1,], info[ix1,])$nbreak
# nbreaks_WithoutSNP2 = countBreaks(windows=c(1,length(ix2)), rils[ix2,], foundergt[ix2,], info[ix2,])$nbreak
# countBreaks(windows=c(1,length(ix3)), rils[ix3,], foundergt[ix3,], info[ix3,])$nbreak
# countBreaks(windows=1:length(ix0), rils[ix0,], foundergt[ix0,], info[ix0,])
# 
# countBreaks(windows=1:2, rils[c(i-1,i+2),], foundergt[c(i-1,i+2),], info[c(i-1,i+2),])
# countBreaks(windows=1:2, rils[c(i-1,i+1),], foundergt[c(i-1,i+1),], info[c(i-1,i+1),])
# countBreaks(windows=1:2, rils[c(i,i+2),], foundergt[c(i,i+2),], info[c(i,i+2),])
# 
# countBreaks(windows=1:length(ix2), rils[ix2,], foundergt[ix2,], info[ix2,])
# 
# #######################################################################
# #######################################################################
# 
# 
# win = 1:51
# 
# fgt = foundergt[win,]
# founderhap = t(unique(t(fgt)))
# gx = rils[win,]
# 
# grouprils = cutree(hclust(dist(t(gx))),h=0)
# npergroup = table(grouprils)
# npergroup = npergroup[npergroup>10]
# npergroup = sort(npergroup, decreasing = T)
# if(length(npergroup)>4) npergroup=npergroup[1:4]
# 
# rilshap = gx[,sapply(as.numeric(names(npergroup)), function(x){which(grouprils==x)[1]})]
# 
# countBreaks(windows=c(1,length(win)), gx, fgt, info[win,])$nbreak
# sum(!(grouprils %in% (as.numeric(names(npergroup)))))
# 
# matchhap = unlist(lapply(1:ncol(founderhap), function(i){
#   x=founderhap[,i]
#   which(apply(abs(x-rilshap), 2, sum, na.rm=T) == 0)
# }))
# 
# if(sum(sort(matchhap) != 1:length(founderhap))!=length(founderhap)){
#   
#   matchhap2 = do.call(rbind,lapply(1:ncol(rilshap), function(i){
#     x=rilshap[,i]
#     x=apply(abs(x-fgt), 2, sum, na.rm=T)
#     c(x[order(x)],order(x))
#   }))
#   
#   matchhap2=do.call(rbind,lapply(1:ncol(fgt), function(i){
#     x=fgt[,i]
#     x=apply(abs(x-rilshap), 2, sum, na.rm=T)
#     c(x[order(x)],order(x))
#   }))
#   
#   ndiff = matchhap2[,1:ncol(rilshap)]
#   matchhap2 =  matchhap2[,(ncol(rilshap)+1):(2*ncol(rilshap))]
#   matchhap2=matchhap2[,1]
#   if(length(unique(matchhap2))<ncol(rilshap)) print("Two distinct hap with equal ndiff")
#   if(max(ndiff[1])>10) print("lot of unmatch")
#   #matchhap2 = order(matchhap2)
#   fgt = rilshap[,matchhap2]
#   
# }
# 
# 
# countBreaks(windows=c(1,length(win)), gx, fgt, info[win,])$nbreak
# sum(!(grouprils %in% (as.numeric(names(npergroup)))))
# 
# 
# 
# 
# ############
# 
# 
# 
# fx
# View(rils[win,])
# 
# cbind(fx, rils[win,1])
# 
# countBreaks(windows=c(1,length(win)), rils[win,], test[win,], info[win,])
# 
# 
# ########################################################################
# ########################################################################
# 
# #win = 1:1000
# #win = (1:nrow(haplotype))[1:nrow(haplotype) %% 10 == 0]
# 
# 
# gpos = info$cM[win]
# #gpos_sub = seq(min(gpos), max(gpos), (max(gpos)-min(gpos))/2000)
# #win = sapply(gpos_sub, function(x){ which.min(abs(gpos-x))})
# #win = win[!duplicated(win)]
# 
# #gpos[win]
# 
# gx = rils[win, ]
# #gx = abs(gx-foundergt[win,2])
# rx = cor(t(gx), use = "pairwise.complete.obs")/cor(t(foundergt[win,]), use = "pairwise.complete.obs")
# rx[lower.tri(rx, diag=T)] = NA
# colnames(rx) = rownames(rx) = win
# rx = reshape2::melt(rx)
# rx = rx[!is.na(rx$value),]
# rx$gpos1 = info$cM[rx$Var1]
# rx$gpos2 = info$cM[rx$Var2]
# 
# ggplot(rx, aes(x=gpos1, y=gpos2, color=value))+
#   theme_minimal()+coord_cartesian(expand=0)+
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
#   geom_point(shape=15, size = 6)+
#   scale_color_gradient2(low="red", high = "blue", mid="white", midpoint = 0)
# 
# 
# rx2 = rx[rx$gpos2>39.5 & rx$gpos2<40.5 & rx$gpos1>39.5 & rx$gpos1<40.5,]
# ggplot(rx2, aes(x=gpos1, y=gpos2, color=value))+
#   theme_minimal()+coord_cartesian(expand=0)+
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
#   geom_point(shape=15, size = 5)+
#   scale_color_gradient2(low="red", high = "blue", mid="white", midpoint = 0)
# 
# 
# #haplotype2=haplotype
# #I:10835139:C-T I     10835139  40.1
# nsnp.win = 10
# i = 12515
# win = haplotype[(i-nsnp.win):(i+nsnp.win), 1]
# gx = gx = rils[ haplotype[win, 1], ]
# gx = abs(gx-haplotype[win,3])
# rownames(gx)=win
# rx = cor(t(gx), use = "pairwise.complete.obs")
# 
# sign(rx[nsnp.win+c(-1,1)])
# rx[1:nsnp.win, 1:nsnp.win]
# 
# 
# 
# cbind(win, abs(founder2[win,]-haplotype[win,3]))
# 
#   
# test = rbind(haplotype[[1]],haplotype[[2]])
# 
# 
# View(rils[349:360,])
# 
# #-21.30448 -21.26806 -21.26804   0.00000  20.50221  20.49523  20.49522  20.01178  20.50221  20.49523
# # 20.49523  20.01180
# 
# 
# cbind(wsnp, lorCis)
# 
# matrix(lorCis, nrow=length(wsnp1))
# 
# 
# phaseshift(i=374,j=375, rils=rils, founder2=founder2)
# 
# View(rils[c(wsnp1,wsnp2),])
# 
# 
# 
# 
# 
# 
# 
# snps[c(11,16),]
# View(rils[8:16,])
# 
# rils[c(445,447),]
# cor(rils[10,], rils[12,], use = "pairwise.complete.obs")
# 
# 
# 
# 
# 
# test = try(phaseshift(i=269,j=271, rils=rils, founder2=founder2))
# grepl("Fixed allele",test[1])
# phaseshift = function(i,j, rils, founder2){
#   
#   snp1 = rils[i,]
#   snp2 = rils[j,]
#   
#   D = calculateD(snp1, snp2)
#   
#   # Expected linkage depending the possible phase in founder 1 & 2
#   gdis = diff(info$cM[c(i,j)])*0.01
#   expectedLDdecrease = (1-gdis)^30
#   
#   # Founder 2 genotype:
#   f2 = founder2[c(i,j),]
#   #f2
#   #If founder 2 is double heterozygous two possible phase (cis or trans), else only one (equivalent to founder 2 genotype)
#   if(sum(f2[,1]!=f2[,2])==2){
#     
#     founder2_PossiblePhase = list(cis =  matrix(c(0,0,1,1),ncol=2),
#                                   trans = matrix(c(1,0,0,1),ncol=2))
#     
#   }else{
#     founder2_PossiblePhase = list(f2)
#   }
#   
#   
#   whichphase = do.call(rbind, lapply(founder2_PossiblePhase, function(f2phase){
#     #f2phase = founder2_PossiblePhase[[1]]
#     
#     # Founder 1 always heterozygous as we focus on heterozygous snp in founder 1
#     # Phase cis:
#     pcis = matrix(c(0,0,1,1),ncol=2)
#     
#     ## Phase trans:
#     ptrans = matrix(c(1,0,0,1),ncol=2)
#     
#     #Genotype at founder crosses stages for different possible phase
#     gcis = cbind(pcis, f2phase)
#     gtrans = cbind(ptrans, f2phase)
#     
#     #rcis = cor(t(gcis))[1,2]
#     #rtrans = cor(t(gtrans))[1,2]
#     
#     Dcis = calculateD(gcis[1,],gcis[2,])
#     Dtrans = calculateD(gtrans[1,],gtrans[2,])
#     
#     testCis = unlist(c(test.D(snp1, snp2, theoDAB=expectedLDdecrease*Dcis)))
#     testTrans = unlist(c(test.D(snp1, snp2, theoDAB=expectedLDdecrease*Dtrans)))
#     
#     
#     OUT = data.frame(D=D,
#                Dcis = Dcis,
#                Dtrans=Dtrans,
#                probObsCis = ifelse(testCis[1]>1e-10,testCis[1], 1e-10),
#                probObsTrans = ifelse(testTrans[1]>1e-10,testTrans[1], 1e-10),
#                pvalDiffCis=testCis[2],
#                pvalDiffTrans=testTrans[2])
#     
#     colnames(OUT)= c("D", "DTheoCis","DTheoTrans","probObsCis","probObsTrans","pvalDiffCis","pvalDiffTrans")
#     
#     OUT$LORcis = log(OUT$probObsCis/OUT$probObsTrans)
#                      
#     OUT
#   }))
#   
#   whichphase = whichphase[which.max(apply(whichphase[,c("probObsCis","probObsTrans")], 1, max)),]
# 
#   whichphase$LORcis
# }
# 
# 
# 
# 
# #founder 2 = Hom
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# hap1 = haplotype[[7]]
# hap2 = haplotype[[9]]
# 
# haplotype[[8]]
# 
# rils[72,]
# 
# test = cor(t(rils[c(hap1[,"wsnp"],hap2[,"wsnp"]),]), use =  "pairwise.complete.obs")
# 
# test = do.call(rbind, lapply(haplotype, function(x){x}))
# 
# 
# 
# founders2 = founders[heterozygous,]
# 
# 
# 
# 
# 
# 
# # i=87
# # 
# # while(i<(nsnp-2)){
# #   
# #   start=i+1
# #   hap = c(0,1,NA)
# #   
# #   i = start
# #   end = F
# #   
# #   while(!end){
# #     snp1 = rils[i,]
# #     snp2 = rils[i+1,]
# #     notmissing = !is.na(snp1) & !is.na(snp2)
# #     r = cor(snp1[notmissing],snp2[notmissing])
# #     if(is.na(r)) end = T
# #     if(!end & (r<th & r> -th)) end = T
# #     if(!end & r > th) hap = rbind(hap, c(hap[nrow(hap),c(1,2)],r))
# #     if(!end & r < -th) hap = rbind(hap, c(hap[nrow(hap),c(2,1)],r))
# #     
# #     if(!end) i=i+1
# #   }
# #   
# #   hap = cbind(start:i, hap)
# #   colnames(hap) = c("wsnp", "g1", "g2", "r_withPrevious")
# #   haplotype = c(haplotype, list(hap))
# # }
# 
# 
# 
# 
# calculateD = function(snp1, snp2){
#   
#   notmissing = !is.na(snp1) & !is.na(snp2)
#   PAB = sum(snp1[notmissing] == 1 & snp2[notmissing] == 1)/sum(notmissing)
#   PA = sum(snp1[notmissing] == 1)/sum(notmissing)
#   PB = sum(snp2[notmissing] == 1)/sum(notmissing)
#   if(PA %in% c(0,1) | PB %in% c(0,1)) stop("Fixed allele")
#   D = PAB - PA*PB
#   
#   D
#   
# }
# 
# test.D = function(snp1, snp2, theoDAB){
#   notmissing = !is.na(snp1) & !is.na(snp2)
#   n = sum(notmissing)
#   PA = sum(snp1[notmissing] == 1)/n
#   PB = sum(snp2[notmissing] == 1)/n
#   if(PA %in% c(0,1) | PB %in% c(0,1)) stop("Fixed allele")
#   theo.PAB = theoDAB + PA*PB # Not exact because allele frequencies change between founders and rils
#   if(theo.PAB<0.01) theo.PAB=0.01
#   
#   obs.nAB =  sum(snp1[notmissing] == 1 & snp2[notmissing] == 1)
#   
#   #P-value test if phase is DIFFERENT FROM the observed PAB count
#   # => The "good" phase have a high P-value
#   prob = dbinom(obs.nAB,size=n,prob=theo.PAB)
#   pval = binom.test(x=obs.nAB,n=n,p=theo.PAB)$p.value
#   pval = as.numeric(pval)
#   data.frame(probObs=prob, pvalDiffFromTheo = pval)
# }
# 
# #phase.likelihood(i=9829,j=9836, rils, founder2, info)
# 
# phase.likelihood = function(i,j, rils, founder2, info){
#   
#   snp1 = rils[i,]
#   snp2 = rils[j,]
#   
#   D = calculateD(snp1, snp2)
#   
#   # Expected linkage depending the possible phase in founder 1 & 2
#   gdis = diff(info$cM[c(i,j)])*0.01
#   expectedLDdecrease = (1-gdis)^10
#   
#   #nmaxBreaks = sum(dbinom(x=0:10,size=10,prob=gdis)[0:10 %% 2 == 1])*ncol(rils)
#   
#   
#   # Founder 2 genotype:
#   f2 = founder2[c(i,j),]
#   #f2
#   #If founder 2 is double heterozygous two possible phase (cis or trans), else only one (equivalent to founder 2 genotype)
#   if(sum(f2[,1]!=f2[,2])==2){
#     
#     founder2_PossiblePhase = list(cis =  matrix(c(0,0,1,1),ncol=2),
#                                   trans = matrix(c(1,0,0,1),ncol=2))
#     
#   }else{
#     founder2_PossiblePhase = list(f2)
#   }
#   
#   
#   stat = do.call(rbind, lapply(founder2_PossiblePhase, function(f2phase){
#     #f2phase = founder2_PossiblePhase[[2]]
#     
#     # Founder 1 always heterozygous as we focus on heterozygous snp in founder 1
#     # Phase cis:
#     pcis = matrix(c(0,0,1,1),ncol=2)
#     
#     ## Phase trans:
#     ptrans = matrix(c(1,0,0,1),ncol=2)
#     
#     #Genotype at founder crosses stages for different possible phase
#     gcis = cbind(pcis, f2phase)
#     gtrans = cbind(ptrans, f2phase)
#     
#     #breaksCis = sum(apply(rils[c(i,j),], 2, function(x){ min(apply(gcis, 2, function(y) sum(x==y))) }),na.rm=T)
#     #breaksTrans = sum(apply(rils[c(i,j),], 2, function(x){ min(apply(gtrans, 2, function(y) sum(x==y))) }),na.rm=T)
#     
#     #rcis = cor(t(gcis))[1,2]
#     #rtrans = cor(t(gtrans))[1,2]
#     
#     Dcis = calculateD(gcis[1,],gcis[2,])
#     Dtrans = calculateD(gtrans[1,],gtrans[2,])
#     
#     testCis = unlist(c(test.D(snp1, snp2, theoDAB=expectedLDdecrease*Dcis)))
#     testTrans = unlist(c(test.D(snp1, snp2, theoDAB=expectedLDdecrease*Dtrans)))
#     
#     
#     OUT = data.frame(D=D,
#                      Dcis = Dcis*expectedLDdecrease,
#                      Dtrans=Dtrans*expectedLDdecrease,
#                      probObsCis = ifelse(testCis[1]>1e-10,testCis[1], 1e-10),
#                      probObsTrans = ifelse(testTrans[1]>1e-10,testTrans[1], 1e-10),
#                      pvalDiffCis=testCis[2],
#                      pvalDiffTrans=testTrans[2])
#     #nmaxBreaks.theo = nmaxBreaks,
#     #breaksCis=breaksCis,
#     #breakstrans=breaksTrans)
#     
#     colnames(OUT)= c("D", "DTheoCis","DTheoTrans","probObsCis","probObsTrans","pvalDiffCis","pvalDiffTrans")
#     
#     OUT$LORcis = log(OUT$probObsCis/OUT$probObsTrans)
#     
#     OUT
#   }))
#   
#   stat = stat[which.max(apply(stat[,c("probObsCis","probObsTrans")], 1, max)),]
#   
#   return(stat)
# }
# 
# 
# #phase.likelihood(i=12209,j=12217, rils=rils, founder2=founder2, info=info)
# 
# 
# #phase.likelihood(i=12212,j=12217, rils=rils, founder2=founder2, info=info)
# 
# phaseshift = function(i,j, rils, founder2, info, LORth=3, Dth=0.02){
#   plh = phase.likelihood(i=i,j=j, rils=rils, founder2=founder2, info=info)
#   shift = plh$LORcis
#   shift = ifelse(abs(shift)<LORth, NA, sign(shift))
#   #if(plh$pvalDiffCis < 1e-4 &  plh$pvalDiffTrans < 1e-4 ) shift = NA
#   if(abs(plh$D) < Dth ) shift = NA
#   return(shift)
# }
# 
# 
# #hap1 = haplotype[[1]]
# #hap2 = haplotype[[2]]
# 
# #hapPhase(hap1=hap1, hap2=hap2, nmax=7, LORth=2)
# 
# hapPhase = function(hap1, hap2, nmax=5, LORth=2){
#   
#   wsnp1 = hap1[,"wsnp"]
#   if(length(wsnp1)>nmax) wsnp1 = wsnp1[(length(wsnp1)-nmax+1):length(wsnp1)]
#   
#   wsnp2 = hap2[,"wsnp"]
#   if(length(wsnp2)>nmax) wsnp2 = wsnp2[1:nmax] 
#   
#   gx = rbind(hap1[hap1[,"wsnp"] %in% wsnp1,], hap2[hap2[,"wsnp"] %in% wsnp2,])
#   
#   wsnp = expand.grid(wsnp1,wsnp2)
#   
#   # 1 indcate that the snp are more likely in the same phase
#   # -1 indicate that the snp are more likely in opposite phase
#   phase = unlist(lapply(1:nrow(wsnp), function(i){
#     #out = try(phaseshift(i=which(rownames(gx)==wsnp[i,1]),
#     #                     j=which(rownames(gx)==wsnp[i,2]),
#     #                     rils=gx,
#     #                     founder2=founder2[c(wsnp1,wsnp2),],
#     #                     info = info[c(wsnp1,wsnp2),]))
#     
#     out = try(phaseshift(i=wsnp[i,1],j=wsnp[i,2], rils=rils, founder2=founder2, info = info, LORth = LORth), silent = T)
#     if(class(out)=='try-error') out = NA
#     
#     pairphase = ifelse(gx[gx[,"wsnp"]== wsnp[i,1], "g2"] ==  gx[gx[,"wsnp"]== wsnp[i,2], "g2"], 1,-1)
#     out = out*pairphase
#     out
#   }))
#   
#   phase[phase==0] = NA
#   phase = mean(phase, na.rm=T)
#   phase[abs(phase)!=1]=NA
#   phase[is.nan(phase)]=NA
#   return(phase)
#   
# }






# phaseHaplotypes = function(founderhaplotypes){
#   phases = do.call(rbind, lapply(1:(length(founderhaplotypes)-1), function(i){
#     #print(i)
#     fi = founderhaplotypes[[i]]$foundergt
#     fj = founderhaplotypes[[i+1]]$foundergt
#     wini = founderhaplotypes[[i]]$win
#     winj = founderhaplotypes[[i+1]]$win
#     midwin = c(wini[floor(length(wini)/2):length(wini)], winj[1:floor(length(winj)/2)])
#     winall = c(wini, winj)
#     winall = c(winall, max(winall)+1) # Just some add on to compensate a dumb thing that I am to lazy to correct nicely right now
#     
#     potentialShift = setNames(expand.grid(c(1,-1), c(1,-1)), c("shiftF1", "shiftF2"))
#     
#     nbreaks = unlist(lapply(1:nrow(potentialShift), function(i){
#       
#       if(potentialShift[i,1]==1){shiftF1=c(1,2)}else{shiftF1=c(2,1)}
#       if(potentialShift[i,2]==1){shiftF2=c(3,4)}else{shiftF2=c(4,3)}
#       shift = c(shiftF1,shiftF2)
#       
#       fused = rbind(fi[wini %in% midwin,], fj[winj %in% midwin,shift])
#       countBreaks(windows=c(1,nrow(fused)), rils=rils[midwin,], foundergt=fused, info=info[midwin,])$nbreak
#     }))
#     
#     shift = as.data.frame(potentialShift[which.min(nbreaks),])
#     shift = c(apply(shift, 2, function(x){if(x==1){c(1,2)}else{c(2,1)}})) + c(0,0,2,2)
#     
#     haploall_phased = rbind(fi, fj[,shift])[1:nrow(haploall),]
#     nbreaksPhased = nbreaks[which.min(nbreaks)]
#     
#     haploall_infered = InferFoundersHaploBlocks(founder1=founder1[winall,],
#                              founder2=founder2[winall,],
#                              rils=rils[winall,],
#                              info=info[winall,],
#                              nsnpWin = length(winall)-1, 
#                              maxSizeCM = 1e6)[[1]] # Just an absuerdely big number so no window subdivision 
#     
#     nbreaksInfered = haploall_infered$nbreaks
#     
#     haploall_infered = haploall_infered$foundergt
#     
#     
#     hapdist = sum(abs(haploall - rbind(fi, fj[,shift])[1:nrow(haploall),]),na.rm=T)
#     
#     #countBreaks(windows=c(1,nrow(haploall)), rils=rils[winall,], foundergt=rbind(fi, fj[,shift])[1:nrow(haploall),], info=info[winall ,])$nbreak
#     #countBreaks(windows=c(1,nrow(haploall)), rils=rils[winall,], foundergt=haploall, info=info[winall ,])$nbreak
#     
#     #nbreaks2 = countBreaks(windows=c(1,nrow(haploall)), rils=rils[winall,], foundergt=haploall, info=info[winall ,])$nbreak
#     
#     #nbreaks2 = InferFoundersHaploBlocks(founder1=founder1[winall,],
#     #                                   founder2=founder2[winall,],
#     #                                   rils=rils[winall,],
#     #                                   info=info[winall,],
#      #                                  nsnpWin = length(winall)-1)[[1]]$nbreaks
#     
#     
#     out = cbind( as.data.frame(potentialShift[which.min(nbreaks),]),
#                  nbreaksPhased = nbreaks[which.min(nbreaks)],
#                  nbreaksHapInfered = nbreaks2,
#                  PhasedInfered_dift = hapdist)
#     
#     return(out) 
#   }))
#   
#   
#   
#     
#     
#   return(phases)
# }



#i=which(founders$ID=="I:2777336:A-G")

# wrongFounderGenotypes = sapply(1:(nrow(rils)-1), function(i){
#   
#   if(i  %% 1000 == 0) print(i)
#   
#   wrong = F
#   r = cor(rils[i,], rils[i+1,], use="pairwise.complete.obs")
#   if(is.na(r)) r = 0
#   nalt1 = sum(founder1[c(i,i+1),])
#   nalt2 = sum(founder2[c(i,i+1),])
#   nhom1 = sum(founder1[c(i,i+1),1] == founder1[c(i,i+1),2])
#   nhom2 = sum(founder2[c(i,i+1),1] == founder2[c(i,i+1),2])
#   
#   # r>0.99 = haplo 0;0 or 1;1 in rils (coupling)
#   # absent of founders in nalt in founders = 1
#   #Same reasoning for repulsion
#   repulsionImposs = (nalt1 %in% c(0,1) & nalt2 %in% c(0,1)) | (nalt1 == 2 & nhom1 ==2 & nalt2==2 & nhom2==2) 
#   couplingImposs = (nalt1 %in% c(0,4) & nalt2 %in% c(0,4))
#   nmiss = sum(is.na(rils[i,]) | is.na(rils[i+1,]))
#   
#   
#   if(r>0.99 & repulsionImposs  & nmiss < 20 ) wrong = T
#   if(r< -0.99 & couplingImposs & nmiss < 20) wrong = T
#   
#   return(wrong)
#   
#   })
# 
# wrongFounderGenotypes = unique(sort(c(which(wrongFounderGenotypes), which(wrongFounderGenotypes)+1)))
# badgeno = founders[wrongFounderGenotypes,]
# 
# rils = rils[-wrongFounderGenotypes,]
# founders = founders[-wrongFounderGenotypes,]
# info = info[-wrongFounderGenotypes,]
# 
# 
# steps = seq(1, nrow(rils), 25)
# winsize = 100
# rth = 0.9
# 
# steps = data.frame(steps=steps, p1=steps-(winsize/2),p2=steps+(winsize/2))
# steps$p1[steps$p1<1]=1
# steps$p2[steps$p2>nrow(rils)]=nrow(rils)
# #info$cM[steps$p2]-info$cM[steps$p1]
# 
# cmx = info$cM[info$ID == "I:12504710:C-T"]
# #which(info$cM[steps$p2]>cmx & info$cM[steps$p1]<cmx)
# 
# #872
# i=873
# 
# 
# #which(info$ID[wrong2] == "I:12504710:C-T")
# 
# wrong2 = lapply(1:nrow(steps), function(i){
#   if(i %% 100 == 0) print(i)
#   step = steps$step[i]
#   p1 = steps$p1[i]
#   p2 = steps$p2[i]
#   win = p1:p2
#   
#   if(length(win)<3) return(NULL)
#   
#   fx = cbind(founder1[win,],founder2[win,])
#   r_ratio = cor(t(rils[win,]),use="pairwise.complete.obs")/cor(t(fx),use="pairwise.complete.obs")
#   colnames(r_ratio) = rownames(r_ratio) = win
#   r_shift = sign(r_ratio+rth)
#   nshift = apply(r_shift, 1, function(x){sum(x==-1, na.rm=T)})
#   
#   while(max(nshift)>0 & ncol(r_shift)>2){
#     i = which.max(nshift)
#     r_shift = r_shift[-i,-i]
#     nshift = apply(r_shift, 1, function(x){sum(x==-1,na.rm=T)})
#   }
#   
#   
#   bad = win[!(win %in% as.numeric(rownames(r_shift)))]
#   if(length(bad)==0) bad = NULL
#   return(bad)
#   
# })
# 
# wrong2 = unique(unlist(wrong2))
# 
# 
# 
# rils = rils[-wrong2,]
# founders = founders[-wrong2,]
# info = info[-wrong2,]
# 
# #ggplot(info, aes(POS,cM))+geom_point()

