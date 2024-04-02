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
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)]],SPLIT[[length(SPLIT)-1]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}



runlsei = function(gx,predictors){
  # gx = vector with single ril genotype with 0/0.5/1 for ref hom / het / alt hom
  # predictors in the founders genotype matrix for the same window
  # The function return a vector of length = ncol(predictors) with founder proportion
  # Example:
  # c(0,0,0,1) => 100% match with the 4th founder
  # c(0,0.2,0.8,0) => likely mean that there is a part of this genomic window corresponding to founder 1 and another to founder 2
  # c(0,0.5,0.5,0) can also mean that founder 2 and 3 are both an exact match (founders are identicals)
  
  # First, SNP with at least missing values
  predictNotMissing = apply(cbind(predictors, gx),1,function(x) sum(is.na(x))==0) 
  if(sum(predictNotMissing)<2){return(rep(NA, ncol(predictors)))}
  predictors = as.matrix(predictors[predictNotMissing,])
  gx = gx[predictNotMissing]
  
  #If a founder is heterozygous at some loci, consider it is the genotype of the ril at this loci
  # at SNP 10, ril is 0 and founder A is 0.5, change founder A SNP 10 in 0
  # It assumes that the founder is truely heterozygous (no genotyping error) 
  # At that whatever is the ril genotype at the herozygous snp should thus be a match
  # Or that we just don't know and still assume a match
  # You may want to question this choice because heterozygous allele is supposedly homozygous founder can be genotyping error
  # So you may just want change to founders[founders==0.5]=NA
  # If you do that a SNP which is heterozygous is one founder will be also considered as SNP with one missing value a be eliminated (above)
  # So you would loose information
  # Anyway, if the typical haplotype blocks contains a lot of different snp between founders it should not matters much
  # (If some founder are only different from a couple of SNP, it likely matters)
  heterozygous = which(predictors==0.5)
  if(length(heterozygous)>0){
    nsnp = nrow(predictors)
    rowhet=sapply(heterozygous, function(h){
      #wcol = floor(h/nsnp)+1
      wrow = (h/nsnp - floor(h/nsnp))*nsnp
      if(wrow == 0) wrow = nsnp#; wcol = wcol-1
      as.integer(round(wrow,digits = 1))
    })
    
    predictors[heterozygous]=gx[rowhet]
  }
  
  predictors = round(predictors)
  predictors[predictors>1]=1 
  
  # First look if there is a perfect match with a founder
  # If yes, the function return c(0,1,0,0) (perfect match with founder 2)
  # Alternatively, it would return something like  c(0.333, 0.333, 0.333, 0) if there is perfect match with multiple founders
  matchgeno = sapply(1:ncol(predictors), function(founder){ sum(predictors[,founder] == gx)/nrow(predictors)})
  matchgeno =  round(matchgeno,digits = 7) == 1
  if(sum(matchgeno)>0){out=ifelse(matchgeno,1/sum(matchgeno), 0); return(out)}
  
  #If not, use lsei to solves a least squares problem under both equality and inequality constraints
  # Basically, we get a vector of founder frequency that would explain well SNP frequencies in the ril
  d = ncol(predictors)
  A = predictors
  E = t(matrix(rep(1,d)))
  Ff = 1
  G = diag(rep(1,d))
  H = matrix(rep(0.000,d))
  Y = gx
  
  out = lsei(A=A,B=Y,E=E,F=Ff,G=G,H=H,verbose=F)
  
  if(out$IsError){out = rep(NA,d)}else{
    out = out$X
    
    # I don't think this part of the code is useful but it does noy harm. I will think about it.
    whichfounder = which(out>0)
    whichpolym = ispolymorphicSNP(subfounders=matrix(predictors[, whichfounder], ncol=length(whichfounder)))
    whichpolym = whichpolym[!(whichpolym %in% which(is.na(gx)))]
    matchgeno = sapply(whichfounder, function(founder){ sum(predictors[whichpolym,founder] == gx[whichpolym])})
    matchgeno = matchgeno==length(whichpolym) & length(whichpolym)>0
    
    if(sum(matchgeno)>0){y=rep(F,ncol(predictors)); y[whichfounder]=matchgeno; out=ifelse(y,1/sum(y), 0)}
  }
  
  return(out)
  
}



pushleft = function(b1,b2,ril, founders, nsnp, rule = 1){
  # This function start from a window between b1 & b2
  # and increase the window range on the left (=decrease b1)
  # to check what how it affects the output of runlsei
  # And to stop once optimize the chosen parameter
  # rule 1 = maximize the max founder weight
  # rule 2 = minimize the number of match founder
  # Basically, this function is useful if we know that window b1=>b2 does match with founder X
  # and that we want to know the left boudary of this haplotype block
  # If at b1 - 10 the output of runlsei still indicate a perfect match with founder X, 
  # it indicates that the founder X haplotype block left boundaries can be pushed by -10
  # if the founder weights output by lsei decrease, it mean that we cannot push the boundary
  # i.e. the breakpoint is before
  
  b1push = b1
  minsize = nsnp + 1
  ni=0
  
  outlsei = runlsei(gx=ril[b1:b2],predictors = founders[b1:b2,])
  maxw = round(max(outlsei), digits = 3)
  nmatch = vtar = length(which(outlsei>0))
  
  if(rule==1){
    vtar = maxw
    condition = (maxw >= vtar & b1push > minsize)
  } 
  
  if(rule==2){
    vtar = nmatch 
    condition = (nmatch <= vtar & b1push > minsize)
  } 
  
  if(is.na(vtar)){stop()}
  
  while(condition==T){
    ni = ni+1
    b1push = b1push-nsnp
    win = b1push:b2
    predictors = founders[win,]
    gx = ril[win]
    
    outlsei = try(runlsei(gx=gx,predictors = predictors),silent = T)
    maxw = round(max(outlsei), digits = 3)
    nmatch = length(which(outlsei>0))
    
    if(rule == 1){
      if(is.na(maxw)) maxw = 0
      condition = (maxw >= vtar & b1push > minsize)
      vtar = maxw
    }
    
    if(rule==2){
      if(nmatch==0) nmatch = 1e6
      condition = (nmatch <= vtar & b1push > minsize)
      vtar = nmatch 
    } 
    
    
  }
  
  b1push = b1push+ifelse(ni>0, nsnp, 0)
  return(b1push)
}



pushright = function(b1,b2,ril, founders, nsnp, rule = 1){
  #rule 1 = maximaize the max weight
  #rule 2 = minimize the number of match founder
  # see pushleft function for details
  
  b2push = b2
  maxsize = length(ril)-nsnp
  ni=0
  
  outlsei = runlsei(gx=ril[b1:b2],predictors = founders[b1:b2,])
  maxw = round(max(outlsei), digits = 3)
  nmatch = vtar = length(which(outlsei>0))
  
  if(rule==1){
    vtar = maxw
    condition = (maxw >= vtar & b2push < maxsize)
  } 
  
  if(rule==2){
    vtar = nmatch 
    condition = (nmatch <= vtar & b2push < maxsize)
  } 
  
  if(is.na(vtar)){stop()}
  
  while(condition==T){
    ni = ni+1
    b2push = b2push+nsnp
    win = b1:b2push
    predictors = founders[win,]
    gx = ril[win]
    
    outlsei = try(runlsei(gx=gx,predictors = predictors),silent = T)
    maxw = round(max(outlsei), digits = 3)
    nmatch = length(which(outlsei>0))
    
    if(rule == 1){
      if(is.na(maxw)) maxw = 0
      condition = (maxw >= vtar & b2push < maxsize)
      vtar = maxw
    }
    
    if(rule==2){
      if(nmatch==0) nmatch = 1e6
      condition = (nmatch <= vtar & b2push < maxsize)
      vtar = nmatch 
    } 
    
    
  }
  
  b2push = b2push-ifelse(ni>0, nsnp, 0)
  return(b2push)
}


# pushright = function(b1,b2,ril, founders, nsnp){
#   
#   maxw = vtar = round(max(runlsei(gx=ril[b1:b2],predictors = founders[b1:b2,])), digits = 3) 
#   if(is.na(vtar)){stop()}
#   b2push = b2
#   
#   maxsize = length(ril)-nsnp
#   ni = 0
#   
#   while(maxw >= vtar & b2push < maxsize){
#     ni = ni+1
#     b2push = b2push+nsnp
#     win = b1:b2push
#     predictors = founders[win,]
#     gx = ril[win]
#     
#     maxw = try(max(runlsei(gx=gx,predictors = predictors)))
#     if(is.na(maxw)){maxw=0}
#     maxw = round(maxw, digits = 3)
#   }
#   
#   b2push = b2push-ifelse(ni>0, nsnp, 0)
#   return(b2push)
# }


ispolymorphicSNP = function(subfounders){
  #return where are the polymorphic snps among founders (generally subset of all founders)
  if(ncol(subfounders)>1){
    mx = apply(subfounders, 1, mean)
    dx = apply(subfounders,2,function(x){abs(x-mx)})
    whichpolym = which(apply(dx,1,sum)>0)
  }else{
    whichpolym=which(F)
  }
  
  return(whichpolym)
}




find.border = function(b1,b2,ril, founders,rule=1){
  # Starting from a initial window between b1 & b2 matching a founder
  # the function expand (push) the window to find the boundaries/breakpoints of the founder haplotype block
  # It use pushleft and pushright function (more comments in pushleft function)
  b1push=pushleft(b1=b1,b2=b2,ril=ril, founders=founders, nsnp=50, rule=rule)
  b1push=pushleft(b1=b1push,b2=b2,ril=ril, founders=founders, nsnp=10, rule=rule)
  b1=pushleft(b1=b1push,b2=b2,ril=ril, founders=founders, nsnp=1, rule=rule)
  
  
  #Search for right border
  b2push = pushright(b1=b1,b2=b2,ril=ril, founders=founders, nsnp=50, rule=rule)
  b2push = pushright(b1=b1,b2=b2push,ril=ril, founders=founders, nsnp=10, rule=rule)
  b2 = pushright(b1=b1,b2=b2push,ril=ril, founders=founders, nsnp=1, rule=rule)
  
  # Return the founder haplotype
  win = b1:b2
  predictors = founders[win,]
  gx = ril[win]
  
  wfounder = runlsei(gx=gx,predictors = predictors)
  
  return(cbind(data.frame(pos1 = b1, pos2=b2), matrix(wfounder, nrow=1)))
  
}

supressEncased = function(HAPLO){
  
  # Search for inferred haplotype encased in another
  # i.e. hap1 go from 1 to 150 & hap2 go from 10 to 140
  # and choose the larger haplo
  
  HAPLO = as.data.frame(HAPLO)
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  h = 1
  while(h < nrow(HAPLO)){
    
    pos1=HAPLO[h,1];pos2=HAPLO[h,2]
    
    #encase =  HAPLO[,1] >= pos1 & HAPLO[,2] <= pos2
    #if(sum(encase)>1)  HAPLO=HAPLO[-which(encase)[which(encase)!=h],]
    
    isEncased = HAPLO[,1] <= pos1 & HAPLO[,2] >= pos2
    if(sum(isEncased)>1){HAPLO=HAPLO[-h,]}else{h=h+1}
  }
  
  return(HAPLO)
}


splithaplo = function(b1,b2,ril, founders){
  
  # when the output of runlsei indicate that between b1 & b2 the haplotype is a combination of several founders
  # This particular function can be used to try to split the window to find the breakpoint between  the different founder haplotypes
  # i.e. within win to 1 to 100 there is a combination of founder A and B
  # the function can return that 1 to 70 is founder B and 71 to 100 is founder A
  
  # First we want to find which are the polymorphic snps between the possible founders
  win = b1:b2
  gx = ril[win]
  predictors = as.matrix(founders[win,])
  
  founderweights = runlsei(gx,predictors)
  whichfounder =  which(founderweights>0)
  
  heterozygous = which(predictors==0.5)
  if(length(heterozygous)>0){
    nsnp = nrow(predictors)
    rowhet=sapply(heterozygous, function(h){
      #wcol = floor(h/nsnp)+1
      wrow = (h/nsnp - floor(h/nsnp))*nsnp
      if(wrow == 0) wrow = nsnp#; wcol = wcol-1
      as.integer(round(wrow,digits = 1))
    })
    
    predictors[heterozygous]=gx[rowhet]
  }
  
  whichpolym = ispolymorphicSNP(subfounders=predictors[, whichfounder])
  whichpolym = whichpolym[!(whichpolym %in% which(is.na(gx)))]
  
  #View(cbind(gx, predictors[,whichfounder]))
  
  # Then, at the polymorphic snps, we look at the that corresponding founder
  # search for indication of breakpoints between the founders
  if(length(whichpolym)>0){
    
    breaks = c(0,whichpolym,length(win)+1)#+b1-1
    breaks = t(sapply(2:(length(breaks)-1), function(i){c(breaks[i-1]+1, breaks[i+1]-1)}))
    breaks = matrix(breaks, ncol=2)
    matchgeno = sapply(whichpolym, function(wsnp){ which(predictors[wsnp,whichfounder] == gx[wsnp])})
    matchgeno = matrix(matchgeno, nrow = length(whichpolym))
    matchgeno = t(sapply(whichpolym, function(wsnp){
      out = rep(NA, length(whichfounder))
      wx = which(predictors[wsnp,whichfounder] == gx[wsnp])
      out[wx]=1
      out[is.na(out)]=0
      out}))
    
    if(nrow(matchgeno)>1){
      matchgeno=apply(matchgeno, 2, function(x){
        wbreak = diff(c(0,which(x==0), length(x)+1))-1
        wbreak = wbreak[wbreak>0]
        wbreak = unlist(sapply(wbreak, function(wb){rep(wb,wb)}))
        x[x==1]=wbreak
        x
      })
    }
    
    
    matchgeno2 = apply(matchgeno, 1, which.max)
    
    tofuse = which(diff(matchgeno2)==0)
    if(length(tofuse)>0){
      for(tf in tofuse){breaks[tf+1,1] = breaks[tf,1]}
      breaks = breaks[-tofuse,]
    }
    
    breaks = matrix(breaks, ncol=2)
    #1:nrow(breaks)
    splitted = do.call(rbind,lapply(1:nrow(breaks), function(i){
      #print(i)
      start = breaks[i,1]
      end = breaks[i,2]
      outf = try(runlsei(gx=gx[start:end],predictors=predictors[start:end,]),silent = T)
      if(class(outf)=="try-error") outf = matrix(rep(NA, ncol(predictors)), nrow=1)
      
      x = cbind(data.frame(pos1 = start+b1-1, pos2=end+b1-1), matrix(outf, nrow=1))
    }))
    
    return(splitted)
    
  }else{
    return(cbind(data.frame(pos1 = b1, pos2=b2), matrix(founderweights, nrow=1)))
  }
  
  
  
}





haplosearch = function(ril, founders, snps){
  # Main function
  
  # 1) We look within many different window for correspondence with founders
  # This is done with runlsei out with return a vector of weight of lenght = ncol(founders)
  # 0,1,0,0 Perfect match with founder 2
  # several numbers > 0 & < 1, several founders possible or combination of founders
  # see runlsei for details
  winsizes = c(300, 200, 150, 125, 100, 75, 50)
  lsout = lapply(winsizes, function(wsize){
    #print(wsize)
    windows = get.win(nrow(founders), winsize=wsize, minsize=50)
    
    out = do.call(cbind, lapply(windows, function(win){
      #print(win[1])
      #which(unlist(lapply(windows, function(x){x[[1]]}))==12001)
      #win = unlist(windows[160])
      predictors = founders[win,]
      gx = ril[win]
      
      p = try(runlsei(gx=gx,predictors = predictors), silent = T)
      if(class(p)=="try error") p = rep(NA,ncol(predictors))
      p=matrix(rep(p,  length(win)), ncol=length(win), nrow=ncol(predictors), byrow = F)
      p
      
    }))
    
    out
  })
  
  #lseout = matrix with row correponding to different window size
  # columns correspoinding to snps
  # and value to the max founder weights 
  # Select the max value
  # = select the window size a each position that maximize the founder weights (ideally a perfect match to one founder = weight of 1)
  maxf = do.call(rbind, lapply(lsout, function(x){apply(x,2,max)}))
  maxf = unlist(apply(maxf, 2, function(x){ if(sum(is.na(x))<nrow(maxf)){which.max(x)}else{NA}}))
  
  
  haploweights = do.call(cbind, lapply(1:length(maxf), function(i){
    if(!is.na(maxf[i])){as.numeric(lsout[[ maxf[i] ]][,i])}else{rep(NA, ncol(founders))}
  }))
  
  haploweights = round(haploweights, digits = 3)

  
  HAPLO = NULL
  
  maxweight = apply(haploweights, 2, max)
  uniqmaxweight = sort(unique(maxweight[maxweight>0]),decreasing = T)
  
  while(max(uniqmaxweight)>0){
    
    vtar = max(uniqmaxweight,na.rm=T)
    #print(vtar)
    
    wvtar = which(maxweight==vtar)
    
    breakpoints = sort(c(wvtar[1], wvtar[which(diff(wvtar)>1)], wvtar[which(diff(wvtar)>1)+1], wvtar[length(wvtar)]))
    breakpoints = matrix(breakpoints, ncol=2, byrow=T)
    
    breakpoints = matrix(breakpoints[breakpoints[,2]-breakpoints[,1] > 0,], ncol=2)
    
    if(nrow(breakpoints)>0){
      
      breakpoints = do.call(rbind, lapply(1:nrow(breakpoints), function(bi){
        b1 = breakpoints[bi,1]
        b2 = breakpoints[bi,2]
        
        wfounder = apply(haploweights[,b1:b2], 2, which.max)
        n = length(unique(wfounder))
        
        if(n>1){
          newbreaks = b1 + which(diff(wfounder)!=0)
          #newbreaks = b1+cumsum(unlist(lapply(split(wfounder,wfounder), length)))
          newbreaks = c(b1,newbreaks,b2+1)
          newbreaks = cbind(newbreaks[1:(length(newbreaks)-1)], newbreaks[2:length(newbreaks)]-1)
          return(newbreaks)
        }else{
          return(breakpoints[bi,])
        }
        
      }))
      
      
      
      haplo=do.call(rbind,lapply(1:nrow(breakpoints), function(bi){
        
        #print(bi)
        b1 = breakpoints[bi,1]
        b2 = breakpoints[bi,2]
        
        #Search for left border
        
        x = try(find.border(b1=b1,b2=b2,ril=ril, founders=founders),silent = T)
        if(class(x)=="try-error") x = matrix(c(b1,b2,rep(NA, ncol(founders))), nrow=1)
        
        wfounder = which(x[,3:ncol(x)]>0)
        
        if(length(wfounder)>1){
          pos1 = x[,1]
          pos2 = x[,2]
          
          x=splithaplo(b1=pos1,b2=pos2,ril=ril, founders=founders)
        }
        colnames(x) = c("pos1", "pos2", 1:ncol(founders))
        return(x)
        
        
      }))
      
      
      for(bi in 1:nrow(haplo)){
        b1 = haplo[bi,1]
        b2 = haplo[bi,2]
        v = ifelse(is.na(haplo[bi,3]), NA, -1)
        fx = matrix(v, nrow = ncol(haplo)-2, ncol = length(b1:b2))
        haploweights[,b1:b2] = fx
      }
      
      colnames(haplo)=c("pos1", "pos2", 1:ncol(founders))
      HAPLO = rbind(HAPLO, haplo)
      
    }else{
      haploweights[,wvtar] = NA
    }
    
    
    maxweight = apply(haploweights, 2, max)
    uniqmaxweight = sort(unique(maxweight[maxweight>0]),decreasing = T)
    
  }
  
  
  
  
  HAPLO = HAPLO[!is.na(HAPLO[,3]),]
  
  HAPLO  = do.call(rbind, lapply(1:nrow(HAPLO), function(i){
    #print(i)
    pos1 = unlist(HAPLO[i,1])
    pos2 = unlist(HAPLO[i,2])
    
    x = find.border(b1=pos1,b2=pos2,ril=ril, founders=founders, rule=1)
    x = round(unlist(x), digit = 5)
    x = matrix(x, nrow=1)
    if(length(which(x[,3:ncol(x)]>0))>1){x=splithaplo(b1=x[,1],b2=x[,2],ril, founders)}
    
    colnames(x) = colnames(HAPLO)
    
    x
  }))
  
  HAPLO = HAPLO[!is.na(HAPLO[,3]),]
  
  
  HAPLO = supressEncased(HAPLO)
  
  ########################################################
  #### FUSE overlaping windows with common founders #####
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  i = 1
  while(i < nrow(HAPLO)){
    j=i+1
    #pos1=HAPLO[i,1];pos2=HAPLO[i,2]
    overlapwithnext = HAPLO[j,1] <= HAPLO[i,2]
    
    x1 = HAPLO[i,3:ncol(HAPLO)]
    x2 = HAPLO[j,3:ncol(HAPLO)]
    x1[x1>0] = round(1, digits=0)
    x2[x2>0] = round(1, digits=0)
    
    commonfounder = which((x1 + x2) > 1)
    
    if(length(commonfounder)>0 & overlapwithnext){
      b1 = HAPLO[i,1]
      b2 =HAPLO[j,2]
      x = round(unlist(runlsei(ril[b1:b2],founders[b1:b2,])), digits = 2)
      
      if(!is.na(x[1])){
        x = c(b1,b2,x)
        HAPLO[i,] = x
        HAPLO = HAPLO[-j,]
      }
    }
    
    i = i+1
    
    #}else{i = i+1}
  }
  
  HAPLO = HAPLO[!is.na(HAPLO[,3]),]
  
  
  ##################################################################
  #### FUSE close haplotypes Separated by only one SNP  ###########
  
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  fused = do.call(rbind, lapply(1:nrow(HAPLO), function(i){
    #print(i)
    whichfounder = which(HAPLO[i,3:ncol(HAPLO)]>0)
    founder = matrix(unlist(round(HAPLO[,whichfounder+2])), nrow = nrow(HAPLO))
    founder[founder>0.0] = 1
    nextcommon =  sort(unlist( apply(founder, 2, function(x){which(x>0)}) ))
    nextcommon = min(nextcommon[nextcommon >i])
    if(nextcommon==Inf){return(NULL)}
    
    j = nextcommon 
    
    overlapwithnext = HAPLO[j,1] <= HAPLO[i,2]
    
    whichfounder = HAPLO[i,3:ncol(HAPLO)]
    whichfounder2 = HAPLO[j,3:ncol(HAPLO)]
    whichfounder[whichfounder>0] = round(1, digits=0)
    whichfounder2[whichfounder2>0] = round(1, digits=0)
    
    commonfounder = which((whichfounder + whichfounder2) > 1)
    
    
    if(!overlapwithnext){
      #print(i)
      
      s1 = HAPLO[i,1]
      e1 = HAPLO[i,2]
      s2 =HAPLO[j,1]
      e2 =HAPLO[j,2]
      #print(s1)
      
      gap = e1:s2
      
      divergent = abs(founders[gap, commonfounder] - ril[gap])
      divergent = matrix(unlist(divergent), nrow = length(gap))
      divergent  = unique(apply(divergent, 2, function(x){which(x==1)}))
      
      gaprelsize = diff(snps$pos[gap[c(1,length(gap))]])/diff(snps$pos[c(b1,b2)])
      
      
      if(gaprelsize < 0.03 & length(divergent)<=1){
        win = c(s1:e1, s2:e2)
        
        new = runlsei(gx=ril[win],predictors=founders[win,])
        new = matrix(c(s1,e2,new), nrow=1)
        colnames(new)=c("pos1", "pos2", 1:ncol(founders))
        return(new)
      }else{return(NULL)}
      
    }else{return(NULL)}
    
  }))
  
  if(is.null(fused)){fused = HAPLO[which(F),]}
  
  ### Fuse the fused haplotype that overlap
  i = 1
  while(i < nrow(fused)){
    j=i+1
    #pos1=HAPLO[i,1];pos2=HAPLO[i,2]
    overlapwithnext = fused[j,1] <= fused[i,2]
    
    x1 = fused[i,3:ncol(fused)]
    x2 = fused[j,3:ncol(fused)]
    x1[x1>0] = round(1, digits=0)
    x2[x2>0] = round(1, digits=0)
    
    commonfounder = x1 == x2
    
    if(sum(!commonfounder)==0 & overlapwithnext){
      b1 = fused[i,1]
      b2 =fused[j,2]
      x = fused[i,3:ncol(fused)]
      x = c(b1,b2,x)
      fused[i,] = x
      fused = fused[-j,]
      fused = matrix(fused, ncol=ncol(HAPLO))
    }else{i = i+1}
  }
  
  # Bind with other haplotypes & supress overlaping
  colnames(fused) = colnames(HAPLO)
  HAPLO = rbind(HAPLO, fused)
  HAPLO=supressEncased(HAPLO)
  
  
  
  #########################################
  ##### When two haplotype at the same place, chose 1
  
  i = 2
  while(i<nrow(HAPLO)){
    
    b1  = HAPLO[i,1]
    b2 = HAPLO[i,2]
    
    overlap1 = min(which(HAPLO[,1] < b1 & HAPLO[,2] > b1))
    overlap2 = max(which(HAPLO[,1] < b2 & HAPLO[,2] > b2))
    
    sameplace = which(HAPLO[,1] > HAPLO[overlap1,1] & HAPLO[,2] < HAPLO[overlap2,2])
    sameplace = sameplace[sameplace>i]
    
    if(length(sameplace)>0){
      keep = which.max(HAPLO[c(i,sameplace), 2] -HAPLO[c(i,sameplace), 1])
      out = c(i,sameplace)[which(!(1:(length(sameplace)+1) %in% keep))]
      HAPLO = HAPLO[-out,]
    }
    #print(i)
    #print(nrow(HAPLO))
    i=i+1
    
    
  }
  
  
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  HAPLO[,1] = snps$pos[HAPLO[,1]]
  HAPLO[,2] = snps$pos[HAPLO[,2]]
  
  for(i in 1:(nrow(HAPLO)-1)){
    j = i + 1
    sizei = HAPLO[i,2]- HAPLO[i,1]
    sizej = HAPLO[j,2] - HAPLO[j,1]
    
    wj = sizei/(sizei+sizej)
    wi = sizej/(sizei+sizej)
    
    newlimit = round(weighted.mean(c(HAPLO[i,2],HAPLO[j,1]), c(wi,wj)))
    
    HAPLO[i,2] = newlimit
    HAPLO[j,1] = newlimit -1
  }
  
  
  
  
  HAPLO = HAPLO[order(HAPLO$pos1,HAPLO$pos2),]
  
  colnames(HAPLO) = c("pos1", "pos2", colnames(founders))
  
  return(HAPLO)
  
}



haploDfToMatrix = function(haplotype, snps){
  
  haplomatrix = lapply(split(haplotype, haplotype$rilname), function(x){
    print(x$rilname[1])
    #x=subset(haplotype,rilname == 101)
    xx = do.call(rbind,lapply(1:nrow(x), function(i){
      #print(i)
      y=x[i,] 
      w1 = which.min(abs(y$pos1-snps$pos))
      w2 = which.min(abs(y$pos2-snps$pos))
      founder = which(y[,3:(ncol(y)-1)]>0)
      if(length(founder)>0){
        return(cbind(w1:w2, rep(founder, w2-w1+1)) )}else{
          return(NULL)
        }
    }))
    
    mm = match(1:max(xx[,1]), xx[,1])
    xx = xx[mm,2]
    xx = matrix(xx, ncol=1)
    colnames(xx) = x$rilname[1]
    xx
    #unlist(xx)
  })
  
  nsnp = max(unlist(lapply(haplomatrix, nrow)))
  
  haplomatrix = do.call(cbind, lapply(haplomatrix,function(x){ x[match(1:nsnp,1:nrow(x))] }))
  
  haplomatrix
}




KeepMostFrequentHaplo = function(haplotype, rils){
  hapfreq = do.call(rbind, lapply(1:nrow(rils), function(i){
    if(i %% 1000 == 0){print(i)}
    pos = snps$pos[i]
    hap = haplotype[haplotype[,1] <= pos & haplotype[,2] > pos,]
    apply(hap[,3:(ncol(hap)-1)], 2, sum)/nrow(hap)
    
  }))
  
  #haplotype haplotype2
  haplotype = do.call(rbind, lapply(1:nrow(haplotype), function(i){
    if(i %% 1000 == 0){print(i)}
    x = haplotype[i,]
    wfounder = x[,3:(ncol(x)-1)]>0
    if(sum(wfounder)>1){
      pos1 = x[,1]
      pos2 = x[,2]
      pos1 = which.min(abs(snps$pos-pos1))
      pos2 = which.min(abs(snps$pos-pos2))
      
      freq = hapfreq[pos1:pos2,]#[wfounder]
      ismax = apply(freq, 2, mean)[wfounder]
      ismax = which.max(ismax)
      wfounder = which(wfounder)[ismax]
      x[,3:(ncol(x)-1)] = 0
      x[,wfounder+2] = 1
    }
    
    return(x)
  }))
  
}



hapstat = function(haplotype, winsize=10000){
  haplostat = do.call(rbind, lapply(seq(1, max(haplotype[,2]), winsize), function(pos){
    
    keep = haplotype[,1] <= pos & haplotype[,2] > pos
    x=haplotype[keep,]
    haplength = x[,2] - x[,1]
    nhap = sum(apply(x[,3:(ncol(x)-1)],2, sum) > 0)
    
    data.frame(pos = pos,
               mean.size = mean(haplength, na.rm=T),
               min.size = min(haplength, na.rm=T),
               max.size = max(haplength, na.rm=T),
               nhap = nhap)
    
  }))
  
  haplostat
}
