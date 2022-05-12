
library(OpenMx)
#source('/Users/annie/Dropbox/Affective_Brain_Lab/update_bias/modeling/src/optimism.R')
#source('/Users/annie/Dropbox/Affective_Brain_Lab/update_bias/modeling/src/optimism_bayes.R')
source('/Users/annie/Dropbox/Affective_Brain_Lab/Twins/Analysis/Kaustubh/src/twinModel.R')
#source('/Users/annie/Dropbox/moderation_mediation.R')
#source('/Users/annie/Dropbox/util.R')

library(gdata) #read.xls seems to be more robust than read.xlsx/2
library(gplots)
library(ggplot2)

read.events <- function(f=NULL) {
  #f: xlsx file
  if(is.null(f)) f <- "C:/Users/Kaustubh/Dropbox/Affective_Brain_Lab/Twins/Data/events.xlsx"
  events <- read.xls(f,1,header=F)
  rownames(events) <- events[,1]
  events <- events[,-1,drop=F]
  
  list1 <- read.xls(f,3,header=F)
  rownames(list1) <- list1[,1]
  list1 <- list1[,-1,drop=F]
  
  list2 <- read.xls(f,4,header=F)
  rownames(list2) <- list2[,1]
  list2 <- list2[,-1,drop=F]
  
  list3 <- read.xls(f,5,header=F)
  rownames(list3) <- list3[,1]
  list3 <- list3[,-1,drop=F]
  
  list4 <- read.xls(f,6,header=F)
  rownames(list4) <- list4[,1]
  list4 <- list4[,-1,drop=F]
  
  return(list(events=events, list1=list1, list2=list2, list3=list3, list4=list4))
}

read.dzmz <- function(f=NULL) {
  #f: xlsx file
  if(is.null(f)) f <- "C:/Users/Kaustubh/Dropbox/Affective_Brain_Lab/Twins/Data/MZ_DZ_gender.xlsx"
  dzmz <- read.xls(f,1,header=F)
  rownames(dzmz) <- dzmz[,1]
  dzmz <- dzmz[,-1,drop=F]
  return(dzmz)
}

twin.bdi <- function(x, bdi="BDI", questions=1:14, sums=F) {
  if(class(x)=="character") x <- read.xls(x,1)
  bdi <- paste(bdi,"_",questions,sep="")
  res <- list()  
  for(i in 1:nrow(x)) {
    res[[i]] <- as.numeric(as.character(x[i,bdi]))-1 #-1 because of qualtrics
  }
  if(sums) {
    res <- lapply(res, sum) #do not use na.rm for sum
    res <- do.call(rbind, res)
    colnames(res) <- "bdi"
  }
  return(res)
}

twin.sanity <- function(x, id="ID",zyg="MZ.DZ",auto.correct=T,max.cnt=Inf) {
  pairs <- twin.pairs(x[,id])
  n <- length(pairs[[1]])
  res <- list()
  res$age <- c()
  res$zyg <- c()
  res$maxcnt <- c()
  res$invcnt <- c()
  E <- paste("E",1:40,sep="")
  ids <- c()
  for(i in 1:n) {
    id1 <- pairs[[1]][i]
    id2 <- pairs[[2]][i]
    ids <- c(ids,id1,id2)
    if(x[id1,"Age"]!=x[id2,"Age"]) {
      res$age <- c(res$age,id1, id2)
      x[id1,"Age"] <- x[id2,"Age"] <- min(x[id1,"Age"],x[id2,"Age"])
    }
    
    if(x[id1,zyg]!=x[id2,zyg])
      res$zyg <- c(res$zyg,id1, id2)
    
    cnt1 <- sort(table(as.numeric(as.character(x[id1,E]))),decreasing=T)[1]
    cnt2 <- sort(table(as.numeric(as.character(x[id2,E]))),decreasing=T)[1]
    res$maxcnt <- c(res$maxcnt, cnt1, cnt2)
    if(cnt1>max.cnt || cnt2>max.cnt)
      res$invcnt <- c(res$invcnt, cnt1, cnt2)
    
  }
  
  inv <- c(res$zyg, res$invcnt) #combine all invalids dont take age
  i <- setdiff(ids, inv)
  
  if(auto.correct)
    res$x <- x[i,] #this removes singletones
  res$ids <- ids
  return(res)
}

twin.stats <- function(x) {
  ids <- rownames(x)
  pairs <- twin.pairs(ids)
  ids <- unlist(pairs)
  show(table(x[ids,"MZ.DZ"]))
  show(table(x[ids,"Gender"]))
}

lotr_score <- function(x,reverse=c(3,7,9), valid=c(1,3,4,7,9,10)) {
  #x: vector with 10 elements
  x <- as.numeric(as.character(x))-1
  x[reverse] <- 4-x[reverse]
  return(sum(x[valid]))
}

twin.lotr <- function(x, lotr="LOTR_Table", questions=1:10, reverse=c(3,7,9), valid=c(1,3,4,7,9,10), sums=F) {
  if(class(x)=="character") x <- read.xls(x,1)
  res <- list()
  for(i in 1:nrow(x)) {
    res[[i]] <- list()
    res[[i]]$lotr <- c()
    for(j in questions) {
      if(!any(valid==j)) next
      b <- paste(lotr,j,sep="_")
      xx <- as.numeric(as.character(x[i,b]))-1
      if(any(reverse==j)) xx <- 4-xx
      res[[i]]$lotr <- c(res[[i]]$lotr, xx)
    }		
  }
  if(sums) {
    res <- lapply(res, function(xx) {return(list(lotr=sum(xx$lotr)))} ) #do not use na.rm for sum
    res <- do.call(rbind, res)
    res <- apply(res, c(1,2), function(xx){return(unlist(xx))})
  }
  return(res)	
}

twin.post <- function(x,evt,lists,means=F) {
  #x: data frame
  #evt: events with rownames = events ids
  #lists: list of events with rownames = event ids ("1", "2" etc.)
  #acc: calculate average absolute difference
  
  #first get base rates for all lists
  res <- list()  
  base <- c()
  lmem <- c()
  lfam <- c()
  lexp <- c()
  lviv <- c()
  lneg <- c()
  i <- 0
  for(myevt in lists) { #go over every list
    i <- i+1 
    #memory
    base <- as.numeric(as.character(evt[rownames(myevt),2]))
    lmem <- c(lmem,paste("List_",i,"_Mem_",1:nrow(myevt),sep=""))    
    #familarity    
    lfam <- c(lfam,paste("List_",i,"_Fam_",1:nrow(myevt),sep=""))    
    #past experience    
    lexp <- c(lexp,paste("List_",i,"_Exp_",1:nrow(myevt),sep=""))    
    #vividness
    lviv <- c(lviv,paste("List_",i,"_Viv_",1:nrow(myevt),sep=""))    
    #negativity    
    lneg <- c(lneg,paste("List_",i,"_Neg_",1:nrow(myevt),sep=""))     
  }
  res$mem <- lapply(1:nrow(x), function(j){abs(base - (2+as.numeric(as.character(x[j,lmem]))))}) #2+ as qualtrics returns positions in dropdown
  res$fam <- lapply(1:nrow(x), function(j){as.numeric(as.character(x[j,lfam]))})
  res$exp <- lapply(1:nrow(x), function(j){as.numeric(as.character(x[j,lexp]))})
  res$viv <- lapply(1:nrow(x), function(j){as.numeric(as.character(x[j,lviv]))})
  res$neg <- lapply(1:nrow(x), function(j){as.numeric(as.character(x[j,lneg]))})
  
  names(res$mem) <- rownames(x)
  names(res$fam) <- rownames(x)
  names(res$exp) <- rownames(x)
  names(res$viv) <- rownames(x)
  names(res$neg) <- rownames(x)
  
  if(means)
    res <- lapply(res, function(xx){sapply(xx, mean, na.rm=T)})
  
  return(res)
}

twin.update <- function(x, evt,est1="E",est2="RE",questions=1:40,means=F,mean.center=F,upd.percent=F) {
  #x: data frame from read.xls, rows of x are subjects
  #events: from read.events
  
  if(class(x)=="character") x <- read.xls(x,1)
  
  events <- evt$events
  post <- twin.post(x,events,evt[2:5])
  
  res <- list()
  for(i in 1:nrow(x)) { #i is for each participant
    res[[i]] <- list()
    res[[i]]$est1 <- c()
    res[[i]]$est2 <- c()
    res[[i]]$base <- c()
    res[[i]]$upd <- c()
    res[[i]]$ee <- c()
    isdesirable <- c()
    idesirable <- c()
    iundesirable <- c()
    for(j in questions) { #j is for each question
      E <- paste(est1,j,sep="")
      RE <- paste(est2,j,sep="")
      
      E <- as.numeric(as.character(x[i,E]))
      RE <- as.numeric(as.character(x[i,RE]))
      e <- as.numeric(as.character(events[j,2]))      
      upd <- E-RE
      ee <- abs(E-e) #estimation error is absolute difference
      if(upd.percent) {
        upd <- 100*upd/E
        ee <- 100*ee/E
      }
      
      if(e<E) { #desirabe
        isdesirable <- c(isdesirable, TRUE)
        idesirable <- c(idesirable, j)
      }
      else if(e>E) {#undesirabe
        upd <- -upd
        isdesirable <- c(isdesirable, FALSE)
        iundesirable <- c(iundesirable, j)
      }
      else
        isdesirable <- c(isdesirable, NA)
      
      res[[i]]$est1 <- c(res[[i]]$est1,E)
      res[[i]]$est2 <- c(res[[i]]$est2,RE)
      res[[i]]$base <- c(res[[i]]$base,e)
      res[[i]]$upd <- c(res[[i]]$upd,upd)
      res[[i]]$ee <- c(res[[i]]$ee,ee)
    } #end of j
    res[[i]]$isdesirable <- isdesirable
    res[[i]]$idesirable <- idesirable
    res[[i]]$iundesirable <- iundesirable
    res[[i]]$mem <- post[["mem"]][[i]]
  }
  names(res) <- rownames(x)
  
  rm(idesirable)
  rm(iundesirable)
  
  if(means) {
    #correlation between esimation error and update
    #learn_desirable <- sapply(res, function(xx){cor(xx$upd[xx$idesirable], xx$ee[xx$idesirable],use="pairwise.complete.obs")})
    #learn_undesirable <- sapply(res, function(xx){cor(xx$upd[xx$iundesirable], xx$ee[xx$iundesirable],use="pairwise.complete.obs")})
    learn_desirable <- sapply(res, function(xx) {cor.my(xx$upd[xx$idesirable], xx$ee[xx$idesirable],use="pairwise.complete.obs")})
    learn_undesirable <- sapply(res, function(xx) {cor.my(xx$upd[xx$iundesirable], xx$ee[xx$iundesirable],use="pairwise.complete.obs")})
    
    ndesir <- sapply(res, function(xx) {length(xx$idesirable)})
    nundesir <- sapply(res, function(xx) {length(xx$iundesirable)})
        
    res <- lapply(res, function(xx) {
      return(list(est1=mean(xx$est1,na.rm=T),est2=mean(xx$est2,na.rm=T),                  
                  base=mean(xx$base,na.rm=T),upd=mean(xx$upd,na.rm=T),ee=mean(xx$ee,na.rm=T),                  
                  est1_desirable=mean(xx$est1[xx$idesirable],na.rm=T),est1_undesirable=mean(xx$est1[xx$iundesirable],na.rm=T),
                  est2_desirable=mean(xx$est2[xx$idesirable],na.rm=T),est2_undesirable=mean(xx$est2[xx$iundesirable],na.rm=T),
                  base_desirable=mean(xx$base[xx$idesirable],na.rm=T),base_undesirable=mean(xx$base[xx$iundesirable],na.rm=T),
                  ee_desirable=mean(xx$ee[xx$idesirable],na.rm=T),ee_undesirable=mean(xx$ee[xx$iundesirable],na.rm=T),
                  desirable=mean(xx$upd[xx$idesirable],na.rm=T),undesirable=mean(xx$upd[xx$iundesirable],na.rm=T)
      ))})
    
    res <- do.call(rbind, res)
    res <- apply(res, c(1,2), function(xx){return(unlist(xx))})
    goodbad <- rep(NA, nrow(res))
    res <- cbind(res,goodbad)
    res[,"goodbad"] <- res[,"desirable"] - res[,"undesirable"]
    
    res <- cbind(res,learn_desirable)
    res <- cbind(res,learn_undesirable)
    
    res <- cbind(res,ndesir)
    res <- cbind(res,nundesir)
    rownames(res) <- rownames(x)
  }
  return(res)
}

twin.Falconer <- function(x,vars,zyg="zyg",twin1="_twin1",twin2="_twin2") {
  fal <- matrix(nrow=length(vars),ncol=7)
  colnames(fal) <- c("h2","c2","e2","rdz","rmz","prdz","prmz")
  rownames(fal) <- vars
  dz <- sort(c(which(x[,zyg]=="dz"), which(x[,zyg]=="DZ")))
  for(v in vars) {    
    t1 <- paste(v,twin1,sep="")
    t2 <- paste(v,twin2,sep="")
    rdz <- cor.test(x[dz,t1], x[dz,t2],use="pairwise.complete.obs")
    rmz <- cor.test(x[-dz,t1], x[-dz,t2],use="pairwise.complete.obs")
    h2 <- 2*(rmz$estimate-rdz$estimate)
    c2 <- rmz$estimate - h2
    e2 <- 1-h2-c2
    fal[v,] <- c(h2,c2,e2,rdz$estimate,rmz$estimate,rdz$p.value,rmz$p.value)
  }
  return(fal)
}

#separate in MZ and DZ
dzmz.separate <- function(x, dzmz, cols=c("ID"), extracols=NULL, twin.sep=F)  {  
  
  if(!is.null(extracols)) x <- cbind(x, extracols)
  ids <- x[,"ID"]
  if(class(ids)=="factor") ids <- levels(ids)[ids]
  row.names(x) <- ids
  
  rdzmz <- as.character(rownames(dzmz))
  if(length(rdzmz)==0) rdzmz <- ids
  
  if(class(dzmz)=="matrix" || class(dzmz)=="data.frame") dzmz <- dzmz[,1]
  else if(class(dzmz)=="character") {
    dzmz <- x[, dzmz]
    names(dzmz) <- ids
  }
  
  if(twin.sep){
    twins <- twin.separate(as.numeric(ids))
    stopifnot(length(twins$twin1)==length(twins$twin2))
    mzs <- list()
    dzs <- list()
    mzs$twin1 <- c()
    mzs$twin2 <- c()
    dzs$twin1 <- c()
    dzs$twin2 <- c()
    
    processed <- c()
    for(i in 1:length(twins$twin1)) {
      id1 <- as.character(twins$twin1[i])
      id2 <- as.character(twins$twin2[i])      
      if(any(processed==id1) || any(processed==id2)) next  
      processed <- c(processed, id1,id2)      
      if(dzmz[id1]=="MZ" && dzmz[id2]=="MZ") { #check both as some might be missing
        mzs$twin1 <- c(mzs$twin1, id1)
        mzs$twin2 <- c(mzs$twin2, id2)
      }
      else if(dzmz[id1]=="DZ" && dzmz[id2]=="DZ") {
        dzs$twin1 <- c(dzs$twin1, id1)
        dzs$twin2 <- c(dzs$twin2, id2)
      }
    }
        
# following doesnt work well for missing/non-interleaved data 
#     xx <- list()
#     xx$twin1 <- ids[mzs & (ids %in% twins$twin1)]
#     xx$twin2 <- ids[mzs & (ids %in% twins$twin2)]
#     mzs <- xx
#     
#     xx <- list()
#     xx$twin1 <- ids[dzs & (ids %in% twins$twin1)]
#     xx$twin2 <- ids[dzs & (ids %in% twins$twin2)]
#     dzs <- xx
  }
  else {
    i <- which(dzmz=="MZ")
    mzs <- ids %in% rdzmz[i]
    i <- which(dzmz=="DZ")
    dzs <- ids %in% rdzmz[i]
  }
  
  if(length(cols)==0) cols <- colnames(x)
  
  res <- list()
  for(cl in cols) {
    res[[cl]] <- list()
    
    if(twin.sep) {
      res[[cl]]$DZ[["twin1"]] <- x[dzs$twin1, cl]
      res[[cl]]$DZ[["twin2"]] <- x[dzs$twin2, cl]
      
      res[[cl]]$MZ[["twin1"]] <- x[mzs$twin1, cl]
      res[[cl]]$MZ[["twin2"]] <- x[mzs$twin2, cl]
      
      #do a sanity check
      #if(cl=="ID")
      #if(as.numeric(res[[cl]]$DZ$twin1[1])!=as.numeric(res[[cl]]$DZ$twin1[1])+1 || as.numeric(res[[cl]]$DZ$twin1[1])!=as.numeric(res[[cl]]$DZ$twin1[1])-1) {
       # cat("dzmz.separate is doing a bad job")
        #stop()
      #}
      
    }
    else {
      res[[cl]]$DZ <- x[dzs, cl]
      res[[cl]]$MZ <- x[mzs, cl]	
    }
  }	
  return(res)
}

twin.pairs <- function(ids) {
  ids <- as.numeric(as.character(ids))
  #get twin pairs
  res <- list()
  res$twin1 <- c()
  res$twin2 <- c()
  for(id in ids) {
    if(is.na(id)) next
    if(any(res$twin1==id)) next
    if(any(res$twin2==id)) next
    j <- c(which(ids==(id+1)), which(ids==(id-1)))
    if(length(j==1)) {
      res$twin1 <- c(res$twin1, id)
      res$twin2 <- c(res$twin2, ids[j])
    }
  }
  res <- lapply(res, as.character)
  return(res)
}

twin.separate <- function(ids) {
  #make sure that ids is not factor
  stopifnot(class(ids)!="factor")
  idsn <- as.numeric(ids)
  twin1 <- c()
  twin2 <- c()
  rm <- c()
  for(i in 1:length(ids)) {
    if( (ids[i] %in% twin1) || (ids[i] %in% twin2) ) next
    #there might be characters at the end
    ii <- idsn[i]
    if(is.na(ii)) {rm <- c(rm, i); next}
    else ii <- ii+1
    j <- which(idsn==ii)
    if(length(j)==0) {
      ii <- idsn[i]
      if(is.na(ii)) {rm <- c(rm, i); next}
      else ii <- ii-1
      j <- which(idsn==ii)
    }
    if(length(j)==0) {
      cat("twin not found: ", ids[i],"\n")
      #stop()
    }
    else {
      twin1 <- c(twin1, ids[i])
      twin2 <- c(twin2, ids[j])
    }
  }
  #show(twin1)
  #cat("-------------------------\n")
  #show(twin2)
  return(list(twin1=twin1, twin2=twin2, rm=rm))
}

twin.df <- function(x) {
  #x: output of dzmz.separate
  #only process lists that have MZ and DZ sub-lists
  res <- NULL
  nn <- c("zyg","mz","twin")
  nr <- 0
  for(n in names(x)) {
    if(class(x[[n]])!="list") next
    if(is.null(x[[n]]$MZ) || is.null(x[[n]]$DZ)) next
    if(is.null(x[[n]]$MZ$twin1) || is.null(x[[n]]$DZ$twin1)) next
    if(is.null(x[[n]]$MZ$twin2) || is.null(x[[n]]$DZ$twin2)) next
    if(is.null(res)) {      
      nr <- sum(sapply(x[[n]]$MZ, length), sapply(x[[n]]$DZ, length))
      res <- data.frame(row.names=1:nr)
      zyg <- c(rep("dz",sum(sapply(x[[n]]$DZ, length))), rep("mz",sum(sapply(x[[n]]$MZ, length))))
      res <- cbind(res, zyg)
      mz <- c(rep(0,sum(sapply(x[[n]]$DZ, length))), rep(1,sum(sapply(x[[n]]$MZ, length))))
      res <- cbind(res, mz)
      
      tw <- c(rep(1,length(x[[n]]$DZ$twin1)), rep(2,length(x[[n]]$DZ$twin2)))
      tw <- c(tw, c(rep(1,length(x[[n]]$MZ$twin1)), rep(2,length(x[[n]]$MZ$twin2))))
      res <- cbind(res, tw)      
    }
    res <- cbind(res, c(x[[n]]$DZ$twin1, x[[n]]$DZ$twin2, x[[n]]$MZ$twin1, x[[n]]$MZ$twin2))
    nn <- c(nn,n)
  }
  colnames(res) <- nn
  return(res)
}

twin.df.twin <- function(x) {
  #x: output of dzmz.separate
  #keep twin1 and twin2 separate
  res <- NULL
  nn <- c("zyg","mz")
  nr <- 0
  ndz <- 0
  nmz <- 0
  for(n in names(x)) {
    if(class(x[[n]])!="list") next
    if(is.null(x[[n]]$MZ) || is.null(x[[n]]$DZ)) next
    if(is.null(x[[n]]$MZ$twin1) || is.null(x[[n]]$DZ$twin1)) next
    if(is.null(x[[n]]$MZ$twin2) || is.null(x[[n]]$DZ$twin2)) next
    if(is.null(res)) {
      ndz <- length(x[[n]]$DZ$twin1)
      nmz <-  length(x[[n]]$MZ$twin2)
      res <- data.frame(row.names=1:(ndz+nmz))
      res <- cbind(res, c(rep("dz",ndz), rep("mz",nmz)))
      res <- cbind(res, c(rep(0,ndz), rep(1,nmz)))
    }
    res <- cbind(res, c(x[[n]]$DZ$twin1, x[[n]]$MZ$twin1))
    res <- cbind(res, c(x[[n]]$DZ$twin2, x[[n]]$MZ$twin2))
    nn <- c(nn,c(paste(n,"_twin1",sep=""), paste(n,"_twin2",sep="")))
  }
  colnames(res) <- nn
  return(res)
}

twin.df.balance <- function(x,id="ID",zyg="zyg") {
  #given data frame x remove rows for which twin doesnt exist
  to.rm <- c()
  for(i in 1:nrow(x)) {
    j <- c(which(x[,id]==(x[i,id]+1)),which(x[,id]==(x[i,id]-1)))
    k <- which(x[,zyg]==x[i,zyg])
    k <- intersect(j,k)
    if(length(k)==0) to.rm <- c(to.rm, i)
  }
  if(length(to.rm)>0) {
    x <- x[-to.rm,]
  }
  return(x)
}

twin.lm <- function(x,indep="mz") {
  #x: data frame from twin.df
  x <- na.omit(x)
  x <- twin.df.balance(x)  
  txt <- paste("lr <- glm(",indep,"~desirable+undesirable+lotr+bdi+mem,family=binomial(logit), data=x)")
  eval(parse(text=txt))    
  return(list(x=x,desir=lm1, undesir=lm2, dz=lr))
}

library(effects)
twin.lm.twin <- function(x,zyg="zyg",vars=c(),scale=F,rand=F) {
  #x: data frame from twin.df.twin
  #dont na.omit here as every var can have different NAs
  zyg2 <- zyg
  cat("\nnumber of pairs")
  show(table(x[,zyg]))
  nn <- setdiff(colnames(x), zyg)
  if(scale) x[,nn] <- scale(x[,nn])
  
  if(rand) {    
    for(n in nn)
      x[,n] <- sample(x[,n], nrow(x))
  }
  
  res <- list()
  #par(mfrow=c(4,2))
  if(length(vars)==0)
    vars <- c("bdi","lotr","mem","desirable","undesirable","learn_desirable","learn_undesirable","goodbad","z","npos","nneg")
  for(v in vars) {
    dep <- paste(v,"_twin1",sep="")
    indep <- paste(v,"_twin2",sep="")
    Y <- x[,dep]
    X <- x[,indep]
    i <- c(attr(na.omit(X),"na.action"), attr(na.omit(Y),"na.action"))
    i <- setdiff(1:nrow(x),i)
    assign(dep,Y[i])
    assign(indep,X[i])
    assign(zyg2, x[i,zyg2])
    txt <- paste("res[[v]] <- glm(",dep,"~",indep,"*",zyg2,")",sep="")
    #cat(txt,"\n")
    eval(parse(text=txt))    
  }  
  return(res)
}


plot.effect <- function(mod,jpg="",...) {
  terms <- attr(mod$terms,"term.labels")
  i <- grep(":",terms)
  if(length(i)==0) return
  if(length(i)>1) return
  term <- terms[i]
  eff <- effect(term=term,mod=mod)
  
  cf <- summary(mod)$coefficients
  n <- rownames(cf)
  sub <- ""
  for(i in 1:nrow(cf)) {
    sub <- paste(sub, n[i], "beta=",round(cf[i,1],2), "p=", round(cf[i,4],2),"|")
  }
  
  #if(jpg!="") jpeg(jpg,width = 800, height = 600)
  plot(eff,multiline=T,sub=sub,lines=c(1,1),...)
  #if(jpg!="") dev.off()
}

#plot MZ-DZ separated
twin.plot <- function(x,prop,zyg="MZ.DZ",title="",se=F,eqn=T,coord.ratio=1) {
  ids <- rownames(x)
  pairs <- twin.pairs(ids)
  df <- data.frame(row.names=1:length(pairs$twin1))
  df <- cbind(df, x[pairs$twin1,prop],x[pairs$twin2,prop], x[pairs$twin1,zyg])
  xvar <- "Twin1"
  yvar <- "Twin2"
  colnames(df) <- c(xvar, yvar,"Zyg")
  
  if(eqn) {
    eqns <- 1:nrow(df)
    for(z in unique(df[,"Zyg"])) {
      txt <- paste("m <- lm(",yvar,"~",xvar,",data=subset(df,Zyg=='",z,"'))",sep="")
      eval(parse(text=txt))
      mydf <- subset(df,Zyg==z)
      r <- cor(mydf[,xvar],mydf[,yvar],use="pairwise.complete.obs")
      r <- format(round(r, 2), nsmall=2)
      label <- paste(z,": ",lm_eqn(m, paste(" r=", r,sep="")),sep="")
      i <- which(df[,"Zyg"]==z)
      eqns[i] <- label
    }
    eqns <- as.factor(eqns)
    df[,"Zyg"] <- eqns
  }
  
  #get interaction
  txt <- paste("m <- lm(",yvar,"~",xvar,"*Zyg,data=df)",sep="")
  eval(parse(text=txt))
  s <- summary(m)
  beta <- round(coef(m)[4],2)
  p <- format(s$coefficients[4,4],digits=4,drop0Trailing=F,format="f")
  title <- paste(title, "(interact: beta =",beta," p =",p,")")
  
  p <- ggplot(df, aes_string(x=xvar, y=yvar, color="Zyg", shape="Zyg")) +
    geom_point() + geom_smooth(method=lm,se=se,fullrange=T) + labs(title=title)
  p <- p + theme(legend.position="bottom",legend.title=element_blank()) 
  
  if(coord.ratio)
    p <- p + coord_fixed(ratio=coord.ratio) + geom_abline(intercept=0, slope=1, linetype="dashed")
  
  return(p)
}

lm_eqn <- function(m,...) {
  s <- summary(m)
  beta <- coef(m)[2]
  p <- s$coefficients[2,4]
  r2 <- s$r.squared
  n <- length(m$residuals)
  
  beta <- format(round(beta, 2), nsmall = 2)
  r2 <- format(round(r2, 2), nsmall = 2)
  p <- formatC(p, digits = 4, format = "f")
  
  intercept <- round(coef(m)[1],2)
  paste("N=",n,"\nbeta=",beta," p=",format(p,drop0Trailing=F),"\nr2=",r2,...,sep="")
}

twin.read.data <- function(f=NULL,id="ID") {
  if(is.null(f)) f <- "C:/Users/Kaustubh/Dropbox/Affective_Brain_Lab/Twins/Data/10_09_2013/Compile_twins_only2.xlsx"
  dat <- read.xls(f,1,header=T)
  cat("read",nrow(dat), "rows and",ncol(dat),"columns\n")
  ids <- as.character(dat[,id])
  rownames(dat) <- ids
  return(dat)
}

twin.analyze <- function(f=NULL,fe=NULL,fz=NULL,id="ID",max.age=Inf,min.age=0,min.anseach=0,max.mem=Inf,max.bdi=100,questions=1:40,ml.bal=F,jags=F,lm.rand=F,ctrl=c(),upd.percent=F,dz.bal=F,balmz=F,sanity=T,auto.correct=T,exclude=c(),onlypairs=F,models=c("ACE","reg","Falconer","med")) {
  #read in data
  dat <- twin.read.data(f)
  #get zygosity
  if(!is.null(fz)) {
    zz <- c()
    fz <- "C:/Users/Kaustubh/Dropbox/Affective_Brain_Lab/Twins/Data/MZ_DZ_gender.xlsx"
    z <- read.xls(fz,2,header=T)    
    row.names(z) <- z[,"Study_No"]
    for(n in row.names(dat)) {
      i <- which(row.names(z)==n)
      zz <- c(zz,z[i,"Zygosity"])
    }
    dat <- cbind.data.frame(dat,list("MZ.DZ"=zz))    
  }
  
  show(table(dat[,"MZ.DZ"]))
  events <- read.events(fe) #this is all spreadsheets
  evt <- events$events
  cat("read",nrow(evt), "events\n")
  
  cat("processing data... ")
  dat <- subset(dat, Age<=max.age)
  dat <- subset(dat, Age>min.age)
  
  dat <- dat[setdiff(rownames(dat),c("431","432")),] #431 didnt understand the task
  dat <- dat[setdiff(rownames(dat),exclude),]
  
  if(sanity) {
    san <- twin.sanity(dat, auto.correct=auto.correct)
    if(auto.correct) {
      dat <- san$x
    }
  }
  show(table(dat[,"MZ.DZ"]))
  
  if(dz.bal) {
    i <- which(dat[,"MZ.DZ"]=="DZ")
    dat <- rbind(dat, dat[i,])
  }
  
  if(balmz) { #we have many more MZ    
    ndz <- as.integer(length(which(dat[,"MZ.DZ"]=="DZ"))/2)
    i <- which(dat[,"MZ.DZ"]=="MZ")
    mzpairs <- twin.pairs(dat[i,id])
    nmz <- length(mzpairs[[1]])
    cat("number of pairs (dz, mz):", ndz, nmz,"\n")
    #select ndz pairs
    if(nmz>ndz) {
      i <- sample(1:nmz,nmz-ndz)
      to.rm <- c(mzpairs[[1]][i], mzpairs[[2]][i])
      cat("removing MZ subjects: ", length(to.rm),"\n")
      dat <- dat[!(ids%in%to.rm),]
      ids <- as.character(dat[,id])
      show(table(dat[,"MZ.DZ"]))
    }
  }
  
  #get some statistics
  lotr <- twin.lotr(dat,sums=T)
  bdi <- twin.bdi(dat,sums=T)
  upd <- twin.update(dat,events,means=T,upd.percent=upd.percent,questions=questions)
  post <- twin.post(dat, evt, events[2:5],means=T)
  
  vars <- c("bdi","lotr","mem","fam","exp","viv","neg","est1","est2","ee","desirable","undesirable","learn_desirable","learn_undesirable","goodbad")
  dat <- cbind(dat, cbind(bdi,lotr,upd))
  dat <- cbind(dat,do.call(cbind,post))
  
  #filters
  nas <- c()
  if(min.anseach>0) {
    minanseach <- pmin(upd[,"ndesir"], upd[,"nundesir"])
    nas <- which(minanseach<min.anseach)
    nas <- unique(nas)
    if(length(nas)>0) {
      cat("excluding",length(nas),"subjects due to not enough trials")
      show(table(dat[nas,"MZ.DZ"]))      
    }
  }
  
  i <- which(bdi>max.bdi)
  if(length(i)>0) {
    cat("BDI filter: excluding",length(i),"subjects with BDI >",max.bdi)
    show(table(dat[i,"MZ.DZ"]))
    nas <- c(nas,i)
  }
  
  i <- which(post$mem>max.mem)
  if(length(i)>0) {
    cat("mem filter: excluding",length(i),"subjects with mem err >",max.mem)
    nas <- c(nas,i)
  } 
  
  nas <- unique(nas)
  if(length(nas)>0) dat <- dat[-nas,]
  if(onlypairs) {
    ids <- rownames(dat)
    ids <- unlist(twin.pairs(ids))
    dat <- dat[ids,]
    show(dim(dat))
  }
  
  cat("done\n")
  
  #ML models
  cat("fitting ML models... ")
  est1n <- paste("E",questions,sep="")
  est2n <- paste("RE",questions,sep="")
  base <- as.numeric(as.character(evt[,2,drop=T]))
  z <- c()
  npos <- c()
  nneg <- c()
  ml.value <- c()
  zj <- c()
  nposj <- c()
  nnegj <- c()
  biases <- list()
  for(i in 1:nrow(dat)) { #per subject
    est1 <- as.numeric(sapply(dat[i,est1n],as.character))
    est2 <- as.numeric(sapply(dat[i,est2n],as.character))    
    bias <- optimism_bias(est1,est2,base,jags=jags,bal=ml.bal)
    z <- c(z, bias$ml$par["z"])
    npos <- c(npos, bias$ml$par["npos"])
    nneg <- c(nneg, bias$ml$par["nneg"])
    
    if(jags) {
      zj <- c(zj, bias$jags$par["z"])
      nposj <- c(nposj, bias$jags$par["npos"])
      nnegj <- c(nnegj, bias$jags$par["nneg"])
    }
    biases[[i]] <- bias
  }
  
  cat("done\n")
  
  vars <- c(vars,"z","npos","nneg")
  dat <- cbind(dat,z,npos,nneg)
  if(jags) {
    vars <- c(vars,"zj","nposj","nnegj")
    dat <- cbind(dat,zj,nposj,nnegj)
  }
  
  #control
  if(length(ctrl)>0) {
    if(any(ctrl %in% vars == T)) {
      cat("controllers are in vars")
      stop()
    }
    ctrl <- paste(ctrl,collapse="+")
    cat("controlling for:",ctrl,"\n")    
    for(v in vars) {
      resid <- c()
      txt <- paste(" resid <- glm(",v,"~",ctrl,",data=dat)$residuals")      
      eval(parse(text=txt))
      #glm doesnt like NAs so replace accordingly
      dat[names(resid),v] <- resid
    }
  }
  
  #plot
  plot.optimism(dat)
  
  #twin analysis
  #get twins separated
  if(length(models))
    twins <- dzmz.separate(dat, dzmz="MZ.DZ",cols=c("ID",vars),twin.sep=T)
  
  #run ACE model
  if("ACE" %in% models) {
    cat("ACE models... ")
    ACE <- list()
    for(v in vars) {
      if(v=="ID") next
      cat(v," ")
      mzData <- cbind(twins[[v]]$MZ$twin1, twins[[v]]$MZ$twin2)
      dzData <- cbind(twins[[v]]$DZ$twin1, twins[[v]]$DZ$twin2)
      ACE[[v]] <- twinModel(dzData, mzData, rm_param=c("c","a","ac"))
    }
    cat("done\n")
  }
  
  #regression models
  if("reg" %in% models) {
    cat("Regression models... ")
    twdf <- twin.df.twin(twins)
    reg <- twin.lm.twin(twdf,vars=vars,rand=lm.rand)
    
    #summary results for regressions
    summs <- lapply(reg,function(xx){summary(xx)$coefficients})
    reg.summary <- matrix(nrow=length(summs),ncol=8)
    for(i in 1:length(summs)) {
      reg.summary[i,1:4] <- summs[[i]][,"Estimate"]
      reg.summary[i,5:8] <- summs[[i]][1:4,4]
    }
    rownames(reg.summary) <- names(summs)
    colnames(reg.summary) <- c("Intercept","X","Zyg","X*Zyg","pIntercept","pX","pZyg","pX*Zyg")
    
    cat("done\n")
  }
  
  #Falconer
  if("Falconer" %in% models) {
    cat("Falconer... ")
    fal <- twin.Falconer(twdf,vars=vars)
    cat("done\n")
  }
  
  i <- intersect(which(is.na(dat[,"desirable"])==F),which(is.na(dat[,"undesirable"])==F))
  j <- intersect(which(is.na(dat[,"learn_desirable"])==F),which(is.na(dat[,"learn_undesirable"])==F))
  cat("\nnumber of subjects with updates and learning data:", length(i), length(j))
  twin.stats(dat[i,])
  twin.stats(dat[j,])
  
  #mediation models
  if("med" %in% models) {
    mvars <- c("bdi","lotr")
    xvars <- c("bdi","lotr")
    yvars <- c("desirable","undesirable")
    medn <- c()
    med <- data.frame()
    for(xv in xvars) {
      for(yv in yvars) {
        if(yv==xv) next
        for(mv in mvars) {
          if(mv==yv) next
          if(mv==xv) next
          mm <- moderation_mediation(dat[,xv], dat[,yv], dat[,mv])$med
          medn <- c(medn, paste(xv,yv,mv,"all",sep="-"))
          med <- rbind.data.frame(med,mm)
          df <- subset(dat,MZ.DZ=="MZ")
          mm <- moderation_mediation(df[,xv], df[,yv], df[,mv])$med
          medn <- c(medn, paste(xv,yv,mv,"MZ",sep="-"))
          med <- rbind.data.frame(med,mm)
          df <- subset(dat,MZ.DZ=="DZ")
          mm <- moderation_mediation(df[,xv], df[,yv], df[,mv])$med
          medn <- c(medn, paste(xv,yv,mv,"DZ",sep="-"))
          med <- rbind.data.frame(med,mm)
        }
      }
    }
    rownames(med) <- medn
  }
  
  res <- list(events=events,dat=dat,min.aneach=min.anseach)
  if(length(models)) res$twins <- twins
  if("ACE" %in% models) res$ACE <- ACE
  
  if("reg" %in% models) {
    res$reg <- reg
    res$reg.summary <- reg.summary
  }
  if("Falconer" %in% models) res$Falconer <- fal
  if("med" %in% models) res$med <- med
  
  res$bias <- biases
  return(res)
}

twin.digin <- function(tw, n) {
  res <- list()
  res$est1 <- tw$bias[[n]]$est1/100
  res$est2 <- tw$bias[[n]]$est2/100
  res$base <- tw$bias[[n]]$base/100
  res$desir <- res$est1>res$base
  res$ll <- optimism_ll_beta_binom(tw$bias[[n]]$ml$par, cbind(res$est1,res$est2,res$base),ret.arr=T)
  res$value <- tw$bias[[n]]$ml$value
  res$par <- tw$bias[[n]]$ml$par
  return(res)
}

# plot.effect(xxd$reg$models$bdi)
# plot.effect(xxd$reg$models$lotr)
# plot.effect(xxd$reg$models$mem)
# plot.effect(xxd$reg$models$desirable)
# plot.effect(xxd$reg$models$undesirable)
# plot.effect(xxd$reg$models$learn_desirable)
# plot.effect(xxd$reg$models$learn_undesirable)