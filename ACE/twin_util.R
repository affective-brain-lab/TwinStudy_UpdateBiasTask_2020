# separate given data in MZ and DZ
# if twin.sep is set to TRUE then the twin pair is also separated
# the first argument is the data matrix
# second argument is a column name of x that contains MZDZ status
twin.dzmz.separate <- function(x, dzmzcol="MZ.DZ", idcol="ID", cols=c(), twin.sep=FALSE) {
  ids <- x[,idcol]
  if(class(ids)=="factor") ids <- levels(ids)[ids]
  row.names(x) <- ids
  
  rdzmz <- as.character(rownames(dzmzcol))
  if(length(rdzmz)==0) rdzmz <- ids
  
  dzmz <- x[,dzmzcol]
  names(dzmz) <- ids
  
  if(twin.sep){
    twins <- twin.separate(as.numeric(ids)) # get the twin ids
    stopifnot(length(twins$twin1)==length(twins$twin2))
    mzs <- list(twin1=c(), twin2=c())
    dzs <- list(twin1=c(), twin2=c())
    
    processed <- c()
    for(i in 1:length(twins$twin1)) {
      id1 <- as.character(twins$twin1[i])
      id2 <- as.character(twins$twin2[i])      
      if(any(processed==id1) || any(processed==id2)) next
      processed <- c(processed, id1,id2)
      if(dzmz[id1]=="MZ" && dzmz[id2]=="MZ") { #check both as some info might be missing
        mzs$twin1 <- c(mzs$twin1, id1)
        mzs$twin2 <- c(mzs$twin2, id2)
      } else if(dzmz[id1]=="DZ" && dzmz[id2]=="DZ") {
        dzs$twin1 <- c(dzs$twin1, id1)
        dzs$twin2 <- c(dzs$twin2, id2)
      }
    }
  }
  else {
    i <- which(dzmz=="MZ")
    mzs <- ids %in% rdzmz[i]
    i <- which(dzmz=="DZ")
    dzs <- ids %in% rdzmz[i]
  }
  
  if(length(cols)==0) cols <- colnames(x)
  
  # create a list of DZ-MZ separated columns and twin separated if asked
  # this can be passed to the ACE model
  res <- list()
  for(cl in cols) {
    res[[cl]] <- list()
    
    if(twin.sep) {
      res[[cl]]$DZ$twin1 <- x[dzs$twin1, cl]
      res[[cl]]$DZ$twin2 <- x[dzs$twin2, cl]
      
      res[[cl]]$MZ$twin1 <- x[mzs$twin1, cl]
      res[[cl]]$MZ$twin2 <- x[mzs$twin2, cl]      
    }
    else {
      res[[cl]]$DZ <- x[dzs, cl]
      res[[cl]]$MZ <- x[mzs, cl]	
    }
  }	
  return(res)
}

# assuming that the twins are numbered consequitevly, find the pairs
twin.separate <- function(ids) {
  #make sure that ids is numeric
  stopifnot(class(ids)=="numeric" || class(ids)=="integer")
  twin1 <- c()
  twin2 <- c()
  torm <- c()
  for(i in 1:length(ids)) {
    if( (ids[i] %in% twin1) || (ids[i] %in% twin2) ) next
    #there might be characters at the end
    if(is.na(ids[i])) {torm <- c(torm, ids[i]); next}
    j <- which(ids==ids[i]+1 | ids==ids[i]-1) # look for id that either higher or lower by 1
    if(length(j)==0) {
      cat("twin not found:", ids[i],"\n")
    } else if(length(j)>1) { # we must find a single twin
      stop("multiple twins found for:",ids[i],"STOPPING!\n")
    } else {
      twin1 <- c(twin1, ids[i])
      twin2 <- c(twin2, ids[j])
    }
  }
  return(list(twin1=twin1, twin2=twin2, rm=torm))
}
