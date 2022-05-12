source('C:/Users/Kaustubh/Dropbox/src/R/plot_util.R')

rs2snp <- function(genfile, rs, rscol=2) {
  #this is really really SLOW!!!
  con <- file(genfile, "r")
  snps <- list()
  n <- 0
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    for(r in rs) {
      if(length(grep(r, line))) {
        snps[[r]] <- line
        break
      }
    }
    n <- n+1
    if((n %% 10000)==0) cat(".")    
    if(length(snps)==length(rs)) break
  }
  close(con)
  return(snps)
}

#get genotypes using triplets  
gen2genotype <- function(x,samplefile="",startcol=6,sampleidcol=1) {
  #x: a line from gen file
  #4th and 5th columns are nucleotides
  
  if(class(x)=="character") x <- unlist(strsplit(x,split="\\s+"))
  n <- length(x)
  AA <- paste(x[4], x[4], sep="")
  AB <- paste(x[4], x[5], sep="")
  BB <- paste(x[5], x[5], sep="")
  genotype <- c()
  nsub <- (n-(startcol-1))/3
  if(nsub != as.integer(nsub)) {
    cat("Non-integer number of subjects:", nsub)
    return(NULL)
  }
  for(i in 1:nsub) {
    gen <- x[startcol:(startcol+2)]
    gt <- NA
    j <- which(gen==max(gen))
    if(length(j)==1) {
      if(j==1) gt <- AA
      else if(j==2) gt <- AB
      else gt <- BB
    }
    genotype <- c(genotype, gt)
    startcol <- startcol+3
  }
  
  #get samplefile if provided
  if(samplefile!="") {
    s <- read.table(samplefile, header=T)
    s <- s[-c(1),] #remove first row
    names(genotype) <- s[,sampleidcol]
  }
  return(genotype)
}

phenogenotype2barplot <- function(p,g,title="",dir=NULL,type="svg",suffix="") {
  #p: phenotype df with each column as different phenotype
  #g: genotype df
  if(is.null(colnames(p))) colnames(p) <- 1:ncol(p)
  if(is.null(colnames(g))) colnames(g) <- 1:ncol(g)
  
  #get common pheno-genotypes
  i <- intersect(rownames(p), rownames(g))
  p <- p[i,,drop=F]
  g <- g[i,,drop=F]
  for(pt in colnames(p)) {
    for(gt in colnames(g)) {
      cat(pt, gt,"\n")
      mygt <- as.vector(g[,gt])
      mypt <- list()
      for(gg in sort(unique(mygt))) {
        i <- which(mygt==gg)
        mypt[[gg]] <- p[i,pt,drop=TRUE]
      }
      #we have all we need for this phenotype
      #plot it
      if(!is.null(dir)) {
        if(type=="svg") svg(paste(dir,"/",pt,"_",gt,suffix,".svg",sep=""))
        else if(type=="png") png(paste(dir,"/",pt,"_",gt,suffix,".png",sep=""))
        else if(type=="jpeg") jpeg(paste(dir,"/",pt,"_",gt,suffix,".jpeg",sep=""))
      }
      pl <- barplot.sem3(mypt,xlab=gt,ylab=pt,title=title)
      print(pl)
      if(!is.null(dir)) dev.off()
    }
  }
}
