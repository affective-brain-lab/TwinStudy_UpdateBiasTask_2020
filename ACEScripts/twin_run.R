# Twin data analysis using ACE and Falconer models
# Kaustubh R. Patil, April 2018

# You need a single CSV file which contains all the information
# It should contain the subject IDs, Zygosity, Age, Gender etc
# as well as the variables which you want to calculate the heritibility of.
# The IDs MUST be numeric.
# Make sure that the data is aligned as per the IDs such that
# each row represents data of exactly the same subject.
# Once you have this file then start R or RStudio and change the session
# to the directory containing this file "twin_run.R".
# Open this file and source it.
# type in the following command and ptress enter:
# results <- twin_run()
# It will ask you to make a series of choices including selecting 
# the csv file with the data, the independent variable, confounds etc.
# After that it will run the models and return the results in the 
# results variable.
# You can look at the detailts using the following commands
# results$Falconer
# results$ACE
# To get more details follow the lists, e.g.
# results$ACE$summary

source("twin_util.R")
source("util.R")
source("twinModel.R")

# run twin analysis
twin_run <- function(dat=NULL, indep=NULL, ctrl=c(), ctrl.factor=rep(FALSE, length(ctrl)), dzmzcol=NULL, idcol=NULL, models=c("ACE", "Falconer"), optimizer="NPSOL") {
  # set optimizer
  # options CSOLNP, NPSOL, 
  if(is.null(optimizer)) {
    optimizers <- c("NPSOL", "CSOLNP", "SLSQP")
    optimizer <- menu(optimizers,graphics=TRUE,title="Select optimizer")
    optimizer <- optimizers[optimizer]
  }
  mxOption(NULL, "Default optimizer", optimizer)
  
  # read in the data
	if(is.null(dat)) dat <- file.choose()
	if(class(dat)=="character") dat <- read.csv(dat, header=TRUE, row.names=NULL)
	stopifnot(!is.null(colnames(dat))) # we must have column names
	
	# get the variables and ncessary cols if not provided
	if(is.null(idcol)) {
		idcol <- menu(colnames(dat),graphics=TRUE,title="Select ID column")
		idcol <- colnames(dat)[idcol]
	}
	
	if(is.null(dzmzcol)) {
		dzmzcol <- menu(colnames(dat),graphics=TRUE,title="Select DZMZ column")
		dzmzcol <- colnames(dat)[dzmzcol]
	}
	
	stopifnot(any(colnames(dat)==idcol))
	stopifnot(any(colnames(dat)==dzmzcol))
	
	if(is.null(indep)) {
	  indep <- menu(colnames(dat),graphics=TRUE,title="Select independent variable")
	  indep <- colnames(dat)[indep]
	}
	
	# get the confounds if not provided
	if(length(ctrl)==0) {
		ctrl <- select.list(colnames(dat),multiple=TRUE,title='Select confounds (press control key)',graphics=TRUE)
		if(length(ctrl)) {
			ff <- select.list(ctrl,multiple=TRUE,title='Select factors (press control key)',graphics=TRUE)
			ctrl.factor <- rep(FALSE, length(ctrl))
			ctrl.factor[ctrl %in% ff] <- TRUE
		}
	}
	
	# control for confounds if asked
	if(length(ctrl)) {
		stopifnot(all(ctrl %in% colnames(dat)))
		res <- control_confounds(dat, indep, ctrl, ctrl.factor)
		ctrl.initial <- paste(sapply(ctrl, function(xx) unlist(strsplit(xx,""))[1]),collapse="")
		indep <- paste(indep,".ctrl.",ctrl.initial,sep="")
		dat[[indep]] <- res
	}
	
	# get the twin separated data
	twins <- twin.dzmz.separate(dat, dzmzcol=dzmzcol, idcol=idcol, cols=indep, twin.sep=TRUE)
	
	res <- list()
	
	#run ACE model
	if("ACE" %in% models) {
		cat("ACE model... ")
		mzData <- cbind(twins[[indep]]$MZ$twin1, twins[[indep]]$MZ$twin2)
		dzData <- cbind(twins[[indep]]$DZ$twin1, twins[[indep]]$DZ$twin2)
		storage.mode(mzData) <- "double"
		storage.mode(dzData) <- "double"
		res$ACE <- twinModel(dzData, mzData, rm_param=c("c","a","ac"))
		cat("done\n")
	}
	
	if("Falconer" %in% models) {
		rdz <- cor.test(twins[[indep]]$DZ$twin1, twins[[indep]]$DZ$twin2,use="pairwise.complete.obs")
		rmz <- cor.test(twins[[indep]]$MZ$twin1, twins[[indep]]$MZ$twin2,use="pairwise.complete.obs")
		h2 <- 2*(rmz$estimate-rdz$estimate)
		c2 <- rmz$estimate - h2
		e2 <- 1-h2-c2
		res$Falconer <- c(h2,c2,e2,rdz$estimate,rmz$estimate,rdz$p.value,rmz$p.value)
		names(res$Falconer) <- c("h2","c2","e2","R.dz","R.mz","Pval.dz","Pval.mz")
	}
	
	res
}

results <- twin_run() 
results
summary(results$ACE$fits$ACE)
