#tutorial: http://www.vipbg.vcu.edu/HGEN619/HGEN619_OpenMx.shtml

if(!require("OpenMx")){
  #install.packages("OpenMx")
  # Kaustubh make sure that we install the version that supports NPSOL
  source('https://vipbg.vcu.edu/vipbg/OpenMx2/software/getOpenMx.R')
  library(OpenMx)
  # change the optimizer to NPSOL instaed of default CSOLNP
  # CSOLNP can end in endless loop even though is should have been fixed in version 2.9.6 https://openmx.ssri.psu.edu/node/4339
  mxOption(NULL,"Default optimizer","NPSOL")
}

source('GenEpiHelperFunctions.R') # for checking the assumptions of the saturated model

twinModel <- function(dzData, mzData, std=FALSE, rm_param=c("c","a","ac"), model="ACE", verbose=FALSE) {
  # mzData & dzData: matrices with two columns
  #std: scale the data
  #rm_param: parameters to remove; 'a', 'c', 'ac'
  
  selVars <- colnames(mzData)
  if(is.null(selVars)[1]) selVars <- c("T1","T2")
  
  colnames(mzData) <- selVars
  colnames(dzData) <- selVars
  
  if(std) {
    mzData <- scale(mzData)
    dzData <- scale(dzData)
  }
  nv <- 1
  
  twinSatModel <- mxModel("sat",
                          mxModel("MZ",
                                  mxMatrix("Full", 1, 2, T, c(0,0), name="expMeanMZ"), 
                                  mxMatrix("Lower", 2, 2, T, .5, name="CholMZ"), 
                                  mxAlgebra(CholMZ %*% t(CholMZ), name="expCovMZ"), 
                                  mxData(mzData, type="raw"),
                                  mxFIMLObjective("expCovMZ", "expMeanMZ", selVars),
                                  mxAlgebra( expMeanMZ[1,1], name="expMeanMZt1"),  		# mean for twin1 from expMeanMZ 
                                  mxAlgebra( expMeanMZ[1,2], name="expMeanMZt2"),			# mean for twin2 from expMeanMZ 
                                  mxAlgebra( expCovMZ[1,1], name="expVarMZt1"),			# variance for twin1 from expCovMZ 
                                  mxAlgebra( expCovMZ[2,2], name="expVarMZt2"),			# variance for twin2 from expCovMZ 
                                  mxAlgebra( expCovMZ[2,1], name="CovMZ")				# covariance between twin1 and twin2 from expCovMZ 
                                  
                                  ),
                          mxModel("DZ",
                                  mxMatrix("Full", 1, 2, T, c(0,0), name="expMeanDZ"), 
                                  mxMatrix("Lower", 2, 2, T, .5, name="CholDZ"), 
                                  mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ"),
                                  mxData(dzData, type="raw"), 
                                  mxFIMLObjective("expCovDZ", "expMeanDZ", selVars),
                                  # Algebra's needed for equality constraints    
                                  mxAlgebra( expMeanDZ[1,1], name="expMeanDZt1"),    	# mean for twin1 from expMeanDZ 
                                  mxAlgebra( expMeanDZ[1,2], name="expMeanDZt2"),			# mean for twin2 from expMeanDZ 
                                  mxAlgebra( expCovDZ[1,1], name="expVarDZt1"),			# variance for twin1 from expCovDZ 
                                  mxAlgebra( expCovDZ[2,2], name="expVarDZt2"),			# variance for twin2 from expCovDZ 
                                  mxAlgebra( expCovDZ[2,1], name="CovDZ")				# covariance between twin1 and twin2 from expCovDZ
                                  ),
                          mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
                          mxAlgebraObjective("twin")
                    ) #mxModel sat
  twinSatFit <- mxRun(twinSatModel, suppressWarnings=TRUE)
  assum <- twinModel.satAssumptions(twinSatFit)
  i <- which(assum[,"p"]<.05)
  if(length(i)) assum <- assum[i,"comparison"]
  else assum <- NULL
  
#   twinSatModelSub1 <- mxModel(twinSatModel, name="satSub1")
#   twinSatModelSub1$MZ$expMeanMZ <- mxMatrix("Full", 1, 2, T, 0, "mMZ", name="expMeanMZ")
#   twinSatModelSub1$DZ$expMeanDZ <- mxMatrix("Full", 1, 2, T, 0, "mDZ", name="expMeanDZ")
#   twinSatFitSub1 <- mxRun(twinSatModelSub1, suppressWarnings=TRUE)
#   
#   twinSatModelSub2 <- mxModel(twinSatModelSub1, name="satSub2")
#   twinSatModelSub2$MZ$expMeanMZ <- mxMatrix("Full", 1, 2, T, 0, "mean", name="expMeanMZ")
#   twinSatModelSub2$DZ$expMeanDZ <- mxMatrix("Full", 1, 2, T, 0, "mean", name="expMeanDZ")
#   twinSatFitSub2 <- mxRun(twinSatModelSub2, suppressWarnings=TRUE)
#   
#   LL_Sat <- mxEval(objective, twinSatFit)
#   LL_Sub1 <- mxEval(objective, twinSatFitSub1)
#   LRT1 <- LL_Sub1 - LL_Sat
#   LL_Sub2 <- mxEval(objective, twinSatFitSub2)
#   LRT2 <- LL_Sub2 - LL_Sat
  
  twinModel <- mxModel(model,
                       mxModel("univACE",
                               # Specify matrices a, c, and e to store a, c, and e path coefficients
                               mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0.6, label="a11", name="a" ), 
                               mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0.6, label="c11", name="c" ), 
                               mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=0.6, label="e11", name="e" ),  
                               # Matrices A, C, and E compute variance components
                               mxAlgebra( expression=a %*% t(a), name="A" ),			# getting the variance by multiplying a and the transpose of a
                               mxAlgebra( expression=c %*% t(c), name="C" ),
                               mxAlgebra( expression=e %*% t(e), name="E" ),
                               # Algebra to compute total variances and standard deviations
                               mxAlgebra( expression=A+C+E, name="V" ),			# creating a vector with the total variance (A, C, and E) called V
                               mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I"),
                               # create an identity matrix without free parameters (a value of 1 for all entries on the diagonal and 0 in all off-diagonal entries)
                               mxAlgebra( expression=solve(sqrt(I*V)), name="iSD"),	
                               # create a vector with the SDs by taking the inverse ("solve") of the square root of the total variance multiplied by the identity matrix
                               # Algebra to compute standardized path estimates and variance components
                               mxAlgebra( expression=a%*%iSD, name="sta"),			# create a vector with the standardised path coefficients for the additive genetic effects
                               mxAlgebra( expression=c%*%iSD, name="stc"),
                               mxAlgebra( expression=e%*%iSD, name="ste"),
                               mxAlgebra( expression=A/V, name="a2"),				
                               # create a vector with the genetic variance (A) devided by the total variance (V) to get standardized estimates of the proportion of variance due to genes (heritability)
                               mxAlgebra( expression=C/V, name="c2"),
                               mxAlgebra( expression=E/V, name="e2"),
                               # Note that the rest of the mxModel statements do not change for bivariate/multivariate case
                               # Matrix & Algebra for expected means vector
                               mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, label="mean", name="Mean" ),
                               # create a Matrix for the expected means
                               mxAlgebra( expression= cbind(Mean,Mean), name="expMean"), # create a matrix containing means for twin 1 and twin 2
                               # Algebra for expected variance/covariance matrix in MZ		
                               # "rbind" joins the expected variances and covariances together vertically, and "cbind" joins the matrices horizontally
                               mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),		
                                                               cbind(A+C   , A+C+E)), name="expCovMZ" ),
                               # Algebra for expected variance/covariance matrix in DZ (using the Kronicke product of 0.5 and A as only half of the genes are shared between the DZ twins
                               mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),
                                                               cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 
                       ),
                       mxModel("MZ",
                               mxData( observed=mzData, type="raw" ),
                               mxFIMLObjective( covariance="univACE.expCovMZ", means="univACE.expMean", dimnames=selVars )
                               # specify that the MZ data should be analyzed using the mean and MZ covariances computed in the ACE mxModel
                       ),
                       mxModel("DZ", 
                               mxData( observed=dzData, type="raw" ),
                               mxFIMLObjective( covariance="univACE.expCovDZ", means="univACE.expMean", dimnames=selVars ) 
                               # specify that the DZ data should be analyzed using the mean and DZ covariances computed in the ACE mxModel
                       ),
                       mxAlgebra( expression=MZ.objective + DZ.objective, name="m2ACEsumll" ),
                       mxAlgebraObjective("m2ACEsumll"),
                       mxCI(c('univACE.a2', 'univACE.c2', 'univACE.e2'))					# get the confidence intervals on the standardized variance components
  )
  
  
  if(model=="AC") {
    txt <- "twinModel$univACE.e <- mxMatrix(type='Full',nrow=1,ncol=1,free=F,values=0,label='e11')"
    eval(parse(text=txt))
  }
  else if(model=="AE") {
    txt <- "twinModel$univACE.c <- mxMatrix(type='Full',nrow=1,ncol=1,free=F,values=0,label='c11')"
    eval(parse(text=txt))
  }
  else if(model=="CE") {
    txt <- "twinModel$univACE.a <- mxMatrix(type='Full',nrow=1,ncol=1,free=F,values=0,label='a11')"
    eval(parse(text=txt))
  }
  
  fits <- list()
  cmp <- list()
  fits[[model]] <- mxRun(twinModel,intervals=T)
  cmp <- mxCompare(twinSatFit, fits[[model]])
  
  for(r in rm_param) {
    myName <- gsub(toupper(r),"","ACE")
    myModel <- mxModel(twinModel, name=myName)
    for(rr in unlist(strsplit(r,""))) {
      #set corresponding matrix to 0 and non-free
      txt <- paste("myModel$univACE.",rr, " <- mxMatrix(type='Full',nrow=1,ncol=1,free=F,values=0,label='",rr,"11')",sep="")
      eval(parse(text=txt))
    }
    fits[[myName]] <- mxRun(myModel,intervals=T)
    cmp <- rbind(cmp,mxCompare(fits[[model]], fits[[myName]])[2,])
  }
  summ <- do.call(rbind,lapply(fits,twinModel.summary))
  fits[["sat"]] <- twinSatFit
  
  npairs <- c(nrow(dzData), nrow(mzData))
  npairs_complete <- c(nrow(dzData[complete.cases(dzData),]), nrow(mzData[complete.cases(mzData),]))
  names(npairs) <- c("DZ", "MZ")
  names(npairs_complete) <- c("DZ", "MZ")
  
  res <- list(fits=fits,summary=list(summary=summ,compare=cmp,violations=assum,npairs=npairs,npairs_complete=npairs_complete))
  res
}

twinModel.summary <- function(x) {

  summs <- summary(x)
  twinModel.summ <- rep(NA,14)
  LL <- mxEval(objective, x)
  
    twinModel.summ[1:3] <- summs$CI[1:3,"estimate"]
    twinModel.summ[4:5] <- summs$CI[1,c("lbound","ubound")]
    twinModel.summ[6:7] <- summs$CI[2,c("lbound","ubound")]
    twinModel.summ[8:9] <- summs$CI[3,c("lbound","ubound")]
    twinModel.summ[10:11] <- summs$information[1:2,"par"] #AIC and BIC
    twinModel.summ[12:13] <- summs$information[1:2,"sample"] #AIC and BIC
    twinModel.summ[14] <- LL
  
  names(twinModel.summ) <- c("A","C","E","lbA","ubA","lbC","ubC","lbE","ubE","AIC.par","BIC.par","AIC.sample","BIC.sample","minus2LL")
  
  twinModel.summ
}

twinModel.satAssumptions <- function(univTwinSatFit) {
  #taken and modified from: http://genepi.qimr.edu.au/staff/sarahMe/files/Assumption_Testing.R
  #details: http://www.vipbg.vcu.edu/HGEN619/OpenMxModelAssumptions.pdf
  #it is incremental testing such that each model adds constraints to the previous
  
  # Assumption Testing - Means - check for mean differences between the different groups 
  # -----------------------------------------------------------------------
  # Check for birth order effects - equate means of twin 1 and twin 2 for MZ and DZ twins
  equateMeansModel1 <- mxModel(univTwinSatFit, name="Mean t1=t2",
                               mxConstraint( MZ.expMeanMZt1 == MZ.expMeanMZt2, name="MeanMZt1t2"), 
                               mxConstraint( DZ.expMeanDZt1 == DZ.expMeanDZt2, name="MeanDZt1t2")
                      )     
  equateMeans1Fit <- mxRun(equateMeansModel1, suppressWarnings=TRUE)  			# run reduced model with means for twin 1 and 2 equated
  
  # Check for DZ and MZ differences - equate means of MZ and DZ twins
  equateMeansModel2 <- mxModel(equateMeansModel1 , name="Mean DZ=MZ",
                               mxConstraint( MZ.expMeanMZt1 == DZ.expMeanDZt1, name="Mean")	# as twin1 and twin2 are equated (above) it is enough to equate twin1 of MZ/DZ group
  )     
  equateMeans2Fit <- mxRun(equateMeansModel2, suppressWarnings=TRUE)
  
  # Assumption Testing - Variances - check for variance differences between the different groups
  # -----------------------------------------------------------------------
  # Check for birth order effects - equate variances of twin 1 and twin 2
  equateVarModel1 <- mxModel(equateMeansModel2, name="Var t1=t2",
                             mxConstraint( MZ.expVarMZt1 == MZ.expVarMZt2, name="VarMZt1t2"),
                             mxConstraint( DZ.expVarDZt1 == DZ.expVarDZt2, name="VarDZt1t2")
  )     
  equateVar1Fit <- mxRun(equateVarModel1, suppressWarnings=TRUE)
  
  # Check for DZ and MZ differences - equate variances of MZ and DZ twins
  equateVarModel2 <- mxModel(equateVarModel1 , name="Var DZ=MZ",
                             mxConstraint( MZ.expVarMZt1 == DZ.expVarDZt1, name="Var")
  )     
  equateVar2Fit <- mxRun(equateVarModel2, suppressWarnings=TRUE)
  
  # Assumption Testing - Covariances
  # -----------------------------------------------------------------------
  # Check for DZ and MZ differences - check for covariance differences between MZ and DZ groups
  equateCovModel1 <- mxModel(equateVarModel2 , name="Cov DZ=MZ",
                             mxConstraint( MZ.CovMZ == DZ.CovDZ, name="OneCov")
  )     
  equateCov1Fit <- mxRun(equateCovModel1, suppressWarnings=TRUE)
  cmp <- mxCompare(univTwinSatFit, list(equateMeans1Fit,equateMeans2Fit,equateVar1Fit,equateVar2Fit,equateCov1Fit  ))
  return(cmp)
}

#write.csv(twinModel.summary)
