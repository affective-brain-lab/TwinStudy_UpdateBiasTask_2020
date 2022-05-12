# Matrix style model - Raw data - Continuous data
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

# Load Libraries & Options
rm(list=ls())
library(OpenMx)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
twindata<-read.table("desc.txt", header=T)

# Select Variables for Analysis
covVars      <- c('a1', 'a2', 'sex1','sex2','md1','md2','ed1','ed2','fd1','fd2','vd1','vd2','nd1','nd2') # list of variables names
nv        <- 1                          # number of variables
ntv       <- nv*2                       # number of total variables
selVars   <- c('Twin1', 'Twin2')

# Select Data for Analysis
mzData    <- subset(twindata, zygosity==1, c(selVars,covVars))
dzData    <- subset(twindata, zygosity==2, c(selVars,covVars))
cov(mzData,use="complete")
cov(dzData,use="complete")

# Set Starting Values
sMu      <- 14                        # start value for means
sVa      <- 40                        # start value for path coefficient
sVe      <- 40                        # start value for path coefficient for e

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean Matrices
intercept     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=sMu, labels="interC", name="intercept" )
betaS         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaS", name="bS" )
betaA         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaA", name="bA" )
betaM         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaM", name="bM" )
betaE         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaE", name="bE" )
betaF         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaF", name="bF" )
betaV         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaV", name="bV" )
betaN         <- mxMatrix( type="Full", nrow=1, ncol=1,   free=TRUE, values=0, labels="betaN", name="bN" )

defSex        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.sex1","data.sex2"), name="Sex" )
defAge        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.a1","data.a2"), name="Age" )
defmem        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.md1","data.md2"), name="mem" )
defexp        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.ed1","data.ed2"), name="exp" )
deffam        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.fd1","data.fd2"), name="fam" )
defviv        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.vd1","data.vd2"), name="viv" )
defneg        <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, labels=c("data.nd1","data.nd2"), name="neg" )
expMean       <- mxAlgebra( expression = intercept + bS*Sex + bA*Age + bM*mem+  bE*exp+  bF*fam+  bV*viv+  bN*neg, name="expMean" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVa, label="VA11", name="VA" ) 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVa, label="VC11", name="VC" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVe, label="VE11", name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ),
                                           cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ),
                                           cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
defs      <- list(defAge, defSex, defmem, defexp, deffam, defviv, defneg)
pars      <- list( intercept, betaS, betaA, betaM, betaE, betaF, betaV, betaN,  covA, covC, covE, covP )
modelMZ   <- mxModel( pars, defs, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, defs, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Variance Components
rowVC     <- rep('VC',nv)
colVC     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
estVC     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="VarC", dimnames=list(rowVC,colVC) )

# Create Confidence Interval Objects
ciACE     <- mxCI( "VarC[1,1:6]" )

# Build Model with Confidence Intervals
modelACE  <- mxModel( "ACEvc", pars,  modelMZ, modelDZ, multi, estVC, ciACE )
mxOption(NULL,"Default optimizer","NPSOL")

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACE Model
fitACE    <- mxRun( modelACE, intervals=T )
sumACE    <- summary( fitACE )
sumACE
# ----------------------------------------------------------------------------------------------------------------------

# RUN SUBMODELS

# Run AE model
modelAE   <- mxModel( fitACE, name="oneAEvc" )
modelAE   <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
fitAE     <- mxRun( modelAE, intervals=T )

# Run CE model
modelCE   <- mxModel( fitACE, name="oneCEvc" )
modelCE   <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
modelCE   <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
fitCE     <- mxRun( modelCE, intervals=T )

# Run E model
modelE    <- mxModel( fitAE, name="oneEvc" )
modelE    <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
fitE      <- mxRun( modelE, intervals=T )

# Print Comparative Fit Statistics
mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) )
round(rbind(fitACE$VarC$result,fitAE$VarC$result,fitCE$VarC$result,fitE$VarC$result),4)

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
