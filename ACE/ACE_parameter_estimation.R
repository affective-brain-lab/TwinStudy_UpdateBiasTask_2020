twinACEModel <- mxModel("twinACE",
    mxModel("ACE",
    # Matrices a, c, and e to store a, c, and e path coefficients
        mxMatrix(
            type="Lower",
            nrow=1,
            ncol=1,
            free=TRUE,
            values=0.6,
            labels="a11",
            name="a"
        ),
        mxMatrix(
            type="Lower",
            nrow=1,
            ncol=1,
            free=TRUE,
            values=0.6,
            labels="c11",
            name="c"
        ),
        mxMatrix(
            type="Lower",
            nrow=1,
            ncol=1,
            free=TRUE,
            values=0.6,
            labels="e11",
            name="e"
        ),
    # Matrices A, C, and E compute variance components
        mxAlgebra(
            expression=a %*% t(a),
            name="A"
        ),
        mxAlgebra(
            expression=c %*% t(c),
            name="C"
        ),
        mxAlgebra(
            expression=e %*% t(e),
            name="E"
        ),
    # Matrix & Algebra for expected means vector
        mxMatrix(
            type="Full",
            nrow=1,
            ncol=1,
            free=TRUE,
            values=20,
            label="mean",
            name="Mean"
        ),
        mxAlgebra(
            expression= cbind(Mean,Mean),
            name="expMean"
        ),
    # Algebra for expected variance/covariance matrix in MZ
        mxAlgebra(
            expression=rbind (cbind(A + C + E , A + C),
                              cbind(A + C     , A + C + E)),
            name="expCovMZ"
        ),
    # Algebra for expected variance/covariance matrix in DZ
        mxAlgebra(
            expression=rbind (cbind(A + C + E     , 0.5 %x% A + C),
                              cbind(0.5 %x% A + C , A + C + E)),
            name="expCovDZ"
        )
    ),
    mxModel("MZ",
        mxData(
            observed=mzData,
            type="raw"
        ),
        mxFIMLObjective(
            covariance="ACE.expCovMZ",
            means="ACE.expMean",
            dimnames=selVars
        )
    ),
    mxModel("DZ",
        mxData(
            observed=dzData,
            type="raw"
        ),
        mxFIMLObjective(
            covariance="ACE.expCovDZ",
            means="ACE.expMean",
            dimnames=selVars
        )
    ),
    mxAlgebra(
        expression=MZ.objective + DZ.objective,
        name="minus2loglikelihood"
    ),
    mxAlgebraObjective("minus2loglikelihood")
 )