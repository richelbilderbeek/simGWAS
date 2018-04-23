##' Estimate the expected variance of beta.  This is approximately expected(1/var(U)).
##'
##' Assumes we have a list, GenoProbList, giving the GenoProb values for each X. 
##' @title Estimate the expected variance of beta
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma1	The odds ratios of effect of the true causal SNPs (not including gamma0, the intercept term)
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbList An list of objects giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
expected_vbeta<-function(N0,N1,snps,W,gamma1,freq,GenoProbList){
                                        #check that we have SNPs X and W in the reference dataset
    fvbeta<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
        ## uses Rcpp file ../src/est_zscore.cpp
        vbeta(N0,N1,Ufactor,powerfactor,GenoProbXW[[2]],GenoProbXW[[3]])
    }

    if (!all(c(snps,W) %in% colnames(freq)))
        stop("SNPs of interest not present in reference dataset.")
    if(length(gamma1)!=length(W))
        stop("length mismatch: gamma1 and W")
    if(length(GenoProbList)!=length(snps))
        stop("GenoProbList should have same length and order as snps")
    g0 <- compute_gamma0(N0=N0,N1=N1,W=W,gamma.CV=gamma1,freq=freq)
    ## compute P(Y=1 | W=w)
    N<-N0+N1
    expeta<-exp(g0+rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma1,`*`)))
                                        #compute the constant factors we will multiply by
    Ufactor<-N0*(N-1)*(N0*expeta-N1)/(N^2)
    powerfactor<-N0*(expeta+1)/N
    sapply(seq_along(snps), function(ii) {
        fvbeta(N0,N1,Ufactor,powerfactor,freq,GenoProbList[[ii]])
    })
}


##' Simulate var(beta)
##'
##' Assumes we have a list, GenoProbList, giving the GenoProb values for each X. 
##' @title Compute a simulated var(beta)
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma1	The odds ratios of effect of the true causal SNPs (not including gamma0, the intercept term)
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbList An list of objects giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @export
##' @author Mary Fortune and Chris Wallace
simulated_vbeta<-function(N0,N1,snps,W,gamma1,freq,GenoProbList,nsim=1){
                                        #check that we have SNPs X and W in the reference dataset
    fab<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
        ## uses Rcpp file ../src/est_zscore.cpp
        vbeta_ab(N0,N1,Ufactor,powerfactor,GenoProbXW[[2]],GenoProbXW[[3]])
    }

    if (!all(c(snps,W) %in% colnames(freq)))
        stop("SNPs of interest not present in reference dataset.")
    if(length(gamma1)!=length(W))
        stop("length mismatch: gamma1 and W")
    if(length(GenoProbList)!=length(snps))
        stop("GenoProbList should have same length and order as snps")
    g0 <- compute_gamma0(N0=N0,N1=N1,W=W,gamma.CV=gamma1,freq=freq)
    ## compute P(Y=1 | W=w)
    N<-N0+N1
    expeta<-exp(g0+rowSums(sweep((hcube(rep(3,length(W)))-1),MARGIN=2,gamma1,`*`)))
                                        #compute the constant factors we will multiply by
    Ufactor<-N0*(N-1)*(N0*expeta-N1)/(N^2)
    powerfactor<-N0*(expeta+1)/N
    AB <- lapply(seq_along(snps), function(ii) {
        fab(N0,N1,Ufactor,powerfactor,freq,GenoProbList[[ii]])
    })
    alpha <- sapply(AB,"[[","alpha")
    beta <- sapply(AB,"[[","beta")
    V <- lapply(1:nsim, function(ab) 1/rgamma(length(alpha), shape=alpha, scale=1/beta) *
                                     (N0*N1/(N)))
}

