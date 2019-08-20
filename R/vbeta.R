##' Estimate the expected variance of beta.  This is approximately expected(1/var(U)).
##'
##' Assumes we have a list, GenoProbList, giving the GenoProb values for each X. 
##' @title Estimate the expected variance of beta
##' @export
##' @inheritParams expected_z_score
##' @inheritParams simulated_z_score
##' @return The expected variance of beta for each SNP X, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
##' @examples
##' freq=fake_freq(nhaps=100,nsnps=5) # fake haplotype frequency data
##' EVB=expected_vbeta(N0=1000,N1=2000,snps=paste0("s",1:5),
##'                     W="s1",gamma.W=log(1.5),freq=freq)
##' EVB # causal variant is SNP 1, with OR 1.5
expected_vbeta<-function(N0,N1,snps,W,gamma.W,freq,
                         GenoProbList=make_GenoProbList(snps=snps,W=W,freq=freq)){
                                        #check that we have SNPs X and W in the reference dataset
    fvbeta<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
        ## uses Rcpp file ../src/est_zscore.cpp
        vbeta(N0,N1,Ufactor,powerfactor,GenoProbXW[[2]],GenoProbXW[[3]])
    }

    if (!all(c(snps,W) %in% colnames(freq)))
        stop("SNPs of interest not present in reference dataset.")
    if(length(gamma.W)!=length(W))
        stop("length mismatch: gamma.W and W")
    if(length(GenoProbList)!=length(snps))
        stop("GenoProbList should have same length and order as snps")
    g0 <- compute_gamma0(N0=N0,N1=N1,W=W,gamma.W=gamma.W,freq=freq)
    ## compute P(Y=1 | W=w)
    N<-N0+N1
    expeta<-exp(g0+rowSums(sweep((combinat::hcube(rep(3,length(W)))-1),MARGIN=2,gamma.W,`*`)))
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
##' @inheritParams expected_z_score
##' @inheritParams simulated_z_score
##' @return A simulated variance of beta for each SNP X, assuming the causal SNPs are W
##' @export
##' @author Mary Fortune and Chris Wallace
##' @examples
##' freq=fake_freq(nhaps=100,nsnps=5) # fake haplotype frequency data
##' VB=simulated_vbeta(N0=1000,N1=2000,snps=paste0("s",1:5),
##'                     W="s1",gamma.W=log(1.5),freq=freq)
##' VB # causal variant is SNP 1, with OR 1.5
simulated_vbeta<-function(N0,N1,snps,W,gamma.W,freq,
                          GenoProbList=make_GenoProbList(snps=snps,W=W,freq=freq),
                          nrep=1){
                                        #check that we have SNPs X and W in the reference dataset
    fab<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
        ## uses Rcpp file ../src/est_zscore.cpp
        vbeta_ab(N0,N1,Ufactor,powerfactor,GenoProbXW[[2]],GenoProbXW[[3]])
    }

    if (!all(c(snps,W) %in% colnames(freq)))
        stop("SNPs of interest not present in reference dataset.")
    if(length(gamma.W)!=length(W))
        stop("length mismatch: gamma.W and W")
    if(length(GenoProbList)!=length(snps))
        stop("GenoProbList should have same length and order as snps")
    g0 <- compute_gamma0(N0=N0,N1=N1,W=W,gamma.W=gamma.W,freq=freq)
    ## compute P(Y=1 | W=w)
    N<-N0+N1
    expeta<-exp(g0+rowSums(sweep((combinat::hcube(rep(3,length(W)))-1),MARGIN=2,gamma.W,`*`)))
                                        #compute the constant factors we will multiply by
    Ufactor<-N0*(N-1)*(N0*expeta-N1)/(N^2)
    powerfactor<-N0*(expeta+1)/N
    AB <- lapply(seq_along(snps), function(ii) {
        fab(N0,N1,Ufactor,powerfactor,freq,GenoProbList[[ii]])
    })
    alpha <- sapply(AB,"[[","alpha")
    beta <- sapply(AB,"[[","beta")
    V <- lapply(1:nrep, function(ab) 1/rgamma(length(alpha), shape=alpha, scale=1/beta) *
                                     (N0*N1/(N)))
    1/do.call("rbind",V)
}

cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}
wcor2 <- function (x, w = rep(1, nrow(x))/nrow(x)) {
    ## normalize 
    w <- w / sum(w)
    ## center 
    x <- sweep(x, 2, colSums(x * w))
    ## compute 
    cs <- colSums(w * x**2)
    crossprod(w*x, x) / sqrt(tcrossprod(cs)) 
}
