
##' Compute vector of expected Z scores
##' 
##' @title Compute vector of expected Z Scores
##' @export
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma.W	The log odds ratios of effect of the true causal SNPs (not including gamma0, the intercept term)
##' @param freq Haplotype frequencies as a data.frame, with column Probability indicating relative frequency in controls.  
##' @param GenoProbList An list of objects giving the probability of seeing each {X,W} genotype vector.  This can be calculated within the function if no value supplied, or you can pass a pre-calculated version
##' @return The expected Z Score for all snps in snps, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
##' @examples
##' freq=fake_freq(nhaps=100,nsnps=5) # fake haplotype frequency data
##' EZ=expected_z_score(N0=1000,N1=2000,snps=paste0("s",1:5),
##'                     W="s1",gamma.W=log(1.5),freq=freq)
##' EZ # causal variant is SNP 1, with OR 1.5
expected_z_score<-function(N0,N1,snps,W,gamma.W,freq,
                            GenoProbList=make_GenoProbList(snps=snps,W=W,freq=freq)){
    if(!(is.data.frame(freq) & "Probability" %in% colnames(freq)))
        stop("freq must be a data.frame, with a column 'Probability' giving relative frequency of each haplotype (row)")
    if(!all(snps %in% colnames(freq)))
        stop("not all snps found in freq")
    if(!all(W %in% colnames(freq)))
        stop("not all W found in freq")
    exp_z_score<-est_statistic(N0,N1,snps,W,gamma.W,freq,GenoProbList)
    return(exp_z_score)
}

##' Compute matrix of simulated Z scores about expected values of 0 -
##' ie under a null of no association at any SNP
##' @title Compute a NULL simulated Z Score
##' @export
##' @inheritParams expected_z_score
##' @param nrep Number of replicates (simulated vectors of Z scores)
##'     under this scenario.  Default=1
##' @author Mary Fortune and Chris Wallace
##' @examples
##' freq=fake_freq(nhaps=100,nsnps=5) # fake haplotype frequency data
##'     Z=simulated_z_null(snps=paste0("s",1:5),freq=freq,nrep=3)
##' Z # no causal variants
simulated_z_null<-function(snps,freq, nrep=1){
    exp_z_score<- rep(0,length(snps))
    LD<-wcor2(as.matrix(freq[,setdiff(colnames(freq),"Probability")]),
              freq$Probability)
    sim_z_score<-rmvnorm(n=nrep,mean=exp_z_score,sigma=LD)
    if(nrep==1)
        return(c(sim_z_score))
    sim_z_score
}

##' Compute matrix of simulated Z scores
##' @title Compute a simulated Z Score
##' @export
##' @inheritParams expected_z_score
##' @param nrep Number of replicates (simulated vectors of Z scores)
##'     under this scenario.  Default=1
##' @author Mary Fortune and Chris Wallace
##' @examples freq=fake_freq(nhaps=100,nsnps=5) # fake haplotype frequency data
##'     Z=simulated_z_score(N0=1000,N1=2000,snps=paste0("s",1:5),
##'                         W="s1",gamma.W=log(1.5),freq=freq,nrep=3)
##'     Z # causal variant is SNP 1, with OR 1.5
simulated_z_score<-function(N0,N1,snps,W,gamma.W,freq,
                            GenoProbList=make_GenoProbList(snps=snps,W=W,freq=freq),
                            nrep=1){
    exp_z_score<- expected_z_score(N0,N1,snps,W,gamma.W,freq,GenoProbList)
    LD<-wcor2(as.matrix(freq[,setdiff(colnames(freq),"Probability")]),
              freq$Probability)
    sim_z_score<-rmvnorm(n=nrep,mean=exp_z_score,sigma=LD)
    if(nrep==1)
        return(c(sim_z_score))
    sim_z_score
}
## simulated_z_score.old<-function(N0,N1,snps,W,gamma.W,freq,df_control,nrep=1){
## 	GenoProbList<-make_GenoProbList(snps,W,freq)
## 	gamma.sim<-c(compute_gamma0(N0,N1,W,gamma.W,freq),gamma.W)
## 	exp_z_score<-est_statistic(N0,N1,snps,W,gamma.sim,freq,GenoProbList)
## 	XX<-new("SnpMatrix", as.matrix(df_control))
## 	LD <- snpStats::ld(XX,XX,stat="R",symmetric=TRUE)
## 	LD<-as.matrix(make.positive.definite(LD))
## 	sim_z_score<-c(rmvnorm(n=nrep,mean=exp_z_score,sigma=LD))
## 	return(sim_z_score)
## }


##' Estimates the expected Z Score for a single SNP
##'
##' Assumes the input CVs, and the relationship, gamma, between them and the trait of interest
##'
##' Assumes we have already generated GenoProbXW for all X
##' @title estimate Z score at a single SNP
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param Ufactor	The constant factor used to compute the expectation of U
##' @param powerfactor	The constant factor used to compute the expectation of the genotype of X to some power
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbXW An object giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
est_zscore<-function(N0,N1,Ufactor,powerfactor,freq,GenoProbXW){
    ## uses Rcpp file ../src/est_zscore.cpp
  ## N<-N0+N1
  #Contribution of each W=w to the Sum
  #P(X=1 AND W=w)
  # PX1W<-GenoProbXW[[2]]
  #P(X=2 AND W=w)
  # PX2W<-GenoProbXW[[3]]
    zscore(N0,N1,Ufactor,powerfactor,GenoProbXW[[2]],GenoProbXW[[3]])
}


##' Wrapper function to run est_zscore for all snps in snps
##'
##' Assumes we have a list, GenoProbList, giving the GenoProb values for each X. 
##' @title estimate Z score at a single SNP
##' @param N0 The number of Y=0
##' @param N1	The number of Y=1
##' @param snps The snps at which we wish to compute the expected Z Score
##' @param W	The true causal SNPs (these need not be in "snps")
##' @param gamma.W	The odds ratios of effect of the true causal SNPs (not including gamma0, the intercept term)
##' @param freq Frequencies of SNP appearances (computed using snphap)
##' @param GenoProbList An list of objects giving the probability of seeing each {X,W} genotype vector
##' @return The expected Z Score for SNP X, assuming the causal SNPs are W
##' @author Mary Fortune and Chris Wallace
est_statistic<-function(N0,N1,snps,W,gamma.W,freq,GenoProbList){
                                        #check that we have SNPs X and W in the reference dataset
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
        est_zscore(N0,N1,Ufactor,powerfactor,freq,GenoProbList[[ii]])
    })
}


#returns P(X=x, W=w)
find_PXaW_MK<-function(x,w,GenoProbXW){
  XW<-paste0(x,paste0(w,collapse=""),collapse="")
  if (XW %in% names(GenoProbXW)){
    return(GenoProbXW[XW])
  }else{
    return(0)
  }	
  
}

