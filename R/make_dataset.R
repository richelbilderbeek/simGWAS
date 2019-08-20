## ##' Simulate a dataset
## ##'
## ##' Simulates a dataset given a set of controls (without
## ##' phenotype) and a causal model.  This is used for testing code.
## ##' 
## ##' @export
## ##' @title Simulate a dataset
## ##' @param df control dataset of genotypes to sample from, without
## ##'     phenotype
## ##' @param N0 number of samples with Y=0
## ##' @param N1 number of samples with Y=1
## ##' @param CV which SNPs are causal?
## ##' @param PWgY0 the distribution of causal snps in controls
## ##' @param PWgY1 the distribution of causal snps in cases
## ##' @return simulated genotype data for N0 controls and N1 cases
## ##' @examples
## ##' ## fake genotype data for 3 snps
## ##' df=fake_df(100,5,maf=0.3)
## ##' data=make_dataset(df,N0=100,N1=100,CV="snp1",PWgY0=0.3,PWgY1=0.4)
## ##' summary(data)
## ##' @author Mary Fortune
## make_dataset<-function(df,N0,N1,CV,PWgY0,PWgY1){
##     m <- length(CV)
## 	#make sure we have genotypes in {0,1,2}
## 	if (isTRUE(all.equal(sort(unique(df[,1])),1:3))){
## 		df<-df-1
## 	}
## 	#compute the weights for sampling the cases
## 	ref_data_W<-df[,CV]
## 	Y1weights<-rep(0,nrow(df))
## 	#which order do the Ws come in the vectors of probabilities above?
## 	orderW<-apply(combinat::hcube(rep(3,m))-1, 1,function(x){paste(x,collapse="")})
## 	for (ii in 1:nrow(df)){
## 		if (length(CV)==1){
## 			Wvalue<-paste(ref_data_W[ii])
## 		}else{
## 			Wvalue<-paste(ref_data_W[ii,],collapse="")
## 		}
## 		whichw<-which(orderW==Wvalue)
## 		Y1weights[ii]<-PWgY1[whichw]/PWgY0[whichw]
## 	}	
## 	rows.index0<-sample(nrow(df), size=N0, replace = TRUE)
## 	rows.index1<-sample(nrow(df), size=N1, replace = TRUE, prob=Y1weights)
## 	dfnew_Y0<-cbind(rep(0,N0),df[rows.index0,])
## 	dfnew_Y1<-cbind(rep(1,N1),df[rows.index1,])
## 	colnames(dfnew_Y0)[1]<-"Y"
## 	colnames(dfnew_Y1)[1]<-"Y"
## 	dfnew<-rbind(dfnew_Y0,dfnew_Y1)
## 	return(list(dfnew,c(rows.index0,rows.index1)))
## }

##' create a fake haplotype frequency dataset
##'
##' no attempt is made at biological realism, this is purely for testing code
##' @title fake haplotype frequencies
##' @param nhaps number of haplotypes
##' @param nsnps number of snps
##' @return data.frame of 1,2, nhaps x nsamples + frequency column
##' @export
##' @author Chris Wallace
##' @examples
##' freq=fake_freq(nhaps=100,nsnps=5)
##' dim(freq)
##' head(freq)
fake_freq=function(nhaps=100,nsnps=10) {
    lag <- min(5,nsnps) # genotypes are correlated between neighbouring variants
    maf <- runif(nsnps+lag,0.05,0.5) # common SNPs
    laghaps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps,1,f)))
    haps <- laghaps[,1:nsnps]
    for(j in 1:lag) 
        haps <- haps + laghaps[,(1:nsnps)+j]
    haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))
    
    snps <- colnames(haps) <- paste0("s",1:nsnps)
    freq <- as.data.frame(haps+1)
    freq$Probability <- 1/nrow(freq)
    freq
}

## ##' @title Simulate a dataset given a set of controls (with phenotype) and a causal model
## ##' @param df dataset of genotypes to sample from, with phenotype
## ##' @param N0 number of samples with Y=0
## ##' @param N1 number of samples with Y=1
## ##' @param CV which SNPs are causal?
## ##' @param PWgY0 the distribution of causal snps in controls
## ##' @param PWgY1 the distribution of causal snps in cases
## ##' @return simulated genotype data for N0 controls and N1 cases
## ##' @author Mary Fortune
## make_dataset_pheno<-function(df,N0,N1,CV,PWgY0,PWgY1){
## 	#make sure we have genotypes in {0,1,2}
## 	if (isTRUE(all.equal(sort(unique(df[,2])),1:3))){
## 		df[,-1]<-df[,-1]-1
## 	}
## 	#compute the weights for sampling the cases
## 	ref_data_W<-df[,CV+1]
## 	Y1weights<-rep(0,nrow(df))
## 	#which order do the Ws come in the vectors of probabilities above?
## 	orderW<-apply(combinat::hcube(rep(3,m))-1, 1,function(x){paste(x,collapse="")})
## 	for (ii in 1:nrow(df)){
## 		if (length(CV)==1){
## 			Wvalue<-paste(ref_data_W[ii])
## 		}else{
## 			Wvalue<-paste(ref_data_W[ii,],collapse="")
## 		}
## 		whichw<-which(orderW==Wvalue)
## 		Y1weights[ii]<-PWgY1[whichw]/PWgY0[whichw]
## 	}	
## 	rows.index0<-sample(nrow(df), size=N0, replace = TRUE)
## 	rows.index1<-sample(nrow(df), size=N1, replace = TRUE, prob=Y1weights)
## 	dfnew_Y0<-cbind(rep(0,N0),df[rows.index0,-1])
## 	dfnew_Y1<-cbind(rep(1,N1),df[rows.index1,-1])
## 	colnames(dfnew_Y0)[1]<-"Y"
## 	colnames(dfnew_Y1)[1]<-"Y"
## 	dfnew<-rbind(dfnew_Y0,dfnew_Y1)
## 	return(list(dfnew,c(rows.index0,rows.index1)))
## }

