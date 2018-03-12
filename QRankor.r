# ---------------------- Quantile Rank-score Based Test (QRBT) ------------------------ # 
#                                                                                       #
#  This code serves as the function to conduct Quantile Rank-score based test (QRank)   #
#  for the continuous traits. It can be used to identify expression quantitative        #
#  trait loci (eQTLs) that are associated with the conditional quantile functions       # 
#  of gene expression.                                                                 #
#  Note: the quantile regression package "quantreg" is required. Users can use the     #
#  following codes to obtain "quantreg" package.        	   			   		       #
#  > install.packages("quantreg")                                                      #
#  > library(quantreg)                  		                        			   #
#                                                                                      #
#  *function* 				  				  					      	    	       		  
#	QRank(gene, snp, cov, tau): obtain the p-value of the association between the gene  #
#	expression and the genetic variants adjusting for the covariates under quantile    #
#	rank-score based inference.                                                        #
#                                                                                      #
#  *input*									    			   			   		       #
#	-gene: the vector of gene expression levels. No requirement for distribution	   #
#	-snp: the vector of genetic variant.  	    					   			 	   #
#	-cov: the vector or matrix of covariates. cov=NULL by default. Work for both cts and# 
#          binary covariates.                                                 #
# -tau: the quantile level to be estimated. Multiple taus can be chosen.             # 	
#                          e.g. tau=c( 0.25, 0.5, 0.75)                     #
#                                                                                      #
#  *output*															                   #
#   - composite.pvalue: provides a signle p-value for across all quantile levels under #
#        consideration, testing H0: the genetic variant and gene expression are not    #
#        assocaited                                                                    #
#   - quantile.specific.pvalue: provides p-values for each quantile level user         #
#         specified, testing H0: the genetic variant and gene expression are           #
#         not assocaited at this quantile level.                                        #
#                                                                                      #
#  *Example:*                                                   	           	       #
#  - create the sample dataset                                                         #                                                                       
# set.seed(123)                                                                        #
# n=300                                                                                #
# x=rbinom(n, 2, 0.2)                                                                  #
# y=rnorm(n, mean=0, sd=1)                                                             #
# z=cbind(rbinom(n, 1, 0.3), rnorm(n, mean=2, sd=2))                                   #
# taus=c( 0.25, 0.5, 0.75)                                                             #
# # - run the proposed QRank approach                                                   #
# QRank(gene=y, snp=x, cov=z, tau=taus)                                                 #
#  - output                                                                          #
# $composite.pvalue
# [1] 0.2241873
# 
# $quantile.specific.pvalue
# 0.25       0.5      0.75 
# 0.5452044 0.1821452 0.5938421                                                       # 
# ------------------------------------------------------------------------------------ # 


QRank = function(gene, snp, cov=NULL, tau) {
  ltau = length(tau)
  x = as.matrix(snp)
  y = as.matrix(gene)
  zz=cbind(rep(1, nrow(y)), cov)

  if (ncol(x)!=1) stop("Error: This program is designed for the analysis of a single SNP at a time.")
  if (nrow(x)!=nrow(y)) stop("Error: 'gene' and 'snp' must have the same length.")

  VN = matrix(0,nrow=ltau, ncol=ltau) 
  for(i in 1:ltau){
  	for(j in 1:ltau)
  	{
  		VN[i,j]= min(tau[i],tau[j])-tau[i]*tau[j]     
  	}
  }
	xstar = lm(x~zz-1)$residual
	SN = NULL
	for(i in 1:ltau) {
	  ranks = suppressWarnings(rq.fit.br(zz, y, tau = tau[i])$dual- (1 - tau[i]))
	  Sn = as.matrix(t(xstar)%*%(ranks))
	  SN = c(SN,Sn)
	} 

	VN2= matrix(outer(VN,t(xstar)%*% xstar,"*"), nrow = ltau)
	pvalue1=pchisq(SN^2/diag(VN2), 1, lower.tail=F);names(pvalue1)= tau
	
	e=solve(chol(VN2))
	SN2=t(e)%*%SN
	pvalue=pchisq(sum(SN2^2), ltau, lower.tail=F)
	
	result = list(composite.pvalue=pvalue, quantile.specific.pvalue=pvalue1,
	              snp = x, gene = y, cov = zz,tau = tau)
	class(result) = "QRank"
	return(result)
}	
	
print.QRank = function(x,...) {
  cat("\nComposite.pvalue:\n")
  print(as.numeric(x$composite.pvalue))
  cat("\nQuantile.specific.pvalue:\n")
  print(x$quantile.specific.pvalue)
}


########### heterogeniety index ###########

heter.QRank = function(object,newtaus = NULL) {
  y = object$gene
  x  = object$snp
  zz  = object$cov
  if (is.null(newtaus)) {newtaus = object$tau} ## default is taus used in QRank function
  else                  {newtaus = newtaus}
  ltau = length(newtaus)
  if (ltau == 1)        stop("Please select multiple taus")
  else
  {
  beta = coef(suppressWarnings(rq(y~x+zz-1,tau = newtaus)))
  heterogeneity=log(sd(beta)/abs(mean(beta)))
  heter.result = list(taus = newtaus, beta = beta, heter.index = heterogeneity)
  class(heter.result) = "QRank.heter"
  heter.result
  }
}



print.QRank.heter = function(x,...) {
  cat("Heterogeneity index:\n")
  print(x$heter.index)  
}





