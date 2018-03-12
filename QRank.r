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



plot_Qrank <- function(formula, testvar, data,tau=seq(0.02, 0.98, by=0.02), linequantile=c(0.25,0.5,0.75), subjectID=NULL, newdata=NULL,main='',ylab='',...){
    # check data
    if(!is.data.frame(data)){
        stop('data needs to be a data frame')
    }
    # check whether testvar is a factor itself
    # check whether fomula contains 'factor' whilc testvar is factor
    if (class(data[,testvar])=='factor'){
        if (!any(grepl(paste0('factor\\(',testvar,'\\)'),as.character(formula)))){
            stop(sprintf('Please put \'factor(%s)\' in the formula part !\n Although %s is a factor itself, you still need to put \'factor(%s)\' in the formula.\n ',testvar,testvar,testvar))
        }
    }
    # remove subjects that contains Na
    data <- data[rowSums(is.na(data))==0,]
    
    # if subjectID exists, do some checking
    if (!is.null(subjectID)){
        # check rownames first
        if (length(data[,1])!=length(unique(rownames(data)))){
            stop('Row names have to be unique when subjectID is assigned !')
        }
        # validate subjectID 
        validid <- sapply(as.character(subjectID),function(x){!any(rownames(data)==x)})
        if (any(validid)){
            cat(sprintf('SubjectID = %s does not exist in row names !\n',subjectID[validid]))
            cat(sprintf('Note: these subjectID might be absent in row names OR these involved subjects might contain Na values\n'))
            warning('Invalid subjectID has been dropped')
            subjectID <- subjectID[!validid]
        }
    }
    
    # generate Dx for regression
    Dx <- model.matrix(formula,data)
    Dx <- as.data.frame(Dx)
    ftm=suppressWarnings(rq(formula ,data=data, tau= tau))
    ft <- ftm$coefficients
    
    # no new data in 
    if (is.null(subjectID)&&is.null(newdata)){
        # treat differently for categorical or continuous test var 
        # categorical test var
        if ((any(grepl(paste0('factor\\(',testvar,'\\)'),as.character(formula))))){
            testindex <- grep(paste0('factor\\(',testvar,'\\)'),colnames(Dx))
            ngroup <- length(testindex)
            coln <- colnames(Dx)
            Dx=colMeans(Dx)
            
            Dx=matrix(rep(Dx,ngroup+1),nrow=ngroup+1,byrow = T)
            colnames(Dx) <- coln
            Dx[,testindex] <- rbind(rep(0,ngroup),diag(1,ngroup))
            
            pred <- Dx %*% ft
            
            step=apply(cbind(pred[,1], pred), 1, function(x) rearrange(stepfun(tau, x)))
            maxy <- max(sapply(step,function(x){ environment(x)$'yright'}))
            miny <- min(sapply(step,function(x){ environment(x)$'yleft'}))
            
            plot(step[[1]], ylab=ylab,ylim=c(miny,maxy), xlab=expression(tau), main=main,lty=1,...)
            
            for (pp in 2:length(step)){
                lines(step[[pp]],  lty=pp,col='grey55')
            }
            lines(step[[1]],lty=1)
            legend("bottomright",  c(paste0("factor(",testvar,")ref"),coln[testindex]), lty= 1:length(step), bty = "n")
            # for continuous test var  
        }else{
            if (length(unique(Dx[,testvar]))<15){
                dquan <- quantile(unique(Dx[,testvar]),probs = linequantile,type = 6) 
            }
            else{
                dquan <- quantile(Dx[,testvar],probs =linequantile) 
            }
            coln <- colnames(Dx)
            Dx=colMeans(Dx)
            lquan <- length(linequantile)
            Dx=matrix(rep(Dx,lquan),nrow=lquan,byrow = T)
            colnames(Dx) <- coln
            Dx[,testvar] <- dquan
            
            pred <- Dx %*% ft
            
            step=apply(cbind(pred[,1], pred), 1, function(x) rearrange(stepfun(tau, x)))
            maxy <- max(sapply(step,function(x){ environment(x)$'yright'}))
            miny <- min(sapply(step,function(x){ environment(x)$'yleft'}))
            
            boldline <- round(median(1:length(step)))
            plot(step[[boldline]], ylab=ylab,ylim=c(miny,maxy), xlab=expression(tau), main=main,lty=1,...)
            ltys <- NULL
            for (pp in (1:length(step))){
                newty <- pp+1
                if (pp>=boldline){
                    if (pp==boldline){
                        newty <- 1
                    }else{
                        newty <- pp
                        lines(step[[pp]],  lty=newty,col='grey55')
                    }
                }else{
                    lines(step[[pp]],  lty=newty,col='grey55')
                }
                ltys <- c(ltys,newty)
            }
            lines(step[[boldline]],lty=1)
            legend("bottomright", paste0(rep(paste0(testvar,' = '),3),format(Dx[,testvar],digits =5)), lty= ltys, bty = "n")
        }
        # with new data in 
    }else{
        if (is.null(newdata)){
            newdata <- sapply(as.character(subjectID),function(x){data[rownames(data)==x,]})
            newdata <- data.matrix(as.data.frame(newdata))
            newdata <- t(newdata)
            mylegend <- paste0('Subject: ',as.character(subjectID))
        }else{
            mylegend <- paste0('Subject: ',1:length(newdata[,1]))
        }
        pred <- predict(ftm,newdata = as.data.frame(newdata))
        step=apply(cbind(pred[,1], pred), 1, function(x) rearrange(stepfun(tau, x)))
        maxy <- max(sapply(step,function(x){ environment(x)$'yright'}))
        miny <- min(sapply(step,function(x){ environment(x)$'yleft'}))
        
        plot(step[[1]], ylab=ylab,ylim=c(miny,maxy), xlab=expression(tau), main=main,lty=1,...)
        if (length(step)>1){
            for (pp in (2:length(step))){
                lines(step[[pp]],  lty=pp)
            }
        }
        legend("bottomright", mylegend, lty= 1:length(step), bty = "n")
    }
}





################

## function
QRank = function(formula, data=NULL, tau=c(0.25,0.5,0.75)) {
    if (is.null(data)||(!is.data.frame(data))){
        stop('Data cannot be empty, and should be a data frame')
    }
    if (max(tau)>1||min(tau)<0){
        warning('tau should be in the interval [0,1]. Use tau=c(0.25,0.5,0.75) instead.')
        tau=c(0.25,0.5,0.75)
    }
    
    ltau = length(tau)
    y <- as.matrix(data[,all.vars(formula)[1]])
    Dx <- model.matrix(formula,data)
    Dx <- as.data.frame(Dx)
    
    Res <- vector('list',length(Dx))
    
    for (kk in 1:length(Dx)){
        if (kk==1&&colnames(Dx)[1]=="(Intercept)"){next}
        testvar <- colnames(Dx)[kk]
        x = as.matrix(Dx[,kk])
        zz = as.matrix(Dx[,-kk])
        
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
                      testvar = x,testvarname=testvar, y = y, cov = zz,tau = tau)
        Res[[kk]] <- result
    }
    if (colnames(Dx)[1]=="(Intercept)"){Res <- Res[-1]}
    class(Res) = "QRank"
    return(Res)
}	

print.QRank = function(x,...) {
    cat(sprintf('Test.Variable  Composite.pvalue  Quantile.Specific.pvalue\n'))
    cat(sprintf('                         '),format(x[[1]]$tau,width = 12),sprintf('\n'))
    for (i in 1:length(x)){
        cat(sprintf('% -15s',x[[i]]$testvarname))
        cat(format(c(x[[i]]$composite.pvalue,x[[i]]$quantile.specific.pvalue),digits=5,width = 12,justify = 'centre'),sprintf('\n'))
    }
    
}

# I made No change for the two functions below!!!!
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





