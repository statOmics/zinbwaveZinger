## got it from https://github.com/Bioconductor-mirror/edgeR/blob/release-3.0/R/glmfit.R


.pvalueAdjustment_kvdb <- function(baseMean, filter, pValue,
                                   theta, alpha=0.05, pAdjustMethod="BH") {
  ## this function has been adapted from the pValueAdjustment function from the DESeq2 Bioconductor package and was originally written by Mike Love.
  require(genefilter)
  # perform independent filtering
  if (missing(filter)) {
    filter <- baseMean
  }
  if (missing(theta)) {
    lowerQuantile <- mean(filter == 0)
    if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
    theta <- seq(lowerQuantile, upperQuantile, length=50)
  }

  # do filtering using genefilter
  stopifnot(length(theta) > 1)
  filtPadj <- filtered_p(filter=filter, test=pValue,
                         theta=theta, method=pAdjustMethod)
  numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
  # prevent over-aggressive filtering when all genes are null,
  # by requiring the max number of rejections is above a fitted curve.
  # If the max number of rejection is not greater than 10, then don't
  # perform independent filtering at all.
  lo.fit <- lowess(numRej ~ theta, f=1/5)
  if (max(numRej) <= 10) {
    j <- 1
  } else {
    residual <- if (all(numRej==0)) {
      0
    } else {
      numRej[numRej > 0] - lo.fit$y[numRej > 0]
    }
    thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
    j <- if (any(numRej > thresh)) {
      which(numRej > thresh)[1]
    } else {
      1
    }
  }
  padj <- filtPadj[, j, drop=TRUE]
  cutoffs <- quantile(filter, theta)
  filterThreshold <- cutoffs[j]
  filterNumRej <- data.frame(theta=theta, numRej=numRej)
  filterTheta <- theta[j]

  return(list(padj=padj, filterThreshold=filterThreshold, filterTheta=filterTheta, filterNumRej = filterNumRej, lo.fit=lo.fit, alpha=alpha))

}


independentFiltering <- function(object, filter, objectType=c("edgeR","limma")){
  if(objectType=="edgeR"){
    hlp <- .pvalueAdjustment_kvdb(filter=filter, pValue=object$table$PValue)
    object$table$padjFilter <- hlp$padj
    return(object)
  } else if(objectType=="limma"){
    hlp <- .pvalueAdjustment_kvdb(filter=filter, pValue=object$P.Value)
    object$padjFilter <- hlp$padj
    return(object)
  } else stop("objectType must be either one of 'edgeR' or 'limma'.")
}

glmWeightedF <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, test="F", ZI=TRUE, independentFiltering=TRUE, filter=NULL)
  ## original function obtained from https://github.com/Bioconductor-mirror/edgeR/blob/release-3.0/R/glmfit.R
  #	Tagwise likelihood ratio tests for DGEGLM
  #	Gordon Smyth, Davis McCarthy and Yunshun Chen.
  #	Created 1 July 2010.  Last modified 22 Nov 2013.
{
  #	Check glmfit
  if(!is(glmfit,"DGEGLM")) {
    if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
      stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
    }
    stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
  }
  if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
  nlibs <- ncol(glmfit)

  #	Check test
  test <- match.arg(test,c("F","f","chisq"))
  if(test=="f") test <- "F"

  #	Check design matrix
  design <- as.matrix(glmfit$design)
  nbeta <- ncol(design)
  if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
  coef.names <- colnames(design)

  #	Evaluate logFC for coef to be tested
  #	Note that contrast takes precedence over coef: if contrast is given
  #	then reform design matrix so that contrast of interest is last column.
  if(is.null(contrast)) {
    if(length(coef) > 1) coef <- unique(coef)
    if(is.character(coef)) {
      check.coef <- coef %in% colnames(design)
      if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
      coef.name <- coef
      coef <- match(coef, colnames(design))
    }
    else
      coef.name <- coef.names[coef]
    logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
  } else {
    contrast <- as.matrix(contrast)
    qrc <- qr(contrast)
    ncontrasts <- qrc$rank
    if(ncontrasts==0) stop("contrasts are all zero")
    coef <- 1:ncontrasts
    if(ncontrasts < ncol(contrast)) contrast <- contrast[,qrc$pivot[coef]]
    logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
    if(ncontrasts>1) {
      coef.name <- paste("LR test of",ncontrasts,"contrasts")
    } else {
      contrast <- drop(contrast)
      i <- contrast!=0
      coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
    }
    Dvec <- rep.int(1,nlibs)
    Dvec[coef] <- diag(qrc$qr)[coef]
    Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
    design <- design %*% Q
  }
  if(length(coef)==1) logFC <- as.vector(logFC)

  #	Null design matrix
  design0 <- design[,-coef,drop=FALSE]

  #	Null fit
  fit.null <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,prior.count=0)

  #	Likelihood ratio statistic
  LR <- fit.null$deviance - glmfit$deviance
  ### ADDED
  if(ZI) fit.null$df.residual <- rowSums(fit.null$weights)-ncol(design0)
  if(ZI) glmfit$df.residual <- rowSums(glmfit$weights)-ncol(design)
  ## END ADDED
  df.test <- fit.null$df.residual - glmfit$df.residual ## okay

  #	Chisquare or F-test
  # LRT.pvalue <- switch(test,
  #                      "F" = {
  #                        phi <- quantile(glmfit$dispersion,p=0.5)
  #                        mu <- quantile(glmfit$fitted.values,p=0.5)
  #                        gamma.prop <- (phi*mu/(1 + phi*mu))^2
  #                        prior.df <- glmfit$prior.df
  #                        if(is.null(prior.df)) prior.df <- 20
  #                        glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
  #                        pf(LR/df.test, df1=df.test, df2=glmfit$df.total, lower.tail = FALSE, log.p = FALSE)
  #                      },
  #                      "chisq" = pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)
  # )
  LRT.pvalue <- {
                         phi <- quantile(glmfit$dispersion,p=0.5)
                         mu <- quantile(glmfit$fitted.values,p=0.5)
                         gamma.prop <- (phi*mu/(1 + phi*mu))^2
                         prior.df <- glmfit$prior.df
                         if(is.null(prior.df)) prior.df <- 20
                         glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
                         pf(LR/df.test, df1=df.test, df2=glmfit$df.total, lower.tail = FALSE, log.p = FALSE)
                       }

  rn <- rownames(glmfit)
  if(is.null(rn))
    rn <- 1:nrow(glmfit)
  else
    rn <- make.unique(rn)
  tab <- data.frame(
    logFC=logFC,
    logCPM=glmfit$AveLogCPM,
    LR=LR,
    PValue=LRT.pvalue,
    row.names=rn
  )
  glmfit$counts <- NULL
  glmfit$table <- tab
  glmfit$comparison <- coef.name
  glmfit$df.test <- df.test
  res <- new("DGELRT",unclass(glmfit))
  if(independentFiltering){
    if(is.null(filter)) filter=rowMeans(glmfit$fitted.values) #aprrox. linear w\ basemean
    res <- independentFiltering(res,filter=filter, objectType="edgeR")
  } else return(res)
}




# glmLRTOld <- function(glmfit,coef=ncol(glmfit$design),contrast=NULL,test="chisq",ZI=TRUE)
# #	Tagwise likelihood ratio tests for DGEGLM
# #	Gordon Smyth, Davis McCarthy and Yunshun Chen.
# #	Created 1 July 2010.  Last modified 22 Nov 2013.
# {
# #	Check glmfit
# 	if(!is(glmfit,"DGEGLM")) {
# 		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
# 			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
# 		}
# 		stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
# 	}
# 	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
# 	nlibs <- ncol(glmfit)
#
# #	Check test
# 	test <- match.arg(test,c("F","f","chisq"))
# 	if(test=="f") test <- "F"
#
# #	Check design matrix
# 	design <- as.matrix(glmfit$design)
# 	nbeta <- ncol(design)
# 	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
# 	coef.names <- colnames(design)
#
# #	Evaluate logFC for coef to be tested
# #	Note that contrast takes precedence over coef: if contrast is given
# #	then reform design matrix so that contrast of interest is last column.
# 	if(is.null(contrast)) {
# 		if(length(coef) > 1) coef <- unique(coef)
# 		if(is.character(coef)) {
# 			check.coef <- coef %in% colnames(design)
# 			if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
# 			coef.name <- coef
# 			coef <- match(coef, colnames(design))
# 		}
# 		else
# 			coef.name <- coef.names[coef]
# 		logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
# 	} else {
# 		contrast <- as.matrix(contrast)
# 		qrc <- qr(contrast)
# 		ncontrasts <- qrc$rank
# 		if(ncontrasts==0) stop("contrasts are all zero")
# 		coef <- 1:ncontrasts
# 		if(ncontrasts < ncol(contrast)) contrast <- contrast[,qrc$pivot[coef]]
# 		logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
# 		if(ncontrasts>1) {
# 			coef.name <- paste("LR test of",ncontrasts,"contrasts")
# 		} else {
# 			contrast <- drop(contrast)
# 			i <- contrast!=0
# 			coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
# 		}
# 		Dvec <- rep.int(1,nlibs)
# 		Dvec[coef] <- diag(qrc$qr)[coef]
# 		Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
# 		design <- design %*% Q
# 	}
# 	if(length(coef)==1) logFC <- as.vector(logFC)
#
# #	Null design matrix
# 	design0 <- design[,-coef,drop=FALSE]
#
# #	Null fit
# 	fit.null <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,prior.count=0)
#
# #	Likelihood ratio statistic
# 	LR <- fit.null$deviance - glmfit$deviance
# 	### ADDED
# 	if(ZI) fit.null$df.residual <- rowSums(fit.null$weights)-ncol(design0)
# 	## END ADDED
# 	df.test <- fit.null$df.residual - glmfit$df.residual ## okay
#
# #	Chisquare or F-test
# 	LRT.pvalue <- switch(test,
# 		"F" = {
# 			phi <- quantile(glmfit$dispersion,p=0.5)
# 			mu <- quantile(glmfit$fitted.values,p=0.5)
# 			gamma.prop <- (phi*mu/(1 + phi*mu))^2
# 			prior.df <- glmfit$prior.df
# 			if(is.null(prior.df)) prior.df <- 20
# 			glmfit$df.total <- glmfit$df.residual + prior.df/gamma.prop
# 			pf(LR/df.test, df1=df.test, df2=glmfit$df.total, lower.tail = FALSE, log.p = FALSE)
# 		},
# 		"chisq" = pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)
# 	)
#
# 	rn <- rownames(glmfit)
# 	if(is.null(rn))
# 		rn <- 1:nrow(glmfit)
# 	else
# 		rn <- make.unique(rn)
# 	tab <- data.frame(
# 		logFC=logFC,
# 		logCPM=glmfit$AveLogCPM,
# 		LR=LR,
# 		PValue=LRT.pvalue,
# 		row.names=rn
# 	)
# 	glmfit$counts <- NULL
# 	glmfit$table <- tab
# 	glmfit$comparison <- coef.name
# 	glmfit$df.test <- df.test
# 	new("DGELRT",unclass(glmfit))
# }

