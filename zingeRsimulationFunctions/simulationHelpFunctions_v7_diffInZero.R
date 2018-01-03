library(edgeR)
library(mgcv)

getExprFraction4 = function(counts, offset){
    #countsModel = counts[counts>0]
    #offsetModel = offset[counts>0]
    countsModel = counts
    offsetModel = offset
  sum(countsModel)/sum(offsetModel)
}

getPhiMoMPositive4 = function(counts, lambda, offset){
    #countsModel = counts[counts>0]
    #offsetModel = offset[counts>0]
    countsModel = counts
    offsetModel = offset
    mu=lambda*offsetModel
    phi = (sum(countsModel^2) - sum(mu^2) - sum(mu)) / sum(mu^2)
    return(phi)
}

reEstimateExprFraction4 = function(counts, offset, lambda, phi){
    countsModel = counts[counts>0]
    offsetModel = offset[counts>0]
    mu=lambda*offsetModel
  sum(countsModel*(1-dnbinom(0,mu=mu,size=1/phi)))/sum(offsetModel)
}

reEstimatePhiMoM4 = function(counts, lambda, offset, phi){
    countsModel = counts[counts>0]
    offsetModel = offset[counts>0]
    mu=lambda*offsetModel
    phi = (sum(countsModel^2 * (1-dnbinom(0,mu=mu,size=1/phi))) - sum(mu^2) - sum(mu)) / sum(mu^2)
    return(phi)
}



getDatasetMoMPositive = function(counts, drop.extreme.dispersion = FALSE, cpm= "AveLogCPM", MoMIter=10){

  #### estimate lambda and overdispersion based on ZTNB.
	d <- DGEList(counts)
	#cp <- cpm(d,normalized.lib.sizes=TRUE)
	dFiltered=d
	dFiltered <- edgeR::calcNormFactors(dFiltered)
  dFiltered$AveLogCPM <- aveLogCPM(dFiltered)
  ## estimate
  lambdaMoM=apply(dFiltered$counts,1,function(x) getExprFraction4(counts=x, offset=colSums(dFiltered$counts)))
  dispMoM = vector(length=nrow(dFiltered$counts))
  for(i in 1:nrow(dFiltered$counts)) dispMoM[i] = getPhiMoMPositive4(counts=dFiltered$counts[i,], offset=colSums(dFiltered$counts), lambda=lambdaMoM[i])
  dispMoM[dispMoM<0] = 1e-3

   ## iterative estimation.
   for(j in 1:MoMIter){
    message(paste0("iteration ",j," in ",MoMIter))
   for(i in 1:nrow(dFiltered$counts)) lambdaMoM[i] = reEstimateExprFraction4(counts=dFiltered$counts[i,], offset=colSums(dFiltered$counts), phi=dispMoM[i], lambda=lambdaMoM[i])
   for(i in 1:nrow(dFiltered$counts)) dispMoM[i] = reEstimatePhiMoM4(counts=dFiltered$counts[i,], offset=colSums(dFiltered$counts), lambda=lambdaMoM[i], phi=dispMoM[i])
   dispMoM[dispMoM<0] = 1e-3 #set negative dispersions very low => Poisson-like.
   }

  ## assume convergence
  params=cbind(dispMoM,lambdaMoM)
	rmRows = which(params[,2]>1) #impossibly high lambda
	rmRows2 = which(params[,2]==0) #zero lambda
	naRows = which(apply(params,1, function(row) any(is.na(row)))) #not fitted
	nonZeroDispRows = which(params[,1]<0 | params[,1]==0) #negative dispersion
	throwRows = c(rmRows,rmRows2,naRows,nonZeroDispRows)
  if(length(throwRows)>0) params = params[-throwRows,]

	### estimate logistic GAM P(zero) ~ s(aveLogCPM)*logLibSize
	### use unfiltered data for this model.
  require(mgcv)
	propZero = colMeans(counts==0)
	propZeroGene = rowMeans(counts==0)
	d <- DGEList(counts)
	d <- edgeR::calcNormFactors(d)
	avCpm <- aveLogCPM(d, normalized.lib.sizes=FALSE)
	cpmHist = hist(avCpm, breaks=150, plot=FALSE)
    	breaks = cpmHist$breaks
    	mids = cpmHist$mids
    	midsHlp=rep(mids,ncol(d$counts))
	logLibSize = log(colSums(counts))
    	logLibHlp=rep(logLibSize,each=length(mids))
	binHlp=sapply(breaks[-length(breaks)],function(x) avCpm>x)
  	binId=apply(binHlp,1,function(x) max(which(x)))
	nonNullCounts = t(sapply(1:length(mids), function(bin){
			    binRows <- binId==bin
			    if(sum(binRows)==0) rep(0,ncol(counts)) else
			    if(sum(binRows)==1) (counts[which(binRows),]!=0)*1 else
				colSums(counts[which(binRows),]!=0)
	    }))
	nullCounts = t(sapply(1:length(mids), function(bin){
		    	binRows <- binId==bin
		    	if(sum(binRows)==0) rep(0,ncol(counts)) else
		    	if(sum(binRows)==1) (counts[which(binRows),]==0)*1 else
			    colSums(counts[which(binRows),]==0)
	    }))
	expectCounts=cbind(c(nullCounts),c(nonNullCounts))
	#zeroFit=mgcv::gam(expectCounts~s(midsHlp)+logLibHlp,family=binomial)
	zeroFit=gam(expectCounts~s(midsHlp,by=logLibHlp),family=binomial)

	### drop extreme dispersions
  dFiltered$AveLogCPM <- aveLogCPM(dFiltered, normalized.lib.sizes=FALSE)
	if(length(throwRows)>0) dFiltered$AveLogCPM <- dFiltered$AveLogCPM[-throwRows]
	if(length(throwRows)>0) propZeroGene = propZeroGene[-throwRows]
	params=data.frame(dispersion=params[,1], lambda=params[,2], aveLogCpm=dFiltered$AveLogCPM, propZeroGene=propZeroGene)
	dispersion <- params$dispersion
	AveLogCPM <- params$aveLogCpm
	lambda <- params$lambda
	propZeroGene <- params$propZeroGene

	if(is.numeric(drop.extreme.dispersion))
	{
		bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE, na.rm=TRUE)
		ids <- dispersion <= bad
		AveLogCPM <- AveLogCPM[ids]
		dispersion <- dispersion[ids]
		lambda <- lambda[ids]
		propZeroGene <- propZeroGene[ids]
		params <- params[ids,]
		dFiltered <- dFiltered[ids,]
	}
	#lambda=lambda/sum(lambda) #make sure they sum to 1
	dataset.AveLogCPM <- AveLogCPM
	dataset.dispersion <- dispersion
	dataset.lambda <- lambda
	dataset.propZeroGene <- propZeroGene
	dataset.lib.size <- d$samples$lib.size
	dataset.nTags <- nrow(d)
	list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags, dataset.propZeroFit=zeroFit, dataset.lambda=lambda, dataset.propZeroGene=propZeroGene, dataset.breaks = breaks, dataset.cpm=cpm)
}


NBsimSingleCell <- function(dataset, group, nTags = 10000, nlibs = length(group), lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1, pUp=.5, foldDiff=3, verbose=TRUE, ind=NULL, params=NULL, cpm="AveLogCPM", max.dispersion=400, min.dispersion=0.1, normalizeLambda=FALSE)
{
  require(edgeR)
  group = as.factor(group)
  expit=function(x) exp(x)/(1+exp(x))
  logit=function(x) log(x/(1-x))

  sample.fun <- function(object)
  {
    nlibs <- object$nlibs
    nTags <- object$nTags
    AveLogCPM <-object$dataset$dataset.AveLogCPM
    dispersion <- object$dataset$dataset.dispersion
    lambda <- object$dataset$dataset.lambda
    #lambda <- (2^AveLogCPM)/1e6
    propZeroGene <- dat$dataset$dataset.propZeroGene
    id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
    object$AveLogCPM <- AveLogCPM[id_r]
    Lambda <- lambda[id_r]
    if(normalizeLambda) Lambda <- Lambda/sum(Lambda) #normalize so they all sum to 1
    Dispersion <- dispersion[id_r]
    Dispersion[Dispersion>max.dispersion] = max.dispersion
    Dispersion[Dispersion<min.dispersion] = min.dispersion
    propZeroGene <- propZeroGene[id_r]
    Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
    object$Lambda <- Lambda
    Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
    object$Dispersion <- Dispersion
    object$propZeroGene <- propZeroGene
    object
  }
  diff.fun <- function(object)
  {
    group <- object$group
    pUp <-  object$pUp
    foldDiff <- object$foldDiff
    Lambda <- object$Lambda
    nTags <- object$nTags
    g <- group == levels(group)[1]
    #AveLogCPM = expandAsMatrix(object$AveLogCPM,dim=c(nTags, nlibs))
    if(length(ind)>0 & !all(foldDiff==1)) {
      fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
      Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
      Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir))
      object$Lambda <- Lambda
      object$indDE <- ind
      object$indNonDE <- (1:nTags)[-ind]
      foldDiff[fcDir==1] <- 1/foldDiff[fcDir==1]
      object$foldDiff <- foldDiff #group2 / group1
    }
    if(all(foldDiff==1)) object$indDE <- NA
    object
  }
  sim.fun <- function(object)
  {
    Lambda <- object$Lambda
    Dispersion <- object$Dispersion
    nTags <- object$nTags
    nlibs <- object$nlibs
    lib.size <- object$lib.size
    zeroFit <- dat$dataset$dataset.propZeroFit
    propZeroGene <- dat$propZeroGene
    propZeroGene[propZeroGene==1] <- 1-1e-4
    propZeroGene[propZeroGene==0] <- 1e-4
    design <- object$design
    avLogCpm <- object$AveLogCPM
    mids <- object$dataset$dataset.mids
    breaks <- object$dataset$dataset.breaks

    ## get matrix of zero probabilities
    libPredict=rep(log(lib.size),each=length(avLogCpm))
    cpmPredict=rep(avLogCpm,length(lib.size))
    zeroProbMatLink = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="link"), byrow=FALSE, ncol=nlibs, nrow=nTags)
    meanDiff = rowMeans(zeroProbMatLink)-logit(propZeroGene)
    zeroProbMat = expit(sweep(zeroProbMatLink,1,meanDiff,"-"))
    #zeroProbMat = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="response"), byrow=FALSE, ncol=nlibs)

    ## simulate negative binomial counts
    # mu=sweep(Lambda,2,lib.size,"*")
    # mu[mu<0.1] = 0.1
    # #adjustment = zeroProbMat*mu
    # #mu=mu+adjustment
    # counts = matrix(rnbinom(n=nTags*nlibs, mu=mu, size=1/Dispersion), nrow=nTags, ncol=nlibs, byrow=FALSE)
    # zeroProbNegBin = matrix(dnbinom(0, mu=mu, size=1/Dispersion), nrow=nTags, ncol=nlibs, byrow=FALSE)
    # expectedZeroProbablityNegBinomial = rowMeans(zeroProbNegBin)

    ## simulate negative binomial counts
    mu=sweep(Lambda,2,lib.size,"*")
    zeroProbNegBin = matrix(dnbinom(0, mu=mu, size=1/Dispersion), nrow=nTags, ncol=nlibs, byrow=FALSE)
    expectedZeroProbablityNegBinomial = rowMeans(zeroProbNegBin)
    dropoutGenes = expectedZeroProbablityNegBinomial < rowMeans(zeroProbMat)
    adjustment = zeroProbMat*mu
    mu[dropoutGenes,]=mu[dropoutGenes,]+adjustment[dropoutGenes,]
    mu[mu<0.1] = 0.1
    counts = matrix(rnbinom(n=nTags*nlibs, mu=mu, size=1/Dispersion), nrow=nTags, ncol=nlibs, byrow=FALSE)

    ## calculate dropouts
    dropoutGenes = expectedZeroProbablityNegBinomial < rowMeans(zeroProbMat)
    message(paste0("Adding extra zeros w.r.t. NB for ",sum(dropoutGenes)," genes"))
    #dropout matrix is 0 for dropout.
    dropoutMatrix = 1-matrix(rbinom(n=nTags*nlibs, size=1, prob=zeroProbMat), nrow=nTags, ncol=nlibs, byrow=FALSE)
    dropoutMatrix[!dropoutGenes,] = 1
    #avoid all dropout genes
    allDropoutId <- which(rowSums(dropoutMatrix)==0)
    while(length(allDropoutId)>0 ){
      dropoutMatrix[allDropoutId,] = 1-matrix(rbinom(n=length(allDropoutId)*nlibs, size=1, prob=zeroProbMat[allDropoutId,]), nrow=length(allDropoutId), ncol=nlibs, byrow=FALSE)
      allDropoutId <- which(rowSums(dropoutMatrix)==0)
    }
    #add dropouts
    dropoutMatrix[counts==0 & dropoutMatrix==0]=1 #if count already zero, it's not a dropout
    counts = counts*dropoutMatrix
    object$dropout = dropoutMatrix

    ## resample positive counts for features with all zero counts
    zeroCountsId <- which(rowSums(counts)==0)
    while(length(zeroCountsId)>0 ){
      counts[zeroCountsId,] = matrix(rnbinom(n=length(zeroCountsId)*nlibs, mu=mu[zeroCountsId,], size=1/Dispersion[zeroCountsId,]), nrow=length(zeroCountsId), ncol=nlibs, byrow=FALSE)
      counts[zeroCountsId,]=counts[zeroCountsId,]*dropoutMatrix[zeroCountsId,]
      zeroCountsId <- which(rowSums(counts)==0)
    }

    ## name features, return object.
    rownames(counts) <- paste("ids", 1:nTags, sep = "")
    colnames(counts) <- paste("sample",1:nlibs,sep="")
    object$counts <- counts
    object
  }

  if(verbose) message("Preparing dataset.\n")
  if(is.null(params)){
    dataset <- getDatasetZTNB(counts = dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
  } else {
    dataset <- params
  }
  dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pUp = pUp, foldDiff = foldDiff))
  if(cpm=="aCpm") dat$dataset$dataset.AveLogCPM = dat$dataset$dataset.aCpm


  if(is.null(dat$lib.size)){
    dat$lib.size <- sample(dataset$dataset.lib.size, nlibs, replace=TRUE)}
  if(is.null(nTags)) dat$nTags <- dat$dataset$dataset.nTags
  if(verbose) message("Sampling.\n")
  dat <- sample.fun(dat)
  if(verbose) message("Calculating differential expression.\n")
  dat <- diff.fun(dat)
  if(verbose) message("Simulating data.\n")
  dat <- sim.fun(dat)
  dat
}


getVersion <-
  function(x)
  {
    ## it is low-level function of pval ##
    x1 <- gsub("(\\_)(\\w+)", "", x)
    v <- unlist(lapply(x1, function(z) {options(warn = -1)
      desp<- packageDescription(z)
      if(length(desp) == 1)
        return("unknown")
      else desp$Version
    }))
    paste0(x,"_", v )
  }

rmVersion <-
  function(x)
  {
    ## it is low-level function of pval ##
    x1 <- strsplit(x, "\\_")
    x1 <- lapply(x1, function(x) x[-length(x)])
    sapply(x1, paste0, collapse = "_")

  }

odd <- function(x)
{
  ## it is low-level function of pval ##
  y <- seq(x)
  idx <- y %% 2 != 0
  x[idx]
}

mainShow <-
  function(count.type, count.name, group, pOutlier)
  {
    ## it is low-level function of pval ##
    pOutlier <- paste(100*pOutlier, "% ", "outliers", sep = "")
    group <- as.factor(group)
    group <- paste0(sum(group == levels(group)[1]), "vs", sum(group == levels(group)[2]))
    if(count.type == "counts")
      paste0("No outliers", "/", count.name, "/", group)
    else
      paste0(pOutlier, "/", count.type, "/",  count.name, "/", group)

  }


resetPar <- function() {
  ## this re-set args of par for plot ##
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}



edgeR.pfun <-
  function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
  {
    ## edgeR standard pipeline ##
    library(edgeR)
    d <- DGEList(counts = counts, group = group )
    d <- edgeR::calcNormFactors(d)
    design = model.matrix(~group)
    d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
    d <- estimateGLMTrendedDisp(d,design=design)
    d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
    f <- glmFit(d, design = design)
    lr <- glmLRT(f, coef=2)
    lfc <- lr$table$logFC
    pval = lr$table$PValue
    padj = p.adjust(pval, "BH")
    out = cbind(pval = pval, padj = padj, lfc = lfc)
    return(out)
  }


edgeRFiltered.pfun <-
  function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
  {
    ## edgeR standard pipeline ##
    library(edgeR) ; library(genefilter)
    d <- DGEList(counts = counts, group = group )
    d <- edgeR::calcNormFactors(d)
    design = model.matrix(~group)
    d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
    d <- estimateGLMTrendedDisp(d,design=design)
    d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
    f <- glmFit(d, design = design)
    lr <- glmLRT(f, coef=2)
    lfc <- lr$table$logFC
    pval = lr$table$PValue
    baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
    hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
    padj <- hlp$padj
    out = cbind(pval = pval, padj = padj, lfc = lfc)
    return(out)
  }

edgeRWeightedOldF.pfun <-
  function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL, weights=matrix(1,nrow=nrow(counts),ncol=ncol(counts)))
  {
    ## edgeR standard pipeline ##
    library(edgeR)
    d <- DGEList(counts = counts, group = group )
    d <- edgeR::calcNormFactors(d)
    design = model.matrix(~group)
    d$weights <- weights
    d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
    d <- estimateGLMTrendedDisp(d,design=design)
    d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
    edger.fit <- glmFit(d, design) #uses weights
    lr <- zinbwave::glmWeightedF(edger.fit,coef=2)
    pval = lr$table$PValue
    padj = p.adjust(pval, "BH")
    out = cbind(pval = pval, padj = padj)
    out[is.na(out)]=1
    return(out)
  }

edgeROldF.pfun <-
  function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
  {
    ## edgeR standard pipeline ##
    library(edgeR)
    d <- DGEList(counts = counts, group = group )
    d <- edgeR::calcNormFactors(d)
    design = model.matrix(~group)
    d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
    d <- estimateGLMTrendedDisp(d,design=design)
    d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
    edger.fit <- glmFit(d, design) #uses weights
    lr <- glmLRTOld(edger.fit,coef=2,test="F", ZI=FALSE)
    pval = lr$table$PValue
    padj = p.adjust(pval, "BH")
    lfc <- lr$table$logFC
    out = cbind(pval = pval, padj = padj, lfc=lfc)
    out[is.na(out)]=1
    return(out)
  }

pvalueAdjustment_kvdb <- function(baseMean, filter, pValue,
                                  theta, alpha=0.05, pAdjustMethod="BH") {
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


zingeREdgeROwnWeights.pfun=function(counts, group, design=NULL, mc.cores=2, niter=NULL, w=NULL){
  library(edgeR) ; library(genefilter)
  d <- DGEList(counts = counts, group = group )
  d <- edgeR::calcNormFactors(d)
  design = model.matrix(~ group)
  d$weights = w
  d=estimateDisp(d,design)
  edger.fit <- glmFit(d, design) #uses weights
  edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  lfc <- edger.lrt$table$logFC
  pval <- edger.lrt$table$PValue
  baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
  hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  padj <- hlp$padj
  out=cbind(pval,padj,lfc)
  return(out)
}


zingeRDESeq2OwnWeights.pfun <-
  function(counts, group, design = NULL, mc.cores = 4, niter=NULL, w=NULL)
  {
    ## implement DESeq2 ##
    library(DESeq2) ; library(genefilter)
    colData <- data.frame(group)
    dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
    colData(dse)$group <- as.factor(colData(dse)$group)
    zeroWeights = w
    dimnames(zeroWeights) = NULL
    assays(dse)[["weights"]] = zeroWeights
    dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
    dse = estimateDispersions(dse)
    dse=nbinomWaldTest(dse, betaPrior=TRUE)
    #dse <- DESeq(dse, betaPrior=TRUE)
    res <- results(dse)
    baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
    pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
    padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
    out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
    out
  }

getDataset <- function(counts, drop.extreme.dispersion = 0.1, drop.low.lambda = TRUE) {
  ## this function generates NB parameters from real dataset ##
  ## it is low-level function of NBsim ##
  d <- DGEList(counts)
  d <- edgeR::calcNormFactors(d)
  cp <- round(cpm(d,normalized.lib.sizes=TRUE),1)
  if(drop.low.lambda){
    d <- d[rowSums(cp>1) >= 2, ]
  }
  d$AveLogCPM <- log2(rowMeans(cpm(d, prior.count = 1e-5)))
  d <- estimateGLMCommonDisp(d)
  d <- estimateGLMTrendedDisp(d)
  d <- estimateGLMTagwiseDisp(d)
  dispersion <- d$tagwise.dispersion
  AveLogCPM <- d$AveLogCPM
  if(is.numeric(drop.extreme.dispersion))
  {
    bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE)
    ids <- dispersion <= bad
    AveLogCPM <- AveLogCPM[ids]
    dispersion <- dispersion[ids]
  }
  dataset.AveLogCPM <- AveLogCPM
  dataset.dispersion <- dispersion
  dataset.lib.size <- d$samples$lib.size
  dataset.nTags <- nrow(d)
  list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags)
}



NBsim <-
  function(dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 10, pDiff=.1, pUp=.5, foldDiff=3, name = NULL, save.file = FALSE, file = NULL, only.add.outlier = FALSE, verbose=TRUE, ind=NULL)

  {
    ## NBsim generate simulated count from the real dataset followed by the NB model ##
    require(edgeR)
    group = as.factor(group)

    sample.fun <- function(object)
    {
      ## it is low-level function of NBsim ##
      ## it samples from the real dataset ##
      nlibs <- object$nlibs
      nTags <- object$nTags
      AveLogCPM <-object$dataset$dataset.AveLogCPM
      dispersion <- object$dataset$dataset.dispersion

      id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
      object$AveLogCPM <- AveLogCPM[id_r] #added by Koen to use for adding zeroes
      Lambda <- 2^(AveLogCPM[id_r])
      Lambda <- Lambda/sum(Lambda)
      Dispersion <- dispersion[id_r]
      id_0<- Lambda == 0
      Lambda <- Lambda[!id_0]
      Dispersion <- Dispersion[!id_0]
      Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
      object$Lambda <- Lambda
      if(!is.na(fix.dispersion))
        Dispersion <- expandAsMatrix(fix.dispersion, dim = c(nTags, nlibs))
      else Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
      object$Dispersion <- Dispersion
      object

    }
    diff.fun <- function(object)
    {

      ## it is low-level function of NBsim ##
      ## it creates diff genes according to foldDiff ##
      group <- object$group
      pDiff <- object$pDiff
      pUp <-  object$pUp
      foldDiff <- object$foldDiff
      Lambda <- object$Lambda
      nTags <- object$nTags
      g <- group == levels(group)[1]
      ## added by Koen to specify DE index yourself and allows to specify foldDiff as matrix
      if(is.null(ind)) ind <- sample(nTags, floor(pDiff*nTags))
      ##
      if(length(ind)>0 & !mean(foldDiff==1)==1 ) {
        fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
        Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
        Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir))
        #Lambda <- t(t(Lambda)/colSums(Lambda))
        object$Lambda <- Lambda
        object$indDE <- ind
        object$indNonDE <- (1:nTags)[-ind]
        object$mask_DEup <- object$mask_DEdown <- object$mask_nonDE <- expandAsMatrix(FALSE, dim = dim(Lambda))
        object$mask_DEup[ind[fcDir == 1], g] <- TRUE
        object$mask_DEup[ind[fcDir == -1], !g] <- TRUE
        object$mask_DEdown[ind[fcDir == -1], g] <- TRUE
        object$mask_DEdown[ind[fcDir == 1], !g] <- TRUE
        object$mask_nonDE[-ind,] <- TRUE
        object$mask_DE <- object$mask_DEup | object$mask_DEdown}
      if(mean(foldDiff==1)==1 | pDiff == 0)
        object$indDE <- NA
      object
    }
    sim.fun <- function(object)
    {
      ## it is low-level function of NBsim ##
      ## it simulate counts using rnbinom ##
      Lambda <- object$Lambda
      Dispersion <- object$Dispersion
      nTags <- object$nTags
      nlibs <- object$nlibs
      lib.size <- object$lib.size
      counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs)
      rownames(counts) <- paste("ids", 1:nTags, sep = "")
      object$counts <- counts
      object
    }

    outlier.fun <- function(object, outlierMech, pOutlier, min.factor = 2, max.factor = 5)
    {
      ## it is low-level function of NBsim ##
      ## it makes outlier ##
      outlierMech <- match.arg(outlierMech, c("S", "M", "R"))
      dim <- dim(object$counts)
      outlier.factor <- function() runif(1, min.factor, max.factor)
      countAddOut <- object$counts
      LambdaAddOut <- object$Lambda
      DispersionAddOut <- object$Dispersion
      switch(outlierMech,
             S = {
               mask_outlier <- expandAsMatrix(FALSE, dim = dim)
               id_r <- which(runif(dim[1]) < pOutlier)
               id_c <- sample(dim[2], length(id_r), replace = TRUE)
               for(i in seq(id_r))
                 mask_outlier[id_r[i], id_c[i]] <- TRUE
               countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
             },
             R = {
               mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
               countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
             },

             M = {
               mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
               LambdaAddOut[mask_outlier] <- sapply(LambdaAddOut[mask_outlier], function(z) z*outlier.factor())
               countAddOut[mask_outlier] <- rnbinom(sum(mask_outlier), mu = t(t(LambdaAddOut)*object$lib.size)[mask_outlier], size = 1/DispersionAddOut[mask_outlier])
             }
      )
      if(!mean(object$foldDiff == 1)==1 & !pDiff == 0)
      {
        indDEupOutlier <- which(apply(object$mask_DEup & mask_outlier, 1, any))
        indDEdownOutlier <- which(apply(object$mask_DEdown & mask_outlier, 1, any))
        indDEnoOutlier <- which(apply((object$mask_DE & !mask_outlier) , 1, all))
        indNonDEOutlier <- which(apply(object$mask_nonDE & mask_outlier, 1, any))
        indNonDEnoOutlier <- which(apply((object$mask_nonDE & !mask_outlier) , 1, all))
        indDEbothOutlier <- NA
        o <- indDEupOutlier %in% indDEdownOutlier
        q <-  indDEdownOutlier %in% indDEupOutlier
        if(any(o))
        {
          indDEupOutlier <- indDEupOutlier[!o]
          indDEbothOutlier <- indDEupOutlier[o]
          indDEdownOutlier <- indDEdownOutlier[!q]
        }
      }
      else
      {
        indDEupOutlier <- indDEdownOutlier <- indDEnoOutlier <- indNonDEOutlier <- indNonDEnoOutlier <- indDEbothOutlier <- NA
      }
      out <- list(countAddOut = countAddOut, outlierMech = outlierMech, pOutlier = pOutlier, mask_outlier = mask_outlier, indDEupOutlier = indDEupOutlier,
                  indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier,
                  indNonDEnoOutlier = indNonDEnoOutlier, LambdaAddOut = LambdaAddOut, DispersionAddOut = DispersionAddOut)

    }

    calProb <- function(x, l) round(1 -(1 - x)^(1/l), 2) ## calculate probability to make sure all the outlierMech produce the same amount of outliers ##
    ##### hlp = as.matrix(islamFilt)

    if(verbose) message("Preparing dataset.\n")
    if(class(dataset) == "DGEList")
    {
      dat <- dataset
      dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
    }
    else if(is.character(dataset))
    {
      load(dataset)
      dat <- get(gsub("(\\.)(\\w+)", "", basename(dataset)))
      dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
    }
    else if(is.matrix(dataset))
    {
      if(is.null(name)) name <- deparse(substitute(dataset))
      dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
      dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))
    }
    else
      dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))

    if(!only.add.outlier)
    {
      if(is.null(lib.size)){
        dat$lib.size <- runif(nlibs, min = 0.7*median(dat$dataset$dataset.lib.size), max = 1.3*median(dat$dataset$dataset.lib.size))
        #propZeroFit=dat$dataset.propZeroFit

      }
      if(is.null(nTags))
        dat$nTags <- dat$dataset$dataset.nTags
      if(verbose) message("Sampling.\n")
      dat <- sample.fun(dat)
      if(verbose) message("Calculating differential expression.\n")
      dat <- diff.fun(dat)
      if(verbose) message("Simulating data.\n")
      dat <- sim.fun(dat)
    }
    if(add.outlier){
      outlierMech <- match.arg(outlierMech,  c("S", "R", "M"), several.ok = TRUE)
      if(length(pOutlier)== 1 & length(outlierMech) > 1 & any(outlierMech == "S"))
      {
        prob <- calProb(pOutlier, length(group))
        pOutlier <- rep(pOutlier, length = length(outlierMech))
        pOutlier[!outlierMech == "S"] <- prob
      }
      else if(!length(pOutlier) == length(outlierMech))
        stop("pOutlier is not equal to outlierMech")
      if(verbose) message("Adding outliers.\n")
      dat[outlierMech] <- mapply(function(x, y) outlier.fun(dat, outlierMech = x, pOutlier = y, min.factor = min.factor, max.factor = max.factor), x = outlierMech, y = pOutlier, SIMPLIFY = FALSE)
      dat$pOutlier <- pOutlier
    }
    if(save.file)
    {

      ## save file for shiny app ##
      if(verbose) message("Saving file.\n")
      if(is.null(file))
      { g <- paste0("g", sum(levels(group)[1] == group), "v", sum(levels(group)[2] == group))
      f <- paste0("f", foldDiff)
      if(add.outlier) o <- paste0("o", sprintf( "%02d",100*pOutlier[1L]))
      else o <- paste0("o", sprintf( "%02d", 0 ))
      file <- paste0(dat$name, "/", g, f, o, ".Rdata")
      dir.create(dat$name, showWarnings = FALSE)
      }
      filenm <- eval(gsub("(\\.)(\\w+)", "", basename(file)))
      assign(filenm, dat)
      save(list = filenm, file = file)
    }
    dat
  }





pval <-
  function(y, ...) ## evaluate DE methods ##
    UseMethod("pval")
pval.default <-
  function(y, group, design = NULL, method = "edgeR", mc.cores = 4, globalEnvir = FALSE, niter=NULL, ...)
  {
    ## evaluate DE methods ##
    ## return to a list of pvalue and runing time ##
    ## pvalue contains pvalue and p-adjust value ##
    gc(FALSE)
    time <- proc.time()
    group <- as.factor(group)
    if(globalEnvir) method <- paste0(method, ".pscript")
    else method <- paste0(method, ".pfun")
    p <- get(method)
    if(globalEnvir)
    {
      L <- list(counts = y, group = group, design = design, mc.cores = mc.cores, p = p)
      e <- list2env(L, envir = .GlobalEnv)
      pvalue <- with(e, eval(p))
      try(rm(list = names(L), envir = e),  silent = TRUE)
      try(rm(pGlobal, envir = e),  silent = TRUE)
    }
    else pvalue <- p(y, group, design, mc.cores, niter, ...)
    pvalue
    new.time <- proc.time()
    output <- list(pvalue = pvalue, time = new.time - time)
  }



pval.DGEList <-
  function(y, method, mc.cores = 4, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), save.file = FALSE, name = deparse(substitute(y)), count.type="counts", niter=NULL)
  {
    ## evaluate DE methods ##
    ## return to a DGEList including pvalue and other necessary indexs for re-analysis and plot ##
    library(parallel)
    counts = y$counts
    pOutlier = mask_outlier = indDEupOutlier = indDEdownOutlier = indDEbothOutlier = indDEnoOutlier = indNonDEOutlier = indNonDEnoOutlier = NA
    names(method) <- method
    group <- y$group
    design <- y$design
    is.parallel <- method %in% parallel.method
    is.globalEnvir <- method %in% globalEnvir.method
    id.re <- !(is.parallel|is.globalEnvir)
    reduced.method <- method[id.re]
    if(any(id.re)) output <-  parallel:::mclapply(reduced.method, function(x) pval(y = counts, group = group, design = design, method = x, niter=niter), mc.cores = mc.cores, mc.preschedule = FALSE) else output <- list()
    if(any(is.parallel))
    {
      for( i in names(method[is.parallel]))
        output[[i]] <- pval(y = counts, group = group, design = design, method = i, mc.cores = mc.cores, niter=niter)
    }
    if(any(is.globalEnvir))
    {
      for( i in names(method[is.globalEnvir]))
        output[[i]] <- pval(y = counts, group = group, design = design, method = i, mc.cores = mc.cores, globalEnvir = TRUE, niter=niter)
    }
    output <- output[method]
    padj <- lapply(output, function(x) x[["pvalue"]][, "padj"])
    pval <- lapply(output, function(x) x[["pvalue"]][, "pval"])
    lfc <- lapply(output, function(x) x[["pvalue"]][, "lfc"])
    time <- lapply(output, function(x) x[["time"]])
    output <- new("DGEList", list(pval = pval, padj = padj, lfc = lfc,  counts = counts, group = group, design = design, indDE = y$indDE, method = names(method), indDEupOutlier = indDEupOutlier, indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, indNonDEnoOutlier = indNonDEnoOutlier, time = time))
    output$main <- mainShow(count.name = y$name, group = group, pOutlier = pOutlier, count.type=count.type)
    output$methodVersion <- getVersion(method)
    output
  }

pval.character <-
  function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), save.file = FALSE, niter=NULL)
  {
    ## for shiny app ##
    fnm <- y
    load(y)
    name <- gsub("(\\.)(\\w+)", "", basename(y))
    y <- get(name)
    pval.DGEList(y = y, method = method, count.type = count.type, mc.cores = mc.cores, parallel.method = parallel.method, globalEnvir.method = globalEnvir.method, save.file = save.file, name = fnm, niter=niter)
  }
pval.FoldList <-
  function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), cut.computing = TRUE)
  {
    ## evaluate DE methods for FoldList ##
    library(parallel)
    count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
    if(count.type == "counts")
    {counts = y$counts
    pOutlier = mask_outlier = indDEupOutlier = indDEdownOutlier = indDEbothOutlier = indDEnoOutlier = indNonDEOutlier = indNonDEnoOutlier = NA}
    else
    {counts = y[[count.type]]$countAddOut
    pOutlier = y[[count.type]]$pOutlier
    mask_outlier = y[[count.type]]$mask_outlier
    indDEupOutlier = y[[count.type]]$indDEupOutlier
    indDEdownOutlier = y[[count.type]]$indDEdownOutlier
    indDEbothOutlier = y[[count.type]]$indDEbothOutlier
    indDEnoOutlier = y[[count.type]]$indDEnoOutlier
    indNonDEnoOutlier = y[[count.type]]$indNonDEnoOutlier
    indNonDEOutlier = y[[count.type]]$indNonDEOutlier}

    names(method) <- method
    group <- y$group
    design <- y$design
    is.parallel <- method %in% parallel.method
    is.globalEnvir <- method %in% globalEnvir.method
    id.re <- !(is.parallel|is.globalEnvir)
    reduced.method <- method[id.re]
    if(any(id.re)) output <- lapply(reduced.method, function(x) parallel:::mclapply(counts, function(w) pval(y = w, method = x, group = group, design = design), mc.cores = mc.cores))
    else output <- list()
    fold_seq <- fold_seq.keep <- y$fold_seq
    if(cut.computing) fold_seq.keep <- odd(fold_seq)
    if(any(is.parallel))
    {
      for(j in fold_seq)
      {
        is.keep <- j %in% fold_seq.keep
        for( i in names(method[is.parallel]))
        {
          if(any(is.keep)) output[[i]][[j]] <- pval(y = counts[[j]], group = group, design = design, method = i, mc.cores = mc.cores)
          else
          {
            output[[i]][[j]][["pavlue"]] <- cbind(pval = NA, padj = NA)
            output[[i]][[j]][["time"]] <- NA
          }
        }
      }
    }
    if(any(is.globalEnvir))
    {
      for(j in fold_seq)
      {
        is.keep <- j %in% fold_seq.keep
        for( i in names(method[is.globalEnvir]))
        {
          if(any(is.keep)) output[[i]][[j]] <- pval(y = counts[[j]], group = group, design = design, method = i, mc.cores = mc.cores, globalEnvir = TRUE)
          else
          {
            output[[i]][[j]][["pavlue"]] <- cbind(pval = NA, padj = NA)
            output[[i]][[j]][["time"]] <- NA
          }
        }
      }
    }
    output <- output[method]
    padj <- try(lapply(output, lapply, function(x) x[["pvalue"]][, "padj"]), silent = TRUE)
    pval <- try(lapply(output, lapply, function(x) x[["pvalue"]][, "pval"]), silent = TRUE)
    time <- try(lapply(output, lapply, function(x) x[["time"]]), silent = TRUE)
    output <- new("FoldList", list(fold_seq = y$fold_seq, pval = pval, padj = padj, counts = counts, count.type = count.type, group = group, design = design, indDE = y$indDE, method = names(method), indDEupOutlier = indDEupOutlier, indDEdownOutlier = indDEdownOutlier,indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, indNonDEnoOutlier = indNonDEnoOutlier, time = time))
    output$main <- mainShow(count.type = count.type, count.name = y$name, group = group, pOutlier = pOutlier)
    output$methodVersion <- getVersion(method)
    output
  }
getPvalVersion <- function(methodVersion, pout = "pval", count.type = "counts", datanm)
{
  ## for shiny app ##
  Type <- switch(count.type, counts = "b", S = "s", M = "m", R = "r")
  filenm <- paste0(pout, "_", Type, "_",  basename(datanm), "_", methodVersion, ".Rdata")
  load(paste0(dirname(datanm),"/", rmVersion(methodVersion),"/", filenm))
  get(methodVersion)
}
getPval <- function(dataset,methodVersion, count.type = c("counts", "S", "R", "M"))
{
  ## for shiny app ##
  load(dataset)
  datanm <- gsub("(\\.)(\\w+)", "", dataset)
  dat <- get(basename(datanm))
  count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
  if(count.type == "counts") Dat <- new("DGEList", dat)
  else
  {
    Dat <- new("DGEList", dat[[count.type]])
    Dat[["counts"]] <- Dat[["countAddOut"]]
  }
  Dat$method <- Dat$methodVersion <- methodVersion
  Dat$group = dat$group
  Dat$indDE = dat$indDE
  Dat$name = dat$name
  Dat$main <- mainShow(count.type = count.type, count.name = Dat$name, group = Dat$group, pOutlier = Dat$pOutlier)
  index <- c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier")
  names(index) <- index
  indDiff <- lapply(index, function(x) Dat[[x]])
  indDiff <- indDiff[!sapply(indDiff, is.null)]
  indDiff <- indDiff[!is.na(indDiff)]
  Dat$index <- names(indDiff)
  names(methodVersion) <- methodVersion
  Dat[["padj"]] <- lapply(methodVersion, getPvalVersion, pout = "padj", count.type = count.type, datanm = datanm)
  Dat[["pval"]] <- lapply(methodVersion, getPvalVersion, pout = "pval", count.type = count.type, datanm = datanm)
  Dat

}


resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

roPlot <-
  function(y, ...)
    UseMethod("roPlot")
## plot ROC curve ##

roPlot.default <-
  function(y, indDiff, returnData=FALSE, plot.max.fpr = 0.4, add = FALSE, cex.axis = 2, threshold = 0.05, col = 1, cex.threshold = 3, plot.max.tpr=1, ...)
  {
    ## plot ROC curve ##
    #old.par <- par(c("mar", "mgp", "cex.axis"))
    #par(mar=c(4,5,3,2))
    #par(mgp = c(2.6, 1, 0))
    #par(cex.axis = cex.axis)
    #on.exit(par(old.par))
    library(ROCR)
    if(any(is.na(y)))
    {
      y[is.na(y)] <- 1
    }
    y = 1 - y
    label <- as.factor(rep("nonDE", length(y)))
    levels(label) <- c("nonDE", "DE")
    label[indDiff] <- "DE"
    pred <- prediction(y, label, label.ordering = c("nonDE", "DE"))
    perf <- performance(pred, "tpr", "fpr")
    if(is.null(plot.max.fpr))
      plot.max.fpr <- 1
    plot(perf, xlim = c(0, plot.max.fpr), ylim=c(0,plot.max.tpr), col = col, add = add, ...)
    if(!is.null(threshold))
    {
      fpr <- approx(y = perf@x.values[[1]], x = perf@alpha.values[[1]], xout = 1- threshold)$y
      tpr <- approx(y = perf@y.values[[1]], x = perf@x.values[[1]], xout = fpr)$y
      points(x = fpr, y = tpr, pch = 4, col = col, cex = cex.threshold, ...)
    }
    #if(returnData) return(perf) #added by Koen Vdb
    #par(resetPar())
  }
roPlot.DGEList <-
  function(y, plot.max.fpr = 0.4, plot.max.tpr=1, pout = "padj", threshold = 0.05, selected.method = NULL, show.legend = TRUE, short.main = FALSE, col = NULL, lty = 1, lwd = 5, box.lwd = 1.5, cex.main = 2.5, cex.axis=2, cex.lab = 2.1, cex = 1.6, cex.threshold = 8)
  {
    ## plot ROC curve for DGEList ##
    library(ROCR)
    main <- y$main
    if(short.main)
    {
      main <- strsplit(main, "/")[[1]]
      main <- main[-c(length(main), length(main)-1)]
      main <- paste0(main, collapse = "/")

    }
    if(is.null(selected.method))
    {
      method <- y$method
      methodVersion <- y$methodVersion
    }
    else
    {
      method <- match.arg(selected.method, y$method, several.ok = TRUE)
      methodVersion <- y$methodVersion[match(method, y$method)]
    }
    pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
    pvalue <- y[[pout]][method]
    pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow", "steelblue", "salmon", "violetred4", "darkred", "skyblue4")
    if(is.null(col)) col <- pre.col[seq(method)]
    roPlot(y = pvalue[[1L]], indDiff = y$indDE, threshold = threshold, plot.max.fpr = plot.max.fpr, plot.max.tpr = plot.max.tpr, col = col[1L], main = main, lty = lty, lwd = lwd, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, cex.threshold = cex.threshold)
    if(is.list(pvalue)) mapply(function(u, v, w) roPlot(y = u, indDiff = y$indDE, threshold = threshold, plot.max.fpr = plot.max.fpr, plot.max.tpr = plot.max.tpr, col = v, lwd = lwd, lty = lty, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, cex.threshold = cex.threshold, add = TRUE), u = pvalue[-1L], v = col[-1L])
    if(show.legend) legend("bottomright", methodVersion, col = col, lty = lty, lwd = lwd, box.lwd = box.lwd, cex = cex)
    if(!is.null(threshold)&show.legend) legend("topleft", paste0("P_adj_value=", threshold), col = "black", pch = 4, lty = NA, lwd = lwd, box.lwd = box.lwd, cex = cex, pt.cex = 1.5*cex)
  }



fdPlot <-
  function(y, ...)
    UseMethod("fdPlot")
## plot False discovery number curve ##

fdPlot.default <-
  function(y, indDiff, add=FALSE, xlab="Number of genes selected",
           ylab="Number of false discoveries", lwd=4, type="l", ... )
  {

    ## plot False discovery number curve ##
    #	old.par <- par(c("mar", "mgp"))
    #	par(mar=c(4,5,3,2))
    #    par(mgp = c(2.6, 1, 0))
    #	on.exit(par(old.par))
    x <- 1:length(indDiff)
    o <- order(y)
    w <- !o[x] %in% indDiff
    y1 <- cumsum(w)
    matplot(x, y1, xlab=xlab, ylab=ylab, lwd=lwd, type=type, add=add, ...)
  }
fdPlot.DGEList <-
  function(y, pout = "padj", selected.method = NULL, short.main = FALSE, show.legend = TRUE, col = NULL, lty = 1, lwd = 5, box.lwd = 1.5, cex.main = 2.5, cex.axis=2, cex.lab = 2.1, cex = 1.3, xlim = NULL)
  {
    ## plot False discovery number curve for DGEList ##
    main <- y$main
    if(short.main)
    {
      main <- strsplit(main, "/")[[1]]
      main <- main[-c(length(main), length(main)-1)]
      main <- paste0(main, collapse = "/")
    }
    if(is.null(selected.method))
    {
      method <- y$method
      methodVersion <- y$methodVersion
    }
    else
    {
      method <- match.arg(selected.method, y$method, several.ok = TRUE)
      methodVersion <- y$methodVersion[match(method, y$method)]
    }
    pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
    pvalue <- y[[pout]][method]

    pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
    if(is.null(col)) col <- pre.col[seq(method)]
    fdPlot(y = pvalue[[1L]], indDiff = y$indDE,  log="y", col = col[1L], main = main, lty = lty, lwd = lwd, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, xlim = xlim)
    if(is.list(pvalue)) mapply(function(u, v, w) fdPlot(y = u, indDiff = y$indDE, log="y", col = v, lwd = lwd, lty = lty, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, add = TRUE), u = pvalue[-1L], v = col[-1L])
    if(show.legend) legend("bottomright", methodVersion, col = col, lty = lty, lwd = lwd, box.lwd = box.lwd, cex = cex)
  }

getPower <-
  function(y, ...)
    UseMethod("getPower")
## plot Power curve ##
getPower.default <-
  function(y, indDiff, threshold)
  {
    ## plot Power curve ##
    if(all(is.na(y)))
      NA
    else if(all(is.na(indDiff)))
      contable(score = y, indDiff = indDiff, threshold = threshold, output = "fpr")
    else
      contable(score = y, indDiff = indDiff, threshold = threshold, output = "power")
  }
getPower.DGEList <-
  function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE, plot = FALSE, col = NULL, cex.main = 2.5, cex.axis=2, cex.sub = 1.5, cex.lab = 2.1, ylim = NULL, ...)
  {
    ## plot Power curve for DGEList ##
    index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
    main <- y$main
    if(short.main)
    {
      main <- strsplit(main, "/")[[1]]
      main <- main[-c(length(main), length(main)-1)]
      main <- paste0(main, collapse = "/")

    }

    if(is.null(selected.method))
    {
      method <- y$method
      methodVersion <- y$methodVersion
    }
    else
    {
      method <- match.arg(selected.method, y$method, several.ok = TRUE)
      methodVersion <- y$methodVersion[match(method, y$method)]
    }
    pout <- match.arg(pout, c("pval", "padj"))
    pvalue <- y[[pout]][method]
    indDiff <- y[[index]]

    if(byAveLogCPM)
    {
      threshold <- rep(threshold, cutCPM)
      diff <- rep(FALSE, nrow(y$counts))
      diff[indDiff] <- TRUE
      d <- DGEList(counts = y$counts, group = y$group)
      d <- edgeR::calcNormFactors(d)
      AveLogCPM <- aveLogCPM(d)
      o <- order(AveLogCPM)
      l <- length(o)/cutCPM
      oo <- split(o, ceiling(seq_along(o)/l))
      AveLogCPM.list <- lapply(oo, function(x) AveLogCPM[x])
      nm <- lapply(AveLogCPM.list, function(x) round(range(x), 2))
      nm <- unlist(lapply(nm, function(x) paste0("(", paste0(x, collapse = ","), "]")))

      power <- list()
      power[["all"]] <- unlist(lapply(pvalue, getPower, indDiff = indDiff, threshold = threshold[1]))
      for( i in 1:length(oo))
      {
        diff_idx <- which(diff[oo[[i]]])
        pvalue_idx <- lapply(pvalue, function(x) x[oo[[i]]])
        power[[nm[i]]] <- unlist(lapply(pvalue_idx, getPower, indDiff = diff_idx, threshold = threshold[i]))
      }
      power <- do.call("rbind", power)
    }
    else
    {
      power <- unlist(lapply(pvalue, getPower, indDiff = indDiff, threshold = threshold))
    }
    pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
    if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$power <- power
    out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
    out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    if(plot)
      if(byAveLogCPM) matPlot(out, output = "power", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    else powBarPlot(out, col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    out
  }



getFDR <-
  function(y, ...)
    UseMethod("getFDR")
## plot Power curve ##
getFDR.default <-
  function(y, indDiff, threshold)
  {
    ## plot Power curve ##
    if(all(is.na(y)))
      NA
    else if(all(is.na(indDiff)))
      NA
    else
      contable(score = y, indDiff = indDiff, threshold = threshold, output = "fdr")
  }
getFDR.DGEList <-
  function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE, plot = FALSE, col = NULL, cex.main = 2.5, cex.axis=2, cex.sub = 1.5, cex.lab = 2.1, ylim = NULL, ...)
  {
    ## plot Power curve for DGEList ##
    index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
    main <- y$main
    if(short.main)
    {
      main <- strsplit(main, "/")[[1]]
      main <- main[-c(length(main), length(main)-1)]
      main <- paste0(main, collapse = "/")

    }

    if(is.null(selected.method))
    {
      method <- y$method
      methodVersion <- y$methodVersion
    }
    else
    {
      method <- match.arg(selected.method, y$method, several.ok = TRUE)
      methodVersion <- y$methodVersion[match(method, y$method)]
    }
    pout <- match.arg(pout, c("pval", "padj"))
    pvalue <- y[[pout]][method]
    indDiff <- y[[index]]

    if(byAveLogCPM)
    {
      threshold <- rep(threshold, cutCPM)
      diff <- rep(FALSE, nrow(y$counts))
      diff[indDiff] <- TRUE
      d <- DGEList(counts = y$counts, group = y$group)
      d <- edgeR::calcNormFactors(d)
      AveLogCPM <- aveLogCPM(d)
      o <- order(AveLogCPM)
      l <- length(o)/cutCPM
      oo <- split(o, ceiling(seq_along(o)/l))
      AveLogCPM.list <- lapply(oo, function(x) AveLogCPM[x])
      nm <- lapply(AveLogCPM.list, function(x) round(range(x), 2))
      nm <- unlist(lapply(nm, function(x) paste0("(", paste0(x, collapse = ","), "]")))

      fdr <- list()
      fdr[["all"]] <- unlist(lapply(pvalue, getFDR, indDiff = indDiff, threshold = threshold[1]))
      for( i in 1:length(oo))
      {
        diff_idx <- which(diff[oo[[i]]])
        pvalue_idx <- lapply(pvalue, function(x) x[oo[[i]]])
        fdr[[nm[i]]] <- unlist(lapply(pvalue_idx, getFDR, indDiff = diff_idx, threshold = threshold[i]))
      }
      fdr <- do.call("rbind", fdr)
    }
    else
    {
      fdr <- unlist(lapply(pvalue, getFDR, indDiff = indDiff, threshold = threshold))
    }
    pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
    if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$fdr <- fdr
    out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
    out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    if(plot)
      if(byAveLogCPM) matPlot(out, output = "fdr", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    else barPlot(out, output = "fdr", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    out
  }
