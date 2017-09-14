NBsimSingleCell <- function(dataset, group, nTags = 10000, nlibs = length(group), lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1, pUp=.5, foldDiff=3, verbose=TRUE, ind=NULL, params=NULL, cpm="AveLogCPM", max.dispersion=400, min.dispersion=0.1)
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
    Lambda <- Lambda/sum(Lambda) #normalize so they all sum to 1
    propZeroGene <- propZeroGene[id_r]    
    Dispersion <- dispersion[id_r]
    Dispersion[Dispersion>max.dispersion] = max.dispersion
    Dispersion[Dispersion<min.dispersion] = min.dispersion
    ####
    nonZeroCounts <- propZeroGene*nlibs
    dispSqueeze <- limma::squeezeVar(Dispersion, df=nonZeroCounts)
    Dispersion <- dispSqueeze$var.post
    ###
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
    zeroProbMatLink = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="link"), byrow=FALSE, ncol=nlibs)
    meanDiff = rowMeans(zeroProbMatLink)-logit(propZeroGene)
    zeroProbMat = expit(sweep(zeroProbMatLink,1,meanDiff,"-"))

    ## simulate negative binomial counts
    mu=sweep(Lambda,2,lib.size,"*")
    mu[mu<0.2] = 1
    adjustment = zeroProbMat*mu
    mu=mu+adjustment
    counts = matrix(rnbinom(n=nTags*nlibs, mu=mu, size=1/Dispersion), nrow=nTags, ncol=nlibs, byrow=FALSE)
    zeroProbNegBin = matrix(dnbinom(0, mu=mu, size=1/Dispersion), nrow=nTags, ncol=nlibs, byrow=FALSE)
    expectedZeroProbablityNegBinomial = rowMeans(zeroProbNegBin)

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



