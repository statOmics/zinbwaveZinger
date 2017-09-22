getParamsZTNB <- function(counts, offset, design=NULL) {
  require(MASS)
  libSize=offset
  #fit a ZTNB model only on positive counts part
  countsModel = counts[counts>0]
  if(length(countsModel)<2) stop("Need at least two positive counts")
  libSizeModel = libSize[counts>0]
  if(is.null(design)){designFit=matrix(1,nrow=length(countsModel),ncol=1)} else {designFit=design[counts>0,]}
  fit=suppressWarnings(try(gamlss(formula=countsModel~-1+designFit+offset(log(libSizeModel)), family="NBIZeroTruncated", control=gamlss.control(trace=FALSE, n.cyc=300)),silent=TRUE))

  if(class(fit)[1]=="try-error") return(c(dispersion=NA, lambda=NA))
  #lambda=exp(mean(fit$mu.coefficients)) #geometric mean
  lambda=mean(countsModel/libSizeModel)
  dispersion=exp(fit$sigma.coefficients)
  return(c(dispersion=dispersion,lambda=lambda))
}



