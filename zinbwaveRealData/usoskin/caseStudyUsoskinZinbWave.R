library(Biobase)
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskin.RData")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")
usoskin=exprs(eset)
keep=rowSums(usoskin>0)>=10
usoskin=usoskin[keep,]
#cellTypeUsoskin=factor(pData(eset)[,7], levels=c("NF", "NP", "PEP", "TH"))
pickingSession=factor(pData(eset)[,2])
cellTypeAll=factor(pData(eset)[,"Level 3"],exclude=TRUE)



## P(zero) ~ library size
expit <- function(x) 1/(1+exp(-x))
## not including batch effect
par(mar=c(5,4.5,4,1))
logLib=log(colSums(usoskin))
pZero=colMeans(usoskin==0)
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m <- glm(pZero~logLib,family="binomial")
a <- coef(m)[1]
b <- coef(m)[2]
lines(x=sort(logLib),y=expit(a+b*sort(logLib)),lwd=2,col="steelblue")
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,pch=1,cex=1.5)

## including main batch effect
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib+pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,lty=1,cex=1.5)


## including batch effect as interaction with logLib
plot(x=logLib,y=pZero, xlab="log library size", ylab="fraction of zeroes", main= "Usoskin et al. 2015", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib*pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]+coef(m2)[5]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
b <- coef(m2)[2]+coef(m2)[6]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)



######### OLD PLOT, still have to adapt code.
## P(zero) ~ library size
png("~/Dropbox/phdKoen/singleCell/figures/supplementary/usoskin.png", width=10,height=5, units="in", res=330)
par(mfrow=c(1,3))
expit <- function(x) 1/(1+exp(-x))
## not including batch effect
par(mar=c(5,4.5,4,1))
logLib=log(colSums(usoskin))
pZero=colMeans(usoskin==0)
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m <- glm(pZero~logLib,family="binomial")
a <- coef(m)[1]
b <- coef(m)[2]
lines(x=sort(logLib),y=expit(a+b*sort(logLib)),lwd=2,col="steelblue")
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,pch=1,cex=1.5)
mtext("a" ,at=4, font=2, cex=4/3)

## including main batch effect
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib+pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,lty=1,cex=1.5)
mtext("b" ,at=4, font=2, cex=4/3)

load("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR/weightsConvergedZingeREdgeRUsoskin.rda")
par(mar=c(5,4.5,4,1))
hist(w[usoskin==0],xlab="Posterior probability",main="",cex.lab=2,cex.axis=2, ylim=c(0,3.5e6))
#load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/weightsUsoskinNoBatchLibSizeDispFast300Iter.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR_noBatch/weightsDESeq2ZeroUsoskin45.RData")
hist(w[usoskin==0],add=TRUE,col=rgb(0.1,0.8,0.1,.2))
legend("topleft",c("picking session + eff. lib. size","eff. lib. size"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n",cex=1.25)
mtext("c" ,at=-.25, font=2, cex=4/3)
dev.off()


####################################################
################# ANALYSES #########################
####################################################
### ZINB-WaVE-edgeR
library(zinbwave) ; library(edgeR)
computeZinbwaveWeights <- function(zinb, counts){
  mu <- getMu(zinb)
  pi <- getPi(zinb)
  theta <- getTheta(zinb)
  theta_mat <- matrix(rep(theta, each = ncol(counts)), ncol = nrow(counts))
  nb_part <- dnbinom(t(counts), size = theta_mat, mu = mu)
  zinb_part <- pi * ( t(counts) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  t(zinbwg)
}

library(doParallel)
registerDoParallel(2)
p=BiocParallel::DoparParam()
design=model.matrix(~cellTypeAll + pickingSession)
zinb_c <- zinbFit(usoskin, X = design, K=2, commondispersion = TRUE, BPPARAM=p)
#save(zinb_c, file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbWaveFitCommonUsoskin.rda")
w = computeZinbwaveWeights(zinb_c, usoskin)
dZinb = DGEList(usoskin)
dZinb = edgeR::calcNormFactors(dZinb)
dZinb$weights = w
dZinb = estimateDisp(dZinb, design)
fitZinb = glmFit(dZinb, design)
L <- matrix(0,nrow=ncol(fitZinb$coefficients),ncol=11)
rownames(L) <- colnames(fitZinb$coefficients)
colnames(L) <- c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH")
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
lrtListZinb=list()
for(i in 1:ncol(L)) lrtListZinb[[i]] <- glmWeightedF(fitZinb,contrast=L[,i],test="F", ZI=TRUE)
padjListZinb=lapply(lrtListZinb, function(x) p.adjust(x$table$PValue,"BH"))
deGenesZinb=unlist(lapply(padjListZinb,function(x) sum(x<.05,na.rm=TRUE)))






### edgeR analysis
design=model.matrix(~cellTypeAll+pickingSession)
dNoWeights=DGEList(usoskin)
dNoWeights=edgeR::calcNormFactors(dNoWeights)
dNoWeights=estimateDisp(dNoWeights,design)
#save(dNoWeights,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dUsoskinNoWeights.rda")
fitNoWeights=glmFit(dNoWeights,design)
L <- matrix(0,nrow=ncol(fitNoWeights$coefficients),ncol=11)
rownames(L) <- colnames(fitNoWeights$coefficients)
colnames(L) <- c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH")
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
lrtListNoWeights=list()
for(i in 1:ncol(L)) lrtListNoWeights[[i]] <- glmWeightedF(fitNoWeights,contrast=L[,i],test="F", ZI=FALSE)
padjListNoWeights=lapply(lrtListNoWeights, function(x) p.adjust(x$table$PValue,"BH"))
deGenesNoWeights=unlist(lapply(padjListNoWeights,function(x) sum(x<.05)))


#### how many genes from the edgeR analysis are also found in the ZINB-WaVE_edgeR analysis and vice versa?
deGenesZinbAll = lapply(lrtListZinb, function(x) which(p.adjust(x$table$PValue,"fdr")<=0.05))
deGenesEdgerAll = lapply(lrtListNoWeights, function(x) which(p.adjust(x$table$PValue,"fdr")<=0.05))
genesInCommon = vector()
for(i in 1:length(deGenesZinbAll)) genesInCommon[i] = mean(deGenesEdgerAll[[i]] %in% deGenesZinbAll[[i]])
for(i in 1:length(deGenesZinbAll)) genesInCommon[i] = mean(deGenesZinbAll[[i]] %in% deGenesEdgerAll[[i]])


### DESeq2 analysis
library(DESeq2)
colData <- data.frame(cellType=cellTypeAll,pickingSession=pickingSession)
dse <- DESeqDataSetFromMatrix(countData = usoskin, colData = colData, design = ~ cellType+pickingSession)
dse <- DESeqDataSetFromMatrix(countData = usoskin, colData = colData, design = ~ cellType+pickingSession)
dse = estimateSizeFactors(dse)
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE)
#save(dse,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dseUsoskin.rda")
resultsNames(dse) #for building contrasts, see ?results
L=matrix(0,nrow=length(resultsNames(dse)),ncol=11)
rownames(L)=resultsNames(dse)
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
resList=list()
for(i in 1:ncol(L)) resList[[i]] = results(dse,contrast=L[,i])
#save(resList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListDESeq2Usoskin.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListDESeq2Usoskin.rda")
lapply(resList,function(x) sum(x$padj<=0.05, na.rm=TRUE))
cbind(deGenesBatch,deGenesNoWeights,deGenesDESeq2)

### DESeq2 poscounts analysis
library(DESeq2)
colData <- data.frame(cellType=cellTypeAll,pickingSession=pickingSession)
dse <- DESeqDataSetFromMatrix(countData = usoskin, colData = colData, design = ~ cellType+pickingSession)
dse = estimateSizeFactors(dse, type="poscounts")
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE)
#save(dse,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dsePoscountsUsoskin.rda")
L=matrix(0,nrow=length(resultsNames(dse)),ncol=11)
rownames(L)=resultsNames(dse)
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
resList=list()
for(i in 1:ncol(L)) resList[[i]] = results(dse,contrast=L[,i])
#save(resList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListPoscountsUsoskin.rda")


### MAST analysis on cpm
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskinCpm.RData")
library(MAST)
tpm=exprs(esetCpm)
sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(cellType=cellTypeAll, pickingSession=pickingSession))
ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
CD <- cData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
cData(sca) <- CD
## differential expression
fit <- zlm.SingleCellAssay(~cngeneson+cellType+pickingSession,sca=sca,method="bayesglm",ebayes=TRUE)
#how many genes have all-zero counts in at least one condition?
mean(apply(fit@coefC[,c(1,3:12)],1,function(row) any(is.na(row))))
L <- matrix(0,nrow=ncol(fit@coefC),ncol=11)
rownames(L) <- colnames(fit@coefC)
colnames(L) <- c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH")
L[3:12,1] <- -1/10 #NF1 vs. others
L[3:12,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[3:12,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[3:12,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[3:12,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[3:12,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[3:12,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[3:12,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[3:12,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[3:12,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[3:12,11] <- c(rep(-1/10,9),1) #TH vs. others
lrList=list()
for(i in 1:ncol(L)) lrList[[i]] = lrTest(fit,hypothesis=matrix(L[,i],ncol=1))
#save(lrList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/lrListMASTUsoskin.rda")
#hurdle model results
padjListHurdle <- lapply(lrList,function(x) p.adjust(x[,'hurdle','Pr(>Chisq)'],"BH"))
unlist(lapply(padjListHurdle,function(x) sum(x<=0.05))) #nr DE genes at 5%
# count component results
padjListCont <- lapply(lrList,function(x) p.adjust(x[,'cont','Pr(>Chisq)'],"BH"))
unlist(lapply(padjListCont,function(x) sum(x<=0.05))) #nr DE genes at 5%





###### JUNK
### slow EM-algorithm
pdf("~/Dropbox/PhD/Research/zeroInflation/singleCell/
   ")
weightsUsoskinBatch = zeroWeightsLibSizeDispFast(counts=usoskin, designZI=model.matrix(~effLogLib+pickingSession), design=model.matrix(~cellTypeAll+pickingSession), maxit=300, plotW=TRUE)
dev.off()



### next 5x10 iterations
w=weightsUsoskinBatch10Iterations
#first I did from 1 to 5 (so 1 is iteration 10-20, 2 is 21-30,...)
#then from 6 to 10
#then from 11 to 15
for(j in 11:15){
    pdf(paste0("~/Dropbox/phdKoen/singleCell/usoskinweightsBatchInMainAndZeroModelAllCellTypes",j,".pdf"))
counts=DGEList(usoskin)
designZI=model.matrix(~effLogLib+pickingSession)
design=model.matrix(~cellTypeAll+pickingSession)
maxit=10
plot=TRUE
plotW=TRUE
    if(plot | plotW) par(mfrow=c(1,plot+plotW))
    zeroId <- counts$counts==0
    #w <- weightsUsoskinBatch10Iterations
    wFinal <- matrix(NA,nrow=nrow(counts),ncol=ncol(counts))
    active <- rowSums(counts$counts==0)>0 #work with genes with at least 1 zero
    #wFinal[!active,]=w[!active,]
    llOld <- matrix(-1e4,nrow=nrow(counts),ncol=ncol(counts))
    likCOld <- matrix(0,nrow=nrow(counts),ncol=ncol(counts))
    maxit=10
    for(i in 1:maxit){
        zeroId <- counts$counts==0
	counts$weights <- w

	### M-step counts
	counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
	counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
	if(plot) plotBCV(counts)
	fit <- glmFit(counts, design)
	if(i>1) likCOld <- likC[!converged,]
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)

	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response")

	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	wOld <- w
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)
	rownames(w) <- rownames(wOld)
	if(plotW) hist(w[zeroId])

	## expected complete data log-likelihood
	if(i>1) llOld <- ll[!converged,]
	ll <- w*log(pi0HatMat) + (1-w)*log(1-pi0HatMat) + (1-w)*log(likC)

	converged=FALSE
	if(any(converged)){
	    wFinal[as.numeric(rownames(w)[converged]),] = w[converged,]
	    w <- matrix(w[!converged,],ncol=ncol(counts), dimnames=list(c(rownames(w)[!converged]),NULL))
	    counts <- counts[!converged,]
	}
	#cat(paste0("mean diff in L: ",round(mean(rowSums(exp(ll)-exp(llOld))),2),". active features: ",nrow(counts),"\n"))
	cat(paste0("iteration ",i))
	#metagenomeSeq for example seems to report mean(ll) instead of difference.
	if(all(converged)) break
	if(i==maxit){
	    wFinal[apply(wFinal,1,function(row) any(is.na(row))),] = w
	    break
	}
    }
    dev.off()
    save(w,file=paste0("/Users/koenvandenberge/PhD_Data/singleCell/usoskinCase/weightsUsoskinBatch",j,".rda"))
}



### ignoring the batch effect
pdf("~/Dropbox/phdKoen/singleCell/usoskinweights100IterationsBatchInMainNotInZeroModelAllCellTypes.pdf")
weightsUsoskinNoBatch50Iterations = zeroWeightsLibSizeFast(counts=usoskin, design=model.matrix(~cellTypeAll+pickingSession), maxit=100, plot=TRUE, plotW=TRUE)
dev.off()
