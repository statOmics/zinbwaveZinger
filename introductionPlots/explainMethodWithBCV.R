source("../zingeRsimulationFunctions/simulationHelpFunctions_v7_diffInZero.R")
source("../method/glmLRTOld.R")
library(zinbwave) ; library(BiocParallel) ; library(doParallel)
## set up
NCORES <- 2
registerDoParallel(NCORES)
register(DoparParam())
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
###################################
######## RNA-seq prep #############
###################################
#### no zero inflation simulation
load("../datasets/bottomly_eset.RData")
bottomly=exprs(bottomly.eset)
bottomly=bottomly[!rowSums(bottomly)==0,]
nSamp <- 10
nTags <- 20e3
set.seed(2)
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(8e6,10e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = bottomly, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = FALSE, lib.size=libSize, drop.low.lambda=TRUE)


data=dataNoZI$counts
d=DGEList(data)
d=calcNormFactors(d)
design=model.matrix(~grp)
d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
dOriginal=d
dOriginalTrend=estimateGLMTrendedDisp(d,design)
plotBCV(d, ylim=c(0,3))

# add zeros
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
set.seed(3)
samp=sample(1:length(zeroId),floor(length(zeroId)*propZeroes))
zeroId[samp]=0
zeroId[dataNoZI$counts==0]=1 #if it already was a zero it is not zero-inflated.
samp=samp[!samp%in%which(dataNoZI$counts==0)] #same
dataZeroes$counts = dataZeroes$counts*zeroId
data=dataZeroes$counts
d=DGEList(data)
d=calcNormFactors(d)
design=model.matrix(~grp)
d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
dZeroes=d
dZeroesTrend=estimateGLMTrendedDisp(dZeroes,design)
plotBCV(d, ylim=c(0,3))

dataZeroEset = SummarizedExperiment(dataZeroes$counts, colData=data.frame(grp=grp))
zinbModelRNAseq = zinbFit(dataZeroEset, X="~grp")
zeroWeightsSim = computeZinbwaveWeights(zinbModelRNAseq, dataZeroes$counts)
hist(zeroWeightsSim[samp],xlab="post. prob. on count component")

## downweighting
d=DGEList(data)
d=calcNormFactors(d)
d$weights=zeroWeightsSim
d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
dWeighted=d
dWeightedTrend=estimateGLMTrendedDisp(dWeighted,design)
plotBCV(d, ylim=c(0,3))

## ROC
pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
falses=which(dataNoZI$counts==0)
tpr=fpr=vector(length=length(pvalSeq))
for(i in 1:length(pvalSeq)){
   excessID <- which(zeroWeightsSim<=pvalSeq[i])
   tpr[i] <- mean(samp%in%excessID)
   fpr[i] <- mean(falses%in%excessID)
}
plot(x=fpr,y=tpr,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")


###################################
######## scRNA-seq prep ###########
###################################
library(scales)
library(GEOquery)
data = read.delim("../datasets/expressionTabAdapted_kvdb.txt")
load("../datasets/seriesMatrix.rda")
countData = data[,8:ncol(data)]
rownames(countData)=data[,1]
countData=countData[,1:92]
countData = countData[-c(1:8),] #rm spike-ins
well = factor(substr(colnames(countData),1,1))
fibroID <- grep(x=pData(seriesMatrix[[1]])$title,pattern="fibroblast")
stemCellID <- grep(x=pData(seriesMatrix[[1]])$title,pattern="stem cell")
colnames(countData)[fibroID] = paste("fibro",1:length(fibroID),sep="_")
colnames(countData)[stemCellID] = paste("stemCell",1:length(stemCellID),sep="_")
cellType=vector(length=ncol(countData))
cellType[fibroID] = "fibro"
cellType[stemCellID] <- "stemCell"
islam = as.matrix(countData)
islam = islam[!rowSums(islam>0)<5,]


### simulate data
paramsIslam = getDatasetMoMPositive(counts = islam)
nSamples=80
grp=as.factor(rep(0:1, each = nSamples/2)) #two-group comparison
nTags=10000 #nr of features
set.seed(11)
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #10% DE
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #fold changes
libSizes=sample(colSums(islam),nSamples,replace=TRUE) #library sizes
simDataIslam <- NBsimSingleCell(foldDiff = fcSim, ind = DEind,
                                dataset = islam, nTags = nTags,
                                group = grp,
                                verbose = TRUE, params = paramsIslam,
                                lib.size = libSizes, normalizeLambda=TRUE)

### make BCV from simulated data
design=model.matrix(~grp)
dIslam=DGEList(simDataIslam$counts)
dIslam=calcNormFactors(dIslam)
dIslam=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dIslam,design, interval=c(0,10)),design,prior.df=0)
dIslamTrend=estimateDisp(dIslam,design,prior.df=0)

## get ZINB-WaVE weights
islamEset = core <- SummarizedExperiment(simDataIslam$counts, colData = data.frame(group = grp))
zinbModel = zinbFit(islamEset, X="~group", epsilon=1e12)
zeroWeights = computeZinbwaveWeights(zinbModel, simDataIslam$counts)
## calculate dispersion
dW=DGEList(simDataIslam$counts)
dW=calcNormFactors(dW)
dW$weights=zeroWeights
dW=estimateDisp(dW,design,prior.df=0)

## ROC simulated data zero excess
pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
falses=which(simDataIslam$dropout==1)
samp=which(simDataIslam$dropout==0)
tprSc=fprSc=vector(length=length(pvalSeq))
for(i in 1:length(pvalSeq)){
   excessID <- which(zeroWeights<=pvalSeq[i])
   tprSc[i] <- mean(samp%in%excessID)
   fprSc[i] <- mean(falses%in%excessID)
}
plot(x=fprSc,y=tprSc,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")



## BCV after weighting





## BCV for islam
# design=model.matrix(~cellType)
# dIslam=DGEList(islam)
# dIslam=calcNormFactors(dIslam)
# dIslam=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dIslam,design, interval=c(0,10)),design,prior.df=0)
# dIslamTrend=estimateDisp(dIslam,design,prior.df=0)

## get ZINB-WaVE weights
##estimate ZINB model
islamEset = core <- SummarizedExperiment(islam, colData = data.frame(cellType = cellType))
zinbModel = zinbFit(islamEset, X="~cellType")
zeroWeights = computeZinbwaveWeights(zinbModel, islam)
## calculate dispersion
dW=DGEList(islam)
dW=calcNormFactors(dW)
dW$weights=zeroWeights
dW=estimateDisp(dW,design,prior.df=0)

## ROC
pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
falses=which(dataNoZI$counts==0)
tpr=fpr=vector(length=length(pvalSeq))
for(i in 1:length(pvalSeq)){
   excessID <- which(zeroWeightsSim<=pvalSeq[i])
   tpr[i] <- mean(samp%in%excessID)
   fpr[i] <- mean(falses%in%excessID)
}
plot(x=fpr,y=tpr,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")


#######################################
### one composite plot, version 2 ################
#######################################
wMean=sapply(1:nrow(zeroWeightsSim), function(i){
	  mean(zeroWeightsSim[i,dataZeroes$counts[i,]==0])
})
wMean[is.na(wMean)]=1
library(Hmisc)
cuts=cut2(wMean,cuts=c(0,0.25,0.75,1))


png("~/Dropbox/phdKoen/singleCell/zinbwaveZingerGithub/zinbwaveZinger/introductionPlots/introBCV_v3.png", width=10,height=5, units="in", res=330)
#pdf("~/Dropbox/phdKoen/singleCell/zinbwaveZingerGithub/zinbwaveZinger/introductionPlots/introBCV_v2.pdf", width=10,height=5)
#dev.new(width=10,height=5)
##### plot for paper: RNA-seq
library(scales)
par(mar=c(4.1,4.25,3,1),bty="l", mfrow=c(2,3), cex.lab=1.5, cex.axis=1.5)

cols=c("red","gold","black")
plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average log CPM", ylab="BCV",col=alpha(cols[as.numeric(cuts)],1/3), type="n")
#sapply(c(4,1,2,3),function(i){
sapply(c(3,1,2),function(i){
	points(x=dZeroes$AveLogCPM[cuts==levels(cuts)[i]],y=sqrt(dZeroes$tagwise.dispersion)[cuts==levels(cuts)[i]],pch=16,cex=.2, col=alpha(cols[i]))
})
o<- order(dZeroes$AveLogCPM)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
#legend("topleft",legend=c("[0,0.25)", "[0.25,0.5)", "[0.5,0.75)", "[0.75,1]" ),col=cols[1:5],lty=1, bty="n",lwd=2, cex=.8,inset=c(0,-0.045))
legend("topleft",legend=c("[0,0.25)", "[0.25,0.75)", "[0.75,1]" ),col=cols,lty=1, bty="n",lwd=2, cex=.9,inset=c(0,-0.045))
mtext("a" ,at=-5, font=2, cex=4/3)

## ROC
plot(x=c(0,fpr),y=c(0,tpr),type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")
legend("bottomright", "Bottomly ZINB-WaVE weights", lty=1, col="steelblue", bty="n")
mtext("b" ,at=-0.22, font=2, cex=4/3)


plot(x=dWeighted$AveLogCPM,y=sqrt(dWeighted$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dWeighted$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
lines(dWeighted$AveLogCPM[o], sqrt(dWeightedTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, NB simul.","NB model, ZINB simul.","ZINB model, ZINB simul."),lty=1,lwd=2,col=c("red","blue","steelblue1"),bty="n")
mtext("c" ,at=-5, font=2, cex=4/3)


### scRNA-seq
plot(dIslam$AveLogCPM,sqrt(dIslam$tagwise.dispersion),pch=16,cex=.2,col=alpha("black",1/3),xlim=c(1.8,12), xlab="Average log CPM", ylab="BCV", ylim=c(0,15))
o <- order(dIslam$AveLogCPM)
lines(dIslam$AveLogCPM[o], sqrt(dIslamTrend$trended.dispersion)[o], col = "red",lwd = 2)
mtext("d" ,at=-0.5, font=2, cex=4/3)

plot(x=fprSc,y=tprSc,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")
legend("bottomright", "Islam ZINB-WaVE weights", lty=1, col="steelblue", bty="n")
mtext("e" ,at=-0.22, font=2, cex=4/3)


plot(dW$AveLogCPM,sqrt(dW$tagwise.dispersion),pch=16,cex=.2,xlab="Average log CPM", ylab="BCV", xlim=c(1.8,12), ylim=c(0,15),col=alpha("black",1/3))
o=order(dIslamTrend$AveLogCPM)
lines(dIslamTrend$AveLogCPM[o], sqrt(dIslamTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dW$AveLogCPM[o], sqrt(dW$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, scRNA-seq","ZINB model, scRNA-seq"),lty=1,col=c("red","steelblue1"),lwd=2,bty="n")
mtext("f" ,at=-0.8, font=2, cex=4/3)

dev.off()
