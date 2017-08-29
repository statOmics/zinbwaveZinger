source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v7_diffInZero.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")
library(DESeq2) ; library(edgeR) ; library(limma) ; library(scde) ; library(MAST) ; library(mgcv) ; library(MultiAssayExperiment) ; library(SummarizedExperiment) ; library(scales) ; library(iCOBRA)

##################################################################
########################### ISLAM ################################
##################################################################
library(GEOquery)
data = read.delim("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/expressionTabAdapted_kvdb.txt")
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/seriesMatrix.rda")
countData = data[,8:ncol(data)]
rownames(countData)=data[,1]
#(pData(seriesMatrix)$title)[93:96] #negative controls added by the authors: can be discarded.
countData=countData[,1:92]
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

## weight distribution on real data
design=model.matrix(~cellType)
weightsIslamIter50=zeroWeightsLibSize(counts=islam, niter=50, design=design, plotW=FALSE)
weightsIslamFast=zeroWeightsLibSizeDispFast(counts=islam, maxit=200, design=design)
par(mfrow=c(1,2))
hist(weightsIslamIter50[islam==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))
hist(weightsIslamFast[islam==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))

#axis(2,at=c(0,2e4,4e4,6e4),labels=c("0","2e4","4e4","6e4"))

## get gene-wise parameters acc to zero-truncated NB distribution
#paramsIslamAllDesignAveLogCPM=getDatasetZTNB(counts=islam, design=model.matrix(~cellType))
#save(paramsIslamAllDesignAveLogCPM,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsIslamAllDesignAveLogCPM.rda")
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsIslamAllDesignAveLogCPM.rda")

####### 40 samples / condition
set.seed(12)
nSamp <- 80
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 15e3
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(11)
libSizes=sample(colSums(islam),nSamp,replace=TRUE)
dataIslamAllAveLogCPM <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = islam, nTags = nTags, group = grp, verbose = TRUE, params=paramsIslamAllDesignAveLogCPM, lib.size=libSizes, cpm="AveLogCPM")

### compare simulated and empirical data
#empirical BCV
condition = cellType
design = model.matrix(~condition)
dIslam=DGEList(islam)
dIslam=calcNormFactors(dIslam)
dIslam=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dIslam,design, interval=c(0,10)),design,prior.df=0)
#simulated BCV
dataIslam=dataIslamAllAveLogCPM
design=model.matrix(~grp)
dSim=DGEList(dataIslam$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)

png("~/Dropbox/phdKoen/singleCell/figures/supplementary/simulationIslam.png", width=7,height=9, units="in", res=380)
par(mfrow=c(2,2), mar=c(5,4,1,1))
plotBCV(dIslam,main="Empirical", ylim=c(0,12), ylab="BCV")
plotBCV(dSim, main="Simulated", ylim=c(0,12), ylab="BCV",xlim=c(1,15))
## zero ~ libSize
plot(x=log(colSums(islam)),y=colMeans(islam==0), xlab="Log library size", ylab="Fraction of zeros", xlim=c(9.5,15.5), ylim=c(0.25,0.95))
points(x=log(colSums(dataIslamAllAveLogCPM$counts)),y=colMeans(dataIslamAllAveLogCPM$counts==0), xlab="Log library size", ylab="Fraction of zeros",col=2)
legend("topright",c("real","simulated"),pch=1,col=1:2, bty="n", cex=.6)
## zero ~ cpm
plot(x=aveLogCPM(islam),y=rowMeans(islam==0), xlab="Average log CPM", ylab="Fraction of zeros", xlim=c(0.5,16),col=alpha(1,1/2),pch=19,cex=.2)
points(x=aveLogCPM(dataIslamAllAveLogCPM$counts),y=rowMeans(dataIslamAllAveLogCPM$counts==0), xlab="aCPM", ylab="Fraction of zeros",col=alpha(2,1/2),pch=19,cex=.2)
legend("topright",c("real","simulated"),pch=19,col=1:2, bty="n", cex=.6)
dev.off()

### performance characteristics
selectedMethods <- c("edgeREMLibSizeDispFastOldFFilteredEdgeR", "MAST", "limma_voomZeroFiltered", "DESeq2Zero_adjustedDf_posCountsNormZeroWeights", "DESeq2", "DESeq2_poscounts", "edgeR", "edgeROldF", "limma_voom", "limma_voomFiltered", "metagenomeSeq","edgeRFiltered", "NODES")
group=grp
#pvalsIslam <- pval(dataIslamAllAveLogCPM, method=selectedMethods, mc.cores=2, niter=200)
#save(pvalsIslam,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamPaper_v7_diffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamPaper_v7_diffZero.rda")
for(k in 1:length(pvalsIslam$padj)) pvalsIslam$padj[[k]][is.na(pvalsIslam$padj[[k]])]=1
for(k in 1:length(pvalsIslam$pval)) pvalsIslam$pval[[k]][is.na(pvalsIslam$pval[[k]])]=1

## scde
#scdePIslam=scde.pfun(dataIslamAllAveLogCPM$counts,group=grp,mc.cores=1) #gives trouble in parallellization
#save(scdePIslam,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslam.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslam.rda")

## ZINB-WaVE: common dispersion
library(zinbwave)
counts=dataIslamAllAveLogCPM$counts
design = model.matrix(~ grp)
library(doParallel)
registerDoParallel(2)
p=BiocParallel::DoparParam()
#zinb = zinbFit(counts, X=design, K=2, BPPARAM=p, commondispersion=TRUE)
#save(zinb,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbIslamCommonDiffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbIslamCommonDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_edgeR_common= zingeREdgeROwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_edgeR_common[is.na(zinbwave_edgeR_common[,"pval"]),"pval"]=1
zinbwave_edgeR_common[is.na(zinbwave_edgeR_common[,"padj"]),"padj"]=1


#zinb = zinbFit(counts, X=design, K=2, BPPARAM=p, commondispersion=FALSE)
#save(zinb,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbIslamGeneDiffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbIslamGeneDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_edgeR_gene= zingeREdgeROwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_edgeR_gene[is.na(zinbwave_edgeR_gene[,"pval"]),"pval"]=1
zinbwave_edgeR_gene[is.na(zinbwave_edgeR_gene[,"padj"]),"padj"]=1

## ZINB-WaVE DESeq2 common dispersion
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbIslamCommonDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_DESeq2_common= zingeRDESeq2OwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_DESeq2_common[is.na(zinbwave_DESeq2_common[,"pval"]),"pval"]=1
zinbwave_DESeq2_common[is.na(zinbwave_DESeq2_common[,"padj"]),"padj"]=1

## ZINB-WaVE DESeq2 gene dispersion
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbIslamGeneDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
notWeightsOkId=c(5351,8967) #because not a single pos count in a group for these genes.
w[notWeightsOkId[1],1]=1 #ad hoc for plot
w[notWeightsOkId[2],1]=1 #ad hoc for plot
zinbwave_DESeq2_gene= zingeRDESeq2OwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_DESeq2_gene[is.na(zinbwave_DESeq2_gene[,"pval"]),"pval"]=1
zinbwave_DESeq2_gene[is.na(zinbwave_DESeq2_gene[,"padj"]),"padj"]=1


# edgeR hurdle
hurdleWeights = 1-(dataIslamAllAveLogCPM$counts==0+0)
hurdleWeights[hurdleWeights==0]=1e-15
edgeRHurdleRes = zingeREdgeROwnWeights.pfun(counts=counts, group=grp, w=hurdleWeights)




#FDR-TPR
truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPM))
truthIslam[dataIslamAllAveLogCPM$indDE,"status"]=1
cobraIslam <- COBRAData(pval =data.frame(
					#zingeR_edgeR=pvalsIslam$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					 #DESeq2=pvalsIslam$pval$DESeq2,
					 DESeq2_poscounts=pvalsIslam$pval$DESeq2_poscounts,
					 zingeR_DESeq2=pvalsIslam$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					 #limma_voom=pvalsIslam$pval$limma_voom,
					 #zingeR_limma_voom=pvalsIslam$pval$limma_voomZero,
					 #edgeR=pvalsIslam$pval$edgeR,
					 #MAST=pvalsIslam$pval$MAST,
					 #metagenomeSeq=pvalsIslam$pval$metagenomeSeq,
					 #scde=scdePIslam[,"pval"],
					 #NODES=pvalsIslam$pval$NODES,
					 #zinbwave_edgeR_common=zinbwave_edgeR_common[,"pval"],
					 #zinbwave_edgeR_gene=zinbwave_edgeR_gene[,"pval"],
					 zinbwave_DESeq2_common=zinbwave_DESeq2_common[,"pval"],
					 zinbwave_DESeq2_gene=zinbwave_DESeq2_gene[,"pval"],
					 #edgeRHurdle=edgeRHurdleRes[,"pval"],
					 row.names = rownames(dataIslamAllAveLogCPM)),
		   padj = data.frame(
		   			 #zingeR_edgeR=pvalsIslam$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
				     #DESeq2=pvalsIslam$padj$DESeq2,
				     DESeq2_poscounts=pvalsIslam$padj$DESeq2_poscounts,
					   zingeR_DESeq2=pvalsIslam$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     #limma_voom=pvalsIslam$padj$limma_voom,
					   #zingeR_limma_voom=pvalsIslam$padj$limma_voomZero,
					   #edgeR=pvalsIslam$padj$edgeR,
					   #MAST=pvalsIslam$padj$MAST,
					   #metagenomeSeq=pvalsIslam$padj$metagenomeSeq,
					   #scde=scdePIslam[,"padj"],
					   #NODES=pvalsIslam$padj$NODES,
					   #zinbwave_edgeR_common=zinbwave_edgeR_common[,"padj"],
					   #zinbwave_edgeR_gene=zinbwave_edgeR_gene[,"padj"],
						 zinbwave_DESeq2_common=zinbwave_DESeq2_common[,"padj"],
						 zinbwave_DESeq2_gene=zinbwave_DESeq2_gene[,"padj"],
						 #edgeRHurdle=edgeRHurdleRes[,"padj"],
				    row.names = rownames(dataIslamAllAveLogCPM)),
                   truth = truthIslam)
cobraperf <- calculate_performance(cobraIslam, binary_truth = "status")
colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black", zingeR_DESeq2_wald="steelblue", MAST2="salmon", edgeRHurdle="dodgerblue", DESeq2_noCook="gold", zinbwave_edgeR_gene="aquamarine3", zinbwave_edgeR_common="gold", zinbwave_DESeq2_gene="aquamarine3", zinbwave_DESeq2_common="gold")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslam.rda")
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamLimma.rda")
plot_fdrtprcurve(cobraplot, pointsize=2)

### check correct identification of dropout and NB zero
#zingeREdgeRWeights=zingeREdgeRWeights.pfun(counts=dataIslamAllAveLogCPM$counts, group=grp)
#save(zingeREdgeRWeights,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zingeREdgeRWeightsIslamDiffZero.rda")
library(rafalib) ; mypar()
par(mfrow=c(1,2))
hist(w[dataIslamAllAveLogCPM$dropout==0], breaks=seq(0,1,0.05), main="ZINB-WaVE weight distribution for dropout on Islam simulation", cex.main=2/3)
hist(zingeREdgeRWeights[dataIslamAllAveLogCPM$dropout==0], breaks=seq(0,1,0.05), main="zingeR weight distribution for dropout on Islam simulation", cex.main=2/3)

par(mfrow=c(1,2))
hist(w[dataIslamAllAveLogCPM$dropout==1], breaks=seq(0,1,0.05), main="ZINB-WaVE weight distribution for NB zero on Islam simulation", cex.main=2/3)
hist(zingeREdgeRWeights[dataIslamAllAveLogCPM$dropout==1], breaks=seq(0,1,0.05), main="zingeR weight distribution for NB zero on Islam simulation", cex.main=2/3)

# #high FC FDR-TPR
# fcSim=(3 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
# set.seed(11)
# libSizes=sample(colSums(islam),nSamp,replace=TRUE)
# dataIslamAllAveLogCPMHighFC <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = islam, nTags = nTags, group = grp, verbose = TRUE, params=paramsIslamAllDesignAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM", min.dispersion=0.001)
# pvalsIslamHighFC <- pval(dataIslamAllAveLogCPMHighFC, method=selectedMethods, mc.cores=2, niter=200)
# #save(pvalsIslamHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamHighFCPaper.rda")
# load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamHighFCPaper.rda")
# for(k in 1:length(pvalsIslamHighFC$padj)) pvalsIslamHighFC$padj[[k]][is.na(pvalsIslamHighFC$padj[[k]])]=1
# #scdePHighFC=scde.pfun(dataIslamAllAveLogCPMHighFC$counts,group=grp)
# #save(scdePHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslamHighFC.rda")
# load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslamHighFC.rda")
# truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPMHighFC))
# truthIslam[dataIslamAllAveLogCPMHighFC$indDE,"status"]=1
# cobraIslam <- COBRAData(pval =data.frame(
# 					zingeR_edgeR=pvalsIslamHighFC$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
# 					 DESeq2=pvalsIslamHighFC$pval$DESeq2,
# 					 DESeq2_poscounts=pvalsIslamHighFC$pval$DESeq2_poscounts,
# 					 zingeR_DESeq2=pvalsIslamHighFC$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
# 					 limma_voom=pvalsIslamHighFC$pval$limma_voom,
# 					 zingeR_limma_voom=pvalsIslamHighFC$pval$limma_voomZero,
# 					 edgeR=pvalsIslamHighFC$pval$edgeR,
# 					 MAST=pvalsIslamHighFC$pval$MAST,
# 					 metagenomeSeq=pvalsIslamHighFC$pval$metagenomeSeq,
# 					 scde=scdePHighFC[,"pval"],
# 					 NODES=pvalsIslamHighFC$pval$NODES,
# 					 row.names = rownames(dataIslamAllAveLogCPMHighFC)),
# 		   padj = data.frame(
# 		   			zingeR_edgeR=pvalsIslamHighFC$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
# 				     DESeq2=pvalsIslamHighFC$padj$DESeq2,
# 				     DESeq2_poscounts=pvalsIslamHighFC$padj$DESeq2_poscounts,
# 					 zingeR_DESeq2=pvalsIslamHighFC$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
# 				     limma_voom=pvalsIslamHighFC$padj$limma_voom,
# 					 zingeR_limma_voom=pvalsIslamHighFC$padj$limma_voomZero,
# 				     edgeR=pvalsIslamHighFC$padj$edgeR,
# 					MAST=pvalsIslamHighFC$padj$MAST,
# 					 metagenomeSeq=pvalsIslamHighFC$padj$metagenomeSeq,
# 					scde=scdePHighFC[,"padj"],
# 					NODES=pvalsIslamHighFC$padj$NODES,
# 				     	row.names = rownames(dataIslamAllAveLogCPMHighFC)),
#                    truth = truthIslam)
# cobraperf <- calculate_performance(cobraIslam, binary_truth = "status")
# colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black")
# colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
# cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
# #save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamHighFC.rda")
# plot_fdrtprcurve(cobraplot, pointsize=2)
#
# ###### DESeq2 comparisons
# selectedMethodsDESeq2 <- c("DESeq2",
# 		     "DESeq2_noFiltering", #no independent filtering step
# 		     "DESeq2_noShrink", #no shrinkage to mean zero prior of coefs
# 				 "DESeq2_noCook", #no Cook's distance cooksCutoff
# 				 "DESeq2_withImputation",
# 				 "DESeq2_poscounts",
# 				 "DESeq2_poscounts_noShrink",
# 				 "DESeq2_poscounts_noFiltering",
# 				 "DESeq2_poscounts_noCook",
# 				 "DESeq2_poscounts_withImputation"
# 				 )
#
# pvalsIslamDESeq2 <- pval(dataIslamAllAveLogCPM, method=selectedMethodsDESeq2, mc.cores=2, niter=200)
# for(k in 1:length(pvalsIslamDESeq2$padj)) pvalsIslamDESeq2$padj[[k]][is.na(pvalsIslamDESeq2$padj[[k]])]=1
# #save(pvalsIslamDESeq2,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamDESeq2Paper.rda")
# load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamDESeq2Paper.rda")
# truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPM))
# truthIslam[dataIslamAllAveLogCPM$indDE,"status"]=1
# cobraIslamDESeq2 <- COBRAData(pval =data.frame(
# 					 DESeq2=pvalsIslamDESeq2$pval$DESeq2,
# 					 #DESeq2_poscounts=pvalsIslamDESeq2$pval$DESeq2_poscounts,
# 					 DESeq2_noFiltering=pvalsIslamDESeq2$pval$DESeq2_noFiltering,
# 					 DESeq2_noCook=pvalsIslamDESeq2$pval$DESeq2_noCook,
# 					 DESeq2_noShrink=pvalsIslamDESeq2$pval$DESeq2_noShrink,
# 					 DESeq2_withImputation=pvalsIslamDESeq2$pval$DESeq2_withImputation,
# 					 #DESeq2_poscounts_noShrink=pvalsIslamDESeq2$pval$DESeq2_poscounts_noShrink,
# 					 #DESeq2_poscounts_noFiltering=pvalsIslamDESeq2$pval$DESeq2_poscounts_noFiltering,
# 					 row.names = rownames(dataIslamAllAveLogCPM)),
# 		   padj = data.frame(
# 				     DESeq2=pvalsIslamDESeq2$padj$DESeq2,
# 				    # DESeq2_poscounts=pvalsIslamDESeq2$padj$DESeq2_poscounts,
# 						 DESeq2_noFiltering=pvalsIslamDESeq2$padj$DESeq2_noFiltering,
# 						 DESeq2_noCook=pvalsIslamDESeq2$padj$DESeq2_noCook,
# 						 DESeq2_noShrink=pvalsIslamDESeq2$padj$DESeq2_noShrink,
# 						 DESeq2_withImputation=pvalsIslamDESeq2$padj$DESeq2_withImputation,
# 						 #DESeq2_poscounts_noShrink=pvalsIslamDESeq2$padj$DESeq2_poscounts_noShrink,
# 						 #DESeq2_poscounts_noFiltering=pvalsIslamDESeq2$padj$DESeq2_poscounts_noFiltering,
# 				     	row.names = rownames(dataIslamAllAveLogCPM)),
# truth = truthIslam)
# cobraperf <- calculate_performance(cobraIslamDESeq2, binary_truth = "status")
# cobraplot <- prepare_data_for_plot(cobraperf)
# plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.3), yaxisrange=c(0,0.5))
#
# cobraIslamDESeq2poscounts <- COBRAData(pval =data.frame(
# 					 #DESeq2=pvalsIslamDESeq2$pval$DESeq2,
# 					 DESeq2_poscounts=pvalsIslamDESeq2$pval$DESeq2_poscounts,
# 					 #DESeq2_noFiltering=pvalsIslamDESeq2$pval$DESeq2_noFiltering,
# 					 #DESeq2_noCook=pvalsIslamDESeq2$pval$DESeq2_noCook,
# 					 #DESeq2_noShrink=pvalsIslamDESeq2$pval$DESeq2_noShrink,
# 					 DESeq2_poscounts_noShrink=pvalsIslamDESeq2$pval$DESeq2_poscounts_noShrink,
# 					 DESeq2_poscounts_noFiltering=pvalsIslamDESeq2$pval$DESeq2_poscounts_noFiltering,
# 					 DESeq2_poscounts_noCook=pvalsIslamDESeq2$pval$DESeq2_poscounts_noCook,
# 					 DESeq2_poscounts_withImputation=pvalsIslamDESeq2$pval$DESeq2_poscounts_withImputation,
# 					 row.names = rownames(dataIslamAllAveLogCPM)),
# 		   padj = data.frame(
# 				     #DESeq2=pvalsIslamDESeq2$padj$DESeq2,
# 				     DESeq2_poscounts=pvalsIslamDESeq2$padj$DESeq2_poscounts,
# 						 #DESeq2_noFiltering=pvalsIslamDESeq2$padj$DESeq2_noFiltering,
# 						 #DESeq2_noCook=pvalsIslamDESeq2$padj$DESeq2_noCook,
# 						 #DESeq2_noShrink=pvalsIslamDESeq2$padj$DESeq2_noShrink,
# 						 DESeq2_poscounts_noShrink=pvalsIslamDESeq2$padj$DESeq2_poscounts_noShrink,
# 						 DESeq2_poscounts_noFiltering=pvalsIslamDESeq2$padj$DESeq2_poscounts_noFiltering,
# 						 DESeq2_poscounts_noCook=pvalsIslamDESeq2$padj$DESeq2_poscounts_noCook,
# 						 DESeq2_poscounts_withImputation=pvalsIslamDESeq2$padj$DESeq2_poscounts_withImputation,
# 				     	row.names = rownames(dataIslamAllAveLogCPM)),
# truth = truthIslam)
# cobraperf <- calculate_performance(cobraIslamDESeq2poscounts, binary_truth = "status")
# cobraplotDESeq2poscounts <- prepare_data_for_plot(cobraperf)
# plot_fdrtprcurve(cobraplotDESeq2poscounts, pointsize=2, xaxisrange=c(0,0.3), yaxisrange=c(0,0.5))
#
# library(cowplot)
# deseq2plot =plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.3), yaxisrange=c(0,0.5), title="Median-of-ratios normalization")
# deseq2posplot =plot_fdrtprcurve(cobraplotDESeq2poscounts, pointsize=2, xaxisrange=c(0,0.3), yaxisrange=c(0,0.5), title="Positive counts normalization")
#
# pow = plot_grid(deseq2plot + theme(legend.position="none")   + xlab("FDP"), deseq2posplot + theme(legend.position="none") + xlab("FDP"), nrow=1, align='vh', hjust=-1, labels=c("a","b"))
# legend_a <- get_legend(deseq2plot + theme(legend.position="bottom"))
# legend_b <- get_legend(deseq2posplot + theme(legend.position="bottom"))
#
# png("~/Dropbox/phdKoen/singleCell/figures/supplementary/DESeq2Variants_islam.png", width=8,height=8, units="in", res=300)
# plot_grid(pow, legend_a, rel_heights=c(1,.2), ncol=1, nrow=2)
# dev.off()

##################################################################
########################### TRAPNELL #############################
##################################################################
trapnellAssay72 <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL11154.rds")
trapnellAssay72 = updateObject(trapnellAssay72)
trapnellAssay <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL16791.rds")
trapnellAssay = updateObject(trapnellAssay)
trapnellAssay48 <- trapnellAssay[,colData(trapnellAssay)[,"characteristics_ch1.1"] == "hour post serum-switch: 48"]
countsTrapnell72 <- round(assay(experiments(trapnellAssay72)$gene,"count"))
id48=colData(trapnellAssay)[,"characteristics_ch1.1"] == "hour post serum-switch: 48"
countsTrapnell48 <- round(assay(experiments(trapnellAssay)$gene[,id48],"count"))
#wells containing debris
debris72 = colData(trapnellAssay72)[,"characteristics_ch1.2"]=="debris: TRUE"
debris48 = colData(trapnellAssay48)[,"characteristics_ch1.2"]=="debris: TRUE"
#wells that did not contain one cell
one72 = colData(trapnellAssay72)[,"characteristics_ch1.4"]!="cells in well: 1"
one48 = colData(trapnellAssay48)[,"characteristics_ch1.4"]!="cells in well: 1"
# remove
countsTrapnell72 = countsTrapnell72[,(!debris72 & !one72)]
countsTrapnell48 = countsTrapnell48[,(!debris48 & !one48)]
countsTrapnell <- cbind(countsTrapnell48,countsTrapnell72)
countsTrapnell <- countsTrapnell[rowSums(countsTrapnell>0)>9,] #expression in at least 10 out of 149 samples. Remains 24,576 genes and 149 samples.
rm(trapnellAssay)
timePoint=factor(c(rep(48,85),rep(72,64)))

## weight distribution on real data
design=model.matrix(~timePoint)
weightsTrapnellIter50=zeroWeightsLibSize(counts=countsTrapnell, niter=50, design=design, plotW=FALSE)
weightsTrapnellFast=zeroWeightsLibSizeDispFast(counts=countsTrapnell, maxit=200, design=design)
par(mfrow=c(1,2))
hist(weightsTrapnellIter50[countsTrapnell==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))
hist(weightsTrapnellFast[countsTrapnell==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))

#paramsTrapnellAllDesignAveLogCPM=getDatasetZTNB(counts=countsTrapnell, design=model.matrix(~timePoint))
#save(paramsTrapnellAllDesignAveLogCPM,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsTrapnellAllDesignAveLogCPM.rda")
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsTrapnellAllDesignAveLogCPM.rda")


## simulate for 80 samples
set.seed(12)
nSamp <- 160
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 15e3
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(11)
libSizes=sample(colSums(countsTrapnell),nSamp,replace=TRUE)
dataTrapnell <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsTrapnell, nTags = nTags, group = grp, verbose = TRUE, params=paramsTrapnellAllDesignAveLogCPM, lib.size=libSizes, cpm="AveLogCPM")

plot(x=log(colSums(countsTrapnell)),y=colMeans(countsTrapnell==0), xlab="Log library size", ylab="Fraction of zeros")
points(x=log(colSums(dataTrapnell$counts)),y=colMeans(dataTrapnell$counts==0), xlab="Log library size", ylab="Fraction of zero's",col=2)
## zero ~ cpm
#par(mar=c(7,4,1,1))
plot(x=aveLogCPM(countsTrapnell),y=rowMeans(countsTrapnell==0), xlab="Average log CPM", ylab="Fraction of zeros", xlim=c(0.5,16),col=alpha(1,1/2),pch=19,cex=.2)
points(x=aveLogCPM(dataTrapnell$counts),y=rowMeans(dataTrapnell$counts==0), xlab="Average Log CPM", ylab="Fraction of zeroes",col=alpha(2,.8),pch=19,cex=.2)

selectedMethods <- c("edgeREMLibSizeDispFastOldFFilteredEdgeR", "limma_voomZeroFiltered", "DESeq2Zero_adjustedDf_posCountsNormZeroWeights", "DESeq2_poscounts", "edgeR", "edgeROldF", "limma_voom", "limma_voomFiltered", "metagenomeSeq","edgeRFiltered", "NODES")
#I have removed DESeq2 since it cannot normalize because all genes have at least one zero...
group=grp
#pvalsTrapnell = pval(dataTrapnell,method=selectedMethods,niter=200, mc.cores=2)
#save(pvalsTrapnell,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellPaper_v7_diffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellPaper_v7_diffZero.rda")
for(k in 1:length(pvalsTrapnell$pval)) pvalsTrapnell$pval[[k]][is.na(pvalsTrapnell$pval[[k]])]=1
for(k in 1:length(pvalsTrapnell$padj)) pvalsTrapnell$padj[[k]][is.na(pvalsTrapnell$padj[[k]])]=1
#scdePTrapnell=scde.pfun(dataTrapnell$counts,group=grp,mc.cores=1) #gives trouble in parallellization
#save(scdePTrapnell,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePTrapnell.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePTrapnell.rda")


## ZINB-WaVE: common dispersion
library(zinbwave)
counts=dataTrapnell$counts
design = model.matrix(~ grp)
library(doParallel)
registerDoParallel(2)
p=BiocParallel::DoparParam()
#zinb = zinbFit(counts, X=design, K=2, BPPARAM=p, commondispersion=TRUE)
#save(zinb,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbTrapnellCommonDiffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbTrapnellCommonDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_edgeR_common= zingeREdgeROwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_edgeR_common[is.na(zinbwave_edgeR_common[,"pval"]),"pval"]=1
zinbwave_edgeR_common[is.na(zinbwave_edgeR_common[,"padj"]),"padj"]=1


## ZINB-WaVE: gene dispersion
library(zinbwave)
counts=dataTrapnell$counts
design = model.matrix(~ grp)
library(doParallel)
registerDoParallel(2)
p=BiocParallel::DoparParam()
#zinb = zinbFit(counts, X=design, K=2, BPPARAM=p, commondispersion=FALSE)
#save(zinb,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbTrapnellGeneDiffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbTrapnellGeneDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_edgeR_gene= zingeREdgeROwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_edgeR_gene[is.na(zinbwave_edgeR_gene[,"pval"]),"pval"]=1
zinbwave_edgeR_gene[is.na(zinbwave_edgeR_gene[,"padj"]),"padj"]=1

## ZINB-WaVE DESeq2 common dispersion
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbTrapnellCommonDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_DESeq2_common= zingeRDESeq2OwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_DESeq2_common[is.na(zinbwave_DESeq2_common[,"pval"]),"pval"]=1
zinbwave_DESeq2_common[is.na(zinbwave_DESeq2_common[,"padj"]),"padj"]=1

## ZINB-WaVE DESeq2 gene dispersion
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zinbTrapnellGeneDiffZero.rda")
w <- calcZinbWaveWeights(zinb, counts)
zinbwave_DESeq2_gene= zingeRDESeq2OwnWeights.pfun(counts=counts, group=grp, w=w)
zinbwave_DESeq2_gene[is.na(zinbwave_DESeq2_gene[,"pval"]),"pval"]=1
zinbwave_DESeq2_gene[is.na(zinbwave_DESeq2_gene[,"padj"]),"padj"]=1



# edgeR hurdle
hurdleWeights = 1-(dataTrapnell$counts==0+0)
hurdleWeights[hurdleWeights==0]=1e-15
edgeRHurdleRes = zingeREdgeROwnWeights.pfun(counts=dataTrapnell$counts, group=grp, w=hurdleWeights)



#FDR-TPR
library(iCOBRA)
truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnell))
truthTrapnell[dataTrapnell$indDE,"status"]=1
cobraTrapnell <- COBRAData(pval =data.frame(
					    #limma_voom=pvalsTrapnell$pval$limma_voom,
					    #zingeR_limma_voom=pvalsTrapnell$pval$limma_voomZero,
					    #edgeR=pvalsTrapnell$pval$edgeR,
					    #zingeR_edgeR=pvalsTrapnell$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					    #DESeq2=pvalsTrapnell$pval$DESeq2,
					    DESeq2_poscounts=pvalsTrapnell$pval$DESeq2_poscounts,
					    zingeR_DESeq2=pvalsTrapnell$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					    #MAST=pvalsTrapnell$pval$MAST,
					    #metagenomeSeq=pvalsTrapnell$pval$metagenomeSeq,
					    #scde=scdePTrapnell[,"pval"],
					    #NODES=pvalsTrapnell$pval$NODES,
							#zinbwave_edgeR_common=zinbwave_edgeR_common[,"pval"],
							#zinbwave_edgeR_gene=zinbwave_edgeR_gene[,"pval"],
							zinbwave_DESeq2_common=zinbwave_DESeq2_common[,"pval"],
							zinbwave_DESeq2_gene=zinbwave_DESeq2_gene[,"pval"],
							#edgeRHurdle=edgeRHurdleRes[,"pval"],
					    row.names = rownames(dataTrapnell)),
		   padj = data.frame(
				     #limma_voom=pvalsTrapnell$padj$limma_voom,
				     #zingeR_limma_voom=pvalsTrapnell$padj$limma_voomZero,
				     #edgeR=pvalsTrapnell$padj$edgeR,
				     #zingeR_edgeR=pvalsTrapnell$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
				     #DESeq2=pvalsTrapnell$padj$DESeq2,
				     DESeq2_poscounts=pvalsTrapnell$padj$DESeq2_poscounts,
					   zingeR_DESeq2=pvalsTrapnell$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     #MAST=pvalsTrapnell$padj$MAST,
				     #metagenomeSeq=pvalsTrapnell$padj$metagenomeSeq,
				     #scde=scdePTrapnell[,"padj"],
				     #NODES=pvalsTrapnell$padj$NODES,
						 #zinbwave_edgeR_common=zinbwave_edgeR_common[,"padj"],
						 #zinbwave_edgeR_gene=zinbwave_edgeR_gene[,"padj"],
						 zinbwave_DESeq2_common=zinbwave_DESeq2_common[,"padj"],
						 zinbwave_DESeq2_gene=zinbwave_DESeq2_gene[,"padj"],
						 #edgeRHurdle=edgeRHurdleRes[,"padj"],
				     row.names = rownames(dataTrapnell)),
                   truth = truthTrapnell)
cobraperf <- calculate_performance(cobraTrapnell, binary_truth = "status")
colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black", zingeR_DESeq2_wald="steelblue", zinbwave_edgeR_common='gold', zinbwave_edgeR_gene="aquamarine3", edgeRHurdle="dodgerblue", zinbwave_DESeq2_common='gold', zinbwave_DESeq2_gene="aquamarine3")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnell.rda")
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellLimma.rda")
plot_fdrtprcurve(cobraplot, pointsize=2)

### check correct identification of dropout and NB zero
#zingeREdgeRWeights=zingeREdgeRWeights.pfun(counts=dataTrapnell$counts, group=grp)
#save(zingeREdgeRWeights,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zingeREdgeRWeightsTrapnellDiffZero.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/zingeREdgeRWeightsTrapnellDiffZero.rda")
library(rafalib) ; mypar()
par(mfrow=c(1,2))
hist(w[dataTrapnell$dropout==0], breaks=seq(0,1,0.05), main="ZINB-WaVE weight distribution for dropout on Trapnell simulation", cex.main=2/3)
hist(zingeREdgeRWeights[dataTrapnell$dropout==0], breaks=seq(0,1,0.05), main="zingeR weight distribution for dropout on Trapnell simulation", cex.main=2/3)

par(mfrow=c(1,2))
hist(w[dataTrapnell$dropout==1], breaks=seq(0,1,0.05), main="ZINB-WaVE weight distribution for NB zero on Trapnell simulation", cex.main=2/3)
hist(zingeREdgeRWeights[dataTrapnell$dropout==1], breaks=seq(0,1,0.05), main="zingeR weight distribution for NB zero on Trapnell simulation", cex.main=2/3)




# #high FC FDR-TPR
# set.seed(12)
# nSamp <- 160
# grp <- as.factor(rep(0:1, each = nSamp/2))
# nTags = 15e3
# DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
# fcSim=(3 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
# set.seed(11)
# libSizes=sample(colSums(countsTrapnell),nSamp,replace=TRUE)
# dataTrapnellHighFC <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsTrapnell, nTags = nTags, group = grp, verbose = TRUE, params=paramsTrapnellAllDesignAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")
# pvalsTrapnellHighFC = pval(dataTrapnellHighFC,method=selectedMethods,niter=200, mc.cores=2)
# #save(pvalsTrapnellHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellHighFCPaper.rda")
# load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellHighFCPaper.rda")
# for(k in 1:length(pvalsTrapnellHighFC$padj)) pvalsTrapnellHighFC$padj[[k]][is.na(pvalsTrapnellHighFC$padj[[k]])]=1
# for(k in 1:length(pvalsTrapnellHighFC$pval)) pvalsTrapnellHighFC$pval[[k]][is.na(pvalsTrapnellHighFC$pval[[k]])]=1
# #scdePHighFC=scde.pfun(dataTrapnellHighFC$counts,group=grp)
# #save(scdePHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePHighFCTrapnell.rda")
# load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePHighFCTrapnell.rda")
#
# truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnellHighFC))
# truthTrapnell[dataTrapnellHighFC$indDE,"status"]=1
# cobraTrapnell <- COBRAData(pval =data.frame(#edgeRRobust=pvalsTrapnell$pval$edgeR_robust,
# 					    limma_voom=pvalsTrapnellHighFC$pval$limma_voom,
# 					    zingeR_limma_voom=pvalsTrapnellHighFC$pval$limma_voomZero,
# 					    edgeR=pvalsTrapnellHighFC$pval$edgeR,
# 					    #edgeRFiltered=pvalsTrapnellHighFC$pval$edgeRFiltered,
# 					    zingeR_edgeR=pvalsTrapnellHighFC$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
# 					    DESeq2=pvalsTrapnellHighFC$pval$DESeq2,
# 					    DESeq2_poscounts=pvalsTrapnellHighFC$pval$DESeq2_poscounts,
# 					    zingeR_DESeq2=pvalsTrapnellHighFC$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
# 					    MAST=pvalsTrapnellHighFC$pval$MAST,
# 					    metagenomeSeq=pvalsTrapnellHighFC$pval$metagenomeSeq,
# 					    scde=scdePHighFC[,"pval"],
# 					    NODES=pvalsTrapnellHighFC$pval$NODES,
# 					    row.names = rownames(dataTrapnellHighFC)),
# 		   padj = data.frame(#edgeRRobust=pvalsTrapnell$padj$edgeR_robust,
# 				     limma_voom=pvalsTrapnellHighFC$padj$limma_voom,
# 				     zingeR_limma_voom=pvalsTrapnellHighFC$padj$limma_voomZero,
# 				     edgeR=pvalsTrapnellHighFC$padj$edgeR,
# 				     #edgeRFiltered=pvalsTrapnellHighFC$padj$edgeRFiltered,
# 				     zingeR_edgeR=pvalsTrapnellHighFC$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
# 				     DESeq2=pvalsTrapnellHighFC$padj$DESeq2,
# 				     DESeq2_poscounts=pvalsTrapnellHighFC$padj$DESeq2_poscounts,
# 				     zingeR_DESeq2=pvalsTrapnellHighFC$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
# 				     MAST=pvalsTrapnellHighFC$padj$MAST,
# 				     metagenomeSeq=pvalsTrapnellHighFC$padj$metagenomeSeq,
# 				     scde=scdePHighFC[,"padj"],
# 				     NODES=pvalsTrapnellHighFC$padj$NODES,
# 				     row.names = rownames(dataTrapnellHighFC)),
#                    truth = truthTrapnell)
# cobraperf <- calculate_performance(cobraTrapnell, binary_truth = "status")
# colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black")
# colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
# cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
# save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellHighFC.rda")
# plot_fdrtprcurve(cobraplot, pointsize=2)
#
#
# #### empirical vs simulated
# #pdf("~/Dropbox/phdKoen/singleCell/figures/supplementary/simulatedDataTrapnellNoNoise.pdf")
# #empirical BCV
# design=model.matrix(~timePoint)
# dEmp=DGEList(countsTrapnell)
# dEmp=calcNormFactors(dEmp)
# dEmp=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dEmp,design,interval=c(0,10)),design,prior.df=0)
# #plotBCVNoLeg(dEmp,col.tagwise=ifelse(rowSums(countsTrapnell>0)>16,1,rowSums(countsTrapnell>0))) #colored
# #simulated BCV
# design=model.matrix(~grp)
# dSim=DGEList(dataTrapnell$counts)
# dSim=calcNormFactors(dSim)
# dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)
#
# png("~/Dropbox/phdKoen/singleCell/figures/supplementary/simulationTrapnell.png", width=7,height=9, units="in", res=380)
# par(mfrow=c(2,2), mar=c(5,4,1,1))
# plotBCV(dEmp,main="Empirical", ylab="BCV", ylim=c(0,10))
# plotBCV(dSim, main="Simulated", ylab="BCV", ylim=c(0,10))
# ## zero ~ libSize
# plot(x=log(colSums(countsTrapnell)),y=colMeans(countsTrapnell==0), xlab="Log library size", ylab="Fraction of zero's")
# points(x=log(colSums(dataTrapnell$counts)),y=colMeans(dataTrapnell$counts==0), xlab="Log library size", ylab="Fraction of zero's",col=2)
# legend("topright",c("real","simulated"),pch=1,col=1:2, bty="n", cex=.6)
# ## zero ~ cpm
# plot(x=aveLogCPM(countsTrapnell),y=rowMeans(countsTrapnell==0), xlab="average log CPM", ylab="Fraction of zero's",col=alpha(1,1/2),pch=19,cex=.2)
# points(x=aveLogCPM(dataTrapnell$counts),y=rowMeans(dataTrapnell$counts==0), col=alpha(2,1/2),pch=19,cex=.2)
# legend("topright",c("real","simulated"),pch=19,col=1:2, bty="n", cex=.6)
# dev.off()
#
#
# ####### combine two-panel plots with same legend
# ## main plot
# library(cowplot)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslam.rda")
# islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.5), yaxisrange=c(0,0.75))
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnell.rda")
# trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.5), yaxisrange=c(0,0.75))
# #plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
# prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
#            trapnellPlot + theme(legend.position="none") + xlab("FDP"),
#            align = 'vh',
#            labels = c("a", "b"),
#            hjust = -1,
#            nrow = 1
#            )
# legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
# p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
# #png("~/Dropbox/phdKoen/singleCell/figures/scSimulation_composite_cutoff.png", width=7,height=8, units="in", res=300)
# pdf("~/Dropbox/phdKoen/singleCell/figures/scSimulation_composite_cutoff.pdf", width=7,height=8)
# p
# dev.off()
#
# ## main plot not cut off
# library(cowplot)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslam.rda")
# islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnell.rda")
# trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
# #plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
# prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
#            trapnellPlot + theme(legend.position="none") + xlab("FDP"),
#            align = 'vh',
#            labels = c("a", "b"),
#            hjust = -1,
#            nrow = 1
#            )
# legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
# p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
# png("~/Dropbox/phdKoen/singleCell/figures/supplementary/scSimulation_composite.png", width=7,height=8, units="in", res=300)
# p
# dev.off()
#
# ## including ZI limma-voom plot
# library(cowplot)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamLimma.rda")
# islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellLimma.rda")
# trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
# #plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
# prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
#            trapnellPlot + theme(legend.position="none") + xlab("FDP"),
#            align = 'vh',
#            labels = c("a", "b"),
#            hjust = -1,
#            nrow = 1
#            )
# legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
# p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
# png("~/Dropbox/phdKoen/singleCell/figures/supplementary/scSimulation_full_composite.png", width=7,height=8, units="in", res=300)
# p
# dev.off()
#
# ## high fold changes
# ## including ZI limma-voom plot
# library(cowplot)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamHighFC.rda")
# islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
# load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellHighFC.rda")
# trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
# #plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
# prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
#            trapnellPlot + theme(legend.position="none") + xlab("FDP"),
#            align = 'vh',
#            labels = c("a", "b"),
#            hjust = -1,
#            nrow = 1
#            )
# legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
# p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
# png("~/Dropbox/phdKoen/singleCell/figures/supplementary/highFCSimulation_composite.png", width=7,height=8, units="in", res=300)
# p
# dev.off()


#### ENGEL
# engelAssay <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE74596.rds")
# cellTypeEngel <- pData(engelAssay)[,"source_name_ch1"]
# countsEngel <- round(assay(experiments(engelAssay)$gene,"count"))
# countsEngelNKT2 <- countsEngel[,cellTypeEngel%in%c("Single_cell_RNA-seq_NKT2","Single_cell_RNA-seq_NKT1")]
# countsEngelNKT2 <- countsEngelNKT2[!rowSums(cpm(countsEngelNKT2)>2)<5,]
# celltypes=cellTypeEngel[cellTypeEngel%in%c("Single_cell_RNA-seq_NKT2","Single_cell_RNA-seq_NKT1")]
# paramsEngel = getDatasetZTNB(counts=countsEngelNKT2, design=model.matrix(~celltypes))
# save(paramsEngel,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsEngel.rda")


# ## simulate
# set.seed(18)
# nSamp <- 114
# grp <- as.factor(rep(0:1, each = nSamp/2))
# nTags = 15e3
# DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
# fcSim=(2 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
# set.seed(11)
# libSizes=sample(colSums(countsEngelNKT2),nSamp,replace=TRUE)
# dataEngel <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsEngelNKT2, nTags = nTags, group = grp, verbose = TRUE, params=paramsEngel, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")
# plot(x=log(colSums(countsEngelNKT2)),y=colMeans(countsEngelNKT2==0))
# points(x=log(colSums(dataEngel$counts)),y=colMeans(dataEngel$counts==0),col=2)

# plot(x=aveLogCPM(countsEngelNKT2),y=rowMeans(countsEngelNKT2==0),pch=16,cex=.2)
# points(x=aveLogCPM(dataEngel$counts),y=rowMeans(dataEngel$counts==0),pch=16,cex=.2, col=2)


# methodshlp=c("edgeREMLibSizeDispFastOldFFilteredEdgeR",  "edgeR")
# group=grp
# pvalsEngel = pval(dataEngel,method=methodshlp,niter=200, mc.cores=2)

# #FDR-TPR
# library(iCOBRA)
# truthEngel=data.frame(status=rep(0,nTags), row.names=rownames(dataEngel))
# truthEngel[dataEngel$indDE,"status"]=1
# cobraEngel <- COBRAData(pval =data.frame(
# 					    #limma_voom=pvalsTrapnell$pval$limma_voom,
# 					    #zingeR_limma_voom=pvalsTrapnell$pval$limma_voomZero,
# 					    edgeR=pvalsEngel$pval$edgeR,
# 					    #edgeRFiltered=pvalsTrapnell$pval$edgeRFiltered,
# 					    zingeR_edgeR=pvalsEngel$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
# 					    #DESeq2=pvalsTrapnell$pval$DESeq2,
# 					    #DESeq2_poscounts=pvalsTrapnell$pval$DESeq2_poscounts,
# 					    #zingeR_DESeq2=pvalsTrapnell$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
# 					    #metagenomeSeq=pvalsTrapnell$pval$metagenomeSeq,
# 					    #scde=scdePTrapnell[,"pval"],
# 					    #NODES=pvalsTrapnell$pval$NODES,
# 					    row.names = rownames(dataEngel)),
# 		   padj = data.frame(#edgeRRobust=pvalsTrapnell$padj$edgeR_robust,
# 				     #limma_voom=pvalsTrapnell$padj$limma_voom,
# 				     #zingeR_limma_voom=pvalsTrapnell$padj$limma_voomZero,
# 				     edgeR=pvalsEngel$padj$edgeR,
# 				     #edgeRFiltered=pvalsTrapnell$padj$edgeRFiltered,
# 				     zingeR_edgeR=pvalsEngel$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
# 				     #DESeq2=pvalsTrapnell$padj$DESeq2,
# 				     #DESeq2_poscounts=pvalsTrapnell$padj$DESeq2_poscounts,
# 					 #zingeR_DESeq2=pvalsTrapnell$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
# 				     #metagenomeSeq=pvalsTrapnell$padj$metagenomeSeq,
# 				     #scde=scdePTrapnell[,"padj"],
# 				     #NODES=pvalsTrapnell$padj$NODES,
# 				     row.names = rownames(dataEngel)),
#                    truth = truthEngel)
# cobraperf <- calculate_performance(cobraEngel, binary_truth = "status")
# colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black", zingeR_DESeq2_wald="steelblue", MASTold='gold')
# colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
# cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
# #save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnell.rda")
# #save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellLimma.rda")
# plot_fdrtprcurve(cobraplot, pointsize=2)
