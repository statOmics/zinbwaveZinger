pathToParentFolder="~/Dropbox/phdKoen/singleCell/zinbwavezingerGitHub/zinbwaveZinger/"
source(paste0(pathToParentFolder,"zingeRsimulationFunctions/simulationHelpFunctions_v7_diffInZero.R"))
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

#########################################################
############# RNA-Seq data simulation ###################
#########################################################
#### no zero inflation simulation
library(Biobase) ; library(genefilter)
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/bottomly_eset.RData")
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


edgeREMLibSizeDispFastOldFFilteredEdgeR.pfun = function(counts, group, design=NULL, mc.cores=2, niter=200){
  library(edgeR) ; library(genefilter)
  d <- DGEList(counts = counts, group = group )
  d <- edgeR::calcNormFactors(d)
  design = model.matrix(~ group)
  zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
  d$weights = zeroWeights
  d=estimateDisp(d,design)
  edger.fit <- glmFit(d, design) #uses weights
  edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  edger.lrt <- glmWeightedF(edger.fit,coef=2,test="F")
  lfc <- edger.lrt$table$logFC
  pval <- edger.lrt$table$PValue
  baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
  hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  padj <- hlp$padj
  out=cbind(pval,padj,lfc)
  return(out)
}





pvalEdger = edgeR.pfun(counts=dataNoZI$counts, group=grp)
pvalDESeq2 = DESeq2.pfun(counts=dataNoZI$counts, group=grp)
pvalZingerEdgeR = edgeREMLibSizeDispFastOldFFilteredEdgeR.pfun(counts=dataNoZI$counts, group=grp, niter=200)

library(zinbwave)
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
zinbwave_edgeR <- function(counts, group, zinb, ylim = NULL, xlim = NULL, main = 'ZINB-WaVE'){
  d=DGEList(counts)
  d=suppressWarnings(calcNormFactors(d))
  design=model.matrix(~group)
  weights <- computeZinbwaveWeights(zinb, d$counts)
  weights[weights==0] = 1e-15 #errors with zero weights
  d$weights <- weights
  d=estimateDisp(d, design)
  plotBCV(d, ylim = ylim, main = main)
  fit=glmFit(d,design)
  lrt=zingeR::glmWeightedF(fit,coef=2, independentFiltering = TRUE)
  cbind(pval = lrt$table$PValue, padj =lrt$table$padjFilter)
}
core <- SummarizedExperiment(dataNoZI$counts,
                             colData = data.frame(grp = grp))
zinb_c <- zinbFit(core, X = '~ grp', commondispersion = TRUE)
pvalZinbwaveEdger = zinbwave_edgeR(dataNoZI$counts, grp, zinb_c, ylim=ylim, main = 'ZINB-WaVE, common dispersion', xlim = xlim)

pvalZingerEdgeR[is.na(pvalZingerEdgeR[,"pval"]),"pval"]=1
pvalZingerEdgeR[is.na(pvalZingerEdgeR[,"padj"]),"padj"]=1
pvalZinbwaveEdger[is.na(pvalZinbwaveEdger[,"pval"]),"pval"]=1
pvalZinbwaveEdger[is.na(pvalZinbwaveEdger[,"padj"]),"padj"]=1



library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(
					#limma_voom=pvalsNoZI$pval$limma_voom,
					edgeR=pvalEdger[,"pval"],
					zingeR_edgeR=pvalZingerEdgeR[,"pval"],
					zinbwave_edgeR=pvalZinbwaveEdger[,"pval"],
					#DESeq2=pvalsNoZI$pval$DESeq2,
					#DESeq2_poscounts=pvalsNoZI$pval$DESeq2_poscounts,
					#zingeR_DESeq2=pvalsNoZI$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					#MAST=pvalsNoZI$pval$MAST,
					#NODES=pvalsNoZI$pval$NODES,
					#scde=scdeP[,"pval"],
					row.names = rownames(dataNoZI)),
		   padj = data.frame(
				     #limma_voom=pvalsNoZI$padj$limma_voom,
						edgeR=pvalEdger[,"padj"],
						zingeR_edgeR=pvalZingerEdgeR[,"padj"],
						zinbwave_edgeR=pvalZinbwaveEdger[,"padj"],
				     #DESeq2=pvalsNoZI$padj$DESeq2,
				     #DESeq2_poscounts=pvalsNoZI$padj$DESeq2_poscounts,
				     #zingeR_DESeq2=pvalsNoZI$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     #MAST=pvalsNoZI$padj$MAST,
				     #NODES=pvalsNoZI$padj$NODES,
				     #scde=scdeP[,"padj"],
				     row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
colors=c(limma_voom="blue", limma_voomZero="steelblue", edgeRTruth="chocolate1", edgeR="red", zingeR_edgeR="salmon", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", MAST="darkturquoise", metagenomeSeq="green", zinbwave_edgeR="grey", NODES="black")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))

# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
set.seed(46)
samp=sample(1:length(zeroId),floor(length(zeroId)*propZeroes))
zeroId[samp]=0
zeroId[dataNoZI$counts==0]=1 #if it already was a zero it is not zero-inflated.
samp=samp[!samp%in%which(dataNoZI$counts==0)] #same
dataZeroes$counts = dataZeroes$counts*zeroId




pvalEdger = edgeR.pfun(counts=dataZeroes$counts, group=grp)
pvalDESeq2 = DESeq2.pfun(counts=dataZeroes$counts, group=grp)
pvalZingerEdgeR = edgeREMLibSizeDispFastOldFFilteredEdgeR.pfun(counts=dataZeroes$counts, group=grp, niter=200)

library(zinbwave)
library(doParallel)
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
zinbwave_edgeR <- function(counts, group, zinb, ylim = NULL, xlim = NULL, main = 'ZINB-WaVE'){
  d=DGEList(counts)
  d=suppressWarnings(calcNormFactors(d))
  design=model.matrix(~group)
  weights <- computeZinbwaveWeights(zinb, d$counts)
  weights[weights==0] = 1e-15 #errors with zero weights
  d$weights <- weights
  d=estimateDisp(d, design)
  plotBCV(d, ylim = ylim, main = main)
  fit=glmFit(d,design)
  lrt=zingeR::glmWeightedF(fit,coef=2, independentFiltering = TRUE)
  cbind(pval = lrt$table$PValue, padj =lrt$table$padjFilter)
}
core <- SummarizedExperiment(dataZeroes$counts,
                             colData = data.frame(grp = grp))
zinb_c <- zinbFit(core, X = '~ grp', commondispersion = TRUE)
pvalZinbwaveEdger = zinbwave_edgeR(dataZeroes$counts, grp, zinb_c, ylim=ylim, main = 'ZINB-WaVE, common dispersion', xlim = xlim)

pvalZingerEdgeR[is.na(pvalZingerEdgeR[,"pval"]),"pval"]=1
pvalZingerEdgeR[is.na(pvalZingerEdgeR[,"padj"]),"padj"]=1
pvalZinbwaveEdger[is.na(pvalZinbwaveEdger[,"pval"]),"pval"]=1
pvalZinbwaveEdger[is.na(pvalZinbwaveEdger[,"padj"]),"padj"]=1



library(iCOBRA)
truthZeroes=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truthZeroes[dataZeroes$indDE,"status"]=1
cobraZero <- COBRAData(pval =data.frame(
					#limma_voom=pvalsNoZI$pval$limma_voom,
					edgeR=pvalEdger[,"pval"],
					zingeR_edgeR=pvalZingerEdgeR[,"pval"],
					zinbwave_edgeR=pvalZinbwaveEdger[,"pval"],
					#DESeq2=pvalsNoZI$pval$DESeq2,
					#DESeq2_poscounts=pvalsNoZI$pval$DESeq2_poscounts,
					#zingeR_DESeq2=pvalsNoZI$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					#MAST=pvalsNoZI$pval$MAST,
					#NODES=pvalsNoZI$pval$NODES,
					#scde=scdeP[,"pval"],
					row.names = rownames(dataZeroes)),
		   padj = data.frame(
				     #limma_voom=pvalsNoZI$padj$limma_voom,
						edgeR=pvalEdger[,"padj"],
						zingeR_edgeR=pvalZingerEdgeR[,"padj"],
						zinbwave_edgeR=pvalZinbwaveEdger[,"padj"],
				     #DESeq2=pvalsNoZI$padj$DESeq2,
				     #DESeq2_poscounts=pvalsNoZI$padj$DESeq2_poscounts,
				     #zingeR_DESeq2=pvalsNoZI$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     #MAST=pvalsNoZI$padj$MAST,
				     #NODES=pvalsNoZI$padj$NODES,
				     #scde=scdeP[,"padj"],
				     row.names = rownames(dataZeroes)),
                   truth = truthZeroes)
cobraperf <- calculate_performance(cobraZero, binary_truth = "status")
colors=c(limma_voom="blue", limma_voomZero="steelblue", edgeRTruth="chocolate1", edgeR="red", zingeR_edgeR="salmon", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", MAST="darkturquoise", metagenomeSeq="green", zinbwave_edgeR="grey", edgernozi="black")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplotZero <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
plot_fdrtprcurve(cobraplotZero, pointsize=2)
#plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))







## two-panel plot
library(cowplot)
prow <- plot_grid( plot_fdrtprcurve(cobraplotZero, pointsize=2) + theme(legend.position="none") + xlab("FDP"),
           plot_fdrtprcurve(cobraplot, pointsize=2) + xlab("FDP") + theme(legend.position="none"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(plot_fdrtprcurve(cobraplotZero, pointsize=2) + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
#png("~/Dropbox/phdKoen/singleCell/figures/supplementary/RNASeq_composite.png", width=7,height=8, units="in", res=300)
pdf("~/Dropbox/phdKoen/singleCell/zinbwaveZingerGithub/updates/rnaseqCompositeDraft.pdf", width=7,height=8)
p
dev.off()


## plot a histogram of the weights for the introduced zeroes
weights=zeroWeightsLibSizeDispFast(counts=dataZeroes$counts, design=model.matrix(~grp), maxit=100, plotW=TRUE)
hist(weights[samp[dataNoZI$counts[samp]!=0]], xlab="weight", xlim=c(0,1), main="")


## ROC curve for identifying excess zeros
weights=zeroWeightsLibSizeDispFast(counts=dataZeroes$counts, design=model.matrix(~grp), maxit=200, plotW=TRUE)
hist(weights[samp], xlab="weight", xlim=c(0,1), main="")

png("~/Dropbox/phdKoen/singleCell/figures/supplementary/rocExcessZerosRNASeqSimulation_composite.png", width=8,height=6, units="in", res=300)
pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
falses=which(dataNoZI$counts==0)
tpr=fpr=vector(length=length(pvalSeq))
for(i in 1:length(pvalSeq)){
    excessID <- which(weights<=pvalSeq[i])
    tpr[i] <- mean(samp%in%excessID)
    fpr[i] <- mean(falses%in%excessID)
}
par(mfrow=c(1,2))
plot(x=fpr,y=tpr,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")
points(x=fpr[pvalSeq==1e-5],y=tpr[pvalSeq==1e-5],col=2,pch=19)
points(x=fpr[pvalSeq==1e-2],y=tpr[pvalSeq==1e-2],col=2,pch=19)
points(x=fpr[519],y=tpr[519],col=2,pch=19) #w=0.05
text(x=fpr[pvalSeq==1e-5]+0.13,y=tpr[pvalSeq==1e-5],labels="w=1e-5")
text(x=fpr[pvalSeq==1e-2]+0.13,y=tpr[pvalSeq==1e-2],labels="w=1e-2")
text(x=fpr[519]+0.13,y=tpr[519],labels="w=5e-2")

## ROC curve stratified by average expression
d=DGEList(dataZeroes$counts)
d=calcNormFactors(d)
acpm=aveLogCPM(d)
library(Hmisc)
cuts=cut2(acpm,g=5)
tpr2=fpr2=list()
for(k in 1:length(levels(cuts))){
    tpr2[[k]]=vector(length=length(pvalSeq))
    fpr2[[k]]=vector(length=length(pvalSeq))
}
for(j in 1:length(levels(cuts))){
    groupID <- cuts==levels(cuts)[j]
    wSub=weights[groupID,]
    falses=which(dataNoZI$counts[cuts==levels(cuts)[j]]==0)
    trues=which(zeroId[cuts==levels(cuts)[j],]==0)
    for(i in 1:length(pvalSeq)){
        excessID <- which(wSub<=pvalSeq[i])
	tpr2[[j]][i] <- mean(trues%in%excessID)
    	fpr2[[j]][i] <- mean(falses%in%excessID)
    }
}

cols=colorRampPalette(c("skyblue","darkblue"))(length(levels(cuts)))
plot(x=fpr2[[1]],y=tpr2[[1]],col="steelblue",type="n", xlab="False positive rate", ylab="True positive rate")
for(k in 1:length(levels(cuts))) lines(x=fpr2[[k]],y=tpr2[[k]],col=cols[k], lty=k)
legend("bottomright",levels(cuts),col=cols,lty=1:length(levels(cuts)))
dev.off()


### investigate Cook's distance cut-off from DESeq2
library(DESeq2)
colData <- data.frame(group)
dse <- DESeqDataSetFromMatrix(countData = dataZeroes$counts, colData = colData, design = ~ group)
colData(dse)$group <- as.factor(colData(dse)$group)
#dse <- DESeq(dse, betaPrior=TRUE)
dse = estimateSizeFactors(dse)
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, betaPrior=TRUE)
res <- results(dse)


noExcessZeroGenes = apply(zeroId,1,function(row) all(row==1))
mean(noExcessZeroGenes) #60% have no introduced zeros

plot(density(log(mcols(dse)$maxCooks[noExcessZeroGenes])), main="RNA-seq: log of maximum Cook's for a gene")
lines(density(log(mcols(dse)$maxCooks[!noExcessZeroGenes])),col=2)
abline(v=log(qf(p=.99,df1=2,df2=8)))
legend("topright",c("no excess gene","excess gene"), col=1:2, lty=1)

excessZeroCooks = assays(dse)[["cooks"]][zeroId==0]
noExcessZeroCooks = assays(dse)[["cooks"]][zeroId==1]
plot(density(log(noExcessZeroCooks)), "Observation level Cook's distance")
lines(density(log(excessZeroCooks)),col=2)
abline(v=log(qf(p=.99,df1=2,df2=8)))
legend("topleft",c("no excess zero","all counts"), col=1:2, lty=1)



### scRNA-seq simulations
group=grp
library(DESeq2)
colData <- data.frame(group)
dse <- DESeqDataSetFromMatrix(countData = dataIslamAllAveLogCPM$counts, colData = colData, design = ~ group)
colData(dse)$group <- as.factor(colData(dse)$group)
#dse <- DESeq(dse, betaPrior=TRUE)
dse = estimateSizeFactors(dse)
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, betaPrior=TRUE)
res <- results(dse, cooksCutoff=Inf)

plot(density(log(mcols(dse)$maxCooks)), main="scRNA-seq: log of maximum Cook's for a gene")
abline(v=log(qf(p=.99,df1=2,df2=78)))
mean(mcols(dse)$maxCooks > qf(p=.99,df1=2,df2=78))



########### OLD CODE

#### no zero inflation simulation
#library(tweeDEseqCountData)
#data(pickrell)
#pickrell <- as.matrix(exprs(pickrell.eset))
nSamp <- 10
nTags <- 20e3
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(5e6,8e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = pickrell, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE)
selectedMethods <- c("edgeR_robust", "limma_voom","edgeR", "DESeq2", "edgeREMLibSizeOldF")
pvalsNoZI <- pval(dataNoZI, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
hlp=edgeREMLibSizeFastOldF.pfun(counts=dataNoZI$counts, group=grp, niter=1e3)
#### performance curves using iCOBRA
#rnaSeqPerformanceNoZeroes.pdf
library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(edgeRRobust=pvalsNoZI$pval$edgeR_robust, limma_voom=pvalsNoZI$pval$limma_voom, edgeR=pvalsNoZI$pval$edgeR,  edgeREMLibSize=pvalsNoZI$pval$edgeREMLibSize, DESeq2=pvalsNoZI$pval$DESeq2, edgeREMFast=hlp[,1], row.names = rownames(dataNoZI)),
		   padj = data.frame(edgeRRobust=pvalsNoZI$padj$edgeR_robust, limma_voom=pvalsNoZI$padj$limma_voom, edgeR=pvalsNoZI$padj$edgeR, edgeREMLibSize=pvalsNoZI$padj$edgeREMLibSize, DESeq2=pvalsNoZI$padj$DESeq2,edgeREMFast=hlp[,2], row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))


# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.1
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
zeroId[sample(1:length(zeroId),floor(length(zeroId)*propZeroes))]=0
dataZeroes$counts = dataZeroes$counts*zeroId
genesWithAddedZero <- which(rowSums(zeroId)<nSamp)
pvalsZeroes = pval(dataZeroes, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
hlp=edgeREMLibSizeFastOldF.pfun(counts=dataZeroes$counts, group=grp, niter=1e3)
##performance curves
#rnaSeqPerformanceWithZeroes.pdf
truthZero=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truthZero[dataZeroes$indDE,"status"]=1
cobraZero <- COBRAData(pval =data.frame(edgeRRobust=pvalsZeroes$pval$edgeR_robust, limma_voom=pvalsZeroes$pval$limma_voom, edgeR=pvalsZeroes$pval$edgeR, edgeREMLibSize=pvalsZeroes$pval$edgeREMLibSize, DESeq2=pvalsZeroes$pval$DESeq2, edgeREMFast=hlp[,1], row.names = rownames(dataZeroes)),
		   padj = data.frame(edgeRRobust=pvalsZeroes$padj$edgeR_robust, limma_voom=pvalsZeroes$padj$limma_voom, edgeR=pvalsZeroes$padj$edgeR, edgeREMLibSize=pvalsZeroes$padj$edgeREMLibSize, DESeq2=pvalsZeroes$padj$DESeq2 , edgeREMFast=hlp[,2] , row.names = rownames(dataZeroes)),
                   truth = truthZero)
cobraperf <- calculate_performance(cobraZero, binary_truth = "status")
cobraplotZero <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplotZero)
plot_roc(cobraplotZero,xaxisrange=c(0,0.1), yaxisrange=c(0,0.6))

#performanceRNAseq.pdf
p1=plot_fdrtprcurve(cobraplot)
p2=plot_fdrtprcurve(cobraplotZero)
library(scater)
multiplot(p1,p2)
