pathToParentFolder="~/Dropbox/phdKoen/singleCell/zinbwavezingerGitHub/zinbwaveZinger/"
source(paste0(pathToParentFolder,"zingeRsimulationFunctions/simulationHelpFunctions_v7_diffInZero.R"))
#source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

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
#libSize = sample(round(seq(4e6,5e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = bottomly, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = FALSE, lib.size=libSize, drop.low.lambda=TRUE)


edgeR <- function(counts, group, ylim = NULL, xlim = NULL){
  d <- DGEList(counts)
  d <- suppressWarnings(edgeR::calcNormFactors(d))
  design <- model.matrix(~group)
  d <- estimateDisp(d, design)
  plotBCV(d, ylim = ylim, main = 'edgeR', xlim = xlim)
  fit <- glmFit(d,design)
  lrt <- glmLRT(fit, coef = 2)
  pval <- lrt$table$PValue
  padj <- p.adjust(pval, "BH")
  cbind(pval = pval, padj = padj)
}

DESeq2 <- function(counts, group, ylim = NULL, xlim = NULL){
  colData <- data.frame(group = group)
  dse <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~group)
  colData(dse)$group <- as.factor(colData(dse)$group)
  dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
  dse <- estimateDispersions(dse)
  dse <- nbinomWaldTest(dse, betaPrior=TRUE)
  rr <- results(dse)
  cbind(pval = rr$pvalue, padj = rr$padj)
}


zinbwave_edgeR <- function(counts, group, zinb, ylim = NULL, xlim = NULL, main = 'ZINB-WaVE'){
  d=DGEList(counts)
  d=suppressWarnings(edgeR::calcNormFactors(d))
  design=model.matrix(~group)
  weights <- computeObservationalWeights(zinb, d$counts)
  d$weights <- weights
  d=estimateDisp(d, design)
  plotBCV(d, ylim = ylim, main = main)
  fit=glmFit(d,design)
  lrt=zingeR::glmWeightedF(fit,coef=2, independentFiltering = TRUE)
  cbind(pval = lrt$table$PValue, padj =lrt$table$padjFilter)
}

zinbwave_DESeq2 <- function(counts, group, zinb){
  colData=data.frame(group=group)
  design=model.matrix(~group)
  dse=DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~group)
  weights <- computeObservationalWeights(zinb, counts(dse))
  dimnames(weights)=NULL
  assays(dse)[["weights"]]=weights
  dse = DESeq2::estimateSizeFactors(dse, type="poscounts")
  dse = estimateDispersions(dse)
  dse = nbinomWaldTest(dse, betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2)
  res = results(dse)
  cbind(pval = res$pvalue, padj = res$padj)
}


library(zinbwave) ; library(DESeq2) ; library(doParallel) ; library(BiocParallel)
NCORES <- 2
registerDoParallel(NCORES)
register(DoparParam())
core <- SummarizedExperiment(dataNoZI$counts,
                             colData = data.frame(grp = grp))
zinbCommonNoZI <- zinbFit(core, X = '~ grp', commondispersion = TRUE, epsilon=1e12)


pvalEdger = edgeR(counts=dataNoZI$counts, group=grp)
pvalDESeq2 = DESeq2(counts=dataNoZI$counts, group=grp)
pvalZinbwaveEdgeR = zinbwave_edgeR(counts=dataNoZI$counts, group=grp, zinb=zinbCommonNoZI)
pvalZinbwaveDESeq2 = zinbwave_DESeq2(counts=dataNoZI$counts, group=grp, zinb=zinbCommonNoZI)

## set NA p-values to 1
pvalDESeq2[is.na(pvalDESeq2[,"pval"]),"pval"]=1
pvalDESeq2[is.na(pvalDESeq2[,"padj"]),"padj"]=1
pvalZinbwaveEdgeR[is.na(pvalZinbwaveEdgeR[,"pval"]),"pval"]=1
pvalZinbwaveEdgeR[is.na(pvalZinbwaveEdgeR[,"padj"]),"padj"]=1
pvalZinbwaveDESeq2[is.na(pvalZinbwaveDESeq2[,"pval"]),"pval"]=1
pvalZinbwaveDESeq2[is.na(pvalZinbwaveDESeq2[,"padj"]),"padj"]=1



library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(
					edgeR=pvalEdger[,"pval"],
					"ZINB-WaVE_edgeR"=pvalZinbwaveEdgeR[,"pval"],
					DESeq2=pvalDESeq2[,"pval"],
					"ZINB-WaVE_DESeq2"=pvalZinbwaveDESeq2[,"pval"],
					row.names = rownames(dataNoZI)),
		   padj = data.frame(
					edgeR=pvalEdger[,"padj"],
					"ZINB-WaVE_edgeR"=pvalZinbwaveEdgeR[,"padj"],
					DESeq2=pvalDESeq2[,"padj"],
					"ZINB-WaVE_DESeq2"=pvalZinbwaveDESeq2[,"padj"],
				     row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
colors=c(limmavoom="blue", "ZINB-WaVE_limmavoom_common"="steelblue", "ZINB-WaVE_limmavoom_genewise"="darkslategray3", edgeR="red", "ZINB-WaVE_edgeR"="salmon", "ZINB-WaVE_edgeR_genewise"="deeppink2",  DESeq2="brown",  "ZINB-WaVE_DESeq2"="darkseagreen", "ZINB-WaVE_DESeq2_genewise"="darkkhaki",  MAST="darkturquoise", metagenomeSeq="forestgreen", scde="grey", NODES="black",  Seurat="dodgerblue")
#iCOBRA converts '-' to '.'. Redo this.
cobraNames = sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)])
cobraNames = gsub(x=cobraNames, pattern=".", fixed=TRUE, replacement="-")
colsCobra=colors[match(cobraNames,names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
plot_fdrtprcurve(cobraplot, pointsize=2)


## compare total variance
#edgeR
counts=dataNoZI$counts
group=grp
d <- DGEList(counts)
d <- suppressWarnings(edgeR::calcNormFactors(d))
design <- model.matrix(~group)
d <- estimateDisp(d, design)
fitEdger=glmFit(d,design)
#zinbwave-edger
d=DGEList(counts)
d=suppressWarnings(edgeR::calcNormFactors(d))
design=model.matrix(~group)
zinb=zinbCommonNoZI
weights=computeObservationalWeights(zinb,d$counts)
d$weights <- weights
d=estimateDisp(d, design)
fitZinb=glmFit(d,design)
#zinger_edger
d=DGEList(counts)
d=suppressWarnings(edgeR::calcNormFactors(d))
design=model.matrix(~group)
wZinger = zingeR::zeroWeightsLS(counts, design, maxit=200)

zinbDefaultEps <- zinbFit(core, X = '~ grp', commondispersion = TRUE)
weightDefaultEps = computeObservationalWeights(zinbDefaultEps,counts)

## dispersion
plot(x=log(fitEdger$dispersion),y=log(fitZinb$dispersion), pch=16, cex=1/2,xlab="NB dispersion", ylab="ZINB dispersion", col=rowSums(counts==0)+1) ; abline(0,1,col=2)
legend("topleft",c("0 zeros","1 zero","2 zeros", "3 zeros"), col=1:4, bty="n", pch=16, cex=1/2)

## ZINB vs NB var
varTotalNB = fitEdger$fitted.values + fitEdger$dispersion+(fitEdger$fitted.values)^2
piZINB = t(getPi(zinbCommonNoZI))
varTotalZINB=(1-piZINB)*fitZinb$fitted.values*(1+fitZinb$fitted.values*(fitZinb$dispersion+piZINB))
library(RColorBrewer)
cols=brewer.pal(n=8,"Dark2")
plot(x=log(varTotalNB[,1]),y=log(varTotalZINB[,1]),pch=16,cex=1/4, xlab="NB total var edgeR", ylab="ZINB total var ZINB-edgeR")
points(x=log(varTotalNB[rowSums(counts==0)>0,1]),y=log(varTotalZINB[rowSums(counts==0)>0,1]),pch=16,cex=1/2,col=cols[rowSums(counts==0)])
abline(0,1,col=2)

legend("topleft",c("0 zeros","1 zero","2 zeros", "3 zeros"), col=1:4, bty="n", pch=16, cex=1/2)

### compare NB var for genes with zeros
varTotalNB2 = fitZinb$fitted.values + fitZinb$dispersion+(fitZinb$fitted.values)^2
zeroGenes = rowSums(counts==0)>0
plot(x=log(varTotalNB[zeroGenes,1]),y=log(varTotalNB2[zeroGenes,1]),pch=16,cex=1/1.2, xlab="log(NB total var edgeR)", ylab="log(NB total var ZINB-edgeR)",col=cols[rowSums(counts[zeroGenes,]==0)])
abline(0,1,col=2)
legend("topleft",c("1 zero","2 zeros", "3 zeros"), col=cols[1:5], bty="n", pch=16, cex=2)



################################################
#### add zeroes: here we see an obvious gain ###
################################################
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
set.seed(46)
samp=sample(1:length(zeroId),floor(length(zeroId)*propZeroes))
zeroId[samp]=0
zeroId[dataNoZI$counts==0]=1 #if it already was a zero it is not zero-inflated.
samp=samp[!samp%in%which(dataNoZI$counts==0)] #same
dataZeroes$counts = dataZeroes$counts*zeroId

core <- SummarizedExperiment(dataZeroes$counts,
                             colData = data.frame(grp = grp))
zinbCommonZero <- zinbFit(core, X = '~ grp', commondispersion = TRUE, epsilon=1e12)


pvalEdgerZero = edgeR(counts=dataZeroes$counts, group=grp)
pvalDESeq2Zero = DESeq2(counts=dataZeroes$counts, group=grp)
pvalZinbwaveEdgeRZero = zinbwave_edgeR(counts=dataZeroes$counts, group=grp, zinb=zinbCommonZero)
pvalZinbwaveDESeq2Zero = zinbwave_DESeq2(counts=dataZeroes$counts, group=grp, zinb=zinbCommonZero)

## set NA p-values to 1
pvalDESeq2Zero[is.na(pvalDESeq2Zero[,"pval"]),"pval"]=1
pvalDESeq2Zero[is.na(pvalDESeq2Zero[,"padj"]),"padj"]=1
pvalZinbwaveEdgeRZero[is.na(pvalZinbwaveEdgeRZero[,"pval"]),"pval"]=1
pvalZinbwaveEdgeRZero[is.na(pvalZinbwaveEdgeRZero[,"padj"]),"padj"]=1
pvalZinbwaveDESeq2Zero[is.na(pvalZinbwaveDESeq2Zero[,"pval"]),"pval"]=1
pvalZinbwaveDESeq2Zero[is.na(pvalZinbwaveDESeq2Zero[,"padj"]),"padj"]=1

library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataZeroes$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(
					edgeR=pvalEdgerZero[,"pval"],
					"ZINB-WaVE_edgeR"=pvalZinbwaveEdgeRZero[,"pval"],
					DESeq2=pvalDESeq2Zero[,"pval"],
					"ZINB-WaVE_DESeq2"=pvalZinbwaveDESeq2Zero[,"pval"],
					row.names = rownames(dataNoZI)),
		   padj = data.frame(
					edgeR=pvalEdgerZero[,"padj"],
					"ZINB-WaVE_edgeR"=pvalZinbwaveEdgeRZero[,"padj"],
					DESeq2=pvalDESeq2Zero[,"padj"],
					"ZINB-WaVE_DESeq2"=pvalZinbwaveDESeq2Zero[,"padj"],
				     row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
colors=c(limmavoom="blue", "ZINB-WaVE_limmavoom_common"="steelblue", "ZINB-WaVE_limmavoom_genewise"="darkslategray3", edgeR="red", "ZINB-WaVE_edgeR"="salmon", "ZINB-WaVE_edgeR_genewise"="deeppink2",  DESeq2="brown",  "ZINB-WaVE_DESeq2"="darkseagreen", "ZINB-WaVE_DESeq2_genewise"="darkkhaki",  MAST="darkturquoise", metagenomeSeq="forestgreen", scde="grey", NODES="black",  Seurat="dodgerblue")
#iCOBRA converts '-' to '.'. Redo this.
cobraNames = sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)])
cobraNames = gsub(x=cobraNames, pattern=".", fixed=TRUE, replacement="-")
colsCobra=colors[match(cobraNames,names(colors))]
cobraplotZero <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
plot_fdrtprcurve(cobraplotZero, pointsize=2)







## two-panel plot
library(cowplot)
prow <- plot_grid( plot_fdrtprcurve(cobraplotZero, pointsize=2) + theme(legend.position="none") + xlab("FDP"),
           plot_fdrtprcurve(cobraplot, pointsize=2) + xlab("FDP") + theme(legend.position="none"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(plot_fdrtprcurve(cobraplot, pointsize=2) + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/rnaseqComposite.png", width=7,height=8, units="in", res=300)
#pdf("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/rnaseqComposite.pdf", width=7,height=8)
p
dev.off()
