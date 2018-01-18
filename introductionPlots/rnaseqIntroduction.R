source("../zingeRsimulationFunctions/simulationHelpFunctions_v7_diffInZero.R")


#### no zero inflation simulation
library(Biobase)
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
#d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
#dOriginal=d
#dOriginalTrend=estimateGLMTrendedDisp(d,design)
d=estimateDisp(d,design,prior.df=0)
dNoZero=d
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
#d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
#d=estimateGLMTagwiseDisp(estimateGLMTrendedDisp(estimateGLMCommonDisp(d,design),design),design,prior.df=0)
#dZeroes=d
#dZeroesTrend=estimateGLMTrendedDisp(dZeroes,design)
d=estimateDisp(d,design,prior.df=0)
dZero=d
plotBCV(d, ylim=c(0,3))



## downweight introduced zeros
dW=DGEList(data)
dW=calcNormFactors(dW)
zeroId[zeroId==0]=1e-6
zeroId[zeroId==1]=.9999
dW$weights=zeroId
dW=estimateDisp(dW,design,prior.df=0)
plotBCV(dW)


## performance of edgeR or downweighted edgeR on ZI RNA-seq data
edgeRPerf = edgeR.pfun(dataZeroes$counts, group=grp, design)
edgeRWeightedPerf = edgeRWeightedOldF.pfun(dataZeroes$counts, group=grp, design, weights=zeroId)
library(iCOBRA)
truth=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truth[dataZeroes$indDE,"status"]=1
cobra = COBRAData(pval=data.frame(edgeR=edgeRPerf[,"pval"],
				  edgeRWeighted=edgeRWeightedPerf[,"pval"],
				  row.names=rownames(dataZeroes)),
		  padj=data.frame(edgeR=p.adjust(edgeRPerf[,"pval"],"fdr"),
				  edgeRWeighted=p.adjust(edgeRWeightedPerf[,"pval"],"fdr"),
				  row.names=rownames(dataZeroes)),
		  truth=truth)
cobraperf = calculate_performance(cobra, binary_truth="status")
colors = c(edgeR="red", edgeRWeighted="chocolate1")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
hlp = plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange = c(0, 0.4))$data


plot(x=hlp$FDR[hlp$method=="edgeR"], y=hlp$TPR[hlp$method=="edgeR"], type="l", xlim=c(0,0.4), col="red", lwd=2, xlab="False discovery rate", ylab="True positive rate")
lines(x=hlp$FDR[hlp$method=="edgeRWeighted"], y=hlp$TPR[hlp$method=="edgeRWeighted"], type="l", xlim=c(0,0.4), col="chocolate1", lwd=2, xlab="False discovery rate", ylab="True positive rate")


plotBCVIk = function (y, xlab = "Average log CPM", ylab = "Biological coefficient of variation",
    pch = 16, cex = 0.2, col.common = "red", col.trend = "blue",
    col.tagwise = "black", ...)
{
    ## copied from plotBCV function from edgeR package.
    if (!is(y, "DGEList"))
        stop("y must be a DGEList.")
    A <- y$AveLogCPM
    if (is.null(A))
        A <- aveLogCPM(y$counts, offset = getOffset(y))
    disp <- getDispersion(y)
    if (is.null(disp))
        stop("No dispersions to plot")
    if (attr(disp, "type") == "common")
        disp <- rep(disp, length = length(A))
    plot(A, sqrt(disp), xlab = xlab, ylab = ylab, type = "n",
        ...)
    labels <- cols <- lty <- pt <- NULL
    if (!is.null(y$tagwise.dispersion)) {
        points(A, sqrt(y$tagwise.dispersion), pch = pch, cex = cex,
            col = col.tagwise)
        labels <- c(labels, "Tagwise")
        cols <- c(cols, col.tagwise)
        lty <- c(lty, -1)
        pt <- c(pt, pch)
    }
    if (!is.null(y$common.dispersion)) {
        abline(h = sqrt(y$common.dispersion), col = col.common,
            lwd = 2)
        labels <- c(labels, "Common")
        cols <- c(cols, col.common)
        lty <- c(lty, 1)
        pt <- c(pt, -1)
    }
    if (!is.null(y$trended.dispersion)) {
        o <- order(A)
        lines(A[o], sqrt(y$trended.dispersion)[o], col = col.trend,
            lwd = 2)
        labels <- c(labels, "Trend")
        cols <- c(cols, col.trend)
        lty <- c(lty, 1)
        pt <- c(pt, -1)
    }
    #legend("topright", legend = labels, lty = lty, pch = pt,
     #   pt.cex = cex, lwd = 2, col = cols)
    invisible()
}


### add 2 real scRNA-seq datasets from conquer
library(MultiAssayExperiment) ; library(edgeR)
files=list.files("/Users/koenvandenberge/PhD_Data/singleCell/conquer/")
rdsFiles=files[grep(x=files,pattern=".rds$")]
# discard trimmed data
rdsFiles=rdsFiles[-grep(rdsFiles,pattern="trimmed.rds")]
# discard small datasets
rdsFiles=rdsFiles[-c(2:4)]
rdsNames=c("Buettner, 2015", "Deng, 2014", "Shalek, 2014", "Shalek, 2014", "Trapnell, 2014", "Trapnell, 2014", "Patel, 2014", "Kumar, 2014", "Kumar, 2014", "Guo, 2015", "Engel, 2016", "Meyer, 2016")
rdsNamesSub=rdsNames[-c(4,6,9)]
bcvList <- list()
for(i in 1:2){ #get first 2 datasets
    cat(i)
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))
    d=DGEList(countData[,1:10])
    d=edgeR::calcNormFactors(d)
    d=estimateDisp(d, prior.df=0)
    bcvList[[i]]=d
}




### composite plot
png("~/Dropbox/phdKoen/singleCell/zinbwaveZingeR/plots2/introBCVRNAseq.png",width=10,height=12,units="in",res=100)
#pdf("~/Dropbox/phdKoen/singleCell/figures/introBCVRNAseq.pdf",width=10,height=8)
#dev.new(width=10,height=12)
par(mfrow=c(3,2), mar=c(5,5,3,1))

plotBCVIk(bcvList[[1]], col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l", ylim=c(0,7), xlim=c(-2.3,15), main="scRNA-seq: Buettner (2015)")
mtext("a",side=3, at=-4.9, cex=4/3,font=2)
plotBCVIk(bcvList[[2]], col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l", ylim=c(0,7), xlim=c(-2.3,15), main="scRNA-seq: Deng (2014)")
mtext("b",side=3, at=-4.9, cex=4/3,font=2)
plotBCVIk(dNoZero, col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l", main="RNA-seq: simulated Bottomly", ylim=c(0,3))
mtext("c",side=3, at=-3.5, cex=4/3,font=2)
plotBCVIk(dZero, col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l", main="ZI RNA-seq: simulated Bottomly", ylim=c(0,3))
mtext("d",side=3, at=-3.5, cex=4/3,font=2)
plotBCVIk(dW, col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l", main="downweighted ZI RNA-seq: simulated Bottomly", ylim=c(0,3))
mtext("e",side=3, at=-3.5, cex=4/3,font=2)
plot(x=hlp$FDR[hlp$method=="edgeR"], y=hlp$TPR[hlp$method=="edgeR"], type="l", xlim=c(0,0.4), col="red", lwd=2, xlab="False discovery proportion", ylab="True positive rate", cex.axis=1.5, cex.lab=1.5, bty="l")
lines(x=hlp$FDR[hlp$method=="edgeRWeighted"], y=hlp$TPR[hlp$method=="edgeRWeighted"], type="l", xlim=c(0,0.4), col="chocolate1", lwd=2, xlab="False discovery rate", ylab="True positive rate")
legend("bottomright",c("edgeR","weighted edgeR"),bty="n",lty=1,lwd=2,col=c("red","chocolate1"), cex=1.33)
mtext("f",side=3,at=-0.05,cex=4/3,font=2)
dev.off()
