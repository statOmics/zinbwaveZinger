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

library(doParallel) ; library(zinbwave) ; library(BiocParallel)
NCORES <- 2
registerDoParallel(NCORES)
register(DoparParam())
design=model.matrix(~cellType)
zinbIslamRealData = zinbFit(islam, X=design, K=0, epsilon=1e12)


## trapnell
load("~/Dropbox/phdKoen/singleCell/zinbwavezingerGitHub/zinbwaveZinger/zinbwaveSimulations/trapnell_sims_fc2/countsTrapnellProcessed.rda")
timePoint=factor(c(rep(48,85),rep(72,64)))
design=model.matrix(~timePoint)
zinbTrapnellRealData = zinbFit(countsTrapnell, X=design, K=0, epsilon=1e12)

## get weights
wIslam = computeObservationalWeights(zinbIslamRealData, islam)
wTrapnell = computeObservationalWeights(zinbTrapnellRealData, countsTrapnell)


png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/postProbRealData.png", width=9,height=7, units="in", res=300)
par(mfrow=c(1,2))
hist(wIslam[islam==0], breaks=seq(0,1,0.05), main="Islam", xlab="Posterior probability", ylim=c(0,5e5))
hist(wTrapnell[countsTrapnell==0], breaks=seq(0,1,0.05), main="Trapnell", xlab="Posterior probability", ylim=c(0,5e5))
dev.off()
