setwd("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/zinbwaveZingerGithub/zinbwaveZinger/suppFigs/")
library(tweeDEseqCountData) ; library(Biobase)
data(pickrell)
load("../datasets/islam.rda")

##### sample comparison RNA-seq vs. scRNA-seq
png("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/zinbwaveZingeR/plots2/varRNAvsScRNASeq.png", width=7, height=5, units="in", res=100)
pickrell <- as.matrix(exprs(pickrell.eset))
par(mfrow=c(1,2),mar=c(5,5,3,1))
islam1 = islam[,1]
islam2 = islam[,3]
plot(x=log((islam1+islam2)*0.5), y=log2(islam1+1)-log2(islam2+1), pch=16, cex=1/3, ylim=c(-9,9), xlab="log(average count)", ylab="log2 fold change", main="scRNA-seq", cex.main=2, cex.lab=4/3, cex.axis=1.2)
abline(h=0,col=2,lwd=2)
pickrell1 = pickrell[,1]
pickrell2 = pickrell[,3]
plot(x=log((pickrell1+pickrell2)*0.5), y=log2(pickrell1+1)-log2(pickrell2+1), pch=16, cex=1/3, , ylim=c(-9,9), xlab="log(average count)", ylab="log2 fold change", main="RNA-seq", cex.main=2, cex.lab=4/3, cex.axis=1.2)
abline(h=0,col=2,lwd=2)
dev.off()


### conquer EDA over several datasets
library(MultiAssayExperiment) ; library(edgeR)
files=list.files("/Users/koenvandenberge/PhD_Data/singleCell/conquer/")
rdsFiles=files[grep(x=files,pattern=".rds$")]
# discard trimmed data
rdsFiles=rdsFiles[-grep(rdsFiles,pattern="trimmed.rds")]
# discard small datasets
rdsFiles=rdsFiles[-c(2:4)]
rdsNames=c("Buettner, 2015", "Deng, 2014", "Shalek, 2014", "Shalek, 2014", "Trapnell, 2014", "Trapnell, 2014", "Patel, 2014", "Kumar, 2014", "Kumar, 2014", "Guo, 2015", "Engel, 2016", "Meyer, 2016")
rdsNamesSub=rdsNames[-c(4,6,9)]

### P(zero) ~ libsize over datasets in conquer tool.
zeroLibsizeList <- list()
for(i in 1:length(rdsFiles)){
    if(i %in% c(4,6,9)) next
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))
    if(i==3){ #Shalek was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==5){ #Trapnell was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==8){ #Kumar was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    dat <- data.frame(logLibSize=log(colSums(countData)), zeroFraction=colMeans(countData==0))
    zeroLibsizeList[[i]]=dat
}
zeroLibsizeList <- zeroLibsizeList[!unlist(lapply(zeroLibsizeList,function(x) is.null(x)))]

png(filename="/Users/koenvandenberge/Dropbox/phdKoen/singleCell/zinbwaveZingeR/plots2/zeroLibSizeConquer.png",width=11,height=11, units="in", res=350)
par(mfrow=c(3,3), mar=c(3,3,1,1), oma=c(3,3,3,3))
for(i in 1:9){
    if(i%in%c(1,4,7)) par(mar=c(5,5,1,1)) else par(mar=c(5,4,1,1))
    plot(x=zeroLibsizeList[[i]][,1], y=zeroLibsizeList[[i]][,2], xlab="", ylab="", bty="l", main=rdsNamesSub[i], pch=16, cex=2/3, xlim=c(4.5,18), ylim=c(0,1))
    mtext("Fraction of zeros",side=2, outer=TRUE, cex=2)
    mtext("Log library size",side=1, outer=TRUE, cex=2)
}
dev.off()

### P(zero) ~ aveLogCPM over datasets in conquer tool.
zeroCpmList <- list()
for(i in 1:length(rdsFiles)){
    if(i %in% c(4,6,9)) next
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))
    if(i==3){ #Shalek was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==5){ #Trapnell was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==8){ #Kumar was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    dat <- data.frame(avCpm=aveLogCPM(countData), zeroFraction=rowMeans(countData==0))
    zeroCpmList[[i]]=dat
}
zeroCpmList <- zeroCpmList[!unlist(lapply(zeroCpmList,function(x) is.null(x)))]

png(filename="/Users/koenvandenberge/Dropbox/phdKoen/singleCell/zinbwaveZingeR/plots2/zeroCpmConquer.png",width=11,height=11, units="in", res=330)
par(mfrow=c(3,3), mar=c(3,3,1,1), oma=c(3,3,3,3))
for(i in 1:9){
    if(i%in%c(1,4,7)) par(mar=c(5,5,1,1)) else par(mar=c(5,4,1,1))
    plot(x=zeroCpmList[[i]][,1], y=zeroCpmList[[i]][,2], xlab="", ylab="", bty="l", main=rdsNamesSub[i], pch=16, cex=1/3, xlim=c(-3,16), ylim=c(0,1))
    mtext("Fraction of zeros",side=2, outer=TRUE, cex=2)
    mtext("Average log CPM",side=1, outer=TRUE, cex=2)
}
dev.off()

### BCV plot based on 10 samples

### BCV over datasets in conquer tool.
bcvList <- list()
for(i in 1:length(rdsFiles)){
    cat(i)
    if(i %in% c(4,6,9)) next
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))
    if(i==3){ #Shalek was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==5){ #Trapnell was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==8){ #Kumar was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    d=DGEList(countData[,1:10])
    d=edgeR::calcNormFactors(d)
    d=estimateDisp(d, prior.df=0)
    bcvList[[i]]=d
}
bcvList <- bcvList[!unlist(lapply(bcvList,function(x) is.null(x)))]

png(filename="/Users/koenvandenberge/Dropbox/phdKoen/singleCell/zinbwaveZingeR/plots2/bcvConquer.png",width=11,height=11, res=330, units="in")
par(mfrow=c(3,3), mar=c(3,3,1,1), oma=c(3,3,3,3))
for(i in 1:9){
    if(i%in%c(1,4,7)) par(mar=c(5,5,1,1)) else par(mar=c(5,4,1,1))
    #plot(x=zeroLibsizeList[[i]][,1], y=zeroLibsizeList[[i]][,2], xlab="", ylab="", bty="l", main=rdsNamesSub[i], pch=16, cex=1/3)
    plotBCV(bcvList[[i]], xlab="", ylab="", main=rdsNamesSub[i], xlim=c(-3,16), ylim=c(0,8.25))
    mtext("Biological coefficient of variation",side=2, outer=TRUE, cex=2)
    mtext("Average log CPM",side=1, outer=TRUE, cex=2)
}
dev.off()
