require(genefilter)
require(edgeR)

group = grp
zinb=zinb_c
d=DGEList(counts)
d=suppressWarnings(calcNormFactors(d))
design=model.matrix(~group)
weights <- computeZinbwaveWeights(zinb, d$counts)
d$weights <- weights
d=estimateDisp(d, design)
fit=glmFit(d,design)
lrt=glmWeightedF(fit, coef=2, independentFiltering = FALSE)

baseMean = filter = rowMedians(fit$fitted.values)
pValue = lrt$table$PValue


plot(rank(filter)/length(filter), -log10(pValue), pch=16, cex=0.45)
abline(h= -log10(0.01))
abline(v=0.02)

filterChoices = data.frame(
  `median` = rowMedians(fit$fitted.values),
  `mu_medians` = rowMedians(t(getMu(zinb)))
)

rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.05, filter=f, test=pValue, theta=theta, method="BH"))

library("RColorBrewer")
myColours = brewer.pal(ncol(filterChoices), "Set1")

matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("bottomleft", legend=colnames(filterChoices), fill=myColours)





alpha=0.05
pAdjustMethod="BH"

lowerQuantile <- mean(filter == 0)
if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
theta <- seq(lowerQuantile, upperQuantile, length=50)

filtPadj <- filtered_p(filter=filter, test=pValue,
                       theta=theta, method=pAdjustMethod)
numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
plot(numRej)
# prevent over-aggressive filtering when all genes are null,
# by requiring the max number of rejections is above a fitted curve.
# If the max number of rejection is not greater than 10, then don't
# perform independent filtering at all.
lo.fit <- lowess(numRej ~ theta, f=1/5)
plot(lo.fit)
if (max(numRej) <= 10) {
  j <- 1
} else {
  residual <- if (all(numRej==0)) {
    0
  } else {
    numRej[numRej > 0] - lo.fit$y[numRej > 0]
  }
  thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
  j <- if (any(numRej > thresh)) {
    which(numRej > thresh)[1]
  } else {
    1
  }
}
padj <- filtPadj[, j, drop=TRUE]
cutoffs <- quantile(filter, theta)
filterThreshold <- cutoffs[j]
filterNumRej <- data.frame(theta=theta, numRej=numRej)
filterTheta <- theta[j]
filterTheta
