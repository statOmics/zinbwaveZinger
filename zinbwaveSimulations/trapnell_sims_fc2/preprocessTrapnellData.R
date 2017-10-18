# RDS files downloaded from conquer repository http://imlspenticton.uzh.ch:3838/conquer/
library(MultiAssayExperiment)


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
