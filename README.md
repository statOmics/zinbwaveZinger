# Data analysis for the ZINB-WaVE / zingeR paper

This repository is designed to allow interested people to reproduce the results and figures of our paper called 'Unlocking RNA-seq tools for zero inflation and single-cell applications using observation weights'.   

For any questions or issues with the code on this repository, please use the "Issues" tab.

## Dependencies

To be able to run the code in this repo, it is required to have `R` (>=3.4) and the following packages.

### R packages

- [zingeR](https://github.com/statOmics/zingeR)  
- Seurat  
- ggplot2  
- RColorBrewer  
- devtools  
- scales  
- gridBase  
- grid  
- Matrix  
- dplyr  
- plyr  
- cowplot  
- Rtsne  
- FNN
- rARPACK 

### Bioconductor packages

- zinbwave  
- limma  
- edgeR  
- DESeq2  
- BiocParallel  
- doParallel 
- Biobase   
- scde  
- SingleCellExperiment  
- MAST
- [xCell](https://github.com/dviraran/xCell)  
- fgsea  
- clusterExperiment  
- GSEABase  


## Getting started

### Simulations

TODO

### Mock comparisons

To generate the plots related to the false positive rate control, run [FPR_mocks_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/tenx/FPR_mocks_tenx.Rmd) for the 10X genomics dataset and [FPR_mocks.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/usoskin/FPR_mocks.Rmd) for the Usoskin dataset. Additionally, to generate the plots related to PCER when the penalization parameter of ZINB-WaVE is varied, run [FPR_mocks_eps_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/epsilon_tenx/FPR_mocks_eps_tenx.Rmd) for the 10X genomics dataset and [FPR_mocks_eps_usoskin.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/epsilon_usoskin/FPR_mocks_eps_usoskin.Rmd) for the Usoskin dataset.

### Real data

To generate the plots related to the analysis of the 10x Genomics PBMC dataset, first run [createDataObject.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/createdata/createDataObject.Rmd) to create the data files. Then, to generate the data when the clustering is done using PCA, run both [de_seurat.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringPCA/de_seurat.Rmd) and [de_othermethods.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringPCA/de_othermethods.Rmd). There are two files to run instead of one unique file because packages `Seurat` and `zinbwave` load both many packages and R complains that there are too many packages loaded. To generate the plots when the clustering is done using ZINB-WaVE, run [dimredZinbwave.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringW/dimredZinbwave.Rmd), [clusterW.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringW/clusterW.Rmd), and then [de.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringW/de.Rmd). Finally, to generate the plots, run [plotPaper.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/plotPaper.Rmd).  

TODO: Usoskin.

### Time benchmarking

For each of the real datasets Islam, Usoskin, 10X genomics, respectively run [benchmark_islam.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/islam/benchmark_islam.Rmd), [benchmark_usoskin.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/usoskin/benchmark_usoskin.Rmd), [benchmark_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/tenx/benchmark_tenx.Rmd). Finally, run [benchmark_all.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/benchmark_all.Rmd).



