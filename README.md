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

### Time benchmarking

For each of the real datasets Islam, Usoskin, 10X genomics, respectively run [benchmark_islam.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/islam/benchmark_islam.Rmd), [benchmark_usoskin.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/usoskin/benchmark_usoskin.Rmd), [benchmark_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/tenx/benchmark_tenx.Rmd). Finally, run [benchmark_all.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/benchmark_all.Rmd).

### Real data

To generate the plots related to the false positive rate control, run [FPR_mocks_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/tenx/FPR_mocks_tenx.Rmd) for the 10X genomics dataset and [FPR_mocks.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/usoskin/FPR_mocks.Rmd) for the Usoskin dataset. Additionally, to generate the plots related to PCER when the penalization parameter of ZINB-WaVE is varied, run [FPR_mocks_eps_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/epsilon_tenx/FPR_mocks_eps_tenx.Rmd) for the 10X genomics dataset and [FPR_mocks_eps_usoskin.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/epsilon_usoskin/FPR_mocks_eps_usoskin.Rmd) for the Usoskin dataset.



### Simulations

TODO

