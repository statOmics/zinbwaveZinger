# Data analysis for the ZINB-WaVE / zingeR paper

This repository is designed to allow interested people to reproduce the results and figures of our paper called 'Observation weights unlock bulk RNA-seq tools for zero inflation and single-cell applications', currently on bioRxiv at https://www.biorxiv.org/content/early/2018/01/18/250126.
All code in the repository is distributed under the GPL-3 license.

For examples on how to use the method, please see the [zinbwave vignette](https://github.com/drisso/zinbwave/blob/master/vignettes/intro.Rmd#differential-expression) or the [zingeR vignette](https://github.com/statOmics/zingeR/blob/master/vignettes/zingeRVignette_v2.Rmd#differential-expression-analysis) on how to use the observation weights in your analysis.

For any questions or issues with the code on this repository, please use the "Issues" tab.

## Dependencies

To be able to run the code in this repo, it is required to have `R` (>=3.4) and the following packages.

### R packages

- [zingeR](https://github.com/statOmics/zingeR)  
- countsimQC
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

The functions used to estimate parameters based on real scRNA-seq data, and to simulate the expression counts can be found in the [simulationHelpFunctions_v7_diffInZero.R](https://github.com/statOmics/zinbwaveZinger/blob/master/zingeRsimulationFunctions/simulationHelpFunctions_v7_diffInZero.R) file in the zingeRsimulationFunctions folder. This framework has been used to simulate all scRNA-seq datasets.
We have also simulated a bulk RNA-seq dataset and code for this simulation can be found in the [rnaseqSim.R](https://github.com/statOmics/zinbwaveZinger/blob/master/zinbwaveSimulations/RNASeq/rnaseqSim.R) file.
The quality of the simulated datasets has been evaluated using the [countsimQC package](https://github.com/csoneson/countsimQC), and code for this evaluation can be found in the [zinbwaveSimulations/evaluateSimulatedData](https://github.com/statOmics/zinbwaveZinger/tree/master/zinbwaveSimulations/evaluateSimulatedData) folder.
The code for the evaluations on the simulated Islam, Trapnell and 10X datasets can be found in the respective [islam_sims_fc2](https://github.com/statOmics/zinbwaveZinger/tree/master/zinbwaveSimulations/islam_sims_fc2), [trapnell_sims_fc2](https://github.com/statOmics/zinbwaveZinger/tree/master/zinbwaveSimulations/trapnell_sims_fc2) and [tenX_sims_fc2](https://github.com/statOmics/zinbwaveZinger/tree/master/zinbwaveSimulations/tenX_sims_fc2) folders. In the respective files, the FDP-TPR plots are saved and the final Figures can be recreated with the [fdrTprPlots.R](https://github.com/statOmics/zinbwaveZinger/blob/master/zinbwaveSimulations/fdrTprPlots.R) file.
We have investigated the effect of the penalty parameter on the simulated Islam and 10x datasets, which can respectively be found in the [islam_sims_fc2_epsilon](https://github.com/statOmics/zinbwaveZinger/tree/master/zinbwaveSimulations/islam_sims_fc2_epsilon) and [tenX_sims_fc2_epsilon](https://github.com/statOmics/zinbwaveZinger/tree/master/zinbwaveSimulations/tenX_sims_fc2_epsilon) folders.

### Mock comparisons

To generate the plots related to the false positive rate control, run [FPR_mocks_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/tenx/FPR_mocks_tenx.Rmd) for the 10X genomics dataset and [FPR_mocks.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/usoskin/FPR_mocks.Rmd) for the Usoskin dataset. Additionally, to generate the plots related to PCER when the penalization parameter of ZINB-WaVE is varied, run [FPR_mocks_eps_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/epsilon_tenx/FPR_mocks_eps_tenx.Rmd) for the 10X genomics dataset and [FPR_mocks_eps_usoskin.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/fpr/epsilon_usoskin/FPR_mocks_eps_usoskin.Rmd) for the Usoskin dataset.

### Real data

To generate the plots related to the analysis of the 10x Genomics PBMC dataset, first run [createDataObject.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/createdata/createDataObject.Rmd) to create the data files. Then, to generate the data when the clustering is done using PCA, run both [de_seurat.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringPCA/de_seurat.Rmd) and [de_othermethods.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringPCA/de_othermethods.Rmd). There are two files to run instead of one unique file because packages `Seurat` and `zinbwave` load both many packages and R complains that there are too many packages loaded. To generate the plots when the clustering is done using ZINB-WaVE, run [dimredZinbwave.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringW/dimredZinbwave.Rmd), [clusterW.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringW/clusterW.Rmd), and then [de.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/clusteringW/de.Rmd). Finally, to generate the plots, run [plotPaper.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/plotPaper.Rmd).  

To generate the results and plots for the differential expression analysis between the cell types identified in the Usoskin dataset, run the [deAnalysis.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/usoskin/deAnalysis.Rmd) file.

### Time benchmarking

For each of the real datasets Islam, Usoskin, 10X genomics, respectively run [benchmark_islam.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/islam/benchmark_islam.Rmd), [benchmark_usoskin.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/usoskin/benchmark_usoskin.Rmd), [benchmark_tenx.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/tenx/benchmark_tenx.Rmd). Finally, run [benchmark_all.Rmd](https://github.com/statOmics/zinbwaveZinger/blob/master/timebenchmark/benchmark_all.Rmd).



