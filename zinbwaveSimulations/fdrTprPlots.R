 setwd("~/Dropbox/phdKoen/singleCell/zinbwavezingerGitHub/zinbwaveZinger/zinbwaveSimulations/")
library(iCOBRA)
library(cowplot)


####### combine two-panel plots with same legend
## islam and Trapnell
load("./islam_sims_fc2/cobraplotIslam.rda")
islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
load("./trapnell_sims_fc2/cobraplotTrapnell.rda")
trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
           trapnellPlot + theme(legend.position="none") + xlab("FDP"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_composite.png", width=7,height=8, units="in", res=300)
p
dev.off()


## islam and Trapnell without limma-voom
load("./islam_sims_fc2/cobraplotIslamNoLimma.rda")
islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.5), yaxisrange=c(0,0.7))
load("./trapnell_sims_fc2/cobraplotTrapnellNoLimma.rda")
trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.5), yaxisrange=c(0,0.7))
#plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
           trapnellPlot + theme(legend.position="none") + xlab("FDP"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_composite_cutoff_noLimma.png", width=7,height=8, units="in", res=300)
p
dev.off()


#10X
load("./tenX_sims_fc2/cobraplot10x.rda")
tenxPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_10x.png", width=7,height=8, units="in", res=300)
tenxPlot
dev.off()

#10X no ZINB-WaVE limma
load("./tenX_sims_fc2/cobraplot10xNoLimma.rda")
tenxPlot=plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.41), yaxisrange=c(0,0.4))
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_10x_cutoff_noLimma.png", width=7,height=8, units="in", res=300)
tenxPlot
dev.off()

######### FDR-TPR plots including genewise ZINB-WaVE methods
## islam and Trapnell
load("./islam_sims_fc2/cobraplotIslamAllMethods.rda")
islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
load("./trapnell_sims_fc2/cobraPlotTrapnellAllMethods.rda")
trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
prow <- plot_grid( islamPlot + theme(legend.position="none") + xlab("FDP"),
           trapnellPlot + theme(legend.position="none") + xlab("FDP"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_composite_allMethods.png", width=9,height=8, units="in", res=300)
p
dev.off()


### 10x
load("./tenX_sims_fc2/cobraplot10x.rda")
tenxPlot=plot_fdrtprcurve(cobraplot, pointsize=2) + xlab("FDP")
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_10x_allMethods.png", width=7,height=8, units="in", res=300)
tenxPlot
dev.off()

load("./tenX_sims_fc2/cobraplot10xNoGenewise.rda")
tenxPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_10x.png", width=7,height=8, units="in", res=300)
tenxPlot
dev.off()

load("./tenX_sims_fc2/cobraplot10xNoLimma.rda")
tenxPlot=plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange=c(0,0.41), yaxisrange=c(0,0.4))
png("~/Dropbox/phdKoen/singleCell/zinbwaveZinger/plots2/scSimulation_10x_cutoff_noLimma.png", width=7,height=8, units="in", res=300)
tenxPlot
dev.off()
