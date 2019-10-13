library(tidyverse)
library(data.table)
library(gridExtra)
library(Matrix)
library(ggbeeswarm)
library(scico)

source("./latentFactors_scripts/classes.R")
source("./src/evaluateModels/utils.R")
######################
#### parse files #####
######################

cogaps <- parseMetricsFiles("./SLE_pipeline/output/evaluate_output/" , "cogaps.*rds")
hpf <- parseMetricsFiles("./SLE_pipeline/output/evaluate_output/" , "hpf.*rds")
scvi <- parseMetricsFiles("./SLE_pipeline/output/evaluate_output/" , "scvi.*rds")
lda <- parseMetricsFiles("./SLE_pipeline/output/evaluate_output/" , "lda.*rds")
cogaps$algo <- "cogaps"
hpf$algo <- "hpf"
scvi$algo <- "scvi"
lda$algo <- "lda"
metrics <- dplyr::bind_rows(cogaps, hpf, scvi, lda)

######################
## SILHOUETTE plots ##
######################
mycolors <- c(scico::scico(2, palette="devon", begin=0.4, end=0.6),
              scico::scico(2, palette="lajolla", begin=0.3, end=0.6))
names(mycolors)<-c("hpf","cogaps","scvi","lda")

metricsPlot <- plotMetrics(metrics, 
                      c("scCoGAPS","scHPF","LDA", "scVI"),
                      mycolors)

labs_facet = c("scCoGAPS","scHPF","LDA", "scVI")
names(labs_facet) = c("cogaps","hpf","lda","scvi")
loss <- plotLoss(metrics, labs_facet)

#toSave <- arrangeGrob(grobs = list(sil, loss), ncol = 1, nrow = 2)
ggsave("./res/figures/metrics_plot_sle.pdf", metricsPlot, device="pdf", width=12, height = 4)
ggsave("./SLE_pipeline/output/figures/lossPlot.pdf", loss, device="pdf", width=6, height = 4)

