library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(clusterProfiler)
library(feather)
library(ggraph)
library(tidygraph)
library(gridExtra)
#source scripts
source("./src/cdbAnalysis/utils.R")
source("./latentFactors_scripts/classes.R")

scset <- readRDS("./res/RA/deviance_scset_markers_large2.rds")

#get genesets and filter by adjusted pvalues
setAllFiltered <- 
  read_feather("./RA_pipeline/output/allSignifCollectionsFiltered.feather") %>%
  dplyr::filter(method == "hpf") %>%
  dplyr::mutate(p.adj = p.adjust(pval, method = "bonferroni")) %>%
  dplyr::filter(p.adj < 0.01)

#get all factors from hpf
factors.df <- getAllAllFactors("./RA_pipeline/output/evaluate_output",
                            "hpf")
rownames(factors.df) <- colnames(scset)

selectedSets <- setAllFiltered %>%
  dplyr::select(index, variable, pval, set) %>%
  dplyr::mutate(factor = gsub("Loading","Factor",index),
                pval = -log10(pval)) %>%
  dplyr::rename(loading = "index") 

allClusters <- clusterSets(selectedSets, factors.df, 10)

allClusters <- allClusters %>%
  as_tibble() %>%
  group_by(ont, factor) %>%
  dplyr::filter(meanSil==max(meanSil) | is.na(meanSil)) %>%
  ungroup() %>%
  dplyr::mutate(clust.memb=ifelse(meanSil<0.5 | is.na(meanSil), 1, clust.memb)) %>%
  dplyr::mutate(variable = ont) %>%
  dplyr::mutate(ont = paste0(ont,"_C",clust.memb))

selectedFactors <- setAllFiltered %>%
  dplyr::filter(variable=="Immune_response_Oncostatin_M_signaling_via_JAK-Stat" | 
                  variable=="Immune_response_Oncostatin_M_signaling_via_MAPK"
                ) %>%
  dplyr::select(index, variable, pval) %>%
  dplyr::mutate(factor = gsub("Loading","Factor",index),
                pval = -log10(pval)) %>%
  dplyr::rename(loading = "index") 

factors.df.filtered <- factors.df[,selectedFactors$factor] 

corrmat <- cor(factors.df.filtered, method = "pearson")
rownames(corrmat) <-  gsub("Immune_response_Oncostatin_M_signaling",
                          "OSMR_signaling",selectedFactors$variable)
colnames(corrmat) <- gsub("Loading.*_","",selectedFactors$loading)

pdf("./RA_pipeline/output/figures/corrplot_OSMRSignaling.pdf",width=6, height=4)
corrplot::corrplot(corrmat,
                   #type = "upper",
                   tl.col = "black",
                   method = "square",
                   is.corr = FALSE,
                   #tl.pos = "d",
                   cl.pos = "b",
                   addrect = 7,
                   order = "hclust",
                   hclust.method="ward.D2")
dev.off()

colnames(corrmat) <- selectedFactors$factor
clustObj <- hclust(as.dist(1-corrmat), method="ward.D2")
cluster1 <-  colnames(corrmat)[clustObj$order[1:3]]
cluster2 <-  colnames(corrmat)[clustObj$order[6:11]]
cluster3 <-  colnames(corrmat)[clustObj$order[12:14]]
cluster4 <-  colnames(corrmat)[clustObj$order[16:17]]

#create factor matrix for scater plotexpression
factors.mat <-  matrix(nrow=nrow(scset),ncol = ncol(scset))
colnames(factors.mat) <- colnames(scset)

factors.mat[1,] <- rowMeans(
  as.matrix(factors.df.filtered[,cluster1])) 
factors.mat[2,] <- rowMeans(
  as.matrix(factors.df.filtered[,cluster2])) 
factors.mat[3,] <- rowMeans(
  as.matrix(factors.df.filtered[,cluster3])) 
factors.mat[4,] <- rowMeans(
  as.matrix(factors.df.filtered[,cluster4])) 

scset@assays$data$weights <- factors.mat

#get colors
mycolors <- readRDS("./res/RA/myColorsCluster2.rds")

#beeswarm plot                                
cells_cluster1 <- scset@colData$type=="Fibroblast"
plot_cluster1 <- scater::plotExpression(scset[,cells_cluster1], 1,
                       x = "graphCluster", 
                       exprs_values = "weights", 
                       colour_by = "graphCluster",
                       one_facet= FALSE)+
  labs(y= "Weight", x=NULL)+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

cells_cluster2 <- scset@colData$type=="Fibroblast"
plot_cluster2 <- scater::plotExpression(scset[,cells_cluster2], 2,
                                        x = "graphCluster", 
                                        exprs_values = "weights", 
                                        colour_by = "graphCluster",
                                        one_facet= FALSE)+
  labs(y= "Weight", x=NULL)+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

cells_cluster3 <- scset@colData$type=="Fibroblast"
plot_cluster3 <- scater::plotExpression(scset[,cells_cluster3], 3,
                                        x = "graphCluster", 
                                        exprs_values = "weights", 
                                        colour_by = "graphCluster",
                                        one_facet= FALSE)+
  labs(y= "Weight", x=NULL)+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

cells_cluster4 <- scset@colData$type=="Fibroblast"
plot_cluster4 <- scater::plotExpression(scset[,cells_cluster4], 4,
                                        x = "graphCluster", 
                                        exprs_values = "weights", 
                                        colour_by = "graphCluster",
                                        one_facet= FALSE)+
  labs(y= "Weight", x=NULL)+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave("./RA_pipeline/output/figures/OSMRSignaling_clusters_beeswar.pdf",grid.arrange(plot_cluster1,
             plot_cluster2,
             plot_cluster3,
             plot_cluster4,
             nrow=2), width=6, height=4, device="pdf")

#compositional plot
factors.tib <- as_tibble(factors.mat[1:4,] %>% t())
colnames(factors.tib) <- c("cluster1","cluster2","cluster3","cluster4") 

cell.comp <- scset@colData[,c("graphCluster","type")] %>%
  as_tibble() %>%
  #dplyr::mutate(tmp = 0) %>%
  dplyr::count(type,graphCluster) %>%
  dplyr::rename(cluster = "graphCluster",
                tot = "n")

factors.tib %>%
  dplyr::mutate(cluster = scset@colData$graphCluster) %>%
  gather(key="clusterPath",
         value="weight",
         -cluster) %>%
  dplyr::filter(weight>0.5) %>%
  dplyr::count(cluster, clusterPath) %>% 
  inner_join(., cell.comp, by="cluster") %>%
  dplyr::mutate(pct = n/tot) %>%
  #dplyr::filter(type == "Fibroblast") %>%
ggplot(.) +
  geom_bar(aes(x=clusterPath, y=n, fill = cluster), 
           position = "fill", stat = "identity")+
  scale_fill_manual(values = mycolors)+
  labs(x="OncostatinM signaling", fill ="cluster", y = "pct")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  guides(fill=guide_legend(ncol=2)) +
  ggsave("./RA_pipeline/output/figures/OSMRSignaling_clusters_compositional.pdf", device="pdf", width=6, height=4)

#heatmap of chemokine signature from IBD, read the full expression object
scset_full <- readRDS("./res/RA/full_scset_large.rds")

genes <- c("CXCL9", "CCL2", "ICAM1",
           "IL6", "ACTA2", "CDH11","IL11","CXCL11","CXCL10", "IL1B", "CXCL1","OSMR"#,
           #"CSF2","FAP"
           #"IL6","PDPN","CXCL2"#,"CXCL3"#,"CXCL5","IL1B","CXCL11","CXCL6",
)

idx <- which(genes %in% rownames(scset_full))

quant1 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster1])), 0.9)
quant2 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster2])), 0.9)
quant3 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster3])), 0.9)
quant4 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster4])), 0.9)

cells1 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster1])) > quant1 & scset@colData$type == "Fibroblast"]
cells2 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster2])) > quant2 & scset@colData$type == "Fibroblast"]
cells3 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster3])) > quant3 & scset@colData$type == "Fibroblast"]
cells4 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster4])) > quant4 & scset@colData$type == "Fibroblast"]

listInput <- list(cells1 = cells1, cells2 = cells2, cells3 = cells3, cells4 = cells4)

pdf("./RA_pipeline/output/figures/upsetplot_OSMRSignaling.pdf",width=6, height=4)
UpSetR::upset(UpSetR::fromList(listInput), order.by = "freq")
dev.off()

cells1_spec <- cells1[which(!(cells1 %in% cells2 | cells1 %in% cells3 | cells1 %in% cells4))]
cells2_spec <- cells2[which(!(cells2 %in% cells1 | cells2 %in% cells3 | cells2 %in% cells4))]
cells3_spec <- cells3[which(!(cells3 %in% cells1 | cells3 %in% cells2 | cells3 %in% cells4))]
cells4_spec <- cells4[which(!(cells4 %in% cells1 | cells4 %in% cells2 | cells4 %in% cells3))]

cells <- colnames(scset[,scset@colData$type=="Fibroblast"])

metaMat <- data.frame(
  graphCluster = case_when(
    cells %in% cells1_spec ~ "OSMR signaling via MAPK",
    cells %in% cells2_spec ~ "OSMR signaling mixed",
    cells %in% cells3_spec ~ "OSMR signaling via JAK-Stat 1",
    cells %in% cells4_spec ~ "OSMR signaling via JAK-Stat 2"
  )
  ,
  cell_name = scset[,cells]@colData$cell_name
) %>%
  na.omit()

matToPlot <- SingleCellExperiment::counts(scset_full[genes[idx],as.character(metaMat$cell_name)])

mycolors <- scico::scico(n=4, palette = "batlow",begin=0.1,end=0.9, direction = 1)
names(mycolors) <- unique(metaMat$graphCluster)
rownames(metaMat) <- colnames(scset[,as.character(metaMat$cell_name)])

source("./src/cdbAnalysis/utils.R")

pdf("./RA_pipeline/output/figures/OSMR_high_markers.pdf", width=5, height=3.2, onefile=FALSE)
plotHeatmapMarkers(matToPlot, metaMat, "OSMR markers", metaMat$graphCluster, mycolors, 1)
dev.off()

cells <- colnames(scset[,scset@colData$type=="Fibroblast"])
matToPlot <- SingleCellExperiment::counts(scset_full[genes[idx],cells])

metaMat <- data.frame(
  graphCluster = scset[,cells]@colData$graphCluster,
  cell_name = scset[,cells]@colData$cell_name
)
rownames(metaMat) <- colnames(scset[,cells])
mycolors <- readRDS("./res/RA/myColorsCluster2.rds")

pdf("./RA_pipeline/output/figures/OSMR_high_markers_clusters.pdf", width=5, height=3.2, onefile=FALSE)
plotHeatmapMarkers(matToPlot, metaMat, "OSMR markers", metaMat$graphCluster, mycolors, 1)
dev.off()

#compositional plot for activity clusters
idx <- which(genes %in% rownames(scset_full))

quant1 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster1])), 0.9)
quant2 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster2])), 0.9)
quant3 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster3])), 0.9)
quant4 <- quantile(rowMeans(as.matrix(factors.df.filtered[,cluster4])), 0.9)

cells1 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster1])) > quant1] #& scset@colData$type == "Fibroblast"]
cells2 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster2])) > quant2] #& scset@colData$type == "Fibroblast"]
cells3 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster3])) > quant3] #& scset@colData$type == "Fibroblast"]
cells4 <- colnames(scset)[rowMeans(as.matrix(factors.df.filtered[,cluster4])) > quant4] #& scset@colData$type == "Fibroblast"]

listInput <- list(cells1 = cells1, cells2 = cells2, cells3 = cells3, cells4 = cells4)


cells1_spec <- cells1[which(!(cells1 %in% cells2 | cells1 %in% cells3 | cells1 %in% cells4))]
cells2_spec <- cells2[which(!(cells2 %in% cells1 | cells2 %in% cells3 | cells2 %in% cells4))]
cells3_spec <- cells3[which(!(cells3 %in% cells1 | cells3 %in% cells2 | cells3 %in% cells4))]
cells4_spec <- cells4[which(!(cells4 %in% cells1 | cells4 %in% cells2 | cells4 %in% cells3))]

cells <- colnames(scset[,])

metaMat <- data.frame(
  graphCluster = case_when(
    cells %in% cells1_spec ~ "OSMR signaling via MAPK",
    cells %in% cells2_spec ~ "OSMR signaling mixed",
    cells %in% cells3_spec ~ "OSMR signaling via JAK-Stat 1",
    cells %in% cells4_spec ~ "OSMR signaling via JAK-Stat 2"
  )
  ,
  cell_name = scset[,cells]@colData$cell_name
) %>%
  na.omit()

mycolors <- readRDS("./res/RA/myColorsCluster2.rds")
metaMat %>%
  dplyr::rename(activity = "graphCluster") %>% 
  inner_join(.,as_tibble(scset@colData[,c("graphCluster","cell_name")])) %>%
  as_tibble() %>%
  #dplyr::select(-one_of(c("cell_name")))
  #dplyr::mutate(tmp = 0) %>%
  dplyr::count(activity,graphCluster) %>%
  ggplot(.) +
  geom_bar(aes(x=activity, y=n, fill = graphCluster), 
           position = "fill", stat = "identity")+
  scale_fill_manual(values = mycolors)+
  labs(x="", fill ="cluster", y = "pct")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  guides(fill=guide_legend(ncol=2))+
  ggsave("./RA_pipeline/output/figures/OSMRSignaling_activity_compositional.pdf", device="pdf", width=6, height=4)

# #compositional plot specific for activity pathways
# cell.comp <- scset@colData[,c("graphCluster","type")] %>%
#   as_tibble() %>%
#   #dplyr::mutate(tmp = 0) %>%
#   dplyr::count(type,graphCluster) %>%
#   dplyr::rename(cluster = "graphCluster",
#                 tot = "n")
# 
# mycolors <- readRDS("./res/RA/myColorsCluster2.rds")
# metaMat %>%
#   as_tibble() %>%
#   dplyr::rename(clusterPath = "graphCluster") %>%
#   inner_join(.,as.data.frame(scset@colData[,c("graphCluster","type","cell_name")]), by="cell_name") %>%
#   dplyr::rename(cluster = "graphCluster") %>%
#   dplyr::select(cluster, clusterPath) %>%
#   dplyr::count(cluster, clusterPath) %>% 
#   inner_join(., cell.comp, by="cluster") %>%
#   dplyr::mutate(pct = n/tot) %>%
#   #dplyr::filter(type == "Fibroblast") %>%
#   ggplot(.) +
#   geom_bar(aes(x=clusterPath, y=n, fill = cluster), 
#            position = "fill", stat = "identity")+
#   scale_fill_manual(values = mycolors)+
#   labs(fill ="cluster", y = "pct")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   guides(fill=guide_legend(ncol=2)) +
#   ggsave("./RA_pipeline/output/figures/OSMRSignaling_clusters_compositional_specific.pdf", device="pdf", width=6, height=4)
# 
# 
# # 
# genes <- c("CXCL9", "CCL2", "ICAM1",
#            "IL6", "ACTA2", "CDH11","IL11","CXCL1","CXCL10", "IL1B"#, "CXCL11"#,
# )
# idx <- which(genes %in% rownames(scset_full))
# cells <- names(factors.mat[1,factors.mat[1,]>0.5])
# cells <- cells[which(cells %in% colnames(scset[,scset@colData$type=="Fibroblast"]))]
# matToPlot <- SingleCellExperiment::logcounts(scset_full[genes[idx],cells])
# 
# metaMat <- data.frame(
#   graphCluster = scset[,cells]@colData$graphCluster,
#   cell_name = scset[,cells]@colData$cell_name
# )
# rownames(metaMat) <- colnames(scset[,cells])
# 
# pdf("./RA_pipeline/output/figures/OSMR_high_markers_cluster1.pdf", width=5, height=3.2, onefile=FALSE)
# plotHeatmapMarkers(matToPlot, metaMat, "OSMR high markers", metaMat$graphCluster, mycolors, 2)
# dev.off()

#THY1 stuff
genes <- c("THY1","PRG4","CLIC5","TSPAN15","IL34","CCL9","TNSF11",
           "MMP3","PDPN"
)
idx <- which(genes %in% rownames(scset_full))
cells <-  colnames(scset[,scset@colData$type=="Fibroblast"])
matToPlot <- SingleCellExperiment::logcounts(scset_full[genes[idx],cells])

metaMat <- data.frame(
  graphCluster = scset[,cells]@colData$graphCluster,
  cell_name = scset[,cells]@colData$cell_name
)
rownames(metaMat) <- colnames(scset[,cells])

pdf("./RA_pipeline/output/figures/THY1_markers.pdf.pdf", width=4, height=3.2, onefile=FALSE)
plotHeatmapMarkers(matToPlot, metaMat, "Fibroblast markers", metaMat$graphCluster, mycolors, 2)
dev.off()

source("./src/Functions.R")
plotUMAP_gene(scset[], "OSMR")+
  ggsave("./RA_pipeline/output/figures/OSMR-expression.pdf", width=6, height=5, device="pdf")

plotUMAP_gene(scset[], "THY1")+
  ggsave("./RA_pipeline/output/figures/THY1-expression.pdf", width=6, height=5, device="pdf")

plotUMAP_gene(scset[], "PRG4")+
  ggsave("./RA_pipeline/output/figures/PRG4-expression.pdf", width=6, height=5, device="pdf")

#plot factors UMAP
scset@colData$factor <- scset@assays$data$weights[1,]
plotUMAP_c(scset, "factor")+
  ggsave("./RA_pipeline/output/figures/OncoMJAKStat-factorUMAP.pdf", 
         width=6, height=5, device="pdf")

#edgeR analysis
# source("./src/Functions.R")
# library(edgeR)
# scset_fibro<-readRDS("./res/RA/scset_fibro_large2.rds")
# 
# combos_fibro <- with(colData(scset_fibro), paste(plate, SUBJECT_ACCESSION, sep="."))
# summed_fibro <- sumCountsAcrossCells(scset_fibro, combos_fibro)
# 
# y_fibro <- DGEList(summed_fibro)
# y_fibro <- y_fibro[aveLogCPM(y_fibro) > 1,]
# y_fibro <- calcNormFactors(y_fibro)
# 
# sum.terms.fibro <- strsplit(colnames(y_fibro), split="\\.")
# sum.plate.fibro <- unlist(lapply(sum.terms.fibro, "[[", i=1))
# sum.donor.fibro <- unlist(lapply(sum.terms.fibro, "[[", i=2))
# sum.disease.fibro <- unique(colData(scset_fibro)[,c("SUBJECT_ACCESSION","disease")])[match(sum.donor.fibro,unique(scset_fibro@colData$SUBJECT_ACCESSION)),"disease"] 
# 
# design.fibro <- model.matrix(~0 + sum.plate.fibro + sum.disease.fibro)# + sum.donor.fibro)# + sum.disease.fibro)
# 
# y_fibro <- estimateDisp(y_fibro, design.fibro)
# fit.fibro <- glmQLFit(y_fibro, design.fibro, robust=TRUE)
# print("Trend dispersion")
# print(summary(y_fibro$trended.dispersion))
# print("Df prior")
# print(summary(fit.fibro$df.prior))
# #res.fibro <- glmQLFTest(fit.fibro)
# res.fibro <- glmLRT(fit.fibro)
# print("Res DE test")
# summary(decideTests(res.fibro))
# 
# #OSMR expression
# scater::plotExpression(scset[,scset@colData$type== "Fibroblast"], c("OSMR"),
#                        x = "graphCluster", exprs_values = "logcounts", colour_by = "graphCluster")+
#   scale_fill_manual(values = mycolors)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   ggsave("./RA_pipeline/output/figures/OSMR-expression-beeswarm.pdf", height=4, width = 5, device="pdf")