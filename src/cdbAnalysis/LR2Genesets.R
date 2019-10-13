library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(clusterProfiler)
library(feather)
library(ggraph)
library(data.table)
library(tidygraph)
library(gridExtra)
library(SingleCellExperiment)
#source scripts
source("./src/cdbAnalysis/utils.R")
source("./latentFactors_scripts/classes.R")

############################
###### GET INPUT DATA ######
############################
#get cdb results adn rest
gmt.df <- read_tsv("./res/notAllGeneSets.tsv") %>% 
  dplyr::mutate(ENTREZID=as.character(ENTREZID))
scset <- readRDS("./res/RA/deviance_scset_markers_large2.rds")
signif_mean <- read_tsv("./RA_pipeline/output/cdb/significant_means.txt")
pval_mean <- read_tsv("./RA_pipeline/output/cdb/pvalues_means.txt")
means <- read_tsv("./RA_pipeline/output/cdb/means.txt")
deconvoluted <- read_tsv("./RA_pipeline/output/cdb/deconvoluted.txt")

#get genesets and filter by adjusted pvalues
setAllFiltered <- 
  read_feather("./RA_pipeline/output/allSignifCollectionsFilteredNotMax.feather") %>%
  dplyr::filter(method == "hpf") %>%
  dplyr::mutate(p.adj = p.adjust(pval, method = "bonferroni")) %>%
  dplyr::filter(p.adj < 0.01)

#process cdb results
gather_mean <- signif_mean %>%
  dplyr::select(-one_of(colnames(signif_mean)[2:10])) %>%
  gather(key="cluster", value="value", -interacting_pair) %>%
  na.omit()

gather_result <- pval_mean %>%
  dplyr::select(-one_of(colnames(pval_mean)[3:9])) %>%
  gather(key="cluster", value="value", -interacting_pair, -id_cp_interaction) %>%
  dplyr::mutate(mean=map_dbl(value, function(x) as.numeric(strsplit(x, " \\| ")[[1]][[1]])),
                pval=map_dbl(value, function(x) as.numeric(strsplit(x, " \\| ")[[1]][[2]]))) %>%
  dplyr::select(-one_of("value")) %>%
  inner_join(.,gather_mean, by=c("interacting_pair", "cluster")) %>%
  dplyr::mutate(pval=ifelse(pval==0,0.0001, pval)) %>%
  dplyr::mutate(log10_pval=-log10(pval),
                log2_mean=log2(mean))
gather_result$unique_id <- paste0("i.",seq(1,nrow(gather_result))) 

# cluster factors by similarity, as done for OSMR
factors.df <- getAllAllFactors("./RA_pipeline/output/evaluate_output",
                               "hpf")
rownames(factors.df) <- colnames(scset)

selectedSets <- setAllFiltered %>%
  dplyr::select(index, variable, pval, set) %>%
  dplyr::mutate(factor = gsub("Loading","Factor",index),
                pval = -log10(pval)) %>%
  dplyr::rename(loading = "index") 

#takes long, parallelize
allClusters <- clusterSets(selectedSets, factors.df, 10)

allClusters <- allClusters %>%
  as_tibble() %>%
  group_by(ont, factor) %>%
  dplyr::filter(meanSil==max(meanSil) | is.na(meanSil)) %>%
  ungroup() %>%
  dplyr::mutate(clust.memb=ifelse(meanSil<0.5 | is.na(meanSil), 1, clust.memb)) %>%
  dplyr::mutate(variable = ont) %>%
  dplyr::mutate(ont = paste0(ont,"_C",clust.memb))

#match genes with clustered sets
matchedGenesFactors <- SingleCellExperiment::rowData(scset) %>%
  as_tibble() %>%
  dplyr::select("ENSEMBL","ENTREZID","symbol") %>%
  dplyr::filter(.,ENSEMBL %in% deconvoluted$ensembl) %>%
  inner_join(.,gmt.df, by="ENTREZID") %>%
  dplyr::rename(variable = "ont") %>%
  inner_join(.,allClusters, by="variable") %>%
  dplyr::select(ENSEMBL,ENTREZID,variable, ont, symbol) %>%
  unique()

##############################
###### START ANALYSIS ########
##############################
#filter log count matrix to include only molecules present in the cdb results
logcount.mat <- scset@assays$data$counts
rownames(logcount.mat) <- SingleCellExperiment::rowData(scset)$ENSEMBL
logcount.mat <- logcount.mat[unique(matchedGenesFactors$ENSEMBL),]


allClustersFiltered <- allClusters %>%
  dplyr::filter(ont %in% matchedGenesFactors$ont)

factors.dt <- factors.df[,unique(allClustersFiltered$factor)] %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() 

factors.dt$factor <- rownames(factors.dt)  

factors.dt.jn <- factors.dt %>%
  inner_join(allClustersFiltered[,c("factor","ont")], by = "factor") %>%
  data.table::as.data.table()

factors.dt.jn[,factor:=NULL]
factors.dt.median <- factors.dt.jn[,lapply(.SD, median,na.rm=TRUE), by=ont]

ont.vec <- factors.dt.median$ont 
factors.dt.median[,ont:=NULL]
factors.dt.median <- data.table::transpose(factors.dt.median)
colnames(factors.dt.median) <- ont.vec

#scale factors and logcount mat for correlations
scaled.factors.dt <- scale(factors.dt.median, center = TRUE, scale =TRUE)
scaled.logcount.mat <- scale(t(logcount.mat), center = TRUE, scale =TRUE)

#get variables for correlations
cell_clusters <- scset@colData$graphCluster

cells <- scset@colData$type == "Monocyte"
corr_df <- cor(scaled.factors.dt[cells ,], scaled.logcount.mat[cells ,], method = "pearson")

MERTKFactors <- dplyr::filter(matchedGenesFactors, ENSEMBL == "ENSG00000153208")
MERTKCorr <- corr_df[MERTKFactors$ont,"ENSG00000153208"][corr_df[MERTKFactors$ont,"ENSG00000153208"]>0.3]

corrmat <- cor(scaled.factors.dt[cells,names(MERTKCorr)], method = "pearson")

pdf("./RA_pipeline/output/figures/corrplot_MERTKSignaling.pdf",width=6, height=6)
corrplot::corrplot(corrmat, 
                   tl.cex=0.5,
                   tl.col = "black",
                   method = "square",
                   is.corr = FALSE,
                   #tl.pos = "d",
                   cl.pos = "b",
                   addrect = 2,
                   order = "hclust",
                   hclust.method="ward.D2")
dev.off()

#function to assign molecules with factors that map to respective gene sets
#visualize correlation plot first
clustObj <- hclust(as.dist(1-corrmat), method="ward.D2")

cluster1 <-  colnames(corrmat)[clustObj$order[1:6]]
cluster2 <-  colnames(corrmat)[clustObj$order[7:15]]

scset@colData$factor <- rowMedians(as.matrix(factors.dt.median[,cluster2, with = FALSE]))
plotUMAP_c(scset, "factor")

quant1 <- quantile(rowMedians(as.matrix(factors.dt.median[,cluster1, with = FALSE])), 0.9)
quant2 <- quantile(rowMedians(as.matrix(factors.dt.median[,cluster2, with = FALSE])), 0.9)

cells1 <- colnames(scset)[rowMedians(as.matrix(factors.dt.median[,cluster1, with = FALSE])) > quant1 & scset@colData$type == "Monocyte"]
cells2 <- colnames(scset)[rowMedians(as.matrix(factors.dt.median[,cluster2, with = FALSE])) > quant2 & scset@colData$type == "Monocyte"]
cells3 <- colnames(scset)[which(
  scset@colData$type == "Monocyte" & !(colnames(scset) %in% union(cells1, cells2))
  )]

listInput <- list(cluster1 = cells1, cluster2 = cells2, no_cluster = cells3)

#upset plot to assess overlap of cells between pathways
pdf("./RA_pipeline/output/figures/upsetplot_MERTKSignaling.pdf",width=6, height=4)
UpSetR::upset(UpSetR::fromList(listInput), order.by = "freq")
dev.off()

#assign cells to specific or shared activities
cells1_spec <- cells1[which(!(cells1 %in% cells2))]
cells2_spec <- cells2[which(!(cells2 %in% cells1))]
cells1_2 <- intersect(cells1, cells2)

#prepare for heatmap
scset_full <- readRDS("./res/RA/full_scset_large.rds")

genes <- c(
           "MERTK","TNFA", "ADMA17","CCL2","TIMD4",
           "CD36","KLF2","KLF4",
           "MARCO","OLR1","P2RY2",
           "USP4","PTPN6"
           #"PDLIM2","SOCS3","HES1",
)

idx <- which(genes %in% rownames(scset_full))

cells <- colnames(scset[,scset@colData$type=="Monocyte"])
matToPlot <- SingleCellExperiment::counts(scset_full[genes[idx],cells])

metaMat <- data.frame(
  graphCluster = case_when(
    cells %in% cells1_spec ~ "cell motility",
    cells %in% cells1_2 ~ "cell motility + endocytosis",
    cells %in% cells2_spec ~ "endocytosis",
    cells %in% cells3 ~ "none"
  )
  ,
  cell_name = scset[,cells]@colData$cell_name
)

mycolors <- scico::scico(n=4, palette = "batlow",begin=0.1,end=0.9, direction = 1)
names(mycolors) <- unique(metaMat$graphCluster)
rownames(metaMat) <- colnames(scset[,cells])

#plot heatmap
source("./src/cdbAnalysis/utils.R")
pdf("./RA_pipeline/output/figures/Efferocytosis_markers_counts.pdf", width=5, height=3.2, onefile=FALSE)
plotHeatmapMarkers(matToPlot, metaMat, "Efferocytosis markers", metaMat$graphCluster, mycolors, 1)
dev.off()

cells <- colnames(scset[,scset@colData$type=="Monocyte"])
matToPlot <- SingleCellExperiment::counts(scset_full[genes[idx],cells])

metaMat <- data.frame(
  graphCluster = scset[,cells]@colData$graphCluster,
  cell_name = scset[,cells]@colData$cell_name
)
rownames(metaMat) <- colnames(scset[,cells])
mycolors <- readRDS("./res/RA/myColorsCluster2.rds")
#plot heatmap for clusters
pdf("./RA_pipeline/output/figures/Efferocytosis_markers_clusters_counts.pdf", width=5, height=3.2, onefile=FALSE)
plotHeatmapMarkers(matToPlot, metaMat, "Efferocytosis markers", metaMat$graphCluster, mycolors, 1)
dev.off()

#compositional analysis for clusters
cells <- colnames(scset[,scset@colData$type=="Monocyte"])

metaMat <- data.frame(
  graphCluster = case_when(
    cells %in% cells1_spec ~ "cell motility",
    cells %in% cells1_2 ~ "cell motility + endocytosis",
    cells %in% cells2_spec ~ "endocytosis",
    cells %in% cells3 ~ "none"
  )
  ,
  cell_name = scset[,cells]@colData$cell_name
)

cell.comp <- scset@colData[,c("graphCluster","type")] %>%
  as_tibble() %>%
  #dplyr::mutate(tmp = 0) %>%
  dplyr::count(type,graphCluster) %>%
  dplyr::rename(cluster = "graphCluster",
                tot = "n")

metaMat %>%
  as_tibble() %>%
  dplyr::rename(ontology = "graphCluster") %>%
  dplyr::inner_join(.,as.data.frame(scset@colData[,c("cell_name","graphCluster")])) %>%
  dplyr::count(ontology, graphCluster) %>%
  dplyr::rename(cluster = "graphCluster") %>%
  inner_join(., cell.comp, by="cluster") %>%
  dplyr::mutate(pct = n/tot) %>%
  ggplot(.) +
  geom_bar(aes(x=ontology, y=n, fill = cluster), 
           position = "fill", stat = "identity")+
  scale_fill_manual(values = mycolors)+
  labs(x=NULL,y = "pct")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  guides(fill=guide_legend(ncol=2)) +
  ggsave("./RA_pipeline/output/figures/MERTKSignaling_clusters_compositional.pdf", device="pdf", width=6, height=4)