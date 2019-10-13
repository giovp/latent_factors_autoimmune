library(tidyverse)
library(feather)
source("./latentFactors_scripts/utils.R")
source("./latentFactors_scripts/classes.R")
source("./src/assignLoadings/utils.R")
source("./src/cdbAnalysis/utils.R")
source("./src/Functions.R")

method = "hpf"
#load files
loss.df <- getBestRun("./SLE_pipeline/output/evaluate_output",
                      "SLE_pipeline",
                      c("hpf","cogaps","scvi","lda"),
                      seq(16,40,2))

path_list <- list.files(path=paste0("./SLE_pipeline/output/evaluate_output"), 
                        pattern = paste0(method,".*.rds"), 
                        full.names = TRUE)

factors_list <- sapply(path_list, function(name){
  
  latent.obj <- readRDS(name)
  
  if(grepl("lda", name)){
    loss.idx <- latent.obj@samp[which(latent.obj@loss==max(latent.obj@loss))[1]]
  }else{
    loss.idx <- latent.obj@samp[which(latent.obj@loss==min(latent.obj@loss))[1]]
  }
  
  cols <- grepl(paste0("_",loss.idx,"_"),colnames(latent.obj@factors.df))
  return(latent.obj@factors.df[,cols, with=FALSE])
})

loadings_list <- sapply(path_list, function(name){
  
  latent.obj <- readRDS(name)
  
  if(grepl("lda", name)){
    loss.idx <- latent.obj@samp[which(latent.obj@loss==max(latent.obj@loss))[1]]
  }else{
    loss.idx <- latent.obj@samp[which(latent.obj@loss==min(latent.obj@loss))[1]]
  }
  
  cols <- grepl(paste0("_",loss.idx,"_"),colnames(latent.obj@loadings.df))
  return(latent.obj@loadings.df[,cols, with=FALSE])
})

scset <- readRDS("./res/SLEcelseq/deviance_scset_clusters_large.rds")
setAllFiltered <- 
  read_feather("./SLE_pipeline/output/allSignifCollectionsFiltered.feather") %>%
  dplyr::filter(method == "hpf") 
setAllFiltered$factor <- gsub("Loading","Factor", setAllFiltered$index)

selectedSets <- setAllFiltered %>%
  dplyr::select(index, variable, pval, set, factor) %>%
  dplyr::mutate(pval = -log10(pval)) %>%
  dplyr::rename(loading = "index") 

factors.df <- getAllAllFactors("./SLE_pipeline/output/evaluate_output",
                               "hpf")
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

factors.dt <- factors.df[,unique(allClusters$factor)] %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() 

factors.dt$factor <- rownames(factors.dt)  

factors.dt.jn <- factors.dt %>%
  inner_join(allClusters[,c("factor","ont")], by = "factor") %>%
  data.table::as.data.table()

factors.dt.jn[,factor:=NULL]
factors.dt.median <- factors.dt.jn[,lapply(.SD, median,na.rm=TRUE), by=ont]

ont.vec <- factors.dt.median$ont 
factors.dt.median[,ont:=NULL]
factors.dt.median <- data.table::transpose(factors.dt.median)
colnames(factors.dt.median) <- ont.vec

factors.dt.median.scaled <- scale(factors.dt.median, center=TRUE, scale = TRUE)

predictor_list <- lapply(seq(1, length(colnames(factors.dt.median.scaled))), function(i){

  response = as.data.frame(factors.dt.median.scaled[,i, drop =FALSE])
  colnames(response) <- "y"
  response$cluster <- scset@colData$graphClusterType
  
  mod <- lm(as.formula("y ~ 0 + cluster"), data = response)
  coeff = mod$coefficients#coef(summary(mod))[,"t value"]#mod$coefficients
  names(coeff) <- gsub("cluster","",names(coeff))
  return(coeff)
})

predictors_all <- dplyr::bind_cols(predictor_list) %>%
  as.data.frame()
colnames(predictors_all)<- colnames(factors.dt.median.scaled)
rownames(predictors_all)<- names(predictor_list[[1]])

collection="GpC7"
#GSE[0-9]*_"
nloadings_response <- allClusters  %>%
  dplyr::filter(set == collection) %>%
  group_by(ont) %>%
  tally() %>%
  tibble::column_to_rownames(.,var = "ont")

response_df <- predictors_all[,] %>% t() %>% as.data.frame()
rownames(response_df) <- colnames(predictors_all)
response_df <- response_df[rownames(nloadings_response ),]
rownames(response_df) <- gsub("GSE[0-9]*_","",rownames(response_df))
rownames(nloadings_response) <- gsub("GSE[0-9]*_","",rownames(nloadings_response))

pdf("./SLE_pipeline/output/figures/C7_cluster_assigned_C1.pdf", width=10, height=14, onefile=FALSE)
response_df_filt <- response_df[grepl("_C1$",rownames(response_df)),] 
nloadings_response_filt <- nloadings_response[grepl("_C1$",rownames(nloadings_response)),, drop=FALSE]
rownames(response_df_filt) <- gsub("_C1$","",rownames(response_df_filt))
rownames(nloadings_response_filt) <- gsub("_C1$","",rownames(nloadings_response_filt))
plotHeatmapLoadingsCluster(response_df_filt,nloadings_response_filt,"C7 collection")
dev.off()

pdf("./SLE_pipeline/output/figures/C7_cluster_assigned_Crest.pdf", width=12, height=8, onefile=FALSE)
response_df_filt <- response_df[!grepl("_C1$",rownames(response_df)),] 
nloadings_response_filt <- nloadings_response[!grepl("_C1$",rownames(nloadings_response)),, drop=FALSE]
plotHeatmapLoadingsCluster(response_df_filt,nloadings_response_filt,"C7 collection")
dev.off()


#clusters <- scset@colData$graphCluster
predictor_list <- sapply(seq(1, length(factors_list)), function(i){
  factors.df <- as.data.frame(scale(factors_list[[i]], center = TRUE, scale = TRUE))
  factors.df$clusters <- scset@colData$graphCluster
  coeff_list <- sapply(seq(1,ncol(factors.df)-1), function(j){
    mod <- lm(as.formula(paste0(colnames(factors.df)[j], "~0 + clusters")), factors.df)
    coeff = mod$coefficients#coef(summary(mod))[,"t value"]#mod$coefficients
    names(coeff) <- gsub("clusters","",names(coeff))
    return(coeff)
  })
  coeff_df <- as.data.frame(coeff_list)
  colnames(coeff_df) <- colnames(factors.df)[-length(colnames(factors.df))]
  return(coeff_df)
})

predictors_all <- dplyr::bind_cols(predictor_list) %>%
  t() %>%
  as.data.frame()
colnames(predictors_all)<- rownames(predictor_list[[1]])

#%WikiPathways.*Homo sapiens
selectedSets <- setAllFiltered %>%
  dplyr::select(index, variable, pval, set) %>%
  dplyr::mutate(factor = gsub("Loading","Factor",index),
                pval = -log10(pval)) %>%
  dplyr::rename(loading = "index") 

factors.df <- getAllAllFactors("./RA_pipeline/output/evaluate_output",
                               "hpf")

allClusters <- clusterSets(selectedSets, factors.df, cores = 4)

allClusters <- allClusters %>%
  as_tibble() %>%
  group_by(ont, factor) %>%
  dplyr::filter(meanSil==max(meanSil) | is.na(meanSil)) %>%
  ungroup() %>% 
  dplyr::mutate(clust.memb=ifelse(meanSil<0.5 | is.na(meanSil), 1, clust.memb)) %>%
  dplyr::mutate(ont = paste0(ont,"_C",clust.memb))

collection="GpC7"
meta_predictors <- allClusters  %>%
  dplyr::filter(factor %in% rownames(predictors_all) & set == collection) %>% 
  #dplyr::select(c("ont","factor","abs_z_score","z_score")) %>% 
  mutate(ont=gsub("GSE[0-9]*_","",ont)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(.,var = "factor")

meta_predictors <- meta_predictors[rownames(predictors_all),]
predictors_all$ont <- meta_predictors$ont

predictors_grouped <- predictors_all %>% 
  group_by(ont) %>%
  summarise_all(list(median)) %>%
  tibble::column_to_rownames(.,var = "ont")

meta_predictors_grouped <- meta_predictors %>%
  group_by(ont) %>%
  tally() %>%
  tibble::column_to_rownames(.,var = "ont")

pdf("./RA_pipeline/output/figures/C7_cluster_assigned_c1.pdf", width=11, height=15, onefile=FALSE)
plotHeatmapLoadingsCluster(predictors_grouped[grepl("_C1$",rownames(predictors_grouped)),], 
                                     meta_predictors_grouped[grepl("_C1$",rownames(meta_predictors_grouped)),, drop=FALSE], "C7 collection")
dev.off()

pdf("./RA_pipeline/output/figures/C7_cluster_assigned_cRest.pdf", width=10, height=5, onefile=FALSE)
plotHeatmapLoadingsCluster(predictors_grouped[!grepl("_C1$",rownames(predictors_grouped)),], 
                           meta_predictors_grouped[!grepl("_C1$",rownames(meta_predictors_grouped)),, drop=FALSE], "C7 collection")
dev.off()

selectedFactors <- allClusters %>%
  dplyr::filter(ont == "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY_C1") %>%
  pull(factor) %>%
  unique()

##########################
###### plot factors ######
##########################

factors.df <- getAllFactors("./RA_pipeline/output/evaluate_output",
                            "hpf",
                            selectedFactors)

rownames(factors.df) <- colnames(scset)

scset@colData$factor <- factors.df[,1]
scset@colData$factor <- rowMedians(as.matrix(factors.df[,]))
plotUMAP_c(scset, "factor")+
  ggsave("./RA_pipeline/output/figures/Bcellreceptor-meanFactors.pdf", width=6, height=5, device="pdf")


#other stuff
plotUMAP_gene(scset, "IFNG")
loadings.df <- getAllLoadings("./RA_pipeline/output/evaluate_output",
                              "hpf",gsub('Factor','Loading', selectedFactors))

loadings_weights <- rowMeans(loadings.df[,])
names(loadings_weights) <- rownames(scset)

cells <- scset@colData$factor > 0.5

scater::plotExpression(scset[,grepl("Bcell", scset@colData$graphClusterType)],
                       names(sort(loadings_weights, decreasing = TRUE)[1:15]),
                       x = "graphClusterType",
                       exprs_values = "logcounts",
                       colour_by = "graphClusterType",
                       one_facet= FALSE)+
  scale_fill_manual(values = mycolors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave("./RA_pipeline/output/figures/TYPEI-IFN-SIGNALING-beeswarm.pdf", width=13, height=13, device="pdf")