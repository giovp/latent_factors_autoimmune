OpenMPController::omp_set_num_threads(snakemake@threads)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(feather)
library(data.table)
reticulate::use_python(snakemake@params$python_vers, required = TRUE)
library(reticulate)

source(paste0(snakemake@scriptdir,"/utils.R"))
source(paste0(snakemake@scriptdir,"/classes.R"))

logfile = gsub(".log",".rds",snakemake@log$evaluate_log)
saveRDS(snakemake, logfile)
#snakemake <- readRDS("./RAsubset_pipeline/output/LOGS_AND_BENCHMARKS/hpf_evaluate_snakemake_16.rds")

kCluster = as.integer(snakemake@wildcards$k)
cores = as.integer(snakemake@threads)
output.rds.obj = snakemake@output$rds_output
output.loadings = snakemake@output$loadings_output
output.factors = snakemake@output$factors_output
output.loadings.best = snakemake@output$loadings_output_best 
output.factors.best = snakemake@output$factors_output_best 
threshold.svcca = as.numeric(snakemake@params$svcca_threshold)

method = case_when(
  grepl("cogaps",snakemake@wildcards$method) ~ "cogaps",
  grepl("lda",snakemake@wildcards$method) ~ "lda",
  grepl("hpf",snakemake@wildcards$method) ~ "hpf",
  grepl("scvi",snakemake@wildcards$method) ~ "scvi"
)

#collect results according to the method used
if (method == "cogaps"){
  
  require(CoGAPS)
  
  file_list <- snakemake@input$file_list_cogaps
  file_list <- file_list[grepl(paste0("_",snakemake@wildcards$k,".rds"),file_list)]
  
  latent.obj <- collectCOGAPSResult(file_list)
  
} else if (method == "lda") {
  
  require(CountClust)
  
  file_list <- snakemake@input$file_list_lda
  file_list <- file_list[grepl(paste0("_",snakemake@wildcards$k,".rds"),file_list)]
  
  latent.obj <- collectLDAResult(file_list)
  
} else if (method == "scvi") {
  
  factors_list <- snakemake@input$factors_list_scvi
  loadings_list <- snakemake@input$loadings_list_scvi
  metrics_list <- snakemake@input$metrics_list_scvi
  
  factors_list <- factors_list[grepl(paste0("_",snakemake@wildcards$k,".feather"),factors_list)]
  loadings_list <- loadings_list[grepl(paste0("_",snakemake@wildcards$k,".feather"),loadings_list)]
  metrics_list <- metrics_list[grepl(paste0("_",snakemake@wildcards$k,".feather"),metrics_list)]
  
  latent.obj <- collectSCVIResult(factors_list, loadings_list, metrics_list)
  
} else if (method == "hpf") {
  
  factors_list <- snakemake@input$factors_list_hpf
  loadings_list <- snakemake@input$loadings_list_hpf
  metrics_list <- snakemake@input$metrics_list_hpf
  
  factors_list <- factors_list[grepl(paste0("_",snakemake@wildcards$k,".feather"),factors_list)]
  loadings_list <- loadings_list[grepl(paste0("_",snakemake@wildcards$k,".feather"),loadings_list)]
  metrics_list <- metrics_list[grepl(paste0("_",snakemake@wildcards$k,".feather"),metrics_list)]
  
  latent.obj <- collectHPFResult(factors_list, loadings_list, metrics_list)
  
}

#compute metrics
latent.obj <- computeAmariDistance(latent.obj)

cca <- reticulate::import_from_path("cca_core", path = snakemake@scriptdir)
latent.obj <- computeSVCCA(latent.obj, cca, threshold.svcca)

#compute metrics for kmedoids clustering
seed=100
#compute kmedoid clustering
kmed <- ClusterR::Cluster_Medoids(
  as.matrix(data.table::transpose(latent.obj@loadings.df)), 
  clusters=kCluster, 
  distance_metric = "euclidean",
  threads = cores, 
  swap_phase = TRUE, 
  verbose = FALSE, 
  seed = seed)

latent.obj@clust.memb <- kmed$clusters
latent.obj@kCluster <-  kCluster

#compute entropy
entropy <- estimateEntropyMixing(latent.obj, 
                                 latent.obj@clust.memb)

#compute silhouette
avg_sil<-sapply(seq(1,length(unique(kmed$silhouette_matrix$clusters))),
                function(i) 
                  median(
                    kmed$silhouette_matrix$silhouette_widths[kmed$silhouette_matrix$clusters==i]))
names(avg_sil)<-unique(kmed$silhouette_matrix$clusters)

latent.obj@sil <- avg_sil
latent.obj@entropy <- entropy
latent.obj@diss <- as.data.table(kmed$dissimilarity_matrix)

#save files
saveRDS(latent.obj, output.rds.obj)
write_feather(latent.obj@loadings.df, output.loadings)
write_feather(latent.obj@factors.df, output.factors)

#save only best iteration in the same file. This is useful for the gene set enrichment

if(snakemake@wildcards$method=="lda"){
  seed_min_iter = latent.obj@samp[which.max(latent.obj@loss)]
}else{
  seed_min_iter = latent.obj@samp[which.min(latent.obj@loss)]
}

min_iter_loadings <- grepl(paste0("_",seed_min_iter,"_"),colnames(latent.obj@loadings.df))
min_iter_factors <- grepl(paste0("_",seed_min_iter,"_"),colnames(latent.obj@factors.df))

write_feather(latent.obj@loadings.df[,min_iter_loadings, with=FALSE], output.loadings.best)
write_feather(latent.obj@factors.df[,min_iter_factors, with=FALSE], output.factors.best)

