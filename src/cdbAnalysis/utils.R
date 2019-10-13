getAllFactors <- function(dir, method, factors){
  
  file.list <- list.files(dir, pattern =paste0(method,".*rds") ,full.names = TRUE)
  
  factors.list <- lapply(file.list, function(file.path){
    latent.obj <- readRDS(file.path)
    selected.factors <- which(colnames(latent.obj@factors.df) %in% factors)
    
    latent.factors <- latent.obj@factors.df[,selected.factors, with=FALSE] %>% as.data.frame()
    return(latent.factors)
  })
  
  factors.df <- dplyr::bind_cols(purrr::compact(factors.list))
  return(factors.df)
}

getAllAllFactors <- function(dir, method){
  
  file.list <- list.files(dir, pattern =paste0(method,".*rds") ,full.names = TRUE)
  
  factors.list <- lapply(file.list, function(file.path){
    latent.obj <- readRDS(file.path)
    
    latent.factors <- latent.obj@factors.df[,] %>% as.data.frame()
    return(latent.factors)
  })
  
  factors.df <- dplyr::bind_cols(purrr::compact(factors.list))
  return(factors.df)
}

getAllLoadings <- function(dir, method, loadings){
  
  file.list <- list.files(dir, pattern =paste0(method,".*rds") ,full.names = TRUE)
  
  loadings.list <- lapply(file.list, function(file.path){
    latent.obj <- readRDS(file.path)
    selected.loadings <- which(colnames(latent.obj@loadings.df) %in% loadings)
    
    latent.loadings <- latent.obj@loadings.df[,selected.loadings, with=FALSE] %>% as.data.frame()
    return(latent.loadings)
  })
  
  loadings.df <- dplyr::bind_cols(purrr::compact(loadings.list))
  return(loadings.df)
}

getGraph <- function(network){
  
  #from https://stackoverflow.com/questions/51851798/using-tidygraph-to-merge-two-edges-from-the-same-two-nodes-into-one/52018050
  pasteCols = function(x, y, sep = ":"){
    stopifnot(length(x) == length(y))
    return(lapply(1:length(x), function(i){paste0(sort(c(x[i], y[i])), collapse = ":")}) %>% unlist())
  }
  
  nodes <- dplyr::bind_rows(
    network %>% dplyr::count(to,unique_id) %>% dplyr::rename(cluster = "to"), 
    network %>% dplyr::count(from,unique_id) %>% dplyr::rename(cluster = "from")
  ) %>%
    unique() %>%
    group_by(cluster) %>%
    summarise(n=sum(n)) %>%
    dplyr::rename(name = "cluster")
  
  edges <- dplyr::bind_rows(
    network,
    network %>% rename(from = "to",to = "from")
  ) %>%
    dplyr::count(from,to) %>%
    dplyr::rename(weight = "n") %>%
    mutate(col_pairs = pasteCols(from, to, sep = ":")) %>% 
    group_by(col_pairs) %>% summarise(sum_weight = sum(weight)) %>% 
    tidyr::separate(col = col_pairs, c("from", "to"), sep = ":") %>%
    dplyr::rename(n="sum_weight")
  
  graph <- tbl_graph(nodes = nodes, 
                     edges = edges, directed = FALSE)
  
  return(graph)
}

plotGraphComparison<-function(graph){
  ggraph(graph, layout = "manual", node.positions = layout[,c("x","y")]) + 
    geom_edge_fan(aes(colour=n, width=n), alpha=0.8)+
    geom_edge_loop(aes(colour=n, width=n), alpha=0.8)+
    scale_edge_width_continuous(
      range = c(0,2),
      limits=c(filt.graph %>% activate(edges) %>% pull(n) %>% min(),
               full.graph %>% activate(edges) %>% pull(n) %>% max()),
      guide = FALSE)+
    scale_edge_colour_distiller(
      palette = "Greys", 
      direction = 1,
      limits=c(filt.graph %>% activate(edges) %>% pull(n) %>% min(),
               full.graph %>% activate(edges) %>% pull(n) %>% max())) +
    scale_edge_alpha(guide = FALSE) +
    geom_node_point(
      aes(size = n, fill = as.factor(name)),shape=21) +
    scale_size(
      limits =c(filt.graph %>% activate(nodes) %>% pull(n) %>% min(),
                full.graph %>% activate(nodes) %>% pull(n) %>% max()),
      range = c(1,10)) +
    geom_node_label(
      aes(label = name), 
      repel = TRUE, 
      label.padding = 0.15)+
    scale_fill_manual(values = mycolors, guide = FALSE)+
    labs(size="degree", edge_colour = "n_interactions")+
    theme_graph(fg_text_colour = 'white', base_family = 'Helvetica')
}

getCorrTable <- function(cell_clusters, all_factors, count_mat){
  
  corr_list <- lapply(unique(cell_clusters), function(clust){
    corr <- cor(count_mat[cell_clusters==clust,], all_factors[cell_clusters==clust,], method = "spearman") %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      as_tibble() %>%
      drop_na()%>%
      gather(key="ont", value = "value", -gene) %>%
      dplyr::filter(value>0.3) %>%
      dplyr::mutate(cluster = clust)
    return(corr)
  })
  return(dplyr::bind_rows(corr_list))
}

#used in next function
computeK <- function(corrmat, methd, k, set){
  
  clustObj<- hclust(as.dist(1 - corrmat), method = methd)
  kList <- lapply(seq(2,k), function(k){
    cutK <- cutree(clustObj, k = k)
    silObj <- cluster::silhouette(cutK, as.dist(1 - corrmat))
    meanSil <- mean(silObj[,"sil_width"])
    df <- data.frame(clust.memb = cutK, ont = rownames(corrmat), set = set, factor = colnames(corrmat))
    df$k<-k 
    df$meanSil<-meanSil
    return(df)
  })
  kDf <- dplyr::bind_rows(kList)
  return(kDf)
}

clusterSets <- function(selectedSets, factors.df, cores){
  
  getOntCorr <- function(idx){
    ont <- unique(selectedSets$variable)[idx]
    selectedFactors <- dplyr::filter(selectedSets, variable == ont)
    
    if(length(selectedFactors$variable)<=2){
      allKCutree <- data.frame(clust.memb = NA, ont = selectedFactors$variable, set = selectedFactors$set, factor = selectedFactors$factor)
      allKCutree$k<-NA 
      allKCutree$meanSil<-NA
    }else{
      factors.df.filtered <- factors.df[,selectedFactors$factor] 
      corrmat <- cor(factors.df.filtered, method = "pearson")
      rownames(corrmat) <- selectedFactors$variable
      colnames(corrmat) <- selectedFactors$factor
      methd = "ward.D2"
      allKCutree <- computeK(corrmat, methd, length(selectedFactors$variable)-1, selectedFactors$set)
    }
    return(allKCutree)
  }
  
  clusterFactors.list <- BiocParallel::bplapply(seq(1,length(unique(selectedSets$variable))), getOntCorr, BPPARAM = BiocParallel::MulticoreParam(cores))
  clusterFactors.df <- dplyr::bind_rows(clusterFactors.list)
  return(clusterFactors.df)
}

plotHeatmapMarkers <- function(matToPlot, metaMat, titl, clusters, mycolors, cuttree){
  
  require(dendsort)
  require(pheatmap)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  matToPlot <- t(matToPlot) %>%
    as_tibble() %>%
    dplyr::mutate(cluster = clusters) %>%
    group_by(cluster) %>%
    summarise_all(list(mean)) %>%
    ungroup() %>%
    tibble::column_to_rownames("cluster")
  
  #matToPlot <- scale(matToPlot, center = TRUE, scale = TRUE)
  
  #mat_cluster_cols <- sort_hclust(hclust(as.dist(1-cor(matToPlot, method = "spearman")), method = "ward.D2"))
  #at_cluster_rows <- sort_hclust(hclust(as.dist(1-cor(t(matToPlot), method = "spearman")), method = "ward.D2"))
  
  mat_cluster_cols <- sort_hclust(hclust(dist(t(matToPlot), method = "euclidean"), method = "ward.D2"))
  mat_cluster_rows <- sort_hclust(hclust(dist(matToPlot, method = "euclidean"), method = "ward.D2"))
  
  # Data frame with column annotations.
  mat_row <- data.frame(cluster = sort(unique(clusters)))
  rownames(mat_row) <- mat_cluster_rows$labels
  
  # List with colors for each annotation.
  mat_colors <- list(cluster = mycolors[which(names(mycolors) %in% mat_row$cluster)])
  
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  mat_breaks <- quantile_breaks(as.matrix(matToPlot), n = 30)
  #mat_breaks <- seq(min(matToPlot),max(matToPlot)+0.01, 0.01)
  #don't use quantile beraks cause your zeros are actually replaced from the matrix
  
  palette <- scico::scico(length(mat_breaks)-1, palette="lajolla", begin = 0.9, end=0.15, direction=-1)
  
  pheatmap(
    mat               = matToPlot,
    color             = palette,
    breaks            = mat_breaks,
    border_color      = NA,
    cluster_cols      = mat_cluster_cols,
    cluster_rows      = mat_cluster_rows,
    cutree_cols = cuttree,
    show_colnames     = TRUE,
    show_rownames     = FALSE,
    #labels_col = metaMat[mat_cluster_cols$labels,]$variable,
    annotation_row    = mat_row,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = titl
  )
}
