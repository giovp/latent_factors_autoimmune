getBestRun <- function(dir,dataset, methods, list_kLatents){
  
  file.list <- list.files(dir, pattern =".*rds" ,full.names = TRUE)
  
  loss.df <- data.frame(hpf=list_kLatents,
                        cogaps=list_kLatents,
                        scvi=list_kLatents,
                        lda=list_kLatents,
                        k=list_kLatents)
  rownames(loss.df) <- list_kLatents
  
  for(method in methods){
    print(method)
    for(k in list_kLatents){
      print(k)
      path <- file.list[grepl(method,file.list)&grepl(k, file.list)]
      latent.obj <- readRDS(path)
      if(grepl("lda", path)){
        loss.idx <- latent.obj@samp[which(latent.obj@loss==max(latent.obj@loss))[1]]
      }else{
        loss.idx <- latent.obj@samp[which(latent.obj@loss==min(latent.obj@loss))[1]]
      }
      loss.df[as.character(k),method]=loss.idx
    }
  }
  
  loss.df <- loss.df %>%
    gather(key = "method", value = "samp", -k) %>%
    mutate(dataset=dataset)
  
  return(loss.df)
}

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

#ont = "Immune_response_Oncostatin_M_signaling_via_MAPK"
clusterSets <- function(selectedSets, factors.df){
  clusterFactors.list <- lapply(unique(selectedSets$variable), function(ont){
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
  })
  clusterFactors.df <- dplyr::bind_rows(clusterFactors.list)
  return(clusterFactors.df)
}

plotHeatmapLoadingsCluster <- function(matToPlot, metaMat, titl){
  
  require(dendsort)
  require(pheatmap)
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  mat_cluster_cols <- sort_hclust(hclust(dist(t(matToPlot))))
  mat_cluster_rows <- sort_hclust(hclust(dist(matToPlot)))
  
  # Data frame with column annotations.
  mat_row <- data.frame(n_factors = as.factor(metaMat[mat_cluster_rows$labels,]))
  rownames(mat_row) <- mat_cluster_rows$labels
  
  # List with colors for each annotation.
  mat_colors <- list(n_factors=scico(length(levels(mat_row$n_factors)), palette="lajolla", direction = -1))
  names(mat_colors$n_factors) <- levels(mat_row$n_factors)
  
  quantile_breaks <- function(mat, n = 30) {
    xs <- mat %>% gather() %>% pull(value)
    breaks <- quantile(xs, probs = seq(0, 0.7, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  #mat_breaks <- quantile_breaks(matToPlot, n = 30)
  #don't use quantile beraks cause your zeros are actually replaced from the matrix
  
  mat_breaks <- c(seq(min(matToPlot), quantile(matToPlot %>% gather() %>% pull(value),0.9), length.out = 15),
                  seq(quantile(matToPlot %>% gather() %>% pull(value),0.90001), max(matToPlot), length.out = 15))
  
  palette <- scico::scico(length(mat_breaks)-1, palette="devon", begin = 0.9, end=0.1, direction=-1)
  
  pheatmap(
    mat               = matToPlot,
    color             = palette,
    breaks            = mat_breaks,
    border_color      = NA,
    cluster_cols      = mat_cluster_cols,
    cluster_rows      = mat_cluster_rows,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    #labels_col = metaMat[mat_cluster_cols$labels,]$variable,
    annotation_row    = mat_row,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = titl
  )
  
}