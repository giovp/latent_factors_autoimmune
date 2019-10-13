
parseMetricsFiles <- function(path_files, pattern){
  
  list_names = list.files(path_files, pattern = pattern,
                          full.names = TRUE)
  df_list <- list()
  for (i in seq(1,length(list_names))){
    file <- readRDS(list_names[[i]])
    
    df <- tibble(
      kCluster = file@kCluster,
      kLatent = list(file@kLatent),
      clust.memb = list(file@clust.memb),
      entropy = list(file@entropy),
      sil = list(file@sil),
      amari = list(file@amari.dist),
      corrToMeanPatt = list(file@corrToMeanPatt),
      loss = list(file@loss),
      mean.svcca = file@mean.svcca
    )
    
    df_list[[i]]<-df
  }
  
  full_df <- dplyr::bind_rows(df_list)
  return(full_df)
}

plotMetrics <- function(tib, labells, mycolors){
  
  plot = tib %>%
    dplyr::mutate(meanKLatent=map_dbl(kLatent, function(x) mean(x))) %>%
    dplyr::mutate(meanAmari=map_dbl(amari, function(x) mean(x)),
                  meanSil=map_dbl(sil, function(x) mean(x))) %>%
    dplyr::select(kCluster,meanKLatent,meanSil, meanAmari, algo, mean.svcca) %>%
    unnest() %>%
    gather(key = "metric", value = "value", -kCluster, -kCluster,  -meanKLatent, -algo) %>%
    dplyr::mutate(
      metric = case_when(
        metric == "meanSil" ~ "mean Silhouette score",
        metric == "meanAmari" ~ "mean Amari distance",
        metric == "mean.svcca" ~ "mean SVCCA score"
      )
    ) %>%
    ggplot() +
    geom_point(aes(x=as.factor(kCluster), y=value, colour=as.factor(algo)))+
    geom_line(aes(x=as.factor(kCluster), y=value, colour=as.factor(algo), group = as.factor(algo)))+
    facet_wrap(.~metric, scales = "free_y")+
    scale_color_manual(values = mycolors, labels = labells)+
    labs(y="value", x="number latent dimensions", colour="Method")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(plot)
}

plotLoss <- function(tib, labs_facet){
  
  plot <- tib %>%
    dplyr::mutate(meanAmari = map_dbl(amari, function(x) mean(x)),
                  meanCorr = map_dbl(corrToMeanPatt, function(x) mean(x)),
                  loss = map(loss, function(x) tail(x, 10))) %>%
    dplyr::select(meanAmari,kCluster,meanCorr, loss, algo) %>% 
    unnest() %>%
    ggplot() +
    geom_beeswarm(aes(x=as.factor(kCluster), 
                      y=loss), 
                  size=1, 
                  alpha=0.5) +
    geom_boxplot(aes(x=as.factor(kCluster), y=loss), 
                 alpha=0.2,
                 lwd=0.4)+
    facet_wrap(.~algo, scales = "free_y", 
               labeller = labeller(algo = labs_facet))+
    labs(x="number latent dimensions", y="loss")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(plot)
} 
