buildLatentNeighborsObjectL2Norm <- function(factors.in, loadings.in, K,S){
  
  out <- new("latentNeighbors", 
             factors.df = factors.in[,lapply(.SD, function(x) x/sqrt(sum(x^2)))],
             loadings.df = loadings.in[,lapply(.SD, function(x) x/sqrt(sum(x^2)))],
             kLatent = as.integer(K),
             samp=as.integer(S))
  
  return(out)
}

buildLatentNeighborsObject <- function(factors.in, 
                                       loadings.in, 
                                       kLatent.in,
                                       corrToMeanPatt.in,
                                       loss.in,
                                       S){
  
  out <- new("latentNeighbors", 
             factors.df = factors.in,
             loadings.df = loadings.in,
             kLatent = kLatent.in,
             loss = loss.in,
             corrToMeanPatt=corrToMeanPatt.in,
             samp=as.integer(S))
  
  return(out)
}

collectCOGAPSResult <- function(file_list){
  
  factors.list = list()
  loadings.list = list()
  kLatent<-c()
  chiSquare <- c()
  corrToMeanPattern <- c()
  samp_list <- c()
  
  for (i in seq(1,length(file_list))){
    
    cogaps <- readRDS(file_list[[i]])
    
    k=cogaps@metadata$params@nPatterns
    samp=cogaps@metadata$params@seed
    
    factors <- as.data.table(cogaps@sampleFactors)
    colnames(factors) <- paste0("Factor",seq(1,dim(factors)[2]),"_",samp,"_",k)
    loadings <- as.data.table(cogaps@featureLoadings)
    colnames(loadings) <- paste0("Loading",seq(1,dim(factors)[2]),"_",samp,"_",k)
    
    factors.list[[i]]<-factors
    loadings.list[[i]]<-loadings
    
    print(dim(cogaps@featureLoadings))
    kLatent <- c(kLatent, dim(cogaps@featureLoadings)[2])
    chiSquare <- c(chiSquare, cogaps@metadata$meanChiSq)
    corrToMeanPattern <- c(corrToMeanPattern,mean(unlist(cogaps@metadata$CorrToMeanPattern)))
    samp_list <- c(samp_list, samp)
  }
  
  #aggregate tables
  factors.full <- setDT(unlist(factors.list, recursive = FALSE), check.names = TRUE)[]
  loadings.full <- setDT(unlist(loadings.list, recursive = FALSE), check.names = TRUE)[]
  
  #build the latentneighbors object. 
  #It also l2norm the matrices column wise if the other function is used
  latent.obj <- buildLatentNeighborsObject(factors.full, 
                                           loadings.full,
                                           kLatent,
                                           corrToMeanPattern,
                                           chiSquare,
                                           samp_list)
  
  return(latent.obj)
}

collectLDAResult <- function(file_list){
  
  factors.list = list()
  loadings.list = list()
  kLatent<-c()
  BF <- c()
  corrToMeanPattern <- c()
  samp_list <- c()
  
  for (i in seq(1,length(file_list))){
    
    lda <- readRDS(file_list[[i]])
    
    k=lda$k
    samp=lda$seed
    
    factors <- as.data.table(lda$fit$omega)
    colnames(factors) <- paste0("Factor",seq(1,dim(factors)[2]),"_",samp,"_",k)
    loadings <- as.data.table(lda$fit$theta)
    colnames(loadings) <- paste0("Loading",seq(1,dim(factors)[2]),"_",samp,"_",k)
    
    factors.list[[i]]<-factors
    loadings.list[[i]]<-loadings
    
    print(dim(lda$fit$theta))
    kLatent <- c(kLatent, dim(lda$fit$theta)[2])
    BF <- c(BF, lda$BF)
    corrToMeanPattern <- c(corrToMeanPattern,0)
    samp_list <- c(samp_list, samp)
  }
  
  #aggregate tables
  factors.full <- setDT(unlist(factors.list, recursive = FALSE), check.names = TRUE)[]
  loadings.full <- setDT(unlist(loadings.list, recursive = FALSE), check.names = TRUE)[]
  
  #build the latentneighbors object. 
  #It also l2norm the matrices column wise if the other function is used
  latent.obj <- buildLatentNeighborsObject(factors.full, 
                                           loadings.full,
                                           kLatent,
                                           corrToMeanPattern,
                                           BF,
                                           samp_list)
  
  return(latent.obj)
}

collectSCVIResult <- function(factors_list, loadings_list, metrics_list){
  
  factors.list = list()
  loadings.list = list()
  kLatent<-c()
  llh <- c()
  corrToMeanPattern <- c()
  samp_list <- c()
  
  for (i in seq(1,length(factors_list))){
    
    metrics <- read_feather(metrics_list[[i]])
    
    k=unique(metrics$n_latent)
    samp=unique(metrics$seed)
    
    factors <- read_feather(factors_list[[i]]) %>% as.data.table()
    colnames(factors) <- paste0("Factor",seq(1,dim(factors)[2]),"_",samp,"_",k)
    loadings <- read_feather(loadings_list[[i]]) %>% as.data.table()
    colnames(loadings) <- paste0("Loading",seq(1,dim(factors)[2]),"_",samp,"_",k)
    
    factors.list[[i]]<-factors
    loadings.list[[i]]<-loadings
    
    print(dim(factors))
    kLatent <- c(kLatent, dim(factors)[2])
    llh <- c(llh, tail(metrics$elbo_test_set)[1])
    corrToMeanPattern <- c(corrToMeanPattern,0)
    samp_list <- c(samp_list, samp)
  }
  
  #aggregate tables
  factors.full <- setDT(unlist(factors.list, recursive = FALSE), check.names = TRUE)[]
  loadings.full <- setDT(unlist(loadings.list, recursive = FALSE), check.names = TRUE)[]
  
  #build the latentneighbors object. 
  #It also l2norm the matrices column wise if the other function is used
  latent.obj <- buildLatentNeighborsObject(factors.full, 
                                           loadings.full,
                                           kLatent,
                                           corrToMeanPattern,
                                           llh,
                                           samp_list)
  return(latent.obj)
  
}

collectHPFResult <- function(factors_list, loadings_list, metrics_list){
  
  factors.list = list()
  loadings.list = list()
  kLatent<-c()
  llh <- c()
  corrToMeanPattern <- c()
  samp_list <- c()
  
  for (i in seq(1,length(factors_list))){
    
    metrics <- read_feather(metrics_list[[i]])
    
    k=unique(metrics$k)
    samp=unique(metrics$seed)
    
    factors <- read_feather(factors_list[[i]]) %>% as.data.table()
    colnames(factors) <- paste0("Factor",seq(1,dim(factors)[2]),"_",samp,"_",k)
    loadings <- read_feather(loadings_list[[i]]) %>% as.data.table()
    colnames(loadings) <- paste0("Loading",seq(1,dim(factors)[2]),"_",samp,"_",k)
    
    factors.list[[i]]<-factors
    loadings.list[[i]]<-loadings
    
    print(dim(factors))
    kLatent <- c(kLatent, dim(factors)[2])
    llh <- c(llh, tail(metrics$loss)[1])
    corrToMeanPattern <- c(corrToMeanPattern,0)
    samp_list <- c(samp_list, samp)
    
  }
  
  #aggregate tables
  factors.full <- setDT(unlist(factors.list, recursive = FALSE), check.names = TRUE)[]
  loadings.full <- setDT(unlist(loadings.list, recursive = FALSE), check.names = TRUE)[]
  
  #build the latentneighbors object. 
  #It also l2norm the matrices column wise if the other function is used
  latent.obj <- buildLatentNeighborsObject(factors.full, 
                                           loadings.full,
                                           kLatent,
                                           corrToMeanPattern,
                                           llh,
                                           samp_list)
  return(latent.obj)
  
}

amariDistance <- function(load1, load2) {
  K <- dim(load1)[2]
  C <- cor(load1, load2)
  return(1 - (sum(apply(C, 1, max)) + sum(apply(C, 2, max))) / (2 * K))
}

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

computeAmariDistance <- function(latent.in){
  
  samp_len=length(latent.in@samp)
  #split_idx <- chunk2(seq(1,ncol(latent.in@loadings.df)),samp)
  
  amari.distances <- lapply(seq(1,samp_len-1), function(i) {
    
    idx.i <- grepl(paste0("_",latent.in@samp[i],"_"),colnames(latent.in@loadings.df))
    load.a <- latent.in@loadings.df[,idx.i, with=FALSE]
    #print(i)
    amari.dist <- sapply(seq(i+1,samp_len), function(j){
      
      idx.j <- grepl(paste0("_",latent.in@samp[j],"_"),colnames(latent.in@loadings.df))
      load.b <- latent.in@loadings.df[,idx.j, with=FALSE]
      #print(j)
      am.dist <- amariDistance(load.a, load.b)
      return(am.dist)
    })
    return(amari.dist)
  })
  
  latent.in@amari.dist <- unlist(amari.distances)
  return(latent.in)
}

computeSVCCA <- function(latent.in, svcca, threshold){
  
  samp_len=length(latent.in@samp)
  #split_idx <- chunk2(seq(1,ncol(latent.in@loadings.df)),samp)
  
  mean.cca.all <- lapply(seq(1,samp_len-1), function(i) {
    
    idx.i <- grepl(paste0("_",latent.in@samp[i],"_"),colnames(latent.in@factors.df))
    load.a <- latent.in@loadings.df[,idx.i, with=FALSE] %>% data.table::transpose()
    #print(i)
    mean.cca <- sapply(seq(i+1,samp_len), function(j){
      
      idx.j <- grepl(paste0("_",latent.in@samp[j],"_"),colnames(latent.in@factors.df))
      load.b <- latent.in@loadings.df[,idx.j, with=FALSE] %>% data.table::transpose()
      #print(j)
      result <- svcca$robust_cca_similarity(
        reticulate::r_to_py(as.matrix(load.a), convert = TRUE),
        reticulate::r_to_py(as.matrix(load.b), convert = TRUE), 
        threshold = threshold)
      
      return(mean(result$cca_coef1))
    })
    return(mean(mean.cca))
  })
  
  latent.in@mean.svcca <- mean(unlist(mean.cca.all))
  return(latent.in)
}

entropyMixingBinary <- function(frequency){
  entropy = -frequency*log(frequency)-(1-frequency)*log(1-frequency)
  if(!is.na(entropy)){
    return(entropy)
  } else {
    #print("Error: na value in calculating entropy")
    return(0)
  }
}

entropyMixing <- function(frequency){
  entropy = -frequency*log(frequency)
  if(!is.na(entropy)){
    return(entropy)
  } else {
    #print("Error: na value in calculating entropy")
    return(0)
  }
}

estimateEntropy <- function(vec, clust){
  ent_list <- sapply(unique(clust),function(cls){
    freq = mean(vec==cls)
    ent = entropyMixingBinary(freq)
    return(ent)
  })
  return(mean(ent_list))
}

estimateEntropyMixing <- function(latent.in, cluster.list.in, samp=10){
  
  factors.list <- sapply(colnames(latent.in@loadings.df),function(name){
    fact.id <- unlist(stringr::str_split(name, "Loading|_"))
    return(fact.id[2])
  })
  
  names(factors.list) <- cluster.list.in
  entropy.vec <- vector()
  for(i in seq(1,length(unique(names(factors.list))))){
    entropy=0
    idx.clust=0
    for(j in seq(1,length(unique(factors.list)))){
      freq = mean(factors.list[names(factors.list)==i]==j)
      entropy = entropy + entropyMixingBinary(freq)
      #print(entropy)
    }
    idx.clust=length(unique(factors.list[names(factors.list)==i]))
    entropy.vec[i]=entropy/idx.clust
  }
  return(entropy.vec)
}
