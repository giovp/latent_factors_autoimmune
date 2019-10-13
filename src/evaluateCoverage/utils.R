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

getCoverage <- function(dir, Methods, Genesets, gmt.df, cores=2){
  
  coverage.list <- lapply(Methods, function(method){
    
    computeCoverage <- function(set){
      file_list <- list.files(path = dir, pattern = paste0(method,"_.*_",set,"_.*"), full.names = TRUE)
      z_score_df<-lapply(file_list,function(file_path)read_feather(file_path))
      z_score_df <- dplyr::bind_rows(z_score_df)
      z_score_df <- z_score_df %>%
        dplyr::mutate(loading=as.integer(gsub("Loading", "", feature))) %>%
        group_by(nIter, kLatent, loading) %>%
        dplyr::mutate(abs_z_score=abs(z_score)) %>%
        dplyr::filter(abs_z_score==max(abs_z_score)) %>%
        ungroup() %>% 
        group_by(nIter, kLatent) %>%
        dplyr::mutate(pval = 2*pnorm(-abs_z_score),
                      bon_alpha= 0.01 / max(loading)) %>%
        #mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
        dplyr::filter(pval<bon_alpha) %>%
        dplyr::select(c("variable", "nIter", "kLatent")) %>%
        unique() %>%
        tally() %>%
        mutate(coverage=n/dplyr::filter(gmt.df, geneSet==set)$size,
               method=method,
               set=set
        )
      print(paste0("Coverage computed for ", method, " method and ", set, " gene set collection"))
      return(z_score_df)
    }
    
    set.list <- BiocParallel::bplapply(Genesets, computeCoverage, BPPARAM = BiocParallel::MulticoreParam(cores))
    set.df <- dplyr::bind_rows(set.list)
    return(set.df)
  })
  
  coverage.df <- dplyr::bind_rows(coverage.list)
  return(coverage.df)
}

getSignifCollections <- function(dir, Methods, Genesets, gmt.df, cores=2){
  
  coverage.list <- lapply(Methods, function(method){
    
    computeCoverage <- function(set){
      file_list <- list.files(path = dir, pattern = paste0(method,"_.*_",set,"_.*"), full.names = TRUE)
      z_score_df<-lapply(file_list,function(file_path)read_feather(file_path))
      z_score_df <- dplyr::bind_rows(z_score_df)
      z_score_df <- z_score_df %>%
        dplyr::mutate(loading=as.integer(gsub("Loading", "", feature))) %>%
        group_by(nIter, kLatent, loading) %>%
        dplyr::mutate(abs_z_score=abs(z_score)) %>%
        dplyr::filter(abs_z_score==max(abs_z_score)) %>%
        ungroup() %>% 
        group_by(nIter, kLatent) %>%
        dplyr::mutate(pval = 2*pnorm(-abs_z_score),
                      bon_alpha= 0.01 / max(loading)) %>%
        #mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
        dplyr::filter(pval<bon_alpha) %>%
        unique() %>%
        mutate(
               method=method,
               set=set
        )
      print(paste0("Coverage computed for ", method, " method and ", set, " gene set collection"))
      return(z_score_df)
    }
    
    set.list <- BiocParallel::bplapply(Genesets, computeCoverage, BPPARAM = BiocParallel::MulticoreParam(cores))
    set.df <- dplyr::bind_rows(set.list)
    return(set.df)
  })
  
  coverage.df <- dplyr::bind_rows(coverage.list)
  return(coverage.df)
}

getSignifCollections_markers <- function(dir, method, Genesets, gmt.df, cores=2){
  
  computeCoverage <- function(set){
    file_list <- list.files(path = dir, pattern = paste0(method,"_",set,"_.*"), full.names = TRUE)
    z_score_df<-read_feather(file_list)
    z_score_df <- z_score_df %>%
      dplyr::mutate(abs_z_score=abs(z_score)) %>%
      #dplyr::filter(abs_z_score==max(abs_z_score)) %>%
      dplyr::mutate(pval = 2*pnorm(-abs_z_score),
                    bon_alpha= 0.01 / length(unique(index))) %>%
      #mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
      dplyr::filter(pval<bon_alpha) %>%
      unique() %>%
      mutate(
        method=method,
        set=set
      )
    print(paste0("Coverage computed for ", method, " method and ", set, " gene set collection"))
    return(z_score_df)
  }
  
  set.list <- BiocParallel::bplapply(Genesets, computeCoverage, BPPARAM = BiocParallel::MulticoreParam(cores))
  set.df <- dplyr::bind_rows(set.list)
  return(set.df)
}

getCoverage_markers <- function(dir, method, Genesets, gmt.df, cores=2){
  
  computeCoverage <- function(set){
    file_list <- list.files(path = dir, pattern = paste0(method,"_",set,"_.*"), full.names = TRUE)
    z_score_df<-read_feather(file_list)
    z_score_df <- z_score_df %>%
      dplyr::mutate(abs_z_score=abs(z_score)) %>%
      group_by(index) %>%
      dplyr::filter(abs_z_score==max(abs_z_score)) %>%
      dplyr::mutate(pval = 2*pnorm(-abs_z_score),
                    bon_alpha= 0.01 / length(unique(index))) %>%
      #mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
      dplyr::filter(pval<bon_alpha) %>%
      dplyr::select(c("variable", "index")) %>%
      unique() %>%
      tally() %>%
      mutate(coverage=n/dplyr::filter(gmt.df, geneSet==set)$size,
             method=method,
             set=set
      )
    print(paste0("Coverage computed for ", method, " method and ", set, " gene set collection"))
    return(z_score_df)
  }
  
  set.list <- BiocParallel::bplapply(Genesets, computeCoverage, BPPARAM = BiocParallel::MulticoreParam(cores))
  set.df <- dplyr::bind_rows(set.list)
  return(set.df)
}