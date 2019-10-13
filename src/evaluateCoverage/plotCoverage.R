library(tidyverse)
library(feather)
library(gridExtra)
library(clusterProfiler)

source("./src/evaluateCoverage/utils.R")
#get gene sets info
gmt.path <- "./dat/BioBombe/3.build-hetnets/data/"
gmt.path.list <- list.files(gmt.path, pattern = ".gmt", full.names = T)

gene.sets <- c(c3.tft.v6.1.entrez.gmt="GpC3TFT",
               c5.bp.v6.1.entrez.gmt="GpC5BP",
               c5.cc.v6.1.entrez.gmt="GpC5CC",
               c5.mf.v6.1.entrez.gmt="GpC5MF",
               c7.all.v6.1.entrez.gmt="GpC7",
               xcell_all_entrez.gmt="GpXCELL",
               c2.cp.reactome.v6.1.entrez.gmt="GpC2CPREACTOME",
               h.all.v6.1.entrez.gmt="GpH",
               c2.cp.biocarta.v6.1.entrez.gmt="GpC2CPBIOCARTA",
               c2.cp.kegg.v6.1.entrez.gmt="GpC2CPKEGG",
               c2.cgp.v6.1.entrez.gmt="GpC2CPG",
               metabaser_maps.gmt="GpMETABASER",
               `wikipathways-20190410-gmt-Homo_sapiens.gmt`="GpWIKIPATH")
gmt.df <- list()
for(i in seq(1,length(gene.sets))){
  print(i)
  set <- names(gene.sets[i])
  gmt.path <- grepl(set, gmt.path.list)
  gmt.obj <- read.gmt(gmt.path.list[gmt.path])
  gmt.df[[i]]=data.frame(geneSet=gene.sets[[i]],
                 size=length(unique(gmt.obj$ont)))
}
gmt.df <- dplyr::bind_rows(gmt.df)

gmt.path <- "./dat/BioBombe/3.build-hetnets/data/"
gmt.path.list <- list.files(gmt.path, pattern = ".gmt", full.names = T)
gene.sets <- c(c3.tft.v6.1.entrez.gmt="GpC3TFT",
               c5.bp.v6.1.entrez.gmt="GpC5BP",
               c5.cc.v6.1.entrez.gmt="GpC5CC",
               c5.mf.v6.1.entrez.gmt="GpC5MF",
               c7.all.v6.1.entrez.gmt="GpC7",
               xcell_all_entrez.gmt="GpXCELL",
               c2.cp.reactome.v6.1.entrez.gmt="GpC2CPREACTOME",
               h.all.v6.1.entrez.gmt="GpH",
               c2.cp.biocarta.v6.1.entrez.gmt="GpC2CPBIOCARTA",
               c2.cp.kegg.v6.1.entrez.gmt="GpC2CPKEGG",
               c2.cgp.v6.1.entrez.gmt="GpC2CPG",
               metabaser_maps.gmt="GpMETABASER",
               `wikipathways-20190410-gmt-Homo_sapiens.gmt`="GpWIKIPATH")
gmt.list <- lapply(names(gene.sets), function(collection){
  gmt.path <- grepl(collection, gmt.path.list)
  gmt.obj <- read.gmt(gmt.path.list[gmt.path])
  gmt.obj$collection <- gene.sets[which(names(gene.sets)==collection)]
  return(gmt.obj)
})

# gmt.df <- dplyr::bind_rows(gmt.list) %>% 
#   as_tibble() %>%
#   dplyr::rename(ENTREZID = "gene")
# write_tsv(gmt.df , "./res/allGeneSets.tsv")

#set some variables

#get full computed coverage dataframe
input.path <- "./SLE_pipeline/output/hetnets_output"
path.to.plot <- "./SLE_pipeline/output/figures/collectionCoverageplot.pdf"
path.to.save <- "./SLE_pipeline/output/allSignifCollections.feather"
path.to.save.filt <- "./SLE_pipeline/output/allSignifCollectionsFiltered.feather"
#takes long
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
coverage <- getCoverage(input.path,
                        c("hpf","cogaps","scvi","lda"), 
                        gene.sets, 
                        gmt.df, 
                        cores = length(gene.sets))

allSignifCollections <- getSignifCollections(input.path,
                        c("hpf","cogaps","scvi","lda"), 
                        gene.sets, 
                        gmt.df, 
                        cores = length(gene.sets))

#get best run for plotting coverage of the best run
loss.df <- getBestRun("./SLE_pipeline/output/evaluate_output",
                      "SLE_pipeline",
                      c("hpf","cogaps","scvi","lda"),
                      seq(16,40,2))

loss.df <- loss.df %>%
  mutate(kLatent=as.character(k),
         nIter = as.character(samp))
##SAVE GENE SETS
write_feather(allSignifCollections, path.to.save)

allSignifCollectionsFiltered <- allSignifCollections %>%
  inner_join(.,loss.df, by=c("kLatent","nIter","method"))

path.to.save.filt = "./SLE_pipeline/output/allSignifCollectionsFilteredNotMax.feather"
write_feather(allSignifCollectionsFiltered, path.to.save.filt)

###PLOT
mycolors <- c(scico::scico(2, palette="devon", begin=0.4, end=0.6),
              scico::scico(2, palette="lajolla", begin=0.3, end=0.6))

coverage_ra <- coverage %>%
  inner_join(.,loss.df, by=c("kLatent","nIter","method"))

coverage_sle <- coverage %>%
  inner_join(.,loss.df, by=c("kLatent","nIter","method"))

names(mycolors)<-c("scHPF","scCoGAPS","scVI","LDA")

coverage$set <- gsub("Gp","",coverage$set)
coverage$set <- gsub("METABASER","METABASE",coverage$set)
coverage$set <- gsub("WIKIPATH","WIKIPAT",coverage$set)
coverage$set <- gsub("H","HALLMARK",coverage$set)
coverage$set <- gsub("WIKIPAT","WIKIPATHWAY",coverage$set)

coverage$method <- gsub("hpf","scHPF",coverage$method)
coverage$method <- gsub("lda","LDA",coverage$method)
coverage$method <- gsub("cogaps","scCoGAPS",coverage$method)
coverage$method <- gsub("scvi","scVI",coverage$method)

total_coverage <- coverage %>%
  group_by(method,set,kLatent) %>%
  summarise(meanCoverage=mean(coverage)) %>%
  ungroup() %>%
  group_by(method,kLatent) %>%
  summarise(meanCoverage=mean(meanCoverage)) %>%
  mutate(set = "TOTAL COLLECTIONS",
         group = "total collections")

#plot <- 
coverage %>%
  group_by(method,set,kLatent) %>%
  summarise(meanCoverage=mean(coverage)) %>%
  dplyr::mutate(group="single collection") %>%
  dplyr::bind_rows(.,total_coverage) %>%
  ungroup() %>%
  dplyr::mutate(set=factor(set, levels=c(unique(coverage$set), unique(total_coverage$set))))%>%
  ggplot(.)+
  geom_point(
    aes(x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method), shape=as.factor(group))
    )+
  geom_line(
    aes(x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method), group=as.factor(method)
        ),stat="identity",show.legend=FALSE)+
  scale_shape_manual(values=c(1,16))+
  facet_wrap(.~set, scales = "free_y", nrow = 2)+
  scale_color_manual(values = mycolors)+
  theme_bw() +
  labs(x="number latent dimensions", color="method", y="mean collection coverage", shape="type")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+#,
        #strip.background = element_blank(), 
        #strip.placement = "outside")
ggsave(path.to.plot,device="pdf", height=4, width=15)


#### MEAN ACROSS COLLECTIONS
input.path <- "./SLE_pipeline/output/hetnets_output"
#takes long
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
coverage_sle <- getCoverage(input.path,
                        c("hpf","cogaps","scvi","lda"), 
                        gene.sets, 
                        gmt.df, 
                        cores = length(gene.sets))

loss_sle <- getBestRun("./SLE_pipeline/output/evaluate_output",
                      "SLE_pipeline",
                      c("hpf","cogaps","scvi","lda"),
                      seq(16,40,2))

loss_sle <- loss_sle %>%
  mutate(kLatent=as.character(k),
         nIter = as.character(samp))

coverage_sle <- coverage_sle %>%
  inner_join(.,loss_sle, by=c("kLatent","nIter","method"))

input.path <- "./RA_pipeline/output/hetnets_output"
#takes long
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
coverage_ra <- getCoverage(input.path,
                            c("hpf","cogaps","scvi","lda"), 
                            gene.sets, 
                            gmt.df, 
                            cores = length(gene.sets))

loss_ra <- getBestRun("./RA_pipeline/output/evaluate_output",
                      "RA_pipeline",
                      c("hpf","cogaps","scvi","lda"),
                      seq(16,40,2))

loss_ra <- loss_ra %>%
  mutate(kLatent=as.character(k),
         nIter = as.character(samp))

coverage_ra <- coverage_ra %>%
  inner_join(.,loss_ra, by=c("kLatent","nIter","method"))

coverage <- dplyr::bind_rows(coverage_sle, coverage_ra) %>%
  dplyr::mutate(dataset=ifelse(dataset == "SLE_pipeline","SLE","RA"))

mycolors <- c(scico::scico(2, palette="devon", begin=0.4, end=0.6),
              scico::scico(2, palette="lajolla", begin=0.3, end=0.6))
names(mycolors)<-c("scHPF","scCoGAPS","scVI","LDA")

coverage$set <- gsub("Gp","",coverage$set)
coverage$set <- gsub("METABASER","METABASE",coverage$set)
coverage$set <- gsub("WIKIPATH","WIKIPAT",coverage$set)
coverage$set <- gsub("H","HALLMARK",coverage$set)
coverage$set <- gsub("WIKIPAT","WIKIPATHWAY",coverage$set)

coverage$method <- gsub("hpf","scHPF",coverage$method)
coverage$method <- gsub("lda","LDA",coverage$method)
coverage$method <- gsub("cogaps","scCoGAPS",coverage$method)
coverage$method <- gsub("scvi","scVI",coverage$method)

coverage %>%
  group_by(method, kLatent, dataset) %>%
  summarise(meanCoverage=mean(coverage)) %>%
  ggplot(.)+
  geom_point(
    aes(x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method))
  )+
  geom_line(
    aes(
      x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method), group=as.factor(method)
    ),stat="identity",show.legend=FALSE)+
  scale_shape_manual(values=c(1,16))+
  facet_wrap(.~dataset)+
  scale_color_manual(values = mycolors)+
  theme_bw() +
  labs(x="K", color="method", y="mean collection coverage", shape="type")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave("./res/figures/coverage_mean_collections.pdf",device="pdf", height=3, width=6)

# shift_legend2 <- function(p) {
#   # ...
#   # to grob
#   gp <- ggplotGrob(p)
#   facet.panels <- grep("^panel", gp[["layout"]][["name"]])
#   empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
#   empty.facet.panels <- facet.panels[empty.facet.panels]
#   
#   # establish name of empty panels
#   empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
#   names <- empty.facet.panels$name
#   # example of names:
#   #[1] "panel-3-2" "panel-3-3"
#   
#   # now we just need a simple call to reposition the legend
#   ret <- lemon::reposition_legend(p, 'center', panel=names)
#   return(ret)
# }
# 
# tosave <- shift_legend2(plot)
# 
# coverage %>%
#   group_by(method,set,kLatent) %>%
#   summarise(meanCoverage=mean(coverage)) %>%
#   ungroup() %>%
#   group_by(method,kLatent) %>%
#   summarise(meanCoverage=mean(meanCoverage)) %>%
#   ggplot(.)+
#   geom_point(aes(x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method)))+
#   geom_line(
#     aes(x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method), group=as.factor(method)),
#     stat="identity")+
#   #facet_wrap(.~set, scales = "free_y", nrow = 2)+
#   scale_color_manual(values = mycolors)+
#   theme_bw() +
#   labs(x="number latent dimensions", color="method", y="mean collection coverage")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
# ggsave(path.to.save.total, device="pdf", height=3, width=4)
