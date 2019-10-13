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
input.path <- "./RA_pipeline/output/hetnets_output"
path.to.plot <- "./RA_pipeline/output/figures/collectionCoverageplot_markers.pdf"
path.to.save <- "./RA_pipeline/output/allSignifCollections_markers.feather"
#takes long
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))

allSignifCollections_markers <- getSignifCollections_markers(input.path,
                        "markers", 
                        gene.sets, 
                        gmt.df, 
                        cores = length(gene.sets))

coverage_markers <- getCoverage_markers(input.path,
                                "markers", 
                                gene.sets, 
                                gmt.df, 
                                cores = length(gene.sets))

coverage <- getCoverage(input.path,
                        c("hpf","cogaps","scvi","lda"), 
                        gene.sets, 
                        gmt.df, 
                        cores = length(gene.sets))
##SAVE GENE SETS
write_feather(allSignifCollections, path.to.save)

###PLOT
mycolors <- c(scico::scico(2, palette="devon", begin=0.4, end=0.6),
              scico::scico(2, palette="lajolla", begin=0.3, end=0.6),
              scico::scico(1, palette="tokyo", begin=0.3, end=0.6))

names(mycolors)<-c("scHPF","scCoGAPS","scVI","LDA","markers")

coverage$set <- gsub("Gp","",coverage$set)
coverage$set <- gsub("METABASER","METABASE",coverage$set)
coverage$set <- gsub("WIKIPATH","WIKIPAT",coverage$set)
coverage$set <- gsub("H","HALLMARK",coverage$set)
coverage$set <- gsub("WIKIPAT","WIKIPATHWAY",coverage$set)

coverage_markers$set <- gsub("Gp","",coverage_markers$set)
coverage_markers$set <- gsub("METABASER","METABASE",coverage_markers$set)
coverage_markers$set <- gsub("WIKIPATH","WIKIPAT",coverage_markers$set)
coverage_markers$set <- gsub("H","HALLMARK",coverage_markers$set)
coverage_markers$set <- gsub("WIKIPAT","WIKIPATHWAY",coverage_markers$set)

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
         group = "total collections") %>%
  dplyr::filter(kLatent <= 18)

total_coverage_markers <- coverage_markers %>%
  group_by(method,set,index) %>%
  summarise(meanCoverage=mean(coverage)) %>%
  ungroup() %>%
  group_by(method) %>%
  summarise(meanCoverage=mean(meanCoverage)) %>%
  mutate(set = "TOTAL COLLECTIONS",
         group = "total collections") %>%
  dplyr::mutate(kLatent = as.character(17))

#plot <- 
total_coverage_markers %>%
  dplyr::bind_rows(.,total_coverage) %>%
  dplyr::select(-one_of("group")) %>%
  #dplyr::mutate(set=factor(set, levels=c(unique(coverage$set), unique(total_coverage$set))))%>%
  ggplot(.)+
  geom_point(
    aes(x=as.factor(kLatent), y=meanCoverage, colour=as.factor(method))
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
ggsave(path.to.plot,device="pdf", height=4, width=4)

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
