library(tidyverse)

dataset="SLEcelseq"

signif_mean <- read_tsv(paste0("./res/",dataset,"/cdb/significant_means.txt"))
pval_mean <- read_tsv(paste0("./res/",dataset,"/cdb/pvalues_means.txt"))
# View(signif_mean)
# View(pval_mean)

gather_mean <- signif_mean %>%
  dplyr::select(-one_of(colnames(signif_mean)[2:10])) %>%
  gather(key="cluster", value="value", -interacting_pair) %>%
  na.omit()

gather_result <- pval_mean %>%
  dplyr::select(-one_of(colnames(pval_mean)[2:9])) %>%
  gather(key="cluster", value="value", -interacting_pair) %>%
  dplyr::mutate(mean=map_dbl(value, function(x) as.numeric(strsplit(x, " \\| ")[[1]][[1]])),
                pval=map_dbl(value, function(x) as.numeric(strsplit(x, " \\| ")[[1]][[2]]))) %>%
  dplyr::select(-one_of("value")) %>%
  inner_join(.,gather_mean, by=c("interacting_pair", "cluster")) %>%
  dplyr::mutate(pval=ifelse(pval==0,0.0001, pval)) %>%
  dplyr::mutate(log10_pval=-log10(pval),
                log2_mean=log2(mean)) 

#plot results of cdb analysis  
cell_type="Tcell"
plot_result <- gather_result[grepl(paste0("^",cell_type,"_"), gather_result$cluster),]
cluster <- cluster::agnes(plot_result[,c("mean","log10_pval")],
                         metric = "euclidean",
                         keep.diss=FALSE,
                         keep.data=FALSE)

plot_result <- plot_result[cluster$order,]
plot_result$interacting_pair <- factor(plot_result$interacting_pair, levels=unique(plot_result$interacting_pair))

plot_result$cluster <- factor(plot_result$cluster, levels=sort(unique(plot_result$cluster)))

plot_result %>%
  dplyr::filter((log2_mean>1 | log2_mean<(-1)) & log10_pval>2) %>%
  ggplot(.)+
  geom_point(aes(x=cluster, y=interacting_pair, size=log10_pval, colour=log2_mean))+
  scale_size(range = c(1, 5))+
  scico::scale_color_scico(palette = "acton", direction = -1)+
  labs(size="-log10(pval)", colour="log2(mean)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))#+
ggsave("./RA_pipeline/output/figures/Tcell_interactions.pdf", device="pdf", width=13, height = 10)





