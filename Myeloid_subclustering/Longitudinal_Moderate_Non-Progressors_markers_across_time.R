library(ggplot2)
library(ggpubr)
library(gplots)
library(gridExtra)
library(RCAv2)
library(dbplyr)
library(gplots)
library(ggrepel)
library(gridExtra)
library(Seurat)
library(Mfuzz)
library(Biobase)

PBMC.integrated.good <- readRDS("PBMC.integrated.Mono.good.rds")

all.df <- FetchData(PBMC.integrated.good, 
                    vars = c("UMAP_1","UMAP_2","cell.annotation","sampleBatch", "severity"))

graph_nnMat <- readRDS("12.2D_comparative_enrichment_plot/graph_nnMat.rds")

########## Moderate non progression different time course

mod.df <- read.delim("../new_RCA_whole/logitudinal_information_all/Moderate_nonprogression.txt",
                       header = TRUE, sep = "\t")
mod.df$sampleBatch <- paste(mod.df$batch, mod.df$sample, sep = "_")

# using rising as reference
early.df <- mod.df[(mod.df$type == "rising systemic inflammation"),]
all.df.early <- all.df[(all.df$sampleBatch %in% early.df$sampleBatch),]

ggplot()+
  geom_point(all.df, 
             mapping = aes(x=UMAP_1, y=UMAP_2), 
             size=1, stroke=0, color="grey")+
  geom_point(all.df.early, 
             mapping = aes(x=UMAP_1, y=UMAP_2), 
             size=1, stroke=0, color="red")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave(paste0("14.Different_stage_of_moderate_nonprogression/0.UMAP_with_moderate_rising.png"),device = "png",
         width = 7, height = 6)

ref_num <- nrow(all.df.early)
# calculate total number of ref cells in each cell's 300 neighbours
graph_nnMat.ref <- graph_nnMat[,rownames(all.df.early)]
graph_nnMat.ref.sum <- Matrix::rowSums(graph_nnMat.ref)

all.stage <- c("peak","late")
all.stage.full <- c("peak COVID-19 hyper-inflammation",
                    "resolution or persistent COVID-19 inflammation")

enrich.matrix.all <- as.data.frame(matrix(nrow = nrow(all.df),
                                          ncol = length(all.stage)))
rownames(enrich.matrix.all) <- rownames(all.df)
colnames(enrich.matrix.all) <- all.stage
for (i in seq_along(all.stage)){
  all.stage.i <- all.stage[i] 
  all.stage.full.i <- all.stage.full[i]
  stage.df <- mod.df[(mod.df$type == all.stage.full.i),]
  
  all.df.i <- all.df[(all.df$sampleBatch %in% stage.df$sampleBatch),]
  
  # calculate number of ref and comp cells in each cell's 300 neighbours
  graph_nnMat.comp <- graph_nnMat[,rownames(all.df.i)]
  graph_nnMat.comp.sum <- Matrix::rowSums(graph_nnMat.comp)
  
  comp.ref.df <- data.frame(comp = graph_nnMat.comp.sum,
                            ref = graph_nnMat.ref.sum)
  
  comp_FC <- (nrow(all.df.i)) / (ref_num)
  
  comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
  comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
  
  comp.ref.df[(comp.ref.df$comp==0 & comp.ref.df$ref==0),"log2Ratio"] <- 0
  comp.ref.df <- comp.ref.df[rownames(all.df),]
  
  enrich.matrix.all[,i] <- comp.ref.df$log2Ratio
  
  #plot UMAP
  all.df$log2Ratio <- comp.ref.df$log2Ratio
  
  all.df[(all.df$log2Ratio > 4),"log2Ratio"] <- 4
  all.df[(all.df$log2Ratio < -4),"log2Ratio"] <- -4
  
  ggplot()+
      geom_point(all.df, 
                 mapping = aes(x=UMAP_1, y=UMAP_2, 
                               color=log2Ratio), 
                 size=0.8, stroke=0)+
      scale_color_gradientn(colors = c("yellow","#666200","black","#420042","magenta" ),
                            values=c(1, .65, .5, .35, 0),
                            limits=c(-4,4))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
      ggsave(paste0("14.Different_stage_of_moderate_nonprogression/",i,".2D_comparative_enrichment_plot_of_",
                    all.stage.i,"_vs_early.png"),device = "png",
             width = 7, height = 6)
    
    ggplot()+
      geom_point(all.df, 
                 mapping = aes(x=UMAP_1, y=UMAP_2), 
                 size=1, stroke=0, color="grey")+
      geom_point(all.df.i, 
                 mapping = aes(x=UMAP_1, y=UMAP_2), 
                 size=1, stroke=0, color="red")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
      ggsave(paste0("14.Different_stage_of_moderate_nonprogression/",i
                    ,".UMAP_with_moderate_",all.stage.i,".png"),device = "png",
             width = 7, height = 6)
  
}

enrich.matrix.all.sub <- enrich.matrix.all[(abs(enrich.matrix.all$peak) > 2 |
                                              abs(enrich.matrix.all$late) > 2 ),]

enrichment.matrix.all.sub.s <- ExpressionSet(assayData=as.matrix(enrich.matrix.all.sub))
m1 <- mestimate(enrichment.matrix.all.sub.s)
cl_wt<-mfuzz(enrichment.matrix.all.sub.s,c=20,m=m1)
saveRDS(cl_wt, "14.Different_stage_of_moderate_nonprogression/cl_wt.rds")
cl_wt <- readRDS("14.Different_stage_of_moderate_nonprogression/cl_wt.rds")

for (i in seq(1,20,1)){
  cluster <- cl_wt$cluster
  cluster.names.i <- names(cluster[cluster==i])
  
  enrich.matrix.all.sub.i <- enrich.matrix.all.sub[cluster.names.i,]
  
  bk = unique(c(seq(-6, 6, length=80)))
  col = colorRampPalette(c("darkblue","white","darkorange" ))(length(bk)-1)
  
  pdf(paste0("14.Different_stage_of_moderate_nonprogression/heatmap_of_enrichment_ratios_for_mfuzz_cluster_",
             i,".pdf"))
  heatmap.2(as.matrix(enrich.matrix.all.sub.i), Rowv = F, Colv = F, col = col, 
            breaks=bk,symm=F, symbreaks=F, scale="none",
            trace="none", cexCol = .5,density.info = "none", labRow = NA,
            main = paste0("cell = ",nrow(enrich.matrix.all.sub.i) ))
  dev.off()
}

increase <- c() # cell states enriched with time
decrease <- c() # cell states depleted with time
for (i in seq(1,20,1)){
  cluster <- cl_wt$cluster
  cluster.names.i <- names(cluster[cluster==i])
  
  enrich.matrix.all.sub.i <- enrich.matrix.all.sub[cluster.names.i,]
  colmean.i <- colMeans(enrich.matrix.all.sub.i)
  colmean.i <- c(0, colmean.i)
  
  cor.i <- cor(colmean.i, seq(1, length(colmean.i),1), method = "pearson")
  
  if (cor.i > 0.5) {
    increase <- c(increase, i)
  } 
  
  if (cor.i < -0.5) {
    decrease <- c(decrease, i)
  }
  
}



cluster <- cl_wt$cluster
cluster.names.c <- names(cluster[cluster %in% increase])
cluster.names.b <- names(cluster[cluster %in% decrease])




all.df$type <- "none"
all.df[cluster.names.c,"type"] <- "increase"
all.df[cluster.names.b,"type"] <- "decrease"

all.df$type <- factor(all.df$type, levels = c("none","increase","decrease"))

ggplot()+
  geom_point(all.df, 
             mapping = aes(x=UMAP_1, y=UMAP_2, color=type),
             size=.5, stroke=0)+
  scale_color_manual(values = c("grey","darkorange2","blue"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("14.Different_stage_of_moderate_nonprogression/integrated_UMAP_highlighting_cells_with_different_enrichment_trends.png",
         device = "png", width = 6, height = 5)



######
###### marker genes in cell states enriched with time
######

all.df.increase <- all.df[(all.df$type == "increase"),]
all.df.increase.tb <- data.frame(table(all.df.increase$cell.annotation))

all.df.increase.tb <- all.df.increase.tb[(all.df.increase.tb$Freq > 1000),]
unique_cell <- as.character(all.df.increase.tb$Var1)

for (cell in unique_cell){
  umap.df.enrich.cytocd8t <- all.df[(all.df$cell.annotation == cell),]
  
  # test sample representation
  umap.df.enrich.cytocd8t.inc <- umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type=="increase"),]
  freq.i <- data.frame(table(umap.df.enrich.cytocd8t.inc$sampleBatch))
  freq.i <- freq.i[order(freq.i$Freq,decreasing = TRUE),]
  freq.i$Var1 <- factor(freq.i$Var1, levels = freq.i$Var1)
  ggplot(freq.i, mapping = aes(x=Var1, y=Freq))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.1*sum(freq.i$Freq),
               color="red")+
    geom_hline(yintercept = 0.2*sum(freq.i$Freq),
               color="red")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0("number of cell = ",sum(freq.i$Freq),"\n",
                   "number of sample = ", nrow(freq.i)))
    ggsave(paste0("14.Different_stage_of_moderate_nonprogression/sample_count_analysis_for_downsampling_in_DEG/2.Barplot_of_sample_counts_in_INCREASE_",
                  gsub(" ","_",cell),".png"),
           device = "png")
  
  
  umap.df.enrich.cytocd8t.inc <- umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type!="increase"),]
  freq.i <- data.frame(table(umap.df.enrich.cytocd8t.inc$sampleBatch))
  freq.i <- freq.i[order(freq.i$Freq,decreasing = TRUE),]
  freq.i$Var1 <- factor(freq.i$Var1, levels = freq.i$Var1)
  ggplot(freq.i, mapping = aes(x=Var1, y=Freq))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.1*sum(freq.i$Freq),
               color="red")+
    geom_hline(yintercept = 0.2*sum(freq.i$Freq),
               color="red")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0("number of cell = ",sum(freq.i$Freq),"\n",
                   "number of sample = ", nrow(freq.i)))
    ggsave(paste0("14.Different_stage_of_moderate_nonprogression/sample_count_analysis_for_downsampling_in_DEG/2.Barplot_of_sample_counts_in_NOT_INCREASE_",
                 gsub(" ","_",cell),".png"),
           device = "png")
  
  ## DEG
  cell.high <- rownames(umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type=="increase"),])
  cell.low <- rownames(umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type!="increase"),])
  
  de.ij <- RCAv2:::ComputePairWiseDE(object = PBMC.integrated.good@assays$RNA@data, 
                                     cells.1 = cell.high,
                                     cells.2 = cell.low,
                                     features = NULL,
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox",
                                     min.pct = 0.1,
                                     min.diff.pct = -Inf,
                                     only.pos = FALSE,
                                     max.cells.per.ident = Inf,
                                     random.seed = 1,
                                     min.cells.group = 3,
                                     pseudocount.use = 1,
                                     MeanExprsThrs = 0.5,
                                     p.adjust.methods = "BH")
  
  de.ij$gene <- rownames(de.ij)
  de.ij <- de.ij[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                      x = de.ij$gene, invert = TRUE),]
  
  de.ij[(de.ij$p_val_adj==0), "p_val_adj"] <- 1e-300
  de.ij$logP <- -log10(de.ij$p_val_adj)
  de.ij.sig <- de.ij[(de.ij$p_val_adj < 0.1 & abs(de.ij$avg_logFC) > 0.25 ),]
  
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC > 0),],
              paste0("14.Different_stage_of_moderate_nonprogression/DEG_analysis/upregulated_genes_in_",
                     gsub(" ","_",cell),"_increase_vs_non-increase.txt"),
              col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC < 0),],
              paste0("14.Different_stage_of_moderate_nonprogression/DEG_analysis/downregulated_genes_in_",
                     gsub(" ","_",cell),"_increase_vs_non-increase.txt"),
              col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  
}



# union of upregulated gene
all.gene <- c()
for (i in seq_along(unique_cell)){
  cell.i <- unique_cell[i]
  deg.file <- read.table(paste0("./14.Different_stage_of_moderate_nonprogression/DEG_analysis/upregulated_genes_in_",
                                gsub(" ","_", cell.i),"_increase_vs_non-increase.txt"),
                         sep="\t", header = TRUE)
  
  all.gene <- c(all.gene,deg.file$gene)
}
all.gene <- unique(all.gene)
write.table(all.gene,
            "14.Different_stage_of_moderate_nonprogression/DEG_analysis/union_of_upregulated_genes_in_increase.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


######
###### marker genes in cell states depleted with time
######

all.df.increase <- all.df[(all.df$type == "decrease"),]
all.df.increase.tb <- data.frame(table(all.df.increase$cell.annotation))

all.df.increase.tb <- all.df.increase.tb[(all.df.increase.tb$Freq > 1000),]
unique_cell <- as.character(all.df.increase.tb$Var1)

for (cell in unique_cell){
  umap.df.enrich.cytocd8t <- all.df[(all.df$cell.annotation == cell),]
  
  # test sample representation
  umap.df.enrich.cytocd8t.inc <- umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type=="decrease"),]
  freq.i <- data.frame(table(umap.df.enrich.cytocd8t.inc$sampleBatch))
  freq.i <- freq.i[order(freq.i$Freq,decreasing = TRUE),]
  freq.i$Var1 <- factor(freq.i$Var1, levels = freq.i$Var1)
  ggplot(freq.i, mapping = aes(x=Var1, y=Freq))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.1*sum(freq.i$Freq),
               color="red")+
    geom_hline(yintercept = 0.2*sum(freq.i$Freq),
               color="red")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0("number of cell = ",sum(freq.i$Freq),"\n",
                   "number of sample = ", nrow(freq.i)))
    ggsave(paste0("14.Different_stage_of_moderate_nonprogression/sample_count_analysis_for_downsampling_in_DEG/2.Barplot_of_sample_counts_in_DECREASE_",
                  gsub(" ","_",cell),".png"),
           device = "png")
  
  
  umap.df.enrich.cytocd8t.inc <- umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type!="decrease"),]
  freq.i <- data.frame(table(umap.df.enrich.cytocd8t.inc$sampleBatch))
  freq.i <- freq.i[order(freq.i$Freq,decreasing = TRUE),]
  freq.i$Var1 <- factor(freq.i$Var1, levels = freq.i$Var1)
  ggplot(freq.i, mapping = aes(x=Var1, y=Freq))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.1*sum(freq.i$Freq),
               color="red")+
    geom_hline(yintercept = 0.2*sum(freq.i$Freq),
               color="red")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0("number of cell = ",sum(freq.i$Freq),"\n",
                   "number of sample = ", nrow(freq.i)))
    ggsave(paste("14.Different_stage_of_moderate_nonprogression/sample_count_analysis_for_downsampling_in_DEG/2.Barplot_of_sample_counts_in_NOT_DECREASE_",
                 gsub(" ","_",cell),".png"),
           device = "png")
  
  ## DEG
  cell.high <- rownames(umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type=="decrease"),])
  cell.low <- rownames(umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type!="decrease"),])
  
  de.ij <- RCAv2:::ComputePairWiseDE(object = PBMC.integrated.good@assays$RNA@data, 
                                     cells.1 = cell.high,
                                     cells.2 = cell.low,
                                     features = NULL,
                                     logfc.threshold = 0.25,
                                     test.use = "wilcox",
                                     min.pct = 0.1,
                                     min.diff.pct = -Inf,
                                     only.pos = FALSE,
                                     max.cells.per.ident = Inf,
                                     random.seed = 1,
                                     min.cells.group = 3,
                                     pseudocount.use = 1,
                                     MeanExprsThrs = 0.5,
                                     p.adjust.methods = "BH")
  
  de.ij$gene <- rownames(de.ij)
  de.ij <- de.ij[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                      x = de.ij$gene, invert = TRUE),]
  
  de.ij[(de.ij$p_val_adj==0), "p_val_adj"] <- 1e-300
  de.ij$logP <- -log10(de.ij$p_val_adj)
  de.ij.sig <- de.ij[(de.ij$p_val_adj < 0.1 & abs(de.ij$avg_logFC) > 0.25 ),]
  
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC > 0),],
              paste0("14.Different_stage_of_moderate_nonprogression/DEG_analysis/upregulated_genes_in_",
                     gsub(" ","_",cell),"_decrease_vs_non-decrease.txt"),
              col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC < 0),],
              paste0("14.Different_stage_of_moderate_nonprogression/DEG_analysis/downregulated_genes_in_",
                     gsub(" ","_",cell),"_decrease_vs_non-decrease.txt"),
              col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  
}

# union of upregulated gene
all.gene <- c()
for (i in seq_along(unique_cell)){
  cell.i <- unique_cell[i]
  deg.file <- read.table(paste0("./14.Different_stage_of_moderate_nonprogression/DEG_analysis/upregulated_genes_in_",
                                gsub(" ","_", cell.i),"_decrease_vs_non-decrease.txt"),
                         sep="\t", header = TRUE)
  
  all.gene <- c(all.gene,deg.file$gene)
}
all.gene <- unique(all.gene)
write.table(all.gene,
            "14.Different_stage_of_moderate_nonprogression/DEG_analysis/union_of_upregulated_genes_in_decrease.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


