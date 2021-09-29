library(Seurat)
library(ggplot2)
library(ggpubr)
library(gplots)
library(gridExtra)
library(RCAv2)
library(dbplyr)
library(ggrepel)
library(gridExtra)
library(Mfuzz)
library(Biobase)

PBMC.integrated.good <- readRDS("PBMC.integrated.Mono.good.rds")


#### comparative density plot across different severities
DefaultAssay(PBMC.integrated.good) <- "integrated"
PBMC.integrated.good <- FindNeighbors(object = PBMC.integrated.good, 
                                      reduction = "pca", dims = 1:23, 
                                      k.param = 300, verbose = TRUE, compute.SNN = FALSE)

umap.df.enrich <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","cell.annotation",
                                                           "sampleBatch","severity"))

# get graph matrix
graph_nnMat <- Matrix::Matrix(PBMC.integrated.good@graphs$integrated_nn, sparse = TRUE)
severity.df <- FetchData(PBMC.integrated.good, vars = c("severity","sampleBatch"))

saveRDS(graph_nnMat,"12.2D_comparative_enrichment_plot/graph_nnMat.rds")
graph_nnMat <- readRDS("12.2D_comparative_enrichment_plot/graph_nnMat.rds")


### average enrichment ratio across 100 times of downsampling
# add infection status 
infect.stat <- read.delim("../infection_status_for_all_sampleBatches.txt",
                          sep = "\t", header = TRUE)
infect.stat.sub <- infect.stat[,c("sampleBatch","cat.x")]

severity.df$barcode <- rownames(severity.df)
severity.infect.df <- merge(severity.df,
                            infect.stat.sub, by = "sampleBatch")


## Mild1
## 19 27 21 18
## Mild2
## 9 32 37 31 
## Moderate
## 3 8 13 13 
## Severe
## 1 4 8 2
## Critical
## 1 6 10 18

comp <- c("Asymptomatic", "Mild2","Moderate","Severe","Critical")
ref <- "Mild1"

ref.down.list <- list(c(9,27,21,18),
                      c(6,16,20,18),
                      c(3,12,21,6),
                      c(3,18,21,18))
comp.down.list <- list(c(9,27,21,18),
                       c(3,8,10,9),
                       c(1,4,7,2),
                       c(1,6,7,6))
cat.list <- c(2,3,4,5)

enrichment.matrix.all <- as.data.frame(matrix(nrow = nrow(umap.df.enrich),
                                              ncol = length(comp)))
rownames(enrichment.matrix.all) <- rownames(umap.df.enrich)
colnames(enrichment.matrix.all) <- comp

for (a in seq_along(comp)){
  comp.a <- comp[a]
  severity.infect.df.comp <- severity.infect.df[(severity.infect.df$severity == comp.a),]
  severity.infect.df.ref <- severity.infect.df[(severity.infect.df$severity == ref),]
  
  if (a == 1){
    ref_num <- nrow(severity.infect.df.ref)
    # calculate total number of ref cells in each cell's 300 neighbours
    graph_nnMat.ref <- graph_nnMat[,(severity.infect.df.ref$barcode)]
    graph_nnMat.ref.sum <- Matrix::rowSums(graph_nnMat.ref)
    
    # calculate number of ref and comp cells in each cell's 300 neighbours
    graph_nnMat.comp <- graph_nnMat[,(severity.infect.df.comp$barcode)]
    graph_nnMat.comp.sum <- Matrix::rowSums(graph_nnMat.comp)
    
    comp.ref.df <- data.frame(comp = graph_nnMat.comp.sum,
                              ref = graph_nnMat.ref.sum)
    
    comp_FC <- (nrow(severity.infect.df.comp)) / (ref_num)
    
    comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
    comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
    
    comp.ref.df[(comp.ref.df$comp==0 & comp.ref.df$ref==0),"log2Ratio"] <- 0
    
    comp.ref.df <- comp.ref.df[rownames(umap.df.enrich),]
    
    enrichment.matrix.all[,a] <- comp.ref.df$log2Ratio
  } else {
    ref.down.list.a <- ref.down.list[[a-1]]
    comp.down.list.a <- comp.down.list[[a-1]]
    
    print(ref.down.list.a)
    print(comp.down.list.a)
    
    # downsampling
    ratio.matrix.all <- as.data.frame(matrix(nrow = nrow(umap.df.enrich),
                                             ncol = 100))
    rownames(ratio.matrix.all) <- rownames(umap.df.enrich)
    for (i in 1:100){
      comp.all.sample <- c()
      ref.all.sample <- c()
      for (b in 1:length(cat.list)){
        cat.b <- cat.list[b]
        
        # get comp sample in cat b
        sample.comp.b <- sample(unique(severity.infect.df.comp[(severity.infect.df.comp$cat.x == cat.b),"sampleBatch"]),
                                comp.down.list.a[b])
        comp.all.sample <- c(comp.all.sample,sample.comp.b)
        
        # get comp sample in cat b
        sample.ref.b <- sample(unique(severity.infect.df.ref[(severity.infect.df.ref$cat.x == cat.b),"sampleBatch"]),
                               ref.down.list.a[b])
        ref.all.sample <- c(ref.all.sample,sample.ref.b)
        
      }
      # get subset of downsampling in comp
      severity.infect.df.comp.i <- severity.infect.df.comp[(severity.infect.df.comp$sampleBatch %in% comp.all.sample),]
      # get subset of downsampling in ref
      severity.infect.df.ref.i <- severity.infect.df.ref[(severity.infect.df.ref$sampleBatch %in% ref.all.sample),]
      
      
      ref_num <- nrow(severity.infect.df.ref.i)
      # calculate total number of ref cells in each cell's 300 neighbours
      graph_nnMat.ref <- graph_nnMat[,(severity.infect.df.ref.i$barcode)]
      graph_nnMat.ref.sum <- Matrix::rowSums(graph_nnMat.ref)
      
      # calculate number of ref and comp cells in each cell's 300 neighbours
      graph_nnMat.comp <- graph_nnMat[,(severity.infect.df.comp.i$barcode)]
      graph_nnMat.comp.sum <- Matrix::rowSums(graph_nnMat.comp)
      
      comp.ref.df <- data.frame(comp = graph_nnMat.comp.sum,
                                ref = graph_nnMat.ref.sum)
      
      comp_FC <- (nrow(severity.infect.df.comp.i)) / (ref_num)
      
      comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
      comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
      
      comp.ref.df[(comp.ref.df$comp==0 & comp.ref.df$ref==0),"log2Ratio"] <- 0
      
      comp.ref.df <- comp.ref.df[rownames(umap.df.enrich),]
      
      ratio.matrix.all[,i] <- comp.ref.df$log2Ratio
      
    }
    
    ratio.matrix.all.mean <- Matrix::rowMeans(ratio.matrix.all)
    
    enrichment.matrix.all[,a] <- ratio.matrix.all.mean
    
  }
  
}

write.table(enrichment.matrix.all, "12.2D_comparative_enrichment_plot_downsampling/table_of_average_enrichmeent_ratio_mild1_reference_across_downsampling.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)



## profile enrichment of each condition against Mild1

comp <- c("Asymptomatic", "Mild2","Moderate","Severe","Critical")

for (i in seq_along(comp)){
  comp.i <- comp[i]
  
  umap.df.enrich$log2Ratio <- enrichment.matrix.all[,comp.i]
  umap.df.enrich[(umap.df.enrich$log2Ratio > 4),"log2Ratio"] <- 4
  umap.df.enrich[(umap.df.enrich$log2Ratio < -4),"log2Ratio"] <- -4
  
  ggplot()+
    geom_point(umap.df.enrich, 
               mapping = aes(x=UMAP_1, y=UMAP_2, 
                             color=log2Ratio), 
               size=0.7, stroke=0)+
    scale_color_gradientn(colors = c("yellow","#666200","black","#420042","magenta" ),
                          values=c(1, .65, .5, .35, 0),
                          limits=c(-4,4))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
    ggsave(paste0("12.2D_comparative_enrichment_plot_downsampling/",i,".2D_comparative_enrichment_plot_of_",
                  comp.i,"_vs_",ref,".png"),device = "png",
           width = 7, height = 6)
  
}



##### using mfuzz to find cell states enriched and depleted
enrichment.matrix.all <- read.table("12.2D_comparative_enrichment_plot_downsampling/table_of_average_enrichmeent_ratio_mild1_reference_across_downsampling.txt",
                                    sep = "\t", header = TRUE, row.names = 1)
enrichment.matrix.all.sub <- enrichment.matrix.all[(abs(enrichment.matrix.all$Mild2) > 2 |
                                                      abs(enrichment.matrix.all$Moderate) > 2 |
                                                      abs(enrichment.matrix.all$Severe) > 2 |
                                                      abs(enrichment.matrix.all$Critical) > 2 ),
                                                   c("Mild2","Moderate","Severe","Critical")]



enrichment.matrix.all.sub.s <- ExpressionSet(assayData=as.matrix(enrichment.matrix.all.sub))
m1 <- mestimate(enrichment.matrix.all.sub.s)
cl_wt<-mfuzz(enrichment.matrix.all.sub.s,c=20,m=m1)
saveRDS(cl_wt, "12.2D_comparative_enrichment_plot_downsampling/cl_wt.rds")
cl_wt <- readRDS("12.2D_comparative_enrichment_plot_downsampling/cl_wt.rds")

for (i in seq(1,20,1)){
  cluster <- cl_wt$cluster
  cluster.names.i <- names(cluster[cluster==i])
  
  enrichment.matrix.all.sub.i <- enrichment.matrix.all.sub[cluster.names.i,
                                                           c("Mild2","Moderate",
                                                             "Severe","Critical")]
  
  bk = unique(c(seq(-6, 6, length=80)))
  col = colorRampPalette(c("darkblue","white","darkorange" ))(length(bk)-1)
  
  pdf(paste0("12.2D_comparative_enrichment_plot_downsampling/heatmap_of_enrichment_ratios_for_mfuzz_cluster_",
             i,".pdf"))
  heatmap.2(as.matrix(enrichment.matrix.all.sub.i), Rowv = F, Colv = F, col = col, 
            breaks=bk,symm=F, symbreaks=F, scale="none",
            trace="none", cexCol = .5,density.info = "none", labRow = NA,
            main = paste0("cell = ",nrow(enrichment.matrix.all.sub.i) ))
  dev.off()
}




increase <- c() # cell states enriched with severity 
decrease <- c() # cell states depleted with severity
for (i in seq(1,20,1)){
  cluster <- cl_wt$cluster
  cluster.names.i <- names(cluster[cluster==i])
  
  enrich.matrix.all.sub.i <- enrichment.matrix.all.sub[cluster.names.i,]
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



umap.df.enrich$type <- "none"
umap.df.enrich[cluster.names.c,"type"] <- "increase"
umap.df.enrich[cluster.names.b,"type"] <- "decrease"

umap.df.enrich$type <- factor(umap.df.enrich$type, 
                              levels = c("none","increase","decrease"))

ggplot()+
  geom_point(umap.df.enrich, 
             mapping = aes(x=UMAP_1, y=UMAP_2, color=type),
             size=.5, stroke=0)+
  scale_color_manual(values = c("grey","darkorange2","blue"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("12.2D_comparative_enrichment_plot_downsampling/integrated_UMAP_highlighting_cells_with_different_enrichment_trends.png",
         device = "png", width = 6, height = 5)


######
###### marker genes in cell states enriched with severity
######
### CD14 
umap.df.enrich.cd4tcm <- umap.df.enrich[(umap.df.enrich$cell.annotation == "CD14 Mono"),]

# test sample representation
umap.df.enrich.cd4tcm.inc <- umap.df.enrich.cd4tcm[(umap.df.enrich.cd4tcm$type=="increase"),]
freq.i <- data.frame(table(umap.df.enrich.cd4tcm.inc$sampleBatch))
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
ggsave("12.2D_comparative_enrichment_plot_downsampling/sample_count_analysis_for_downsampling_in_DEG/1.Barplot_of_sample_counts_in_INCREASE_CD14_Mono.png",
         device = "png")


umap.df.enrich.cd4tcm.inc <- umap.df.enrich.cd4tcm[(umap.df.enrich.cd4tcm$type!="increase"),]
freq.i <- data.frame(table(umap.df.enrich.cd4tcm.inc$sampleBatch))
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
ggsave("12.2D_comparative_enrichment_plot_downsampling/sample_count_analysis_for_downsampling_in_DEG/1.Barplot_of_sample_counts_in_NOT_INCREASE_CD14_Mono.png",
         device = "png")
## no need downsampling
## DEG
cell.high <- rownames(umap.df.enrich.cd4tcm[(umap.df.enrich.cd4tcm$type=="increase"),])
cell.low <- rownames(umap.df.enrich.cd4tcm[(umap.df.enrich.cd4tcm$type!="increase"),])

de.ij <- RCAv2:::ComputePairWiseDE(object = PBMC.integrated.good@assays$RNA@data, 
                                   cells.1 = cell.high,
                                   cells.2 = cell.low,
                                   features = NULL,
                                   logfc.threshold = 0,
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

write.table(de.ij.sig[which(de.ij.sig$avg_logFC > 0),"gene"],
            "12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/upregulated_genes_in_CD14_Mono_increase_vs_non-increase.txt",
            col.names = FALSE, quote = FALSE, row.names = FALSE)
write.table(de.ij.sig[which(de.ij.sig$avg_logFC < 0),"gene"],
            "12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/downregulated_genes_in_CD14_Mono_increase_vs_non-increase.txt",
            col.names = FALSE, quote = FALSE, row.names = FALSE)
write.table(de.ij,
            "12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/DEG_in_CD14_Mono_increase_vs_non-increase.txt",
            col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")


ggplot()+
  geom_point(de.ij, mapping = aes(x = avg_logFC, y = logP), color="grey")+
  geom_point(de.ij.sig[which(de.ij.sig$avg_logFC>0),],
             mapping = aes(x = avg_logFC, y = logP), color="red")+
  geom_point(de.ij.sig[which(de.ij.sig$avg_logFC<0),],
             mapping = aes(x = avg_logFC, y = logP), color="blue")+
  geom_label_repel(data=de.ij.sig[which(de.ij.sig$avg_logFC> 0.5),], 
                   aes(x = avg_logFC, y = logP,label=gene),
                   box.padding   = 1, 
                   point.padding = 0.6,
                   segment.color = 'red',
                   size=6, color="red",max.overlaps =100)+
  ylim(c(0,650))+
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))
ggsave("12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/volcano_plot_of_DEGs_in_CD14_Mono_increase_vs_non-increase.pdf",
         device = "pdf",width = 5, height = 3)





### CD16
umap.df.enrich.cd4tcm <- umap.df.enrich[(umap.df.enrich$cell.annotation == "CD16 Mono"),]

# test sample representation
umap.df.enrich.cd4tcm.inc <- umap.df.enrich.cd4tcm[(umap.df.enrich.cd4tcm$type=="increase"),]
freq.i <- data.frame(table(umap.df.enrich.cd4tcm.inc$sampleBatch))
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
ggsave("12.2D_comparative_enrichment_plot_downsampling/sample_count_analysis_for_downsampling_in_DEG/1.Barplot_of_sample_counts_in_INCREASE_CD16_Mono.png",
         device = "png")


umap.df.enrich.cd4tcm.inc <- umap.df.enrich.cd4tcm[(umap.df.enrich.cd4tcm$type!="increase"),]
freq.i <- data.frame(table(umap.df.enrich.cd4tcm.inc$sampleBatch))
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
ggsave("12.2D_comparative_enrichment_plot_downsampling/sample_count_analysis_for_downsampling_in_DEG/1.Barplot_of_sample_counts_in_NOT_INCREASE_CD16_Mono.png",
         device = "png")
## too few cells after downsampling


######
###### marker genes in cell states depleted with severity
######

umap.df.enrich.decrease <- umap.df.enrich[(umap.df.enrich$type == "decrease"),]
umap.df.enrich.decrease.tb <- data.frame(table(umap.df.enrich.decrease$cell.annotation))

umap.df.enrich.decrease.tb <- umap.df.enrich.decrease.tb[(umap.df.enrich.decrease.tb$Freq > 100),]
unique_cell <- as.character(umap.df.enrich.decrease.tb$Var1)

for (cell in unique_cell){
  umap.df.enrich.cytocd8t <- umap.df.enrich[(umap.df.enrich$cell.annotation == cell),]
  
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
    ggsave(paste0("12.2D_comparative_enrichment_plot_downsampling/sample_count_analysis_for_downsampling_in_DEG/2.Barplot_of_sample_counts_in_BRIDGE_",
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
    ggsave(paste("12.2D_comparative_enrichment_plot_downsampling/sample_count_analysis_for_downsampling_in_DEG/2.Barplot_of_sample_counts_in_NOT_BRIDGE_",
                 gsub(" ","_",cell),".png"),
           device = "png")
  
  ## DEG
  cell.high <- rownames(umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type=="decrease"),])
  cell.low <- rownames(umap.df.enrich.cytocd8t[(umap.df.enrich.cytocd8t$type!="decrease"),])
  
  de.ij <- RCAv2:::ComputePairWiseDE(object = PBMC.integrated.good@assays$RNA@data, 
                                     cells.1 = cell.high,
                                     cells.2 = cell.low,
                                     features = NULL,
                                     logfc.threshold = 0,
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
  
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC > 0),"gene"],
              paste0("12.2D_comparative_enrichment_plot_downsampling//DEG_analysis/upregulated_genes_in_",
                     gsub(" ","_",cell),"_decrease_vs_non-bridge.txt"),
              col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC < 0),"gene"],
              paste0("12.2D_comparative_enrichment_plot_downsampling//DEG_analysis/downregulated_genes_in_",
                     gsub(" ","_",cell),"_decrease_vs_non-decrease.txt"),
              col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(de.ij,
              paste0("12.2D_comparative_enrichment_plot_downsampling//DEG_analysis/DEG_in_",
                     gsub(" ","_",cell),"_decrease_vs_non-decrease.txt"),
              col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  ggplot()+
    geom_point(de.ij, mapping = aes(x = avg_logFC, y = logP), color="grey")+
    geom_point(de.ij.sig[which(de.ij.sig$avg_logFC>0),],
               mapping = aes(x = avg_logFC, y = logP), color="red")+
    geom_point(de.ij.sig[which(de.ij.sig$avg_logFC<0),],
               mapping = aes(x = avg_logFC, y = logP), color="blue")+
    geom_label_repel(data=de.ij.sig[which(de.ij.sig$avg_logFC> 0.4),], 
                     aes(x = avg_logFC, y = logP,label=gene),
                     box.padding   = .5, 
                     point.padding = 0.3,
                     segment.color = 'red',
                     size=6, color="red",max.overlaps =100)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+ylim(c(0,100))
    ggsave(paste0("12.2D_comparative_enrichment_plot_downsampling//DEG_analysis/volcano_plot_of_DEGs_in_",
                  gsub(" ","_",cell),"_decrease_vs_non-decrease.pdf"),
           device = "pdf",width = 5, height = 3)
  
  
  
  
  
}


### UMAP heatmap showing metagene expression

# marker genes in the cell states enriched with severity in CD14 Mono
go <- read.delim("12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/GO_of_upregulated_genes_in_CD14_Mono_increase_vs_non-increase/CUSTOM140587348375024.human.enrichr.reports.txt",
                 header = TRUE, sep = "\t")
go <- go[order(go$Adjusted.P.value, decreasing = FALSE),]
go <- go[1:3,]
gene.1 <- unlist(strsplit(go$Genes[1], split = ";"))
gene.2 <- unlist(strsplit(go$Genes[2], split = ";"))
gene.3 <- unlist(strsplit(go$Genes[3], split = ";"))
gene.c <- unique(c(gene.1, gene.2, gene.3))

gene.exp <- as.data.frame(as.matrix(PBMC.integrated.good@assays$RNA@data[gene.c,]))
gene.exp.scale <- as.data.frame(t(scale(t(gene.exp), center = T, scale = T)))
gene.exp.mean <- Matrix::colMeans(gene.exp.scale)

umap.df <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","cell.annotation"))
gene.exp.mean <- gene.exp.mean[rownames(umap.df)]
umap.df$gene <- gene.exp.mean

umap.df[which(umap.df$gene > 2),"gene"] <- 2
umap.df[which(umap.df$gene < -2),"gene"] <- -2

umap.df.sub <- umap.df[(umap.df$cell.annotation == "CD14 Mono"),]

ggplot()+
  geom_point(data = umap.df,
             aes(x=UMAP_1, y=UMAP_2), 
             size=.5,stroke=0, color="black",alpha=.2)+
  geom_point(data = umap.df.sub[(abs(umap.df.sub$gene) < 1 ),],
             aes(x=UMAP_1, y=UMAP_2, color=gene), 
             size=.5,stroke=0)+
  geom_point(data = umap.df.sub[(abs(umap.df.sub$gene) >= 1),],
             aes(x=UMAP_1, y=UMAP_2, color=gene), 
             size=.5,stroke=0)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_colour_gradient2(low = "blue",high = "red",mid = "grey", limits=c(-2,2))
ggsave("12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/UMAP_of_heatmap_showing_neutrophil_genes_in_Mono.png",
         width = 5, height = 4)


# marker genes in the cell states depleted with severity in CD14 Mono
go <- read.delim("12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/GO_of_upregulated_genes_in_CD14_Mono_decrease_vs_non-bridge/CUSTOM140587348375104.human.enrichr.reports.txt",
                 header = TRUE, sep = "\t")
go <- go[order(go$Adjusted.P.value, decreasing = FALSE),]
go <- go[1:2,]
gene.1 <- unlist(strsplit(go$Genes[1], split = ";"))
gene.2 <- unlist(strsplit(go$Genes[2], split = ";"))
gene.c <- unique(c(gene.1, gene.2))


gene.exp <- as.data.frame(as.matrix(PBMC.integrated.good@assays$RNA@data[gene.c,]))
gene.exp.scale <- as.data.frame(t(scale(t(gene.exp), center = T, scale = T)))
gene.exp.mean <- Matrix::colMeans(gene.exp.scale)

umap.df <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","cell.annotation"))
gene.exp.mean <- gene.exp.mean[rownames(umap.df)]
umap.df$gene <- gene.exp.mean

umap.df[which(umap.df$gene > 2),"gene"] <- 2
umap.df[which(umap.df$gene < -2),"gene"] <- -2

umap.df.sub <- umap.df[(umap.df$cell.annotation == "CD14 Mono"),]

ggplot()+
  geom_point(data = umap.df,
             aes(x=UMAP_1, y=UMAP_2), 
             size=.5,stroke=0, color="black", alpha=.2)+
  geom_point(data = umap.df.sub[(abs(umap.df.sub$gene) < 1 ),],
             aes(x=UMAP_1, y=UMAP_2, color=gene), 
             size=.5,stroke=0)+
  geom_point(data = umap.df.sub[(abs(umap.df.sub$gene) >= 1),],
             aes(x=UMAP_1, y=UMAP_2, color=gene), 
             size=.5,stroke=0)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_colour_gradient2(low = "blue",high = "red",mid = "grey", limits=c(-2,2))
ggsave("12.2D_comparative_enrichment_plot_downsampling/DEG_analysis/UMAP_of_heatmap_showing_antigen_presentation_genes_in_Mono.png",
         width = 5, height = 4)


