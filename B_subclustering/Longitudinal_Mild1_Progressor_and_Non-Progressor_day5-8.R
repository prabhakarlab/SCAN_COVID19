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

PBMC.integrated.good <- readRDS("PBMC.integrated.B.good.rds")


############
########### enrichment plot between progression and non-progression
###########
graph_nnMat <- readRDS("12.2D_comparative_enrichment_plot/graph_nnMat.rds")

all.df <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","sampleBatch","batch","severity",
                                                   "cell.annotation","age"))

### progression
prog <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression.txt",
                   header = TRUE, sep="\t")
prog$sampleBatch <- paste(prog$batch,prog$sample, sep="_")

prog <- prog[(prog$sampleBatch %in% all.df$sampleBatch),]
prog <- prog[(prog$type == "rising systemic inflammation"),]

all.df.prog <- all.df[(all.df$sampleBatch %in% prog$sampleBatch),]

### non-progression
nonprog <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_nonprogression_rising_stage.txt",
                      header = TRUE, sep="\t")
nonprog <- nonprog[(nonprog$sampleBatch %in% all.df$sampleBatch),]
# remove redundant samples for each patient
nonprog <- nonprog[!(nonprog$sampleBatch %in% c("2_CG1","12_CG30","3_CG47","8_CG55")),]

all.df.nonprog <- all.df[(all.df$sampleBatch %in% nonprog$sampleBatch),]

### plot non-progression as ref
ref_num <- nrow(all.df.nonprog)
# calculate total number of ref cells in each cell's 300 neighbours
graph_nnMat.ref <- graph_nnMat[,rownames(all.df.nonprog)]
graph_nnMat.ref.sum <- Matrix::rowSums(graph_nnMat.ref)

# calculate number of ref and comp cells in each cell's 300 neighbours
graph_nnMat.comp <- graph_nnMat[,rownames(all.df.prog)]
graph_nnMat.comp.sum <- Matrix::rowSums(graph_nnMat.comp)

comp.ref.df <- data.frame(comp = graph_nnMat.comp.sum,
                          ref = graph_nnMat.ref.sum)

comp_FC <- (nrow(all.df.prog)) / (ref_num)

comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)

comp.ref.df[(comp.ref.df$comp<=5 & comp.ref.df$ref<=5),"log2Ratio"] <- 0

comp.ref.df <- comp.ref.df[rownames(all.df),]
all.df$log2Ratio <- comp.ref.df$log2Ratio

all.df[(all.df$log2Ratio > 4),"log2Ratio"] <- 4
all.df[(all.df$log2Ratio < -4),"log2Ratio"] <- -4

ggplot()+
  geom_point(all.df, 
             mapping = aes(x=UMAP_1, y=UMAP_2, 
                           color=log2Ratio), 
             size=1, stroke=0)+
  scale_color_gradientn(colors = c("yellow","#666200","black","#420042","magenta" ),
                        values=c(1, .65, .5, .35, 0),
                        limits=c(-4,4))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave(paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/UMAP_of_enrichment_ratio_between_progression_and_nonprogression.png")
         ,device = "png",
         width = 7, height = 6)


#############
#### DEG between progressor and non-progressor
############

all.df <- FetchData(PBMC.integrated.good, vars = c("sampleBatch","batch","severity",
                                                   "cell.annotation","age"))

# progressor
logi.df <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression.txt", header = TRUE, sep="\t")
logi.df$sampleBatch <- paste(logi.df$batch,logi.df$sample, sep = "_")
logi.df <- logi.df[(logi.df$sampleBatch %in% all.df$sampleBatch),]
logi.df <- logi.df[(logi.df$type == "rising systemic inflammation"),]
dev.df <- all.df[(all.df$sampleBatch %in% logi.df$sampleBatch),]


# non-progressor
mild2.nop <- read.delim("../new_RCA_whole/logitudinal_information_all/Mild1_nonprogression_rising_stage.txt", 
                        sep="\t",header = TRUE)

mild2.nop <- mild2.nop[(mild2.nop$sampleBatch %in% all.df$sampleBatch),]
# remove redundant samples for each patient
mild2.nop <- mild2.nop[!(mild2.nop$sampleBatch %in% c("2_CG1","12_CG30","3_CG47","8_CG55")),]

all.df.mild2 <- all.df[(all.df$sampleBatch %in% mild2.nop$sampleBatch),]


unique_celltype <- unique(all.df$cell.annotation)

## profile cell count for each sample in progression and non-progression 
for (cell.a in  seq_along(unique_celltype)){
  cell <- unique_celltype[cell.a]
  
  
  all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == cell),]
  dev.df.sub <- dev.df[(dev.df$cell.annotation == cell),]
  
  # barplot for non-progression 
  freq.n <- data.frame(table(all.df.mild2.sub$sampleBatch))
  freq.n <- freq.n[order(freq.n$Freq,decreasing = TRUE),]
  freq.n$Var1 <- factor(freq.n$Var1, levels = freq.n$Var1)
  ggplot(freq.n, mapping = aes(x=Var1, y=Freq))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.1*sum(freq.n$Freq),
               color="red")+
    geom_hline(yintercept = 0.2*sum(freq.n$Freq),
               color="red")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0("number of cell = ",sum(freq.n$Freq),"\n",
                   "number of sample = ", nrow(freq.n)))
    ggsave(paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/1.Barplot_of_sample_counts_in_NonProgression_Mild1_",
                  "in_",gsub(" ","_",cell),".png"), device = "png")
  
  # barplot for progression 
  freq.p <- data.frame(table(dev.df.sub$sampleBatch))
  freq.p <- freq.p[order(freq.p$Freq,decreasing = TRUE),]
  freq.p$Var1 <- factor(freq.p$Var1, levels = freq.p$Var1)
  ggplot(freq.p, mapping = aes(x=Var1, y=Freq))+
    geom_bar(stat = "identity")+
    geom_hline(yintercept = 0.1*sum(freq.p$Freq),
               color="red")+
    geom_hline(yintercept = 0.2*sum(freq.p$Freq),
               color="red")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ggtitle(paste0("number of cell = ",sum(freq.p$Freq),"\n",
                   "number of sample = ", nrow(freq.p)))
    ggsave(paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/1.Barplot_of_sample_counts_in_Progression_Mild1_",
                  "in_",gsub(" ","_",cell),".png"), device = "png")
  
}


## downsampling DEG
# Memory B
dev.df.sub <- dev.df[(dev.df$cell.annotation == "Memory B"),]
all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == "Memory B"),]
dev.df.sub.table <- data.frame(table(dev.df.sub$sampleBatch))
sample.down <- as.character(dev.df.sub.table[(dev.df.sub.table$Freq > 60),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  dev.df.sub.i <- dev.df.sub[(dev.df.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(dev.df.sub.i), 60)
  dev.df.sub <- rbind(dev.df.sub[(dev.df.sub$sampleBatch != sample.down.i),],
                      dev.df.sub.i[down.cell.i,])
}

freq.p <- data.frame(table(dev.df.sub$sampleBatch))
freq.p <- freq.p[order(freq.p$Freq,decreasing = TRUE),]
freq.p$Var1 <- factor(freq.p$Var1, levels = freq.p$Var1)
ggplot(freq.p, mapping = aes(x=Var1, y=Freq))+
  geom_bar(stat = "identity")+
  geom_hline(yintercept = 0.1*sum(freq.p$Freq),
             color="red")+
  geom_hline(yintercept = 0.2*sum(freq.p$Freq),
             color="red")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle(paste0("number of cell = ",sum(freq.p$Freq),"\n",
                 "number of sample = ", nrow(freq.p)))
ggsave("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/1.Barplot_of_sample_counts_in_Progression_Mild1_Memory_B_after_downsampling.png", device = "png")


all.df.mild2.sub.table <- data.frame(table(all.df.mild2.sub$sampleBatch))
sample.down <- as.character(all.df.mild2.sub.table[(all.df.mild2.sub.table$Freq > 50),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  all.df.mild2.sub.i <- all.df.mild2.sub[(all.df.mild2.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(all.df.mild2.sub.i),50)
  all.df.mild2.sub <- rbind(all.df.mild2.sub[(all.df.mild2.sub$sampleBatch != sample.down.i),],
                            all.df.mild2.sub[down.cell.i,])
}

freq.p <- data.frame(table(all.df.mild2.sub$sampleBatch))
freq.p <- freq.p[order(freq.p$Freq,decreasing = TRUE),]
freq.p$Var1 <- factor(freq.p$Var1, levels = freq.p$Var1)
ggplot(freq.p, mapping = aes(x=Var1, y=Freq))+
  geom_bar(stat = "identity")+
  geom_hline(yintercept = 0.1*sum(freq.p$Freq),
             color="red")+
  geom_hline(yintercept = 0.2*sum(freq.p$Freq),
             color="red")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle(paste0("number of cell = ",sum(freq.p$Freq),"\n",
                 "number of sample = ", nrow(freq.p)))
ggsave("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/1.Barplot_of_sample_counts_in_NonProgression_Mild1_Memory_B_after_downsampling.png", device = "png")


de.ij <- RCAv2:::ComputePairWiseDE(object = PBMC.integrated.good@assays$RNA@data, 
                                   cells.1 = rownames(dev.df.sub),
                                   cells.2 = rownames(all.df.mild2.sub),
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
write.table(de.ij, paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//14.Single_cell_DEG_between_Mild1_progression_and_nonprogression_in_Memory_B.txt"), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)


# Naive B
dev.df.sub <- dev.df[(dev.df$cell.annotation == "Naive B"),]
all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == "Naive B"),]
dev.df.sub.table <- data.frame(table(dev.df.sub$sampleBatch))
sample.down <- as.character(dev.df.sub.table[(dev.df.sub.table$Freq > 100),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  dev.df.sub.i <- dev.df.sub[(dev.df.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(dev.df.sub.i), 100)
  dev.df.sub <- rbind(dev.df.sub[(dev.df.sub$sampleBatch != sample.down.i),],
                      dev.df.sub.i[down.cell.i,])
}

freq.p <- data.frame(table(dev.df.sub$sampleBatch))
freq.p <- freq.p[order(freq.p$Freq,decreasing = TRUE),]
freq.p$Var1 <- factor(freq.p$Var1, levels = freq.p$Var1)
ggplot(freq.p, mapping = aes(x=Var1, y=Freq))+
  geom_bar(stat = "identity")+
  geom_hline(yintercept = 0.1*sum(freq.p$Freq),
             color="red")+
  geom_hline(yintercept = 0.2*sum(freq.p$Freq),
             color="red")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle(paste0("number of cell = ",sum(freq.p$Freq),"\n",
                 "number of sample = ", nrow(freq.p)))
ggsave("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/1.Barplot_of_sample_counts_in_Progression_Mild1_Naive_B_after_downsampling.png", device = "png")


all.df.mild2.sub.table <- data.frame(table(all.df.mild2.sub$sampleBatch))
sample.down <- as.character(all.df.mild2.sub.table[(all.df.mild2.sub.table$Freq > 200),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  all.df.mild2.sub.i <- all.df.mild2.sub[(all.df.mild2.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(all.df.mild2.sub.i),200)
  all.df.mild2.sub <- rbind(all.df.mild2.sub[(all.df.mild2.sub$sampleBatch != sample.down.i),],
                            all.df.mild2.sub[down.cell.i,])
}

freq.p <- data.frame(table(all.df.mild2.sub$sampleBatch))
freq.p <- freq.p[order(freq.p$Freq,decreasing = TRUE),]
freq.p$Var1 <- factor(freq.p$Var1, levels = freq.p$Var1)
ggplot(freq.p, mapping = aes(x=Var1, y=Freq))+
  geom_bar(stat = "identity")+
  geom_hline(yintercept = 0.1*sum(freq.p$Freq),
             color="red")+
  geom_hline(yintercept = 0.2*sum(freq.p$Freq),
             color="red")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ggtitle(paste0("number of cell = ",sum(freq.p$Freq),"\n",
                 "number of sample = ", nrow(freq.p)))
ggsave("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/1.Barplot_of_sample_counts_in_NonProgression_Mild1_Naive_B_after_downsampling.png", device = "png")


de.ij <- RCAv2:::ComputePairWiseDE(object = PBMC.integrated.good@assays$RNA@data, 
                                   cells.1 = rownames(dev.df.sub),
                                   cells.2 = rownames(all.df.mild2.sub),
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
write.table(de.ij, paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//14.Single_cell_DEG_between_Mild1_progression_and_nonprogression_in_Naive_B.txt"), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)





## volcano plot
cell <- c("Memory B","Naive B")
for (i in seq_along(cell)){
  cell.i <- cell[i]
  
  deg.file <- read.table(paste0("./14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//14.Single_cell_DEG_between_Mild1_progression_and_nonprogression_in_",
                                gsub(" ","_", cell.i),".txt"),
                         sep="\t", header = TRUE)
  
  
  deg.file <- deg.file[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                            x = deg.file$gene, invert = TRUE),]
  
  deg.file[(deg.file$p_val_adj==0), "p_val_adj"] <- 1e-300
  deg.file$logP <- -log10(deg.file$p_val_adj)
  de.ij.sig <- deg.file[(deg.file$p_val_adj < 0.1 & abs(deg.file$avg_logFC) > 0.25 ),]
  
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC > 0),"gene"],
              paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//upregulated_genes_in_",
                     gsub(" ","_",cell.i),"_mild1_progression_vs_non_progression.txt"),
              col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC < 0),"gene"],
              paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//downregulated_genes_in_",
                     gsub(" ","_",cell.i),"_mild1_progression_vs_non_progression.txt"),
              col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  
  
  ggplot()+
    geom_point(deg.file, mapping = aes(x = avg_logFC, y = logP), color="grey")+
    geom_point(de.ij.sig[which(de.ij.sig$avg_logFC>0),],
               mapping = aes(x = avg_logFC, y = logP), color="red")+
    geom_point(de.ij.sig[which(de.ij.sig$avg_logFC<0),],
               mapping = aes(x = avg_logFC, y = logP), color="blue")+
    geom_label_repel(data=de.ij.sig[which(de.ij.sig$avg_logFC> 0.3),], 
                     aes(x = avg_logFC, y = logP,label=gene),
                     box.padding   = .4, 
                     point.padding = 0.3,
                     segment.color = 'red',
                     size=8, color="red",max.overlaps =100)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
    ggsave(paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//volcano_plot_in_",
                  gsub(" ","_",cell.i),"_mild1_progression_vs_non_progression.pdf"),
           width = 11, height = 6)
  
}


################# gene heatmap for upregulated IFN genes in progression compared to non pregression

# upregulated genes
gene.c <- c()
cell <- c("Memory B","Naive B")
for (i in seq_along(cell)){
  cell.i <- cell[i]
  
  gene.i <- read.table(paste0("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage//upregulated_genes_in_",
                              gsub(" ","_",cell.i),"_mild1_progression_vs_non_progression.txt"), header = FALSE,stringsAsFactors = FALSE)
  gene.c <- c(gene.c, gene.i$V1)
}

gene.c <- unique(gene.c)
write.table(gene.c,"14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/union_of_upregulated_genes_across_different_cells.txt",
            col.names = FALSE, row.names = FALSE,quote = FALSE)
go <- read.delim("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/GO_of_union_of_upregulated_genes_across_different_cells.txt",
                 sep = "\t", header = TRUE)
go <- go[c(1,2),]
gene.1 <- unlist(strsplit(go$Genes[1], split = ";"))
gene.2 <- unlist(strsplit(go$Genes[2], split = ";"))
gene.c <- unique(c(gene.1, gene.2))

# start form gene expression matrix
all.df <- FetchData(PBMC.integrated.good, vars = c("sampleBatch","batch","severity",
                                                   "cell.annotation","age"))

# progressor
logi.df <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression.txt", 
                      header = TRUE, sep="\t")
logi.df$sampleBatch <- paste(logi.df$batch,logi.df$sample, sep = "_")
logi.df <- logi.df[(logi.df$type == "rising systemic inflammation"),]

dev.df <- all.df[(all.df$sampleBatch %in% logi.df$sampleBatch),]

mild2.pro.sample <- unique(dev.df$sampleBatch)

mild2.pro.gene <- data.frame(matrix(nrow = length(gene.c),
                                    ncol = length(mild2.pro.sample)))
for (a in seq_along(mild2.pro.sample)){
  mild2.pro.sample.a <- mild2.pro.sample[a]
  dev.df.a <- dev.df[(dev.df$sampleBatch == mild2.pro.sample.a),]
  gene.a <- PBMC.integrated.good@assays$RNA@data[gene.c, rownames(dev.df.a)]
  gene.a.mean <- Matrix::rowMeans(gene.a)
  gene.a.mean <- gene.a.mean[gene.c]
  mild2.pro.gene[,a] <- gene.a.mean
}
colnames(mild2.pro.gene) <- paste("progression",seq(1,length(mild2.pro.sample),1),sep = "" )
rownames(mild2.pro.gene) <- gene.c


# non-progressor
mild2.nop2 <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_nonprogression_rising_stage.txt", 
                         header = TRUE, sep="\t")
# remove redundant samples for each patient
mild2.nop2 <- mild2.nop2[!(mild2.nop2$sampleBatch %in% c("2_CG1","12_CG30","3_CG47","8_CG55")),]

all.df.mild2.2 <- all.df[(all.df$sampleBatch %in% mild2.nop2$sampleBatch),]

all.df.mild2.2.table <- as.data.frame(table(all.df.mild2.2$sampleBatch))

mild2.nonpro.sample2 <- as.character(all.df.mild2.2.table[(all.df.mild2.2.table$Freq > 10),"Var1"])

mild2.nonpro.gene2 <- data.frame(matrix(nrow = length(gene.c),
                                        ncol = length(mild2.nonpro.sample2)))
for (a in seq_along(mild2.nonpro.sample2)){
  mild2.nonpro.sample2.a <- mild2.nonpro.sample2[a]
  all.df.mild2.a <- all.df.mild2.2[(all.df.mild2.2$sampleBatch == mild2.nonpro.sample2.a),]
  gene.a <- PBMC.integrated.good@assays$RNA@data[gene.c, rownames(all.df.mild2.a)]
  gene.a.mean <- Matrix::rowMeans(gene.a)
  gene.a.mean <- gene.a.mean[gene.c]
  mild2.nonpro.gene2[,a] <- gene.a.mean
}
colnames(mild2.nonpro.gene2) <- paste("nonprogression_rising_",seq(1,length(mild2.nonpro.sample2),1),sep = "" )
rownames(mild2.nonpro.gene2) <- gene.c

mild2.gene.all <- cbind(mild2.pro.gene, mild2.nonpro.gene2)

mild2.gene.all.zcore =  t(scale(t(mild2.gene.all), center = T, scale = T))

bk = unique(c(seq(-3, 3, length=80)))
col = colorRampPalette(c("magenta","black","yellow"))(length(bk)-1)

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="euclidean")
pdf("14.sample_counts_subanalysis_for_DEG_analysis_mild1_progression_and_nonprogression_single_cell_rising_stage/heatmap_of_upregulated_genes_in_mild1_progression_and_non_progression_different_timepoints.pdf")
heatmap.2(as.matrix(mild2.gene.all.zcore), Rowv = T, cexRow = .3,cexCol  = .3,
          Colv = F, col = col, breaks=bk,symm=F, symbreaks=F, scale="none",
          trace="none",density.info = "none",hclustfun=hclustfunc, distfun=distfunc)
dev.off()



