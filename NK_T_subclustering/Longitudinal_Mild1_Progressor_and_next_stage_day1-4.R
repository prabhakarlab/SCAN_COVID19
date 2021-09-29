library(Seurat)
library(ggplot2)
library(ggpubr)
library(gplots)
library(gridExtra)
library(RCAv2)
library(dbplyr)
library(ggrepel)

PBMC.integrated.good <- readRDS("PBMC.integrated.T.NK.good.rds")

all.df <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","sampleBatch","batch","severity",
                                                   "cell.annotation","age"))

# mild1 progressor
logi.df <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression.txt", header = TRUE, sep="\t")
logi.df$sampleBatch <- paste(logi.df$batch,logi.df$sample, sep = "_")
logi.df <- logi.df[(logi.df$sampleBatch %in% all.df$sampleBatch),]
# get early infection
logi.df <- logi.df[(logi.df$type == "early infection"),]
# get only one sample for C-T014
logi.df <- logi.df[1:4,]

dev.df <- all.df[(all.df$sampleBatch %in% logi.df$sampleBatch),]

## next stage of mild1
mild2.next <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression_next_stage.txt",
                         header = TRUE, sep = "\t")
mild2.next$sampleBatch <- paste(mild2.next$batch, mild2.next$sample, sep = "_")
mild2.next <- mild2.next[(mild2.next$sampleBatch %in% all.df$sampleBatch),]

mild2.next <- mild2.next[(mild2.next$sample %in% logi.df$sample),]
# keep one sample per patient
mild2.next <- mild2.next[!(mild2.next$sampleBatch %in% c("13_CG51","6_CG60","18_CG75","14_CG79","19_CG43")),]

next.df <- all.df[(all.df$sampleBatch %in% mild2.next$sampleBatch),]

ggplot()+
  geom_point(all.df, mapping = aes(x=UMAP_1, y=UMAP_2), color="grey", size=.5, stroke=0)+
  geom_point(dev.df, mapping = aes(x=UMAP_1, y=UMAP_2), color="#00AE46", size=.5, stroke=0)+
  geom_point(next.df[(next.df$severity == "Mild2"),],
             mapping = aes(x=UMAP_1, y=UMAP_2, color=severity), 
             size=.5, stroke=0, color = "#00997F")+
  geom_point(next.df[(next.df$severity == "Moderate"),],
             mapping = aes(x=UMAP_1, y=UMAP_2, color=severity), 
             size=.5, stroke=0, color = "#1313A9")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave("14.Next_stage_of_mild1_progression_early/UMAP_of_Mild1_progression_and_Next_stage.png",
         device = "png",width = 5, height = 4)



##########
##### 
##### 
##### next stage of mild1 progression
#####
##########
##### 
all.df <- FetchData(PBMC.integrated.good, vars = c("sampleBatch","batch","severity",
                                                   "cell.annotation","age"))
all.df.mild2 <- next.df

unique_celltype <- unique(all.df$cell.annotation)

## profile cell count for each sample to ensure no sample overrepresented
for (cell.a in  seq_along(unique_celltype)){
  cell <- unique_celltype[cell.a]
  
  
  all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == cell),]
  dev.df.sub <- dev.df[(dev.df$cell.annotation == cell),]
  
  # barplot for non-progression mild2
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
    ggsave(paste0("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_NEXT_",
                  "in_",gsub(" ","_",cell),".png"), device = "png")
  
  # barplot for progression mild2
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
    ggsave(paste0("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_",
                  "in_",gsub(" ","_",cell),".png"), device = "png")
  
}


## downsampling DEG
# CD4 Tcm
dev.df.sub <- dev.df[(dev.df$cell.annotation == "CD4 Tcm"),]
all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == "CD4 Tcm"),]

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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_CD4_Tcm_after_downsampling.png", device = "png")


all.df.mild2.sub.table <- data.frame(table(all.df.mild2.sub$sampleBatch))
sample.down <- as.character(all.df.mild2.sub.table[(all.df.mild2.sub.table$Freq > 150),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  all.df.mild2.sub.i <- all.df.mild2.sub[(all.df.mild2.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(all.df.mild2.sub.i), 150)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_NEXT_CD4_Tcm_after_downsampling.png", device = "png")


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
write.table(de.ij, paste0("14.Next_stage_of_mild1_progression_early/14.Single_cell_DEG_between_Mild1_progression_and_NEXT_in_CD4_Tcm.txt"), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)




# CD4 Tem
dev.df.sub <- dev.df[(dev.df$cell.annotation == "CD4 Tem"),]
all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == "CD4 Tem"),]

dev.df.sub.table <- data.frame(table(dev.df.sub$sampleBatch))
sample.down <- as.character(dev.df.sub.table[(dev.df.sub.table$Freq > 150),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  dev.df.sub.i <- dev.df.sub[(dev.df.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(dev.df.sub.i), 150)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_CD4_Tem_after_downsampling.png", device = "png")


all.df.mild2.sub.table <- data.frame(table(all.df.mild2.sub$sampleBatch))
sample.down <- as.character(all.df.mild2.sub.table[(all.df.mild2.sub.table$Freq > 200),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  all.df.mild2.sub.i <- all.df.mild2.sub[(all.df.mild2.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(all.df.mild2.sub.i), 200)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_NEXT_CD4_Tem_after_downsampling.png", device = "png")


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
write.table(de.ij, paste0("14.Next_stage_of_mild1_progression_early/14.Single_cell_DEG_between_Mild1_progression_and_NEXT_in_CD4_Tem.txt"), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)






# CD16bright NK
dev.df.sub <- dev.df[(dev.df$cell.annotation == "CD16bright NK"),]
all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == "CD16bright NK"),]
dev.df.sub.table <- data.frame(table(dev.df.sub$sampleBatch))

sample.down <- as.character(dev.df.sub.table[(dev.df.sub.table$Freq > 600),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  dev.df.sub.i <- dev.df.sub[(dev.df.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(dev.df.sub.i), 150)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_CD16bright_NK_after_downsampling.png", device = "png")


all.df.mild2.sub.table <- data.frame(table(all.df.mild2.sub$sampleBatch))
sample.down <- as.character(all.df.mild2.sub.table[(all.df.mild2.sub.table$Freq > 300),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  all.df.mild2.sub.i <- all.df.mild2.sub[(all.df.mild2.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(all.df.mild2.sub.i), 300)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_NEXT_CD16bright_NK_after_downsampling.png", device = "png")


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
write.table(de.ij, paste0("14.Next_stage_of_mild1_progression_early/14.Single_cell_DEG_between_Mild1_progression_and_NEXT_in_CD16bright_NK.txt"), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)







# Cytotoxic CD8 T
dev.df.sub <- dev.df[(dev.df$cell.annotation == "Cytotoxic CD8 T"),]
all.df.mild2.sub <- all.df.mild2[(all.df.mild2$cell.annotation == "Cytotoxic CD8 T"),]

dev.df.sub.table <- data.frame(table(dev.df.sub$sampleBatch))
sample.down <- as.character(dev.df.sub.table[(dev.df.sub.table$Freq > 300),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  dev.df.sub.i <- dev.df.sub[(dev.df.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(dev.df.sub.i), 300)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_Cytotoxic_CD8T_after_downsampling.png", device = "png")



all.df.mild2.sub.table <- data.frame(table(all.df.mild2.sub$sampleBatch))
sample.down <- as.character(all.df.mild2.sub.table[(all.df.mild2.sub.table$Freq > 600),"Var1"])
for (i in seq_along(sample.down)){
  sample.down.i <- sample.down[i]
  all.df.mild2.sub.i <- all.df.mild2.sub[(all.df.mild2.sub$sampleBatch == sample.down.i),]
  down.cell.i <- sample(rownames(all.df.mild2.sub.i), 500)
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
ggsave("14.Next_stage_of_mild1_progression_early/1.Barplot_of_sample_counts_in_Progression_Mild1_NEXT_Cytotoxic_CD8_T_after_downsampling.png", device = "png")


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
write.table(de.ij, paste0("14.Next_stage_of_mild1_progression_early/14.Single_cell_DEG_between_Mild1_progression_and_NEXT_in_Cytotoxic_CD8_T.txt"), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)


## volcano plot
cell <- c("CD4 Tcm","CD4 Tem","CD16bright NK","Cytotoxic CD8 T")
for (i in seq_along(cell)){
  cell.i <- cell[i]
  
  deg.file <- read.table(paste0("./14.Next_stage_of_mild1_progression_early/14.Single_cell_DEG_between_Mild1_progression_and_NEXT_in_",
                                gsub(" ","_", cell.i),".txt"),
                         sep="\t", header = TRUE)
  
  
  deg.file <- deg.file[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                            x = deg.file$gene, invert = TRUE),]
  
  deg.file[(deg.file$p_val_adj==0), "p_val_adj"] <- 1e-300
  deg.file$logP <- -log10(deg.file$p_val_adj)
  de.ij.sig <- deg.file[(deg.file$p_val_adj < 0.1 & abs(deg.file$avg_logFC) > 0.25 ),]
  
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC > 0),"gene"],
              paste0("14.Next_stage_of_mild1_progression_early/upregulated_genes_in_",
                     gsub(" ","_",cell.i),"_mild1_progression_vs_NEXT.txt"),
              col.names = FALSE, quote = FALSE, row.names = FALSE)
  write.table(de.ij.sig[which(de.ij.sig$avg_logFC < 0),"gene"],
              paste0("14.Next_stage_of_mild1_progression_early/downregulated_genes_in_",
                     gsub(" ","_",cell.i),"_mild1_progression_vs_NEXT.txt"),
              col.names = FALSE, quote = FALSE, row.names = FALSE)
  
  
  
  ggplot()+
    geom_point(deg.file, mapping = aes(x = avg_logFC, y = logP), color="grey")+
    geom_point(de.ij.sig[which(de.ij.sig$avg_logFC>0),],
               mapping = aes(x = avg_logFC, y = logP), color="red")+
    geom_point(de.ij.sig[which(de.ij.sig$avg_logFC<0),],
               mapping = aes(x = avg_logFC, y = logP), color="blue")+
    geom_label_repel(data=de.ij.sig[which(de.ij.sig$avg_logFC> 0.6),], 
                     aes(x = avg_logFC, y = logP,label=gene),
                     box.padding   = .4, 
                     point.padding = 0.3,
                     segment.color = 'red',
                     size=8, color="red",max.overlaps =100)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
    ggsave(paste0("14.Next_stage_of_mild1_progression_early/volcano_plot_in_",
                  gsub(" ","_",cell.i),"_mild1_progression_vs_NEXT.pdf"),
           width = 11, height = 6)
  
}





############
########### enrichment plot between progression and next
###########
graph_nnMat <- readRDS("12.2D_comparative_enrichment_plot/graph_nnMat.rds")

all.df <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","sampleBatch","batch","severity",
                                                   "cell.annotation","age"))

### progression
prog <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression.txt",
                   header = TRUE, sep="\t")
prog$sampleBatch <- paste(prog$batch,prog$sample, sep="_")

prog <- prog[(prog$sampleBatch %in% all.df$sampleBatch),]
# get early infection
prog <- prog[(prog$type == "early infection"),]
# get only one sample for C-T014
prog <- prog[1:4,]

all.df.prog <- all.df[(all.df$sampleBatch %in% prog$sampleBatch),]

### next stage
nonprog <- read.table("../new_RCA_whole/logitudinal_information_all/Mild1_progression_next_stage.txt",
                      header = TRUE, sep="\t")

nonprog$sampleBatch <- paste(nonprog$batch, nonprog$sample, sep = "_")
nonprog <- nonprog[(nonprog$sample %in% prog$sample),]

nonprog <- nonprog[(nonprog$sampleBatch %in% PBMC.integrated.good$sampleBatch),]
# keep one sample per patient
nonprog <- nonprog[!(nonprog$sampleBatch %in% c("13_CG51","6_CG60","18_CG75","14_CG79","19_CG43")),]

all.df.nonprog <- all.df[(all.df$sampleBatch %in% nonprog$sampleBatch),]

### plot prog as ref
ref_num <- nrow(all.df.prog)
# calculate total number of ref cells in each cell's 300 neighbours
graph_nnMat.ref <- graph_nnMat[,rownames(all.df.prog)]
graph_nnMat.ref.sum <- Matrix::rowSums(graph_nnMat.ref)

# calculate number of ref and comp cells in each cell's 300 neighbours
graph_nnMat.comp <- graph_nnMat[,rownames(all.df.nonprog)]
graph_nnMat.comp.sum <- Matrix::rowSums(graph_nnMat.comp)

comp.ref.df <- data.frame(comp = graph_nnMat.comp.sum,
                          ref = graph_nnMat.ref.sum)

comp_FC <- (nrow(all.df.nonprog)) / (ref_num)

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
             size=0.5, stroke=0)+
  scale_color_gradientn(colors = c("yellow","#666200","black","#420042","magenta" ),
                        values=c(1, .65, .5, .35, 0),
                        limits=c(-4,4))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave(paste0("14.Next_stage_of_mild1_progression_early/UMAP_of_enrichment_ratio_between_progression_and_NEXT.png")
         ,device = "png",
         width = 7, height = 6)



