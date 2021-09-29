library(RCAv2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)

# meta data
metadata <- read.table("../metadata_for_integrated_data_after_removing_bad_cells_and_QC.txt",
                       sep = "\t", header = TRUE)
rownames(metadata) <- metadata$barcode

# NK T (removed bad clusters)
nk.t <- readRDS("../NK_T_subclustering/PBMC.integrated.T.NK.good.rds")
metadata.nk <- metadata[which(metadata$cell %in% c("NK cells","T cells")),]
metadata.nk.good <- metadata.nk[colnames(nk.t),]
metadata.nk.good$cell.detailed.annotation <- nk.t$cell.annotation

# monocytes (removed bad clusters)
mono <- readRDS("../Monocytes_mDC/PBMC.integrated.Mono.good.rds")
metadata.mono <- metadata[which(metadata$cell %in% c("Monocytes","mDC")),]
metadata.mono.good <- metadata.mono[colnames(mono),]
metadata.mono.good$cell.detailed.annotation <- mono$cell.annotation

# B (removed bad clusters)
B <- readRDS("../B_subset/PBMC.integrated.B.good.rds")
metadata.b <- metadata[which(metadata$cell %in% c("B cells")),]
metadata.b.good <- metadata.b[colnames(B),]
metadata.b.good$cell.detailed.annotation <- B$cell.annotation

# platelets
metadata.p <-  metadata[which(metadata$cell %in% c("Platelets")),]
metadata.p$cell.detailed.annotation <- "Platelets"

# Plasma B
metadata.pb <-  metadata[which(metadata$cell %in% c("Plasma B")),]
metadata.pb$cell.detailed.annotation <- "Plasma B"

#combined
metadata.new <- rbind(metadata.nk.good,metadata.mono.good,
                      metadata.b.good,metadata.p,metadata.pb)

RCA.old <- readRDS("../RCA.combined.all.after.removing.bad.cells.again.and.QC.rds")

good.barcode <- metadata.new$barcode
good.barcode.index <- which(colnames(RCA.old$raw.data) %in% good.barcode)
RCA.old$raw.data <- RCA.old$raw.data[, good.barcode.index]
RCA.old$data <- RCA.old$data[, good.barcode.index]
RCA.old$projection.data <- RCA.old$projection.data[, good.barcode.index]
RCA.old$cell.Type.Estimate <- RCA.old$cell.Type.Estimate[good.barcode.index]
RCA.old$clustering.out$dynamicColorsList[[1]] <- RCA.old$clustering.out$dynamicColorsList[[1]][good.barcode.index]

RCA.old <- computeUMAP(RCA.old)

saveRDS(RCA.old,"RCA.combined.all.after.subclustering.rds")

RCA.old <- readRDS("RCA.combined.all.after.subclustering.rds")

umap.df <- RCA.old$umap.coordinates
umap.df$barcode <- rownames(umap.df)
metadata.new.sub <- metadata.new[,c("sampleBatch","barcode","cell","TCR","BCR",
                                    "sample","lib","batch","severity","age","sex",
                                    "ethnicity","cell.detailed.annotation")]
umap.df <- merge(x = umap.df, y = metadata.new.sub, by = "barcode")

umap.df$cell.general.annotation <- umap.df$cell.detailed.annotation
umap.df$cell.general.annotation <- as.character(umap.df$cell.general.annotation)
umap.df[(umap.df$cell.detailed.annotation %in% c("CD16bright NK","CD16dim NK")),
        "cell.general.annotation"] <- "NK cell"
umap.df[(umap.df$cell.detailed.annotation %in% c("CD4 Tcm","CD4 Tem","Cytotoxic CD4 T",
                                                 "Treg","Naive CD4 T")),
        "cell.general.annotation"] <- "T cell"
umap.df[(umap.df$cell.detailed.annotation %in% c("CD8 Tem","Cytotoxic CD8 T",
                                                 "Naive CD8 T")),
        "cell.general.annotation"] <- "T cell"
umap.df[(umap.df$cell.detailed.annotation %in% c("CD14 Mono","IM Mono","CD16 Mono")),
        "cell.general.annotation"] <- "Monocyte"
umap.df[(umap.df$cell.detailed.annotation %in% c("Naive B","Memory B")),
        "cell.general.annotation"] <- "B cell"

write.table(umap.df,"metadata.txt",
            sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


pdf("1.integrated_UMAP_with_cell_annotation.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cell.general.annotation)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=c("lightgreen","brown","darkgreen",
                              "darkorange","magenta","red","steelblue","turquoise",
                              "black","darkmagenta","pink")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()


pdf("2.integrated_UMAP_with_TCR.pdf",width=7,height=7)
ggplot() + 
  geom_point(data = umap.df[which(umap.df$TCR==0),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$TCR==1),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="red") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_blank())
dev.off()


######### BCR
pdf("2.integrated_UMAP_with_BCR.pdf",width=7,height=7)
ggplot() + 
  geom_point(data = umap.df[which(umap.df$BCR==0),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$BCR==1),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="blue") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_blank())
dev.off()


##### plot QC
rawdata <- RCA.old$raw.data
# NODG
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
##### QC on each cluster
QC_list <- list()
unique_celltype <- unique(umap.df$cell.general.annotation)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- umap.df[which(umap.df$cell.general.annotation==cell.i),"barcode"]
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.2, stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  QC_list[[i]] <- p.i
}

p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))
pdf("2.QC_plot_of_each_cell_types.pdf", 
    width = 3, height = 3*length(unique(umap.df$cell.general.annotation)))
grid.arrange(p.all,ncol=1)
dev.off()


##### plot severity
umap.severity <- umap.df[, c("UMAP1","UMAP2","severity")]
umap.severity <- umap.severity[(umap.severity$severity != "Undefined"),]
umap.severity$severity <- factor(umap.severity$severity,
                                 levels = c("Asymptomatic","Mild1",
                                            "Mild2","Moderate","Severe","Critical"))


ggplot() + 
  geom_point(data = umap.df, 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0, color="grey")+
  geom_point(data = umap.severity, 
             mapping = aes(x = UMAP1, y = UMAP2, color=severity),
             size = .1, stroke=0) + 
  scale_color_manual(values=c("#E1BC00","#A5D700","#00997F",
                              "#1313A9","#B20086","#F20000")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("4.integrated_UMAP_with_clinical_severity.pdf",
         width=9,height=7)

##########
########## differential abundance 
##########
metadata <- read.table("metadata.txt",
                       sep = "\t", header = TRUE)


abundance.df <- metadata[,c("cell.general.annotation",
                            "sampleBatch",
                            "severity")]

abundance.df$count <- 1

abundance.df.final <- aggregate(count ~ ., abundance.df, FUN = sum)

unique_sampleBatch <- unique(abundance.df.final$sampleBatch)
unique_celltype <- unique(abundance.df.final$cell.general.annotation)

abundance.df.final_new <- data.frame()
# calculate the proportion for each sample
for (i in seq_along(unique_sampleBatch)){
  sampleBatch.i <- unique_sampleBatch[i]
  abundance.df.final.i <- abundance.df.final[which(abundance.df.final$sampleBatch==sampleBatch.i),]
  
  # add row without corresponding cell types
  unique_celltype.no <- unique_celltype[!(unique_celltype %in% abundance.df.final.i$cell.general.annotation)]
  if(length(unique_celltype.no) > 0){
    no.df <- data.frame(cell.general.annotation = unique_celltype.no,
                        sampleBatch = rep(sampleBatch.i, length(unique_celltype.no)),
                        severity= rep(abundance.df.final.i$severity[1], 
                                      length(unique_celltype.no)),
                        count = rep(0,length(unique_celltype.no)))
    abundance.df.final.i <- rbind(abundance.df.final.i, no.df)
  }
  
  
  
  sum.i <- sum(abundance.df.final.i$count)
  abundance.df.final.i$proportion <- 100*abundance.df.final.i$count/sum.i
  abundance.df.final_new <- rbind(abundance.df.final_new,abundance.df.final.i)
}

infect.stat <- read.delim("../infection_status_for_all_sampleBatches.txt",
                          sep = "\t", header = TRUE)
infect.stat.sub <- infect.stat[,c("sampleBatch","cat.x","cat")]

abundance.df.final_new.cat <- merge(abundance.df.final_new, infect.stat.sub,
                                    by = "sampleBatch")
# remove undefined severity
abundance.df.final_new.cat <- abundance.df.final_new.cat[(abundance.df.final_new.cat$severity != "Undefined"),]

## day 1-8
abundance.df.final_new.cat1 <- abundance.df.final_new.cat[(abundance.df.final_new.cat$cat.x %in% c(1,2,3)),]

unique_celltype <- c("B cell","mDC","Monocyte",
                     "NK cell","Plasma B","Platelets","T cell")
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new1 <- abundance.df.final_new.cat1[which(abundance.df.final_new.cat1$cell.general.annotation==unique_celltype.i),]
  
  abundance.df.final_new1$severity <- factor(abundance.df.final_new1$severity,
                                             levels = c("Asymptomatic","Mild1","Mild2",
                                                        "Moderate","Severe","Critical"))
  
  
  stat.test <- wilcox_test(data=abundance.df.final_new1,
                           proportion ~ severity,
                           p.adjust.method = "BH")
  
  stat.test <- stat.test %>% add_xy_position(x = "severity")
  
  stat.test <- stat.test[(stat.test$p.adj<0.05),]

  
  kruskal.i <- kruskal.test(proportion~severity, data=abundance.df.final_new1)
  kruskal.i.res <- kruskal.i$p.value
  
  ggboxplot(abundance.df.final_new1, 
            x = "severity", y = "proportion", fill = "severity",
            outlier.size = 0, alpha=0.7,
            outlier.shape = NA) + 
    geom_jitter(size=1.5, stroke=0, width = .2 )+
    stat_pvalue_manual(stat.test, label = "p.adj")+
    ylab("proportion (%)")+
    xlab("")+
    scale_fill_manual(values = c("#E1BC00","#A5D700","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ggtitle(paste0("Kruskal p=",kruskal.i.res))
    ggsave(paste0("6.differential_abundance_across_samples_in_",
                  gsub(" ", "_", unique_celltype.i),"_early_and_rising_adjusted_pvalue.pdf"),width = 6, height = 4)
  
}


## > day 8
abundance.df.final_new.cat2 <- abundance.df.final_new.cat[(abundance.df.final_new.cat$cat.x %in% c(1,4,5)),]


unique_celltype <- c("B cell","mDC","Monocyte",
                     "NK cell","Plasma B","Platelets","T cell")

## adjusted p-value heatmap
severity.c <- c("Asymptomatic","Mild1","Mild2",
                "Moderate","Severe","Critical")
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new1 <- abundance.df.final_new.cat2[which(abundance.df.final_new.cat2$cell.general.annotation==unique_celltype.i),]
  
  abundance.df.final_new1$severity <- factor(abundance.df.final_new1$severity,
                                             levels = severity.c)
  
  stat.test <- wilcox_test(data=abundance.df.final_new1,
                           proportion ~ severity,
                           p.adjust.method = "BH")
  
  stat.test$log10padj <- as.numeric(format(round(-log10(stat.test$p.adj), 3), nsmall = 3))
  
  heatmap.matrx <- as.data.frame(matrix(ncol=6,
                                        nrow=6))
  rownames(heatmap.matrx) <- severity.c
  colnames(heatmap.matrx) <- severity.c
  for (a in 1:5){
    severity.a <- severity.c[a]
    for (b in (a+1):6){
      severity.b <- severity.c[b]
      heatmap.matrx[severity.a,severity.b] <- stat.test[(stat.test$group1 == severity.a &
                                                           stat.test$group2 == severity.b), "log10padj"]
    }
  }
  
  bk = unique(c(seq(0, 1.5, length=5), seq(1.5, 4, length=20)))
  col = colorRampPalette(c("white","#B30106"))(length(bk)-1)
  
  pdf(paste0("7.differential_abundance_across_samples_in_",
            gsub(" ", "_", unique_celltype.i),"_peak_and_late_adjusted_pvalue_heatmap.pdf"))
  heatmap.2(as.matrix(heatmap.matrx), Rowv = F, Colv = F, col = col,
            breaks=bk,symm=F, symbreaks=F, scale="none" ,
            trace="none",density.info = "none",
            cellnote=as.matrix(heatmap.matrx),
            notecex=1.2,
            notecol="black")
  dev.off()
}


for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new1 <- abundance.df.final_new.cat2[which(abundance.df.final_new.cat2$cell.general.annotation==unique_celltype.i),]
  
  abundance.df.final_new1$severity <- factor(abundance.df.final_new1$severity,
                                             levels = c("Asymptomatic","Mild1","Mild2",
                                                        "Moderate","Severe","Critical"))
  
  kruskal.i <- kruskal.test(proportion~severity, data=abundance.df.final_new1)
  kruskal.i.res <- kruskal.i$p.value
  
  
  ggplot(abundance.df.final_new1, aes(x=severity, y=proportion,fill=severity)) + 
    geom_boxplot(outlier.size = 0, alpha=0.7,
                 outlier.shape = NA) + geom_jitter(size=.8, stroke=0, width = .2 )+
    theme(axis.text.x  = element_text(angle=90, vjust=0.5))+ 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    ylab("proportion (%)")+
    xlab("")+
    scale_fill_manual(values = c("#E1BC00","#A5D700","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ggtitle(paste0("P-value = ",kruskal.i.res))
    ggsave(paste0("7.differential_abundance_across_samples_in_",
                  gsub(" ", "_", unique_celltype.i),"_peak_and_late.pdf"),width = 6, height = 4)
  
}

