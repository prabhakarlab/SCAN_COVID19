library(Seurat)
library(ggplot2)
library(ggpubr)
library(gplots)
library(Matrix)
library(gridExtra)
library(rstatix)

metadata <- read.table("../metadata_for_integrated_data_after_removing_bad_cells_and_QC.txt",
                       sep = "\t", header = TRUE)

metadata.b <- metadata[(metadata$cell %in% c("B cells")),]
b.barcode <- metadata.b$barcode

### get T NK cells matrix list
PBMC.list <- list()
for (i in c(2:8,10:13,15:20)){
  
  rca.path.1 <- paste0("../../../SCBAT",i,"/SCBAT",i,"_1/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
  rca.path.2 <- paste0("../../../SCBAT",i,"/SCBAT",i,"_2/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
  rca.1 <- readRDS(rca.path.1)
  rca.2 <- readRDS(rca.path.2)
  
  rca1.raw <- rca.1$raw.data
  colnames(rca1.raw) <- paste(colnames(rca1.raw),"_",i,"1",sep = "")
  
  rca2.raw <- rca.2$raw.data
  colnames(rca2.raw) <- paste(colnames(rca2.raw),"_",i,"2",sep = "")
  
  rca.all <- Seurat:::RowMergeSparseMatrices(rca1.raw, rca2.raw)
  rca.all.nk.t <- rca.all[,which(colnames(rca.all) %in% b.barcode)]
  #Generate a Seurat object
  seu.i <- CreateSeuratObject(counts = rca.all.nk.t, 
                              min.cells = 0, 
                              min.features = 0)
  name.i <- paste0("batch",i)
  
  PBMC.list[[name.i]] <- seu.i
}

# batch 1 library 1
rca11 <- readRDS("../../../SCBAT1/SCBAT1_1/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
rca11.raw <- rca11$raw.data
colnames(rca11.raw) <- paste(colnames(rca11.raw),"_11",sep = "")
rca11.raw.b <- rca11.raw[,which(colnames(rca11.raw) %in% b.barcode)]
seu.11 <- CreateSeuratObject(counts = rca11.raw.b, 
                             min.cells = 0, 
                             min.features = 0)
PBMC.list[["batch1"]] <- seu.11

# batch 9 library 2
rca92 <- readRDS("../../../SCBAT9/SCBAT9_2/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
rca92.raw <- rca92$raw.data
colnames(rca92.raw) <- paste(colnames(rca92.raw),"_92",sep = "")
rca92.raw.b <- rca92.raw[,which(colnames(rca92.raw) %in% b.barcode)]
seu.92 <- CreateSeuratObject(counts = rca92.raw.b, 
                             min.cells = 0, 
                             min.features = 0)
PBMC.list[["batch9"]] <- seu.92

# batch 14 library 2
rca142 <- readRDS("../../../SCBAT14/SCBAT14_2/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
rca142.raw <- rca142$raw.data
colnames(rca142.raw) <- paste(colnames(rca142.raw),"_142",sep = "")
rca142.raw.b <- rca142.raw[,which(colnames(rca142.raw) %in% b.barcode)]
seu.142 <- CreateSeuratObject(counts = rca142.raw.b, 
                              min.cells = 0, 
                              min.features = 0)
PBMC.list[["batch14"]] <- seu.142

saveRDS(PBMC.list,"PBMC.list.rds")

for (i in seq(1,20,1)){
  PBMC.list.i <- PBMC.list[[i]]
  message( paste0(names(PBMC.list)[i], ": ",ncol(PBMC.list.i)) )
}

#normalize and select features
PBMC.list <- lapply(X = PBMC.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})


# select features for downstream integration, and run PCA on each object in the list
features <- SelectIntegrationFeatures(object.list = PBMC.list)
PBMC.list <- lapply(X = PBMC.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


# using batch 5 as reference
reference_dataset <- which(names(PBMC.list) == "batch5")
anchors <- FindIntegrationAnchors(object.list = PBMC.list, 
                                  reference = reference_dataset, 
                                  reduction = "rpca", 
                                  dims = 1:30)



saveRDS(anchors,"anchors.for.reciprocal.pca.B.rds")

# integrate data
PBMC.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(PBMC.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
PBMC.integrated <- ScaleData(PBMC.integrated)
PBMC.integrated <- RunPCA(PBMC.integrated, npcs = 30)

pdf("1.Elbowplot.pdf")
ElbowPlot(PBMC.integrated, ndims =30)
dev.off()

#UMAP
PBMC.integrated <- RunUMAP(PBMC.integrated, reduction = "pca", dims = 1:15)

# clustering
PBMC.integrated <- FindNeighbors(PBMC.integrated, dims = 1:15)
PBMC.integrated <- FindClusters(PBMC.integrated, resolution = 1)

saveRDS(PBMC.integrated,"PBMC.integrated.B.rds")
PBMC.integrated <- readRDS("PBMC.integrated.B.rds")


pdf("1.integrated_UMAP_without_annotation.pdf")
DimPlot(PBMC.integrated, reduction = "umap", label = TRUE)
dev.off()


## Plot NODG and pMito for each cluster
rawdata <- PBMC.integrated@assays$RNA@counts
# NODG
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)

umap.df$NODG <- nGeneVec[umap.df$barcode]
umap.df$pMito <- pMitoVec[umap.df$barcode]
umap.df$seurat_clusters <- factor(umap.df$seurat_clusters,
                                  levels = seq(0, length(unique(umap.df$seurat_clusters))-1, 1))

pdf("3.NODG_comparisons_across_clusters.pdf",
    width = 6, height = 4)

ggplot(umap.df, aes(x=seurat_clusters, y=NODG)) + 
  geom_boxplot(outlier.size = 0, alpha=0.7,
               outlier.shape = NA) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ylab("NODG")+
  xlab("")+
  geom_hline(yintercept = 1000,color="red")
dev.off()

pdf("3.pMito_comparisons_across_clusters.pdf",
    width = 6, height = 4)

ggplot(umap.df, aes(x=seurat_clusters, y=pMito)) + 
  geom_boxplot(outlier.size = 0, alpha=0.7,
               outlier.shape = NA) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ylab("pMito")+
  xlab("")
dev.off()


######## profile NODG  in UMAP
umap.df2 <- FetchData(PBMC.integrated, vars = c("UMAP_1","UMAP_2"))

# NODG
rawdata <- PBMC.integrated@assays$RNA@counts
nGeneVec <- Matrix::colSums(rawdata>0)
nGeneVec <- nGeneVec[(rownames(umap.df2))]
umap.df2$NODG <- nGeneVec

## all cells UMAP highlighting NODG
umap.df2[which(umap.df2$NODG>4000),"NODG"] <- 4000
pdf("4.integrated_UMAP_highlighting_NODG_before_QC.pdf", width = 9, height = 8)
ggplot()+geom_point(data = umap.df2,
                    aes(x=UMAP_1, y=UMAP_2, color=NODG), size=.5,stroke=0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_gradient(low = "blue",high = "orange", limits=c(200,4000))
dev.off()


############## remove bad clusters (median NODG < 1000)
good.barcode <- umap.df[!(umap.df$seurat_clusters %in% c(4,14,16)),"barcode"]

PBMC.integrated.good <- subset(PBMC.integrated, cells=good.barcode)

PBMC.integrated.good <- RunPCA(PBMC.integrated.good, npcs = 30)

pdf("7.Elbowplot_after_remove_bad_clustrs.pdf")
ElbowPlot(PBMC.integrated.good, ndims =30)
dev.off()


#UMAP
PBMC.integrated.good <- RunUMAP(PBMC.integrated.good, reduction = "pca", dims = 1:15)

# clustering
PBMC.integrated.good <- FindNeighbors(PBMC.integrated.good, dims = 1:15)
PBMC.integrated.good <- FindClusters(PBMC.integrated.good, resolution = 1)


umap.df <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","seurat_clusters"))

### annotated based cluster.final results and marker genes
umap.df$barcode <- rownames(umap.df)



## Plot NODG and pMito for each cluster
rawdata <- PBMC.integrated.good@assays$RNA@counts
# NODG
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)

umap.df$NODG <- nGeneVec[umap.df$barcode]
umap.df$pMito <- pMitoVec[umap.df$barcode]
umap.df$seurat_clusters <- factor(umap.df$seurat_clusters,
                                  levels = seq(0, length(unique(umap.df$seurat_clusters))-1, 1))

pdf("7.NODG_comparisons_across_clusters_after_removing_bad.pdf",
    width = 6, height = 4)

ggplot(umap.df, aes(x=seurat_clusters, y=NODG)) + 
  geom_boxplot(outlier.size = 0, alpha=0.7,
               outlier.shape = NA) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ylab("NODG")+
  xlab("")
dev.off()

pdf("7.pMito_comparisons_across_clusters_after_removing_bad.pdf",
    width = 6, height = 4)

ggplot(umap.df, aes(x=seurat_clusters, y=pMito)) + 
  geom_boxplot(outlier.size = 0, alpha=0.7,
               outlier.shape = NA) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  ylab("pMito")+
  xlab("")
dev.off()

######## profile NODG in UMAP
umap.df2 <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2"))

# NODG
rawdata <- PBMC.integrated.good@assays$RNA@counts
nGeneVec <- Matrix::colSums(rawdata>0)
nGeneVec <- nGeneVec[(rownames(umap.df2))]
umap.df2$NODG <- nGeneVec

## all cells UMAP highlighting NODG
umap.df2[which(umap.df2$NODG>4000),"NODG"] <- 4000
pdf("7.integrated_UMAP_highlighting_NODG_after_removing_bad.pdf", width = 9, height = 8)
ggplot()+geom_point(data = umap.df2,
                    aes(x=UMAP_1, y=UMAP_2, color=NODG), size=.5,stroke=0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_gradient(low = "blue",high = "orange", limits=c(200,4000))
dev.off()



pdf("7.integrated_UMAP_without_annotation_after_removing_bad.pdf")
DimPlot(PBMC.integrated.good, label = TRUE, reduction = "umap")
dev.off()


# profile B marker genes
marker_genes <- c("MS4A1","IL4R","CD19",
                  "CD27","CD38","TNFRSF17",
                  "IL3RA")


unique.cluster <- seq(0,length(unique(umap.df$seurat_clusters))-1,1)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(PBMC.integrated.good@assays$RNA@data)) {
    marker.i <- PBMC.integrated.good@assays$RNA@data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$seurat_clusters==cluster.j),"barcode"]
      marker.i.j <- marker.i[,barcode.j]
      exp_per.i.j <- sum(marker.i.j>0)*100/length(marker.i.j)
      exp_mean.i.j <- mean(marker.i.j)
      exp_per.i <- c(exp_per.i,exp_per.i.j)
      exp_mean.i <- c(exp_mean.i,exp_mean.i.j)
    }
    names(exp_per.i) <- unique.cluster
    
    exp_mean.i<- as.vector(scale(exp_mean.i,center = TRUE, scale = TRUE))
    names(exp_mean.i) <- unique.cluster
    marker.df.tmp <- data.frame(gene=rep(i,length(exp_mean.i)),
                                cluster = seq(1,length(exp_mean.i),1),
                                exp_per=exp_per.i,
                                exp_mean = exp_mean.i)
    
    marker.df <- rbind(marker.df,marker.df.tmp)
  }
  
}

marker.df$exp_per <- marker.df$exp_per/20
colnames(marker.df) <- c("gene","cluster","expresson_per","scaled_expression")

pdf("8.Marker_genes_across_different_clusters_B.pdf", width = 6, height = 3)
ggplot()+geom_point(data = marker.df, 
                    mapping = aes(x=cluster,y=gene,
                                  size=expresson_per, color=scaled_expression))+
  scale_colour_gradient(low = "white",high = "red")+
  scale_x_continuous(breaks=seq(1,length(unique.cluster),1),labels=unique.cluster)+
  scale_y_continuous(breaks=seq(1,length(marker_genes),1),labels=marker_genes)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()





### annotated based on marker genes
umap.df$cell <- NA
umap.df[which(umap.df$seurat_clusters==0),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==1),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==2),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==3),"cell"] <- "Memory B"
umap.df[which(umap.df$seurat_clusters==4),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==5),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==6),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==7),"cell"] <- "Memory B"
umap.df[which(umap.df$seurat_clusters==8),"cell"] <- "Memory B"
umap.df[which(umap.df$seurat_clusters==9),"cell"] <- "Memory B"
umap.df[which(umap.df$seurat_clusters==10),"cell"] <- "Memory B"
umap.df[which(umap.df$seurat_clusters==11),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==12),"cell"] <- "Naive B"
umap.df[which(umap.df$seurat_clusters==13),"cell"] <- "Memory B"


umap.df$cell <- factor(umap.df$cell, levels=c("Naive B","Memory B"))
pdf("9.integrated_UMAP_with_cell_annnotations.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP_1, y = UMAP_2, color=cell)) + 
  geom_point(size = 1, stroke=0) + 
  scale_color_manual(values=c("steelblue","orange"))  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

# marker genes
marker_genes <- c("MS4A1","IL4R","CD19",
                  "CD27","CD38","TNFRSF17",
                  "IL3RA")


unique.cluster <- c("Naive B","Memory B")
# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(PBMC.integrated.good@assays$RNA@data)) {
    marker.i <- PBMC.integrated.good@assays$RNA@data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cell==cluster.j),"barcode"]
      marker.i.j <- marker.i[,barcode.j]
      exp_per.i.j <- sum(marker.i.j>0)*100/length(marker.i.j)
      exp_mean.i.j <- mean(marker.i.j)
      exp_per.i <- c(exp_per.i,exp_per.i.j)
      exp_mean.i <- c(exp_mean.i,exp_mean.i.j)
    }
    names(exp_per.i) <- unique.cluster

    names(exp_mean.i) <- unique.cluster
    marker.df.tmp <- data.frame(gene=rep(i,length(exp_mean.i)),
                                cluster = seq(1,length(exp_mean.i),1),
                                exp_per=exp_per.i,
                                exp_mean = exp_mean.i)
    
    marker.df <- rbind(marker.df,marker.df.tmp)
  }
  
}

marker.df$exp_per <- marker.df$exp_per/20
colnames(marker.df) <- c("gene","cluster","expresson_per","scaled_expression")

pdf("8.Marker_genes_across_different_cell_types_B.pdf", width = 5, height = 4)
ggplot()+geom_point(data = marker.df, 
                    mapping = aes(x=cluster,y=gene,
                                  size=expresson_per, color=scaled_expression))+
  scale_colour_gradient(low = "white",high = "red")+
  scale_x_continuous(breaks=seq(1,length(unique.cluster),1),labels=unique.cluster)+
  scale_y_continuous(breaks=seq(1,length(marker_genes),1),labels=marker_genes)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()

umap.df$lib <- unlist(lapply(umap.df$barcode, function(x) unlist(strsplit(x, split = "_"))[2] ))
umap.df$batch <- substr(umap.df$lib,1,nchar(umap.df$lib)-1)

# add cell annotation, lib and bacth
PBMC.integrated.good$cell.annotation <- umap.df[colnames(PBMC.integrated.good), "cell"]
PBMC.integrated.good$lib <- umap.df[colnames(PBMC.integrated.good), "lib"]
PBMC.integrated.good$batch <- umap.df[colnames(PBMC.integrated.good), "batch"]

# distribution of batches across different clusters
cluster.sample <- FetchData(PBMC.integrated.good, vars = c("seurat_clusters","batch"))
cluster.sample$count <- 1
cluster.sample.final <- aggregate(count ~ ., cluster.sample, FUN = sum)
cluster.sample.final$batch <- factor(cluster.sample.final$batch,
                                     levels = seq(1,20,1))

colors37 = c("#466791","#60BF37","#953ADA","#4FBE6C","#CE49D3","#A7B43D","#5A51DC",
             "#D49F36","#552095","#507F2D","#DB37AA","#84B67C","#A06FDA","#DF462A",
             "#5B83DB","#C76C2D","#4F49A3","#82702D","#DD6BBB","#334C22","#D83979",
             "#55BAAD","#DC4555","#62AAD3","#8C3025","#417D61","#862977","#BBA672",
             "#403367","#DA8A6D","#A79CD4","#71482C","#C689D0","#6B2940","#D593A7",
             "#895C8B","#BD5975")


pdf("13.distribution_of_batches_across_different_clusters_percentge.pdf")
ggplot(cluster.sample.final, aes(fill=batch, y=count, x=seurat_clusters)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = colors37[1:length(unique(cluster.sample.final$batch))])+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()

cluster.sample.final$cluster <- factor(cluster.sample.final$seurat_clusters,
                                       levels = unique(cluster.sample.final$seurat_clusters))
pdf("13.distribution_of_clusters_across_different_batch.pdf")
ggplot(cluster.sample.final, aes(fill=seurat_clusters, y=count, x=batch)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = colors37[1:length(unique(cluster.sample.final$seurat_clusters))])+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()



## add sample clinical data
all.meatadata <- read.table("../metadata_for_integrated_data_after_removing_bad_cells_and_QC.txt",
                            sep = "\t", header = TRUE)
all.meatadata.sub <- all.meatadata[which(all.meatadata$barcode %in% colnames(PBMC.integrated.good)),]
# reorder sub metadata according to PBMC.integrated.good barcode order
rownames(all.meatadata.sub) <- all.meatadata.sub$barcode
all.meatadata.sub <- all.meatadata.sub[colnames(PBMC.integrated.good),]

PBMC.integrated.good$sampleBatch <- all.meatadata.sub$sampleBatch
PBMC.integrated.good$severity <- all.meatadata.sub$severity
PBMC.integrated.good$sex <- all.meatadata.sub$sex
PBMC.integrated.good$age <- all.meatadata.sub$age

##### differential abundance in B
abundance.df <- FetchData(PBMC.integrated.good, c("cell.annotation",
                                                  "sampleBatch",
                                                  "severity"))

# remove undefined severity
abundance.df <- abundance.df[(abundance.df$severity != "Undefined"),]

abundance.df$count <- 1

abundance.df.final <- aggregate(count ~ ., abundance.df, FUN = sum)

unique_sampleBatch <- unique(abundance.df.final$sampleBatch)
unique_celltype <- unique(abundance.df.final$cell.annotation)

abundance.df.final_new <- data.frame()
# calculate the proportion for each sample
for (i in seq_along(unique_sampleBatch)){
  sampleBatch.i <- unique_sampleBatch[i]
  abundance.df.final.i <- abundance.df.final[which(abundance.df.final$sampleBatch==sampleBatch.i),]
  
  # add row without corresponding cell types
  unique_celltype.no <- unique_celltype[!(unique_celltype %in% abundance.df.final.i$cell.annotation)]
  if(length(unique_celltype.no) > 0){
    no.df <- data.frame(cell.annotation = unique_celltype.no,
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

#######################################  day 1-8
abundance.df.final_new.cat1 <- abundance.df.final_new.cat[(abundance.df.final_new.cat$cat.x %in% c(1,2,3)),]


unique_celltype <- c("Memory B","Naive B")

### calculate kruskal
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new1 <- abundance.df.final_new.cat1[which(abundance.df.final_new.cat1$cell.annotation==unique_celltype.i),]
  
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
    ylab("proportion in B (%)")+
    xlab("")+
    scale_fill_manual(values = c("#80DD00","#00AE46","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ggtitle(paste0("Kruskal P-value = ",kruskal.i.res))
  ggsave(paste0("16.1.differential_abundance_across_samples_in_",
                  gsub(" ", "_", unique_celltype.i),"_early_and_rising.pdf"),width = 6, height = 4)
  
}


## adjusted p-value heatmap
severity.c <- c("Asymptomatic","Mild1","Mild2",
                "Moderate","Severe","Critical")
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new.down1 <- abundance.df.final_new.cat1[(abundance.df.final_new.cat1$cell.annotation==unique_celltype.i),]
  
  abundance.df.final_new.down1$severity <- factor(abundance.df.final_new.down1$severity,
                                                  levels = severity.c)
  
  stat.test <- wilcox_test(data=abundance.df.final_new.down1,
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
  
  pdf(paste0("16.1.differential_abundance_across_samples_in_",
             gsub(" ", "_", unique_celltype.i),"_early_and_rising_adjusted_pvalue_heatmap.pdf"))
  heatmap.2(as.matrix(heatmap.matrx), Rowv = F, Colv = F, col = col,
            breaks=bk,symm=F, symbreaks=F, scale="none" ,
            trace="none",density.info = "none",
            cellnote=as.matrix(heatmap.matrx),
            notecex=1.2,
            notecol="black")
  dev.off()
}


####################################### day >8
abundance.df.final_new.cat2 <- abundance.df.final_new.cat[(abundance.df.final_new.cat$cat.x %in% c(1,4,5)),]

unique_celltype <- c("Memory B","Naive B")

### calculate kruskal
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new1 <- abundance.df.final_new.cat2[which(abundance.df.final_new.cat2$cell.annotation==unique_celltype.i),]
  
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
    ylab("proportion in B (%)")+
    xlab("")+
    scale_fill_manual(values = c("#80DD00","#00AE46","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ggtitle(paste0("Kruskal P-value = ",kruskal.i.res))
  ggsave(paste0("16.2.differential_abundance_across_samples_in_",
                  gsub(" ", "_", unique_celltype.i),"_peak_and_late.pdf"),width = 6, height = 4)
  
}


## adjusted p-value heatmap
severity.c <- c("Asymptomatic","Mild1","Mild2",
                "Moderate","Severe","Critical")
for (i in seq_along(unique_celltype)){
  unique_celltype.i <- unique_celltype[i]
  abundance.df.final_new.down1 <- abundance.df.final_new.cat2[(abundance.df.final_new.cat2$cell.annotation==unique_celltype.i),]
  
  abundance.df.final_new.down1$severity <- factor(abundance.df.final_new.down1$severity,
                                                  levels = severity.c)
  
  stat.test <- wilcox_test(data=abundance.df.final_new.down1,
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
  
  pdf(paste0("16.3.differential_abundance_across_samples_in_",
             gsub(" ", "_", unique_celltype.i),"_peak_and_late_adjusted_pvalue_heatmap.pdf"))
  heatmap.2(as.matrix(heatmap.matrx), Rowv = F, Colv = F, col = col,
            breaks=bk,symm=F, symbreaks=F, scale="none" ,
            trace="none",density.info = "none",
            cellnote=as.matrix(heatmap.matrx),
            notecex=1.2,
            notecol="black")
  dev.off()
}


#############
##### BCR clonotype
#############
PBMC.integrated.good <- readRDS("PBMC.integrated.B.good.rds")
source("../utils/quy_BCR_utils.R")
clonotype <- read.table("BCR_barcodes_across_all_batch_clonotypes.txt", sep = "\t", header = TRUE)

clonotype$lib <- paste0(ceiling(clonotype$batch/2), 
                        (2 - ceiling(clonotype$batch/2)*2 + clonotype$batch))
# remove bad libraries 1-2, 9-1, 14-1
clonotype <- clonotype[!(clonotype$lib %in% c("12","91","141")),]

clonotype$barcode <- paste(clonotype$barcode, clonotype$lib, sep = "_")

celltype <- FetchData(PBMC.integrated.good, vars = c("UMAP_1","UMAP_2","cell.annotation","severity",
                                                     "sampleBatch"))
celltype$barcode <- rownames(celltype)
clonotype.sub <- clonotype[which(clonotype$barcode %in% celltype$barcode),c("barcode","BCRgene")]
rownames(clonotype.sub) <- clonotype.sub$barcode

celltype$BCRgene <- NA
celltype$Frequency <- NA
celltype$Proportion <- NA
celltype[clonotype.sub$barcode, "BCRgene"] <- clonotype.sub$BCRgene



## check distribution of deteced BCR numbers of across sample
celltype.t <- celltype[!(is.na(celltype$BCRgene)),]
celltype.t.table <- as.data.frame(table(celltype.t$sampleBatch))
pdf("18.distribution_of_detected_BCR_number_per_sample.pdf")
hist(celltype.t.table$Freq, breaks = seq(0,750,10),
     xlab = "TCR number", main = "", xlim = c(0,800), ylim = c(0,30))
abline(v = 30,col="red")
abline(v = median(celltype.t.table$Freq),col="blue")
dev.off()
## end check distribution of deteced BCR numbers of across sample


celltype.sub <- celltype[(!is.na(celltype$BCRgene)),]
number_of_BCR <- c()
unique_sample <- unique(celltype.sub$sampleBatch)
for (i in seq_along(unique_sample)){
  sample.i <- unique_sample[i]
  celltype.sub.i <- celltype.sub[which(celltype.sub$sampleBatch==sample.i),]
  number_of_BCR <- c(number_of_BCR,nrow(celltype.sub.i))
  if (nrow(celltype.sub.i) < 30){
    celltype[celltype.sub.i$barcode,"Frequency"] <- 0
    celltype[celltype.sub.i$barcode,"Proportion"] <- 0
  } else {
    celltype.sub.i.freq <- calClonalFreq(celltype.sub.i[,c("barcode","BCRgene")], 
                                         cloneType="BCRgene")
    celltype[celltype.sub.i.freq$barcode,"Frequency"] <- celltype.sub.i.freq$Frequency
    celltype[celltype.sub.i.freq$barcode,"Proportion"] <- celltype.sub.i.freq$Proportion
  }
  
}


ggplot()+
  geom_point(celltype[(is.na(celltype$Proportion)),], 
             mapping=aes(x=UMAP_1,y=UMAP_2), size=1, color="grey")+
  geom_point(celltype[((!is.na(celltype$Proportion)) & 
                         celltype$Proportion <= 0.05),], 
             mapping=aes(x=UMAP_1,y=UMAP_2, color=Proportion), size=1)+
  geom_point(celltype[((!is.na(celltype$Proportion)) & 
                         celltype$Proportion > 0.05 & celltype$Proportion <= 0.1),], 
             mapping=aes(x=UMAP_1,y=UMAP_2, color=Proportion), size=1)+
  geom_point(celltype[((!is.na(celltype$Proportion)) & 
                         celltype$Proportion > 0.1 & celltype$Proportion <= 0.2),], 
             mapping=aes(x=UMAP_1,y=UMAP_2, color=Proportion), size=1)+
  geom_point(celltype[((!is.na(celltype$Proportion)) & 
                         celltype$Proportion > 0.2),], 
             mapping=aes(x=UMAP_1,y=UMAP_2, color=Proportion), size=1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  scale_color_gradient2(low = "darkblue",mid="orange",high = "yellow", limits=c(0,0.2), midpoint = 0.1)
ggsave(paste0("18.integrated_UMAP_with_TCR_clonatype_proportion.png"),device = "png",width = 6, height = 5)


############# calculate top 5 expanded BCR ratio across different severity (top30)
celltype.bk.sub <- celltype[(!is.na(celltype$BCRgene)),]
# remove undefined severity
celltype.bk.sub <- celltype.bk.sub[(celltype.bk.sub$severity != "Undefined"),]



#### check distribution of unique BCR numbers across sample
uniq.bcr <- data.frame()
unique_sample <- unique(celltype.bk.sub$sampleBatch)
for (i in seq_along(unique_sample)){
  unique_sample.i <- unique_sample[i]
  celltype.bk.sub.i <- celltype.bk.sub[(celltype.bk.sub$sampleBatch == unique_sample.i),
                                       c("sampleBatch","BCRgene")]
  celltype.bk.sub.i <- unique(celltype.bk.sub.i)
  new.add <- data.frame(sampleBatch = unique_sample.i,
                        BCR_num = nrow(celltype.bk.sub.i))
  uniq.bcr <- rbind(uniq.bcr,new.add)
}
pdf("18.distribution_of_detected_unique_BCR_number_per_sample.pdf")
hist(uniq.bcr$BCR_num, breaks = seq(0,710,10),xlab = "Unique BCR number", main = "")
abline(v = 30,col="red")
dev.off()
#### end check distribution of unique BCR numbers across sample



# get infection status
infect.stat <- read.delim("../infection_status_for_all_sampleBatches.txt",
                          sep = "\t", header = TRUE)
infect.stat.sub <- infect.stat[,c("sampleBatch","cat.x","cat")]
celltype.bk.cat.sub <- merge(celltype.bk.sub, infect.stat.sub,
                             by = "sampleBatch")


# day 1-8
celltype.bk.cat.sub1 <- celltype.bk.cat.sub[(celltype.bk.cat.sub$cat.x %in% c(1,2,3)),]
unique_sampleBatch <- unique(celltype.bk.cat.sub1$sampleBatch)

for (a in seq(3,20,1)){
  top.proportion.df <- data.frame()
  for (i in seq_along(unique_sampleBatch)){
    samplebatch.i <- unique_sampleBatch[i]
    celltype.bk.sub.i <- celltype.bk.cat.sub1[(celltype.bk.cat.sub1$sampleBatch == samplebatch.i),]
    if (length(unique(celltype.bk.sub.i$BCRgene)) < 30){
      next
    }
    BCRgene.freq <- data.frame(table(celltype.bk.sub.i$BCRgene))
    BCRgene.freq <- BCRgene.freq[order(BCRgene.freq$Freq, decreasing = TRUE),]
    
    new.df <- data.frame(sampleBatch = samplebatch.i,
                         severity = celltype.bk.sub.i$severity[1],
                         top_proportion = 100*sum(BCRgene.freq[1:a,"Freq"])/sum(BCRgene.freq[1:30,"Freq"]))
    top.proportion.df <- rbind(top.proportion.df,new.df)
  }
  
  top.proportion.df$severity <- factor(top.proportion.df$severity,
                                       levels = c("Asymptomatic","Mild1","Mild2",
                                                  "Moderate","Severe","Critical"))
  
  stat.test <- wilcox_test(data=top.proportion.df,
                           top_proportion ~ severity,
                           p.adjust.method = "BH")
  
  stat.test <- stat.test %>% add_xy_position(x = "severity")

  kruskal.a <- kruskal.test(top_proportion~severity, data=top.proportion.df)
  kruskal.a.res <- kruskal.a$p.value
  
  ggviolin(top.proportion.df, 
           x = "severity", y = "top_proportion", fill = "severity",
           alpha=0.7, trim = TRUE) + 
    geom_boxplot(fill="white", width=.05, outlier.shape = NA)+
    geom_jitter(size=.8, stroke=0, width = .2,height=0 )+
    scale_fill_manual(values = c("#E1BC00","#A5D700","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ylab(paste0("top",a," proportion (%)"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    stat_pvalue_manual(stat.test, label = "p.adj")+
    ggtitle(paste0("Kruskal p = ",kruskal.a.res))
  ggsave(paste0("18.1.top",a,"_BCR_proportion_comparsion_across_severities_early_rising.pdf"))
  
}

# day > 8
celltype.bk.cat.sub1 <- celltype.bk.cat.sub[(celltype.bk.cat.sub$cat.x %in% c(1,4,5)),]
unique_sampleBatch <- unique(celltype.bk.cat.sub1$sampleBatch)

for (a in seq(3,20,1)){
  top.proportion.df <- data.frame()
  for (i in seq_along(unique_sampleBatch)){
    samplebatch.i <- unique_sampleBatch[i]
    celltype.bk.sub.i <- celltype.bk.cat.sub1[(celltype.bk.cat.sub1$sampleBatch == samplebatch.i),]
    if (length(unique(celltype.bk.sub.i$BCRgene)) < 30){
      next
    }
    BCRgene.freq <- data.frame(table(celltype.bk.sub.i$BCRgene))
    BCRgene.freq <- BCRgene.freq[order(BCRgene.freq$Freq, decreasing = TRUE),]
    
    new.df <- data.frame(sampleBatch = samplebatch.i,
                         severity = celltype.bk.sub.i$severity[1],
                         top_proportion = 100*sum(BCRgene.freq[1:a,"Freq"])/sum(BCRgene.freq[1:30,"Freq"]))
    top.proportion.df <- rbind(top.proportion.df,new.df)
  }
  
  top.proportion.df$severity <- factor(top.proportion.df$severity,
                                       levels = c("Asymptomatic","Mild1","Mild2",
                                                  "Moderate","Severe","Critical"))
  
  stat.test <- wilcox_test(data=top.proportion.df,
                           top_proportion ~ severity,
                           p.adjust.method = "BH")
  
  stat.test <- stat.test %>% add_xy_position(x = "severity")
  
  kruskal.a <- kruskal.test(top_proportion~severity, data=top.proportion.df)
  kruskal.a.res <- kruskal.a$p.value
  
  
  ggviolin(top.proportion.df, 
           x = "severity", y = "top_proportion", fill = "severity",
           alpha=0.7, trim = TRUE) + 
    geom_boxplot(fill="white", width=.05, outlier.shape = NA)+
    geom_jitter(size=.8, stroke=0, width = .2,height=0 )+
    scale_fill_manual(values = c("#E1BC00","#A5D700","#00997F",
                                 "#1313A9","#B20086","#F20000"))+
    ylab(paste0("top",a," proportion (%)"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    stat_pvalue_manual(stat.test, label = "p.adj")+
    ggtitle(paste0("Kruskal p = ",kruskal.a.res))
  ggsave(paste0("18.2.top",a,"_BCR_proportion_comparsion_across_severities_peak_late.pdf"))
  
}



