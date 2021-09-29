library(RCAv2)
library(dplyr)
library(ggplot2)
library(gridExtra)


################################################################################
########################### read singlets rds #################################
################################################################################

# read Batch 1-1
rca <- readRDS("../../SCBAT1/SCBAT1_1/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
rca <- dataLogNormalise(rca)
norm.data.all <- rca$data
colnames(norm.data.all) <- paste(colnames(norm.data.all), "_11", sep = "")
raw.data.all <- rca$raw.data
colnames(raw.data.all) <- paste(colnames(raw.data.all), "_11", sep = "")

# merge gene-cell matrices from 20 batches
for (i in seq(2,20,1)){
  for (j in c(1,2)){
    #skip batch 9-1 due to low quality
    if (i==9 & j==1){
      next
    }
    #skip batch 14-1 due to low quality
    if (i==14 & j==1){
      next
    }
    message(paste(i,j, sep = "_"))
    
    rca.path <- paste0("../../SCBAT",i,"/SCBAT",i,"_",j,"/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
    rca <- readRDS(rca.path)
    rca <- dataLogNormalise(rca)
    
    norm.data <- rca$data
    colnames(norm.data) <- paste(colnames(norm.data), "_",i,j, sep = "")
    
    raw.data <- rca$raw.data
    colnames(raw.data) <- paste(colnames(raw.data), "_",i,j, sep = "")
    
    raw.data.all <- Seurat:::RowMergeSparseMatrices(raw.data.all, raw.data)
    norm.data.all <- Seurat:::RowMergeSparseMatrices(norm.data.all, norm.data)
  }
}

raw.data.all <- raw.data.all[rownames(norm.data.all),]

################################################################################
########################### create RCA object #################################
################################################################################

RCA.all <- createRCAObject(rawData = raw.data.all, normData = norm.data.all)
############ project to all immune panels ################

# remove gene "^MT-|^ERCC|^RPS|^RPL|^HSP"
gene.row <- rownames(RCA.all$raw.data)
gene.row.good <- gene.row[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                               x = gene.row, invert = TRUE)]

RCA.all$data <- RCA.all$data[gene.row.good,]


# get projection results against global panel
RCA.all <- dataProject(RCA.all, method = "GlobalPanel",
                      corMeth = "pearson", scale = TRUE)
global.proj <- as.data.frame(as.matrix(RCA.all$projection.data))
global.proj.immune <- read.table("./rownames_of_glocal_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
RCA.all <- dataProjectMultiPanel(RCA.all,method = list("NovershternPanel", 
                                                     "MonacoPanel"),
                                scale = TRUE,
                                corMeth = "pearson")
two.proj <- as.data.frame(as.matrix(RCA.all$projection.data))

# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")

# Assign projection result to RCA object
RCA.all$projection.data <- proj.all


RCAv2:::elbowPlot(RCA.all, nPCs = 50, approx = TRUE)

RCA.all <- dataSClust(RCA.all, res=2, nPCs = 13, approx = TRUE)

#UMAP
RCA.all <- computeUMAP(RCA.all)

#Estimate the most probable cell type label for each cell
RCA.all <- estimateCellTypeFromProjection(RCA.all,confidence = NULL)

# saveRDS(RCA.all, "RCA.combined.all.rds")
# RCA.all <- readRDS("RCA.combined.all.rds")

### annotate each cluster
color.cluster <- data.frame(cell=unlist(RCA.all$cell.Type.Estimate),
                            color=RCA.all$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)

# Annotate cells based on cluster.final
# clusters according to the major cell type annotations

umap.df <- RCA.all$umap.coordinates
umap.df$barcode <- rownames(umap.df)
umap.df$cluster <- RCA.all$clustering.out$dynamicColorsList[[1]]
umap.df$cell <- NA
umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "BAD"
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "Monocytes"

umap.df[which(umap.df$cluster=="green"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "BAD"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "NK cells"

umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "BAD"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "BAD"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "B cells"


pdf("2.integrated_UMAP_with_RCA_annotation.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cell)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=c("blue","grey","brown","darkgreen",
                              "darkorange","red","pink","turquoise")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

unique.cluster <- unique(umap.df$cluster)
umap.df$cluster <- factor(umap.df$cluster, levels = unique.cluster)

pdf("2.integrated_UMAP_without_annotation.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cluster)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=unique.cluster) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

#### double check with cell annotation with PBMC marker gene expression levels
#### annotate the clusters lacking no major PBMC marker genes as "BAD" above
# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")

unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all$data)) {
    marker.i <- RCA.all$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("3.Marker_genes_across_different_clusters.pdf", width = 12, height = 6)
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



### plot QC metrics for each cluster
rawdata <- RCA.all$raw.data
# get NODG
nGeneVec <- Matrix::colSums(rawdata>0)
# get pMito
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
##### QC on each cluster
QC_list <- list()
unique_celltype <- unique(umap.df$cluster)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- rownames(umap.df[which(umap.df$cluster==cell.i),])
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.5, stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  QC_list[[i]] <- p.i
}

p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))
pdf("3.QC_plot_of_each_clusters.pdf", 
    width = 3, height = 3*length(unique(umap.df$cluster)))
grid.arrange(p.all,ncol=1)
dev.off()


################################################################################
########################### remove bad cells #################################
################################################################################

good.barcode <- umap.df[which(umap.df$cell != "BAD"),"barcode"]
good.barcode.index <- which(colnames(RCA.all$raw.data) %in% good.barcode)
RCA.all$raw.data <- RCA.all$raw.data[, good.barcode.index]
RCA.all$data <- RCA.all$data[, good.barcode.index]
RCA.all$projection.data <- RCA.all$projection.data[, good.barcode.index]
RCA.all$cell.Type.Estimate <- RCA.all$cell.Type.Estimate[good.barcode.index]

RCA.all$clustering.out$dynamicColorsList[[1]] <- RCA.all$clustering.out$dynamicColorsList[[1]][good.barcode.index]
saveRDS(RCA.all, "RCA.combined.all.after.removing.bad.cells.rds")


RCAv2:::elbowPlot(RCA.all, nPCs = 50, approx = TRUE,filename="elbow.plot.after.removing.bad.cells.png")
RCA.all <- dataSClust(RCA.all, res=2, nPCs = 14, approx = TRUE)

saveRDS(RCA.all,"RCA.combined.all.after.removing.bad.cells.rds")

RCA.all <- readRDS("RCA.combined.all.after.removing.bad.cells.rds")

########### profile PBMC marker genes across different cluster after 1st round of removing BAD clusters and reclustering
# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")

umap.df <- data.frame(cluster=RCA.all$clustering.out$dynamicColorsList[[1]])
umap.df$barcode <- colnames(RCA.all$raw.data)


unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all$data)) {
    marker.i <- RCA.all$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("5.Marker_genes_across_different_clusters_after_removing_bad_cells.pdf", width = 12, height = 6)
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


#### plot QC metrics for each cluster
rawdata <- RCA.all$raw.data
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
##### QC on each cluster
QC_list <- list()
unique_celltype <- unique(umap.df$cluster)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- umap.df[which(umap.df$cluster==cell.i),"barcode"]
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.5, stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  QC_list[[i]] <- p.i
}

p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))
pdf("5.QC_plot_of_each_clusters_after_removing_bad_cells.pdf", 
    width = 3, height = 3*length(unique(umap.df$cluster)))
grid.arrange(p.all,ncol=1)
dev.off()

# remove cluster skyblue3
good.barcode <- umap.df[!(umap.df$cluster %in% c("mediumpurple3","lightsteelblue1")) ,"barcode"]

################################################################################
########################### remove bad cells again #############################
########################### cluster mediumpurple3 and lightsteelblue1 ##########
################################################################################

good.barcode.index <- which(colnames(RCA.all$raw.data) %in% good.barcode)
RCA.all$raw.data <- RCA.all$raw.data[, good.barcode.index]
RCA.all$data <- RCA.all$data[, good.barcode.index]
RCA.all$projection.data <- RCA.all$projection.data[, good.barcode.index]
RCA.all$cell.Type.Estimate <- RCA.all$cell.Type.Estimate[good.barcode.index]

RCA.all$clustering.out$dynamicColorsList[[1]] <- RCA.all$clustering.out$dynamicColorsList[[1]][good.barcode.index]

RCAv2:::elbowPlot(RCA.all, nPCs = 50, approx = TRUE,filename="elbow.plot.after.removing.bad.cells.again.png")
RCA.all <- dataSClust(RCA.all, res=2, nPCs = 14, approx = TRUE)


RCA.all <- readRDS("RCA.combined.all.after.removing.bad.cells.again.rds")
# compute umap
RCA.all <- computeUMAP(RCA.all)

saveRDS(RCA.all, "RCA.combined.all.after.removing.bad.cells.again.rds")
# RCA.all <- readRDS("RCA.combined.all.after.removing.bad.cells.again.rds")


# profile marker genes across each cluster
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")

umap.df <- data.frame(cluster=RCA.all$clustering.out$dynamicColorsList[[1]])
umap.df$barcode <- colnames(RCA.all$raw.data)


unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(RCA.all$data)) {
    marker.i <- RCA.all$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
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

pdf("6.Marker_genes_across_different_clusters_after_removing_bad_cells_again.pdf", width = 12, height = 6)
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



### annotate each cluster
color.cluster <- data.frame(cell=unlist(RCA.all$cell.Type.Estimate),
                            color=RCA.all$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)

# Based on cluster.final and marker gene profiling results, annotate each cluster

umap.df <- RCA.all$umap.coordinates
umap.df$barcode <- rownames(umap.df)
umap.df$cluster <- RCA.all$clustering.out$dynamicColorsList[[1]]
umap.df$cell <- NA
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="green"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "T cells"

umap.df[which(umap.df$cluster=="red"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "T cells"


umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "B cells"

umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="plum2"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "Monocytes"

umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "NK cells"


pdf("6.integrated_UMAP_with_RCA_annotation_after_removing_bad_cells_again.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cell)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=c("blue","brown","darkgreen",
                              "darkorange","red","pink","turquoise")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

unique.cluster <- unique(umap.df$cluster)
umap.df$cluster <- factor(umap.df$cluster, levels = unique.cluster)

pdf("6.integrated_UMAP_without_annotation_after_removing_bad_cells_again.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cluster)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=unique.cluster) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()


## profile TCR BCR presence
######### TCR
tcr.all <- read.table("TCR_barcodes_across_all_batch.txt", stringsAsFactors = FALSE,
                      header = TRUE)
tcr.all$lib <- paste0(ceiling(tcr.all$batch/2), 
                      (2 - ceiling(tcr.all$batch/2)*2 + tcr.all$batch))
# remove batch 1-2, 9-1, 14-1
tcr.all <- tcr.all[!(tcr.all$lib %in% c("12","91","141")),]

tcr.all$barcode <- paste(tcr.all$barcode, tcr.all$lib, sep = "_")

# add TCR info in umap
umap.df$TCR <- 0
umap.df[which(umap.df$barcode %in% tcr.all$barcode),"TCR"] <- 1
umap.df$TCR <- factor(umap.df$TCR, levels = c(0,1))

pdf("6.integrated_UMAP_with_TCR_presence_after_removing_bad_cells_again.pdf",width=7,height=7)
ggplot() + 
  geom_point(data = umap.df[which(umap.df$TCR==0),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$TCR==1),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="red") +
  theme_bw(10) 
dev.off()


######### BCR
bcr.all <- read.table("BCR_barcodes_across_all_batch.txt", stringsAsFactors = FALSE,
                      header = TRUE)
bcr.all$lib <- paste0(ceiling(bcr.all$batch/2), 
                      (2 - ceiling(bcr.all$batch/2)*2 + bcr.all$batch))
# remove batch 1-2, 9-1, 14-1
bcr.all <- bcr.all[!(bcr.all$lib %in% c("12","91","141")),]

bcr.all$barcode <- paste(bcr.all$barcode, bcr.all$lib, sep = "_")

# add BCR info in umap
umap.df$BCR <- 0
umap.df[which(umap.df$barcode %in% bcr.all$barcode),"BCR"] <- 1
umap.df$BCR <- factor(umap.df$BCR, levels = c(0,1))

pdf("5.integrated_UMAP_with_BCR_presence_after_removing_bad_cells_again.pdf",width=7,height=7)
ggplot() + 
  geom_point(data = umap.df[which(umap.df$BCR==0),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$BCR==1),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="blue") +
  theme_bw(10) 
dev.off()

################################################################################
########################### QC on each cluster #################################
################################################################################

# plot QC
rawdata <- RCA.all$raw.data
nGeneVec <- Matrix::colSums(rawdata>0)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(rawdata), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(rawdata[mito.genes, , drop = FALSE])/Matrix::colSums(rawdata)
##### QC on each cluster
QC_list <- list()
unique_celltype <- unique(umap.df$cell)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- rownames(umap.df[which(umap.df$cell==cell.i),])
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.5, stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  QC_list[[i]] <- p.i
}

p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))
pdf("7.QC_plot_of_each_cell_type_after_removing_bad_cells.pdf", 
    width = 3, height = 3*length(unique(umap.df$cell)))
grid.arrange(p.all,ncol=1)
dev.off()

# QC on each cluster
##### QC filtering
nodg.low <- c(500,500,500,500,100,1000,500)
nodg.up <- c(2700,2700,3000,3400,1300,4000,3500)
pmt.low <- c(0.001,0.001,0.001,0.001,0.001,0.001,0.001)
pmt.up <- c(0.06,0.08,0.09,0.09,0.15,0.1,0.1)
good.barcode <- c()
QC_list <- list()
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- rownames(umap.df[which(umap.df$cell==cell.i),])
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  rownames(qc.i) <- barcode.i
  nodg.low.i <- nodg.low[i]
  nodg.up.i <- nodg.up[i]
  pmt.low.i <- pmt.low[i]
  pmt.up.i <- pmt.up[i]
  qc.good.i <- qc.i[which(qc.i$NODG >= nodg.low.i & qc.i$NODG <= nodg.up.i &
                            qc.i$pMt >= pmt.low.i & qc.i$pMt <= pmt.up.i),]
  good.barcode <- c(good.barcode, rownames(qc.good.i))
  qc.i$quality <- "bad"
  qc.i[rownames(qc.good.i),"quality"] <- "good"
  qc.i$quality <- factor(qc.i$quality, levels = c("good","bad"))
  p.i <- ggplot2::ggplot() +
    ggplot2::geom_point(data = qc.i, ggplot2::aes(x = NODG, y = pMt,color=quality),size = .2,stroke=0) + 
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d(data = qc.i, mapping=aes(x = NODG, y = pMt)) +
    ggplot2::scale_color_manual(values = c("black","grey"))
  QC_list[[i]] <- p.i
}
p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(QC_list))
pdf("7.QC_plot_of_each_cell_type_after_removing_bad_cells_with_target_areas.pdf", 
    width = 4, height = 3*length(unique(umap.df$cell)))
grid.arrange(p.all,ncol=1)
dev.off()

# after QC

umap.df.good <- umap.df[good.barcode,]
QC_list <- list()
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- rownames(umap.df.good[which(umap.df.good$cell==cell.i),])
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = 0.5, stroke=0) + 
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  QC_list[[i]] <- p.i
}
p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))
pdf("8.QC_plot_of_each_cell_type_after_removing_bad_cells_and_after_QC.pdf", 
    width = 3, height = 3*length(unique(umap.df.good$cell)))
grid.arrange(p.all,ncol=1)
dev.off()

#### add original sample information
umap.df.good$sample <- NA

for (i in seq(1,20,1)) {
  for (j in c(1,2)){
    #skip batch 1-2
    if (i==1 & j==2){
      next
    }
    #skip batch 9-1
    if (i==9 & j==1){
      next
    }
    #skip batch 14-1
    if (i==14 & j==1){
      next
    }
    
    demuxlet.path <- paste0("/mnt/quy_1TB/COVID_project/SCAN_COVID/demuxlet_analysis/SCBAT",i,
                            "_",j,"/SCBAT",i,"_",j,"_demuxlet.best")
    demuxlet.file <- read.table(demuxlet.path, sep = '\t', header = TRUE)
    demuxlet.file$BARCODE <- paste(demuxlet.file$BARCODE,"_",i,j, sep = "")
    demuxlet.file <- demuxlet.file[which(demuxlet.file$BARCODE %in% umap.df.good$barcode),]
    umap.df.good[demuxlet.file$BARCODE,"sample"] <- demuxlet.file$SNG.1ST
  }
}

umap.df.good$lib <- unlist(lapply(umap.df.good$barcode, function(x) unlist(strsplit(x, split = "_"))[2] ))
umap.df.good$batch <- substr(umap.df.good$lib,1,nchar(umap.df.good$lib)-1)

umap.df.good$sampleBatch <- paste(umap.df.good$batch, umap.df.good$sample, sep = "_")

## plot sample
table.df <- data.frame(table(umap.df.good$sampleBatch))
table.df <- table.df[order(table.df$Freq),]
pdf("9.barplot_of_samples.pdf", width = 15, height = 5)
ggplot(table.df, aes(x=Var1, y=Freq))+
  geom_bar(stat="identity",fill="steelblue")+theme_minimal()+
  xlab("sample")+ylab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# remove samplebatch less than 200 cells and unknown severity

sampleBatch.good <- table.df[which(table.df$Freq>=200),"Var1"]
sampleBatch.good <- sampleBatch.good[!(sampleBatch.good %in% c("12_CG61","4_CG17","4_CG18","10_CG54"))]
umap.df.good.filtered <- umap.df.good[which(umap.df.good$sampleBatch %in% sampleBatch.good),]

# post-QC
good.barcode <- umap.df.good.filtered$barcode
good.barcode.index <- which(colnames(RCA.all$raw.data) %in% good.barcode)
RCA.all$raw.data <- RCA.all$raw.data[, good.barcode.index]
RCA.all$data <- RCA.all$data[, good.barcode.index]
RCA.all$projection.data <- RCA.all$projection.data[, good.barcode.index]
RCA.all$cell.Type.Estimate <- RCA.all$cell.Type.Estimate[good.barcode.index]

RCA.all$clustering.out$dynamicColorsList[[1]] <- RCA.all$clustering.out$dynamicColorsList[[1]][good.barcode.index]

saveRDS(RCA.all, "RCA.combined.all.after.removing.bad.cells.again.and.QC.rds")
RCA.all <- computeUMAP(RCA.all)
RCA.all <- readRDS("RCA.combined.all.after.removing.bad.cells.again.and.QC.rds")



umap.df <- RCA.all$umap.coordinates
umap.df$barcode <- rownames(umap.df)
umap.df$cluster <- RCA.all$clustering.out$dynamicColorsList[[1]]
umap.df$cell <- NA
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="green"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="saddlebrown"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="orange"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "T cells"

umap.df[which(umap.df$cluster=="red"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="paleturquoise"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="bisque4"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="steelblue"),"cell"] <- "T cells"


umap.df[which(umap.df$cluster=="darkmagenta"),"cell"] <- "B cells"

umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="sienna3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="black"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="lightsteelblue1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="plum2"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkorange"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="plum1"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown4"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="skyblue3"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="mediumpurple3"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="darkorange2"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellowgreen"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="lightcyan1"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="skyblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkslateblue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="violet"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="orangered4"),"cell"] <- "Monocytes"

umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="white"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="ivory"),"cell"] <- "Platelets"
umap.df[which(umap.df$cluster=="floralwhite"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="darkolivegreen"),"cell"] <- "NK cells"



pdf("9.integrated_UMAP_with_RCA_annotation_after_removing_bad_cells_again_and_QC.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cell)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=c("blue","brown","darkgreen",
                              "darkorange","red","pink","turquoise")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

unique.cluster <- unique(umap.df$cluster)
umap.df$cluster <- factor(umap.df$cluster, levels = unique.cluster)

pdf("9.integrated_UMAP_without_annotation_after_removing_bad_cells_again_and_QC.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cluster)) + 
  geom_point(size = .1, stroke=0) + 
  scale_color_manual(values=unique.cluster) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()


## profile TCR BCR presence
######### TCR
tcr.all <- read.table("TCR_barcodes_across_all_batch.txt", stringsAsFactors = FALSE,
                      header = TRUE)
tcr.all$lib <- paste0(ceiling(tcr.all$batch/2), 
                      (2 - ceiling(tcr.all$batch/2)*2 + tcr.all$batch))
# remove bad libraries 1-2, 9-1, 14-1
tcr.all <- tcr.all[!(tcr.all$lib %in% c("12","91","141")),]

tcr.all$barcode <- paste(tcr.all$barcode, tcr.all$lib, sep = "_")

# add TCR info in umap
umap.df$TCR <- 0
umap.df[which(umap.df$barcode %in% tcr.all$barcode),"TCR"] <- 1

pdf("9.integrated_UMAP_with_TCR_presence_after_removing_bad_cells_again_and_QC.pdf",width=7,height=7)
ggplot() + 
  geom_point(data = umap.df[which(umap.df$TCR==0),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$TCR==1),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="red") +
  theme_bw(10) 
dev.off()


######### BCR
bcr.all <- read.table("BCR_barcodes_across_all_batch.txt", stringsAsFactors = FALSE,
                      header = TRUE)
bcr.all$lib <- paste0(ceiling(bcr.all$batch/2), 
                      (2 - ceiling(bcr.all$batch/2)*2 + bcr.all$batch))
# remove bad libraries 1-2, 9-1, 14-1
bcr.all <- bcr.all[!(bcr.all$lib %in% c("12","91","141")),]

bcr.all$barcode <- paste(bcr.all$barcode, bcr.all$lib, sep = "_")

# add BCR info in umap
umap.df$BCR <- 0
umap.df[which(umap.df$barcode %in% bcr.all$barcode),"BCR"] <- 1

pdf("9.integrated_UMAP_with_BCR_presence_after_removing_bad_cells_again_and_QC.pdf",width=7,height=7)
ggplot() + 
  geom_point(data = umap.df[which(umap.df$BCR==0),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="grey") + 
  geom_point(data = umap.df[which(umap.df$BCR==1),], 
             mapping = aes(x = UMAP1, y = UMAP2),
             size = .1, stroke=0,col="blue") +
  theme_bw(10) 
dev.off()



################################################################################
########################### add sample information #############################
################################################################################
umap.df$sample <- NA

for (i in seq(1,20,1)) {
  for (j in c(1,2)){
    #skip library 1-2
    if (i==1 & j==2){
      next
    }
    #skip library 9-1
    if (i==9 & j==1){
      next
    }
    #skip library 14-1
    if (i==14 & j==1){
      next
    }
    
    demuxlet.path <- paste0("../demuxlet_analysis/SCBAT",i,
                            "_",j,"/SCBAT",i,"_",j,"_demuxlet.best")
    demuxlet.file <- read.table(demuxlet.path, sep = '\t', header = TRUE)
    demuxlet.file$BARCODE <- paste(demuxlet.file$BARCODE,"_",i,j, sep = "")
    demuxlet.file <- demuxlet.file[which(demuxlet.file$BARCODE %in% umap.df$barcode),]
    umap.df[demuxlet.file$BARCODE,"sample"] <- demuxlet.file$SNG.1ST
  }
}

umap.df$lib <- unlist(lapply(umap.df$barcode, function(x) unlist(strsplit(x, split = "_"))[2] ))
umap.df$batch <- substr(umap.df$lib,1,nchar(umap.df$lib)-1)

umap.df$sampleBatch <- paste(umap.df$batch, umap.df$sample, sep = "_")

# sample clinical data
sample.cli <- read.table("sample_metadata.txt", sep = "\t", header = TRUE)

sample.cli$sampleBatch <- paste(sample.cli$batch, sample.cli$sample, sep = "_")
sample.cli <- sample.cli[,c("sampleBatch","severity","age","sex","ethnicity")]
sample.cli <- sample.cli[which(sample.cli$sampleBatch %in% umap.df$sampleBatch), ]

sample.cli <- unique(sample.cli)
umap.df.merge <- merge(x=umap.df, y = sample.cli, by = "sampleBatch")

##### mark some samples as undefined in severity due to the ambiguous severity states
umap.df.merge[(umap.df.merge$sampleBatch %in% c("1_CG66","7_CG63","10_CG91",
                                                "15_CG63","16_CG63")), "severity"] <- "Undefined"
write.table(umap.df.merge, "metadata_for_integrated_data_after_removing_bad_cells_and_QC.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





