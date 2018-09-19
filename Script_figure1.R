library(Seurat)
library(dplyr)
library(RColorBrewer)
source("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/Script_functions.R")

## Get random colors from a ramp based onf Spectral from RColorBrewer
colfunc <- colorRampPalette(brewer.pal(11, "Spectral"))
random_spectral <- sample(colfunc(14))
superman_icecream_palette <- c("#E55748", "#3C7AB6", "#FDD884", "#FDB164", "#5E4FA2", "#FEF3AB", "#4DA7B0", "#F6814C", "#DCF199", "#F5FBAF", "#7BCAA4", "#AFDEA3", "#C82F4C", "#9E0142")

#D145tra.data <- Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2330_Miller_HT189-Tracheal-Epithelium/Sample_HT-189-Tracheal-Epi/outs/HT189-TrachealEpi-filtered_gene_bc_matrices/hg19/")
#colnames(D145tra.data) <- paste("D145tra",colnames(D145tra.data),sep = "_")
#D145tra <- CreateSeuratObject(raw.data = D145tra.data, min.cells = 3, min.genes = 200)

#vivo <- readRDS("/home/qianhui_yu/Work/Lung/Data/Lung_data_master.RDS")
#Master <- MergeSeurat(object1 = vivo, object2 = D145tra)


# pre-process
#mito.genes <- grep(pattern = "^MT-", x = rownames(x = Master@data), value = TRUE)
#percent.mito <- Matrix::colSums(Master@raw.data[mito.genes, ])/Matrix::colSums(Master@raw.data)
#mito.low.cells <- names(which(percent.mito<0.1))
#non.d59.cells <- rownames(Master@meta.data)[which(!grepl("D59", rownames(Master@meta.data)))]
#selected.cells <- intersect(mito.low.cells, non.d59.cells)

#raw.data <- Master@raw.data
#genes <- read.table("~/Work/Annotation/Table_gencode.v22.annotation.protein_coding_gene_name.txt",stringsAsFactors=F,head=T)[,1]
#raw.data <- raw.data[rownames(raw.data)%in%genes, selected.cells]
#Master <- CreateSeuratObject(raw.data=raw.data)
#Master@meta.data <- cbind(Master@meta.data, percent.mito[selected.cells])
#colnames(Master@meta.data)[ncol(Master@meta.data)] <- "percent.mito"
#Master.bk <- Master

#raw.data <- Master@raw.data
genes <- read.table("/mnt/black/scRNA-seq/gencode.v22.annotation.protein_coding_gene_name.txt",stringsAsFactors=F,head=T)
#raw.data <- raw.data[rownames(raw.data)%in%genes, selected.cells]
#Master <- CreateSeuratObject(raw.data=raw.data)
#Master@meta.data <- cbind(Master@meta.data, percent.mito[selected.cells])
#colnames(Master@meta.data)[ncol(Master@meta.data)] <- "percent.mito"
#Master.bk <- Master
#Master <- SubsetData(Master, subset.name = "nGene", accept.high = 8000, accept.low = 1500)
#Master <- NormalizeData(object = Master)

#Master <- ScaleData(object=Master, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
#Master <- FindVariableGenes(object=Master, do.plot=FALSE)
#saveRDS(Master, file="Res_Master_1500.rds")

Master <- readRDS("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/Res_Master_1500.rds")

Master <- SetAllIdent(Master, id = "orig.ident")

setwd("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal")

# identify variable genes in each sample separately
samples <- as.character(unique(Master@meta.data$orig.ident))
seu.by.sample <- list()
for(idx in samples){
  cat(paste(idx, "start\n"))
  work.dir <- paste0(idx, "_2/")
  dir.create(file.path(work.dir))
  
  seu.obj <- SubsetData(Master, ident.use = idx, do.clean=T)
  seu.obj <- ScaleData(object=seu.obj, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
  seu.obj <- FindVariableGenes(object=seu.obj, do.plot=FALSE)
  
  
  seu.obj <- RunPCA(object=seu.obj, pc.genes=seu.obj@var.genes, do.print=FALSE)
  seu.obj <- FindClusters(object=seu.obj, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE)
  seu.obj <- RunTSNE(object=seu.obj, dims.use=1:20, do.fast=TRUE)
  
  p2 <- TSNEPlot(object=seu.obj, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
  p3 <- TSNEPlot(object=seu.obj, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
  pdf(paste0(work.dir,"Plot_tSNE_res0.8-2.pdf"),width=24, height=10)
  print(plot_grid(p2, p3, ncol=2, nrow=1))
  dev.off()
  
  # plot canonical cell type marker expression on tSNE
  markers <- c("PDGFRB", "PDGFRA", "DES", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "FOXF1", "VIM", "ACTA2", "COL1A1", 
               "EPCAM", "CDH1", "SOX2", "SOX9", "TP63", "KRT5", "KRT15", "F3", "EGFR", "IL33", "CHGA", "ELAVL4", 
               "NGFR", "SCGB1A1", "SCGB3A2", "BPIFA1", "MUC5AC", "MUC5B", "FOXJ1", "CFTR", "SFTPC", "SFTPB", 
               "ABCA3", "PDPN", "HOPX", "AGER")
  markers <- markers[markers%in%rownames(seu.obj@data)]
  pdf(paste0(work.dir,"Plot_marker_expression_on_tSNE_hvg.pdf"), width=12, height=(ceiling(length(markers) / 4) * 3))
  FeaturePlot(object=seu.obj, features.plot=markers, cols.use=c("navy", "darkorange1"), no.legend=TRUE, nCol=4, pt.size = 0.25)
  dev.off()
  
  saveRDS(seu.obj, file=paste0(work.dir,"Res_",idx,".rds"))
  
  seu.by.sample[[idx]] <- seu.obj
}

hvg.gene.idx <- matrix(F, nrow=nrow(Master@data), ncol=length(samples))
rownames(hvg.gene.idx) <- rownames(Master@data)
colnames(hvg.gene.idx) <- samples
for(j in seq(length(samples))){
  hvg.gene.idx[which(rownames(hvg.gene.idx)%in%seu.by.sample[[j]]@var.genes),j] <- T
}
# selected highly variable genes should be hvg in at least 2 samples
selected.hvg <- rownames(hvg.gene.idx)[apply(hvg.gene.idx, 1, sum)>1]
save(hvg.gene.idx, selected.hvg, file="Res_hvg_gene_by_sample.rdata")

Master <- RunPCA(object=Master, pc.genes=selected.hvg, do.print=FALSE)
Master <- FindClusters(object=Master, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE)
Master <- RunTSNE(object=Master, dims.use=1:20, do.fast=TRUE)

p2 <- TSNEPlot(object=Master, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=Master, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE_res0.8.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(Master, file="Res_Master_1500.rds")

Master.hvg <- Master
Master.hvg <- RunPCA(object=Master.hvg, pc.genes=Master.hvg@var.genes, do.print=FALSE)
Master.hvg <- FindClusters(object=Master.hvg, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=T)
Master.hvg <- RunTSNE(object=Master.hvg, dims.use=1:20, do.fast=TRUE)

p2 <- TSNEPlot(object=Master.hvg, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=Master.hvg, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE_hvg_res0.8.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(Master.hvg, file="Res_Master_hvg.rds")

# plot canonical cell type marker expression on tSNE
markers <- c("PDGFRB", "PDGFRA", "DES", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "FOXF1", "VIM", "ACTA2", "COL1A1", 
             "EPCAM", "CDH1", "SOX2", "SOX9", "TP63", "KRT5", "KRT15", "F3", "EGFR", "IL33", "CHGA", "ELAVL4", 
             "NGFR", "SCGB1A1", "SCGB3A2", "BPIFA1", "MUC5AC", "MUC5B", "FOXJ1", "CFTR", "SFTPC", "SFTPB", 
             "ABCA3", "PDPN", "HOPX", "AGER")
markers <- markers[markers%in%rownames(Master.hvg@data)]
g1 <- markers
png("Plot_marker_expression_on_tSNE_hvg.png", width=1600, height=4800)
FeaturePlot(object=Master.hvg, features.plot=g1, cols.use=c("navy", "darkorange1"), no.legend=F)
dev.off()
png("Plot_marker_expression_on_tSNE.png", width=1600, height=4800)
FeaturePlot(object=Master, features.plot=g1, cols.use=c("navy", "darkorange1"), no.legend=F)
dev.off()

#####################################################
####### choose Master for downstream analysis #######
#####################################################

# extract the EpCAM+ clusters and do the subclustering
epi.cluster <- c(19, 8, 16, 15, 23, 13, 24, 2, 22, 3)
epi.idx <- rep(F, nrow(Master@meta.data))
epi.idx[which(Master@meta.data$res.0.8%in%epi.cluster)] <- T

epi <- SubsetData(Master, cells.use=Master@cell.names[which(epi.idx)])
epi <- FindVariableGenes(object=epi, do.plot=F)
epi <- RunPCA(object=epi, pc.genes=epi@var.genes, do.print=FALSE)
epi <- FindClusters(object=epi, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
epi <- RunTSNE(object=epi, dims.use=1:20, do.fast=TRUE)

p2 <- TSNEPlot(object=epi, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=epi, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE_epi_clusters_res0.8.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(epi, file="Res_epi.rds")

cluster.marker <- FindAllMarkers(object=epi, only.pos=T, min.pct=0.25, thresh.use=0.25)
saveRDS(cluster.marker, file="Res_Epi_cluster_marker.rds")

library(dplyr)
top <- cluster.marker %>% group_by(cluster) %>% top_n(5, avg_logFC)
png("Plot_EpCAM_positive_subclustering_top5_cluster_marker.png",height=1000,width=2000)
DoHeatmap(object = epi, genes.use = top$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

png("Plot_marker_expression_on_epi_tSNE.png", width=1600, height=4800)
FeaturePlot(object=epi, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
dev.off()

vivo.cl.ave.expr <- foreach(cl.idx=unique(cluster.marker$cluster), .combine='cbind', .multicombine=T)%dopar%{
  cat(paste("Cluster", cl.idx, "start\n"))
  cl.sample.idx <- which(epi@meta.data$res.0.8==cl.idx)
  e <- apply(epi@data[,cl.sample.idx],1,mean)
  return(e)
}
colnames(vivo.cl.ave.expr) <- paste0("Cluster", unique(cluster.marker$cluster))
saveRDS(vivo.cl.ave.expr, file="Dat_InVivo_cluster_ave_expr.rds")

# get highly variable genes in each sample separately
samples <- as.character(unique(epi@meta.data$orig.ident))
dir.create(file.path("epi_by_sample/"))
epi.by.sample <- list()
for(idx in samples){
  cat(paste(idx, "start\n"))
  work.dir <- paste0("epi_by_sample/", idx, "/")
  dir.create(file.path(work.dir))
  
  seu.obj <- SubsetData(epi, cells.use = rownames(epi@meta.data)[which(epi@meta.data$orig.ident==idx)], do.clean=T)
  seu.obj <- ScaleData(object=seu.obj, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
  seu.obj <- FindVariableGenes(object=seu.obj, do.plot=FALSE)
  
  
  seu.obj <- RunPCA(object=seu.obj, pc.genes=seu.obj@var.genes, do.print=FALSE)
  seu.obj <- FindClusters(object=seu.obj, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE)
  seu.obj <- RunTSNE(object=seu.obj, dims.use=1:20, do.fast=TRUE)
  
  p2 <- TSNEPlot(object=seu.obj, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
  p3 <- TSNEPlot(object=seu.obj, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
  pdf(paste0(work.dir,"Plot_tSNE_res0.8-2.pdf"),width=24, height=10)
  print(plot_grid(p2, p3, ncol=2, nrow=1))
  dev.off()
  
  # plot canonical cell type marker expression on tSNE
  markers <- c("PDGFRB", "PDGFRA", "DES", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "FOXF1", "VIM", "ACTA2", "COL1A1", 
               "EPCAM", "CDH1", "SOX2", "SOX9", "TP63", "KRT5", "KRT15", "F3", "EGFR", "IL33", "CHGA", "ELAVL4", 
               "NGFR", "SCGB1A1", "SCGB3A2", "BPIFA1", "MUC5AC", "MUC5B", "FOXJ1", "CFTR", "SFTPC", "SFTPB", 
               "ABCA3", "PDPN", "HOPX", "AGER")
  markers <- markers[markers%in%rownames(seu.obj@data)]
  g1 <- markers
  png(paste0(work.dir,"Plot_marker_expression_on_tSNE_hvg.png"), width=1600, height=4800)
  FeaturePlot(object=seu.obj, features.plot=g1, cols.use=c("navy", "darkorange1"), no.legend=F)
  dev.off()
  
  saveRDS(seu.obj, file=paste0(work.dir,"Res_",idx,".rds"))
  
  cluster.marker <- FindAllMarkers(object=seu.obj, only.pos=T, min.pct=0.25, thresh.use=0.25)
  saveRDS(cluster.marker, file=paste0(work.dir, "Res_Epi_cluster_marker.rds"))
  
  epi.by.sample[[idx]] <- seu.obj
}


##################################
## only focus on distal samples ##
##################################

setwd("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/distal")

d103dis <- readRDS("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/epi_by_sample/D103dis/Res_D103dis.rds")
d125dis <- readRDS("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/epi_by_sample/D125dis/Res_D125dis.rds")
distal <- MergeSeurat(object1 = d103dis, object2 = d125dis)
# perform scaling
distal <- ScaleData(object=distal, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
selected.hvg <- intersect(d103dis@var.genes, d125dis@var.genes)
distal <- RunPCA(object=distal, pc.genes=selected.hvg, do.print=FALSE)
distal <- FindClusters(object=distal, reduction.type="pca", dims.use=1:20, resolution=2, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
distal <- RunTSNE(object=distal, dims.use=1:20, do.fast=TRUE)

p2 <- TSNEPlot(object=distal, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=distal, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()

# do CCA
cca <- RunCCA(object=d103dis, object2=d125dis, num.cc=50, genes.use=selected.hvg)
pdf("Plot_CC_num.pdf", width=10, height=8)
MetageneBicorPlot(cca, grouping.var = "orig.ident", dims.eval = 1:50, display.progress = TRUE)
dev.off()
cca <- AlignSubspace(object=cca, reduction.type="cca", grouping.var="orig.ident", dims.align=1:15)
cca <- FindClusters(object=cca, reduction.type="cca.aligned", dims.use=1:15, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
cca <- RunTSNE(object=cca, reduction.use="cca.aligned", dims.use=1:15, do.fast=TRUE)

#colfunc <- colorRampPalette(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"))

#colfunc <- colorRampPalette(c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2"))

#colfunc <- colorRampPalette(c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"))

colfunc <- colorRampPalette(brewer.pal(9, "Set1"))
p1 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Set1")
colfunc <- colorRampPalette(brewer.pal(8, "Set2"))
p2 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Set2")
colfunc <- colorRampPalette(brewer.pal(12, "Set3"))
p3 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Set3")
colfunc <- colorRampPalette(brewer.pal(12, "Paired"))
p4 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Paired")
colfunc <- colorRampPalette(brewer.pal(8, "Dark2"))
p5 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Dark2")
colfunc <- colorRampPalette(brewer.pal(8, "Accent"))
p6 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Accent")
colfunc <- colorRampPalette(brewer.pal(11, "Spectral"))
p7 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, clabel.size = 6, ols.use=colfunc(nlevels(cca@ident)), plot.title = "Spectral")

leg1 <- get_legend(p1)
leg2 <- get_legend(p2)
leg3 <- get_legend(p3)
leg4 <- get_legend(p4)
leg5 <- get_legend(p5)
leg6 <- get_legend(p6)
leg7 <- get_legend(p7)

pdf("Plot_tSNE_res0.8_afterCCA_distal_palCompare.pdf", width=10, height=70)
plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), leg1, leg2, p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), leg3, leg4, p5 + theme(legend.position="none"), p6 + theme(legend.position="none"), leg5, leg6, p7 + theme(legend.position="none"), leg7, ncol=2, nrow=14)
dev.off()


p2 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=superman_icecream_palette)
p3 <- DimPlot(object = cca, reduction.use = "tsne", group.by = "orig.ident", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, cols.use=superman_icecream_palette)
leg2 <- get_legend(p2)
leg3 <- get_legend(p3)

pdf("Plot_tSNE_res0.8_afterCCA_distal_spectral.pdf", width=10, height=10)
plot_grid(p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), leg2, leg3, ncol=2, nrow=2)
dev.off()

#p2 <- TSNEPlot(object=cca, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
#p3 <- TSNEPlot(object=cca, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
#png("Plot_tSNE_res0.8_afterCCA.png",width=1000, height=500)
#plot_grid(p2, p3, ncol=2, nrow=1)
#dev.off()
#saveRDS(cca, file="Dat_cca.rds")

# plot the marker gene expression on tSNE
source("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/Script_functions.R")
#g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g1 <- c("SFTPC", "ID2","SOX9","SOX11","ETV4","E2F8","HOPX","AQP5","AGER", "PDPN","SFTPB","ABCA3","SCGB3A2","SCGB1A1","FOXJ1","TP63","KRT15","MUC5AC","MUC5B","CHGA","EPCAM","CDH1","HMGA2","RFX6")
g2 <- intersect(g1, rownames(cca@data))
tsne.coor <- cca@dr$tsne@cell.embeddings
pdf("Plot_tSNE_cellTypeMarkers-2_distal.pdf", height=20, width=12)
par(mfrow=c(ceiling(length(g2) / 4), 4))
for(x in g2){
  plotFeature2(coor=tsne.coor, values=cca@data[x,], nCols=c("navy", "darkorange1"), main=x, xlab="", ylab="", bty="n", xaxt="n", yaxt="n",zeroAsGray=T, cex= 0.5)
}
dev.off()


cluster.vec <- paste0("Cluster", as.numeric(cca@meta.data$res.0.8))
sample.cols <- c("#e7298a","#66a61e")
names(sample.cols) <- c("D103dis", "D125dis")
sample.vec <- as.character(cca@meta.data$orig.ident)
pdf("Plot_tSNE_samples.pdf", height=7,width=7)
plotFeature2(coor=tsne.coor, values=sample.vec, gCols=sample.cols, main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
legend("bottomright", legend=names(sample.cols), text.col=sample.cols, bty="n")
dev.off()
pdf("Plot_tSNE_cluster.pdf", height=7, width=7)
TSNEPlot(object=cca, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
dev.off()

# plot the relative expression levels of cluster 0
cm.0 <- FindMarkers(cca, ident.1=0,  min.pct=0.25, logfc.threshold=0.25)
cm.0.h <- rownames(cm.0)[tail(order(cm.0$avg_logFC),25)]
cm.0.l <- rownames(cm.0)[head(order(cm.0$avg_logFC),6)]

# get the average expression levels in each cluster
cluster.num <- max(as.numeric(unique(cca@meta.data$res.0.8)))
cluster.res <- cca@meta.data$res.0.8
cl.ave.expr <- sapply(0:cluster.num, function(i){ 
  cl.sample.idx <- which(cluster.res==i)
  e <- apply(cca@data[c(cm.0.h, cm.0.l),cl.sample.idx],1,mean)
  return(e)
})
colnames(cl.ave.expr) <- paste0("Cluster", 0:cluster.num)

colors <- colorRampPalette(c("midnightblue","white", "darkorange1"))(n=299)
library(gplots)
pdf("Plot_heatmap_cm0-6.pdf")
heatmap.2(cl.ave.expr,trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE)
dev.off()

distal.markers <- FindAllMarkers(object = cca, return.thresh = 0.01, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, thresh.use = 0.25, nthreads = 5)
write.csv(distal.markers, file = "All-markers_resolution-0.8_distal_markersByCluster.csv")




###################################
## only focus on trachea samples ##
###################################

setwd("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/trachea")

# HT-189
sample1.epi <- readRDS("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/epi_by_sample/D145tra/Res_D145tra.rds")
# HT-187
sample2.epi <- readRDS("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/epi_by_sample/D103tra/Res_D103tra.rds")
# HT-182
sample3.epi  <- readRDS("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/epi_by_sample/D125prox/Res_D125prox.rds")

sample1.epi@meta.data[, "BloodType"] <- "RBC"
sample2.epi@meta.data[, "BloodType"] <- "RBC"
sample3.epi@meta.data[, "BloodType"] <- "Bloodless"

## Combine epithelial only seurat object subsets and regress out batch effects
sample1.2.epi <- MergeSeurat(object1 = sample1.epi, object2 = sample2.epi)
sample1.2.epi <- ScaleData(object = sample1.2.epi, vars.to.regress = "orig.ident")

sample1.2.3.epi <- MergeSeurat(object = sample1.2.epi, object2 = sample3.epi)
sample1.2.3.epi <- ScaleData(object = sample1.2.3.epi, vars.to.regress = "BloodType")

## Scale the bloodless dataset in prep for CCA
sample3.epi <- ScaleData(object = sample3.epi)

## Find highest variance genes in the RBS and Bloodless objects 
sample1.2.epi <-  FindVariableGenes(object = sample1.2.epi, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.8)
sample3.epi <-  FindVariableGenes(object = sample3.epi, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.8)

hvg.sample1.2.epi <- rownames(x = head(x = sample1.2.epi@hvg.info, n = 5000))
hvg.sample3.epi <- rownames(x = head(x = sample3.epi@hvg.info, n = 5000))
hvg.union <- union(x = hvg.sample1.2.epi, y = hvg.sample3.epi)

## Get a vector of all expressed genes shared across all seurat objects to remove non-interescting genes
allGenesIntersected <- intersect(x = rownames(sample1.2.epi@data), y = rownames(sample3.epi@data))
corrected.hvg.union <- intersect(x = allGenesIntersected, y = hvg.union)

## Run CCA and combine RBC and Bloodless objects
combTrachEpi <- RunCCA(object = sample1.2.epi, object2 = sample3.epi, genes.use = corrected.hvg.union, num.cc = 50)

## Find the best number of CCAs to use for alignment based on shared sources of variance by CCA
### After 24 CCs the shared variance sources drop off
pdf("MetageneBicorPlot_plot_combIntest.pdf", width=10, height=8)
MetageneBicorPlot(combTrachEpi, grouping.var = "BloodType", dims.eval = 1:50, display.progress = TRUE)
dev.off()

## Align samples by BloodType to remove effects of RBC mRNA contamination as a source of meaningful variance
combTrachEpi <- AlignSubspace(object = combTrachEpi, reduction.type = "cca", grouping.var = "BloodType", dims.align = 1:24)
combTrachEpi <- ScaleData(object = combTrachEpi, vars.to.regress = c("orig.ident"))

## Identify the most variable genes between cells
combTrachEpi <- FindVariableGenes(object = combTrachEpi,
                                  mean.function = ExpMean,
                                  dispersion.function = LogVMR,
                                  x.low.cutoff = 0.0125,
                                  x.high.cutoff = 3,
                                  y.cutoff = 0.8)
## Compute PCs in the nrely scaled object
combTrachEpi <- RunPCA(object = combTrachEpi, 
                       pc.genes = combTrachEpi@var.genes, 
                       do.print = FALSE,
                       pcs.compute = 100)
## Find cluster assignments for all cells based on cca-aligned expression data
combTrachEpi <- FindClusters(object = combTrachEpi, 
                             reduction.type = "cca.aligned", 
                             dims.use = 1:24, 
                             resolution = 0.6, 
                             save.SNN = TRUE,
                             force.recalc = TRUE)
## Run dimensional reduction on cca-aligned expression data
combTrachEpi <- RunTSNE(object = combTrachEpi, 
                        reduction.use = "cca.aligned", 
                        dims.use = 1:24, 
                        do.fast = TRUE, 
                        dim.embed = 2, 
                        max_iter = 2000, 
                        nthreads = 5, 
                        overwrite = TRUE)
## Reorder the cluster numbers so that the most similar clusters have sequencial numbers
pdf("BuildClusterMarker_plot_combTrachEpi.pdf", width=10, height=8)
combTrachEpi <- BuildClusterTree(combTrachEpi, do.reorder = TRUE, reorder.numeric = TRUE, pcs.use = 1:75)
dev.off()


p2 <- TSNEPlot(object=sample1.2.3.epi, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=sample1.2.3.epi, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()

# do CCA
#cca <- RunCCA(object=d103dis, object2=d125dis, num.cc=50, genes.use=selected.hvg)
#pdf("Plot_CC_num.pdf", width=10, height=8)
#MetageneBicorPlot(cca, grouping.var = "orig.ident", dims.eval = 1:50, display.progress = TRUE)
#dev.off()
#cca <- AlignSubspace(object=cca, reduction.type="cca", grouping.var="orig.ident", dims.align=1:15)
#cca <- FindClusters(object=cca, reduction.type="cca.aligned", dims.use=1:15, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
#cca <- RunTSNE(object=cca, reduction.use="cca.aligned", dims.use=1:15, do.fast=TRUE)

cca <- combTrachEpi

colfunc <- colorRampPalette(brewer.pal(11, "Spectral"))

colfunc <- colorRampPalette(brewer.pal(9, "Set1"))
p1 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Set1")
colfunc <- colorRampPalette(brewer.pal(8, "Set2"))
p2 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Set2")
colfunc <- colorRampPalette(brewer.pal(12, "Set3"))
p3 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Set3")
colfunc <- colorRampPalette(brewer.pal(12, "Paired"))
p4 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Paired")
colfunc <- colorRampPalette(brewer.pal(8, "Dark2"))
p5 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Dark2")
colfunc <- colorRampPalette(brewer.pal(8, "Accent"))
p6 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Accent")
colfunc <- colorRampPalette(brewer.pal(11, "Spectral"))
p7 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=colfunc(nlevels(cca@ident)), plot.title = "Spectral")

leg1 <- get_legend(p1)
leg2 <- get_legend(p2)
leg3 <- get_legend(p3)
leg4 <- get_legend(p4)
leg5 <- get_legend(p5)
leg6 <- get_legend(p6)
leg7 <- get_legend(p7)

pdf("Plot_tSNE_res0.8_afterCCA_trachea_palCompare.pdf", width=10, height=70)
plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), leg1, leg2, p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), leg3, leg4, p5 + theme(legend.position="none"), p6 + theme(legend.position="none"), leg5, leg6, p7 + theme(legend.position="none"), leg7, ncol=2, nrow=14)
dev.off()

p2 <- DimPlot(object = cca, reduction.use = "tsne", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, do.label = T, label.size = 6, cols.use=superman_icecream_palette)
p3 <- DimPlot(object = cca, reduction.use = "tsne", group.by = "orig.ident", do.return = TRUE, pt.size = 0.75, no.legend = F, no.axes = T, cols.use=superman_icecream_palette)
leg2 <- get_legend(p2)
leg3 <- get_legend(p3)

pdf("Plot_tSNE_res0.8_afterCCA_trachea_spectral.pdf", width=10, height=10)
plot_grid(p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), leg2, leg3, ncol=2, nrow=2)
dev.off()

saveRDS(cca, file="Dat_cca.rds")



# plot the marker gene expression on tSNE
source("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/Script_functions.R")
#g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g1 <- c("SFTPC", "ID2","SOX9","SOX11","ETV4","E2F8","HOPX","AQP5","AGER", "PDPN","SFTPB","ABCA3","SCGB3A2","SCGB1A1","FOXJ1","TP63","KRT15","MUC5AC","MUC5B","CHGA","EPCAM","CDH1","HMGA2","RFX6")
g2 <- intersect(g1, rownames(cca@data))
tsne.coor <- cca@dr$tsne@cell.embeddings
pdf("Plot_tSNE_cellTypeMarkers-2_trachea.pdf", height=20, width=12)
par(mfrow=c(ceiling(length(g2) / 4), 4))
for(x in g2){
  plotFeature2(coor=tsne.coor, values=cca@data[x,], nCols=c("navy", "darkorange1"), main=x, xlab="", ylab="", bty="n", xaxt="n", yaxt="n",zeroAsGray=T, cex= 0.5)
}
dev.off()


# plot the marker gene expression on tSNE
source("/mnt/black/scRNA-seq/newSeurat/1500_Proximal_distal/Script_functions.R")
#g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g1 <- c("EPCAM", "CDH1", "SOX2", "SOX9", "TP63", "KRT5", "KRT15", "F3", "EGFR", "IL33", "PDPN", "S100A2", "SFTPC", "HOPX", "AQP5", "AGER", "ABCA3", "SFTPB", "SCGB3A2", "CFTR", "SCGB1A1", "MUC5AC", "MUC5B", "FOXJ1", "CHGA")
g2 <- intersect(g1, rownames(cca@data))
tsne.coor <- cca@dr$tsne@cell.embeddings
pdf("Plot_tSNE_cellTypeMarkers-2_trachea-2.pdf", height=20, width=12)
par(mfrow=c(ceiling(length(g2) / 4), 4))
for(x in g2){
  plotFeature2(coor=tsne.coor, values=cca@data[x,], nCols=c("navy", "darkorange1"), main=x, xlab="", ylab="", bty="n", xaxt="n", yaxt="n",zeroAsGray=T, cex= 0.5)
}
dev.off()


g1 <- c("EPCAM", "CDH1", "SOX2", "SOX9", "TP63", "KRT5", "KRT15", "F3", "EGFR", "IL33", "PDPN", "S100A2", "SFTPC", "HOPX", "AQP5", "AGER", "ABCA3", "SFTPB", "SCGB3A2", "CFTR", "SCGB1A1", "MUC5AC", "MUC5B", "FOXJ1", "CHGA", "TOP2A", "MKI67", "CFTR")
g2 <- intersect(g1, rownames(combTrachEpi@data))
tsne.coor <- combTrachEpi@dr$tsne@cell.embeddings
pdf("Plot_tSNE_cellTypeMarkers-2_trachea-3.pdf", height=20, width=12)
par(mfrow=c(ceiling(length(g2) / 4), 4))
for(x in g2){
  plotFeature2(coor=tsne.coor, values=combTrachEpi@data[x,], nCols=c("navy", "darkorange1"), main=x, xlab="", ylab="", bty="n", xaxt="n", yaxt="n",zeroAsGray=T, cex= 0.5)
}
dev.off()




cluster.vec <- paste0("Cluster", as.numeric(cca@meta.data$res.0.8))
sample.cols <- c("#e7298a","#66a61e")
names(sample.cols) <- c("D145tra", "D103tra", "D125prox")
sample.vec <- as.character(cca@meta.data$orig.ident)
pdf("Plot_tSNE_samples.pdf", height=7,width=7)
plotFeature2(coor=tsne.coor, values=sample.vec, gCols=sample.cols, main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
legend("bottomright", legend=names(sample.cols), text.col=sample.cols, bty="n")
dev.off()
pdf("Plot_tSNE_cluster.pdf", height=7, width=7)
TSNEPlot(object=cca, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
dev.off()

# plot the relative expression levels of cluster 0
cm.13 <- FindMarkers(cca, ident.1=13,  min.pct=0.25, logfc.threshold=0.25)
cm.13.h <- rownames(cm.13)[tail(order(cm.13$avg_logFC),25)]
cm.13.l <- rownames(cm.13)[head(order(cm.13$avg_logFC),6)]

combTrachEpi.aveExp <- AverageExpression(cca, genes.use = c(cm.13.h, cm.13.l))

colnames(combTrachEpi.aveExp) <- paste0("Cluster", colnames(combTrachEpi.aveExp))

colors <- colorRampPalette(c("midnightblue","white", "darkorange1"))(n=299)
library(gplots)
pdf("Plot_heatmap_cm0-6.pdf")
heatmap.2(as.matrix(combTrachEpi.aveExp),trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE)
dev.off()

combTrachEpi.markers <- FindAllMarkers(object = combTrachEpi, return.thresh = 0.01, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, thresh.use = 0.25, nthreads = 5)
write.csv(combTrachEpi.markers, file = "All-markers_resolution-0.6_trachea_markersByCluster.csv")
