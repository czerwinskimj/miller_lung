library(Seurat)
library(dplyr)
source("~/Work/commonScript/Script_functions.R")
D145tra.data <- Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2330_Miller_HT189-Tracheal-Epithelium/Sample_HT-189-Tracheal-Epi/outs/HT189-TrachealEpi-filtered_gene_bc_matrices/hg19/")
colnames(D145tra.data) <- paste("D145tra",colnames(D145tra.data),sep = "_")
D145tra <- CreateSeuratObject(raw.data = D145tra.data, min.cells = 3, min.genes = 200)

vivo <- readRDS("/home/qianhui_yu/Work/Lung/Data/Lung_data_master.RDS")
Master <- MergeSeurat(object1 = vivo, object2 = D145tra)


# pre-process
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Master@data), value = TRUE)
percent.mito <- Matrix::colSums(Master@raw.data[mito.genes, ])/Matrix::colSums(Master@raw.data)
mito.low.cells <- names(which(percent.mito<0.1))
non.d59.cells <- rownames(Master@meta.data)[which(!grepl("D59", rownames(Master@meta.data)))]
selected.cells <- intersect(mito.low.cells, non.d59.cells)

raw.data <- Master@raw.data
genes <- read.table("~/Work/Annotation/Table_gencode.v22.annotation.protein_coding_gene_name.txt",stringsAsFactors=F,head=T)[,1]
raw.data <- raw.data[rownames(raw.data)%in%genes, selected.cells]
Master <- CreateSeuratObject(raw.data=raw.data)
Master@meta.data <- cbind(Master@meta.data, percent.mito[selected.cells])
colnames(Master@meta.data)[ncol(Master@meta.data)] <- "percent.mito"
Master.bk <- Master

raw.data <- Master@raw.data
genes <- read.table("~/Work/Annotation/Table_gencode.v22.annotation.protein_coding_gene_name.txt",stringsAsFactors=F,head=T)[,1]
raw.data <- raw.data[rownames(raw.data)%in%genes, selected.cells]
Master <- CreateSeuratObject(raw.data=raw.data)
Master@meta.data <- cbind(Master@meta.data, percent.mito[selected.cells])
colnames(Master@meta.data)[ncol(Master@meta.data)] <- "percent.mito"
Master.bk <- Master
Master <- SubsetData(Master, subset.name = "nGene", accept.high = 8000, accept.low = 1500)
Master <- NormalizeData(object = Master)

Master <- ScaleData(object=Master, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
Master <- FindVariableGenes(object=Master, do.plot=FALSE)
saveRDS(Master, file="Res_Master_1500.rds")

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
	png(paste0(work.dir,"Plot_tSNE_res0.8-2.png"),width=1000, height=500)
	print(plot_grid(p2, p3, ncol=2, nrow=1))
	dev.off()
		
	# plot canonical cell type marker expression on tSNE
	markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_major_cell_type_markers_from_literature_search.txt",sep="\t",stringsAsFactors=F)
	markers <- markers[markers[,2]%in%rownames(seu.obj@data),]
	g1 <- markers[,2]
	png(paste0(work.dir,"Plot_marker_expression_on_tSNE_hvg.png"), width=1600, height=4800)
	FeaturePlot(object=seu.obj, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
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
markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_major_cell_type_markers_from_literature_search.txt",sep="\t",stringsAsFactors=F)
markers <- markers[markers[,2]%in%rownames(Master.hvg@data),]
g1 <- markers[,2]
png("Plot_marker_expression_on_tSNE_hvg.png", width=1600, height=4800)
FeaturePlot(object=Master.hvg, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
dev.off()
png("Plot_marker_expression_on_tSNE.png", width=1600, height=4800)
FeaturePlot(object=Master, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
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
for(idx in samples[-1]){
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
	png(paste0(work.dir,"Plot_tSNE_res0.8.png"),width=1000, height=500)
	print(plot_grid(p2, p3, ncol=2, nrow=1))
	dev.off()
		
	# plot canonical cell type marker expression on tSNE
	markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_major_cell_type_markers_from_literature_search.txt",sep="\t",stringsAsFactors=F)
	markers <- markers[markers[,2]%in%rownames(seu.obj@data),]
	g1 <- markers[,2]
	png(paste0(work.dir,"Plot_marker_expression_on_tSNE_hvg.png"), width=1600, height=4800)
	FeaturePlot(object=seu.obj, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
	dev.off()

	saveRDS(seu.obj, file=paste0(work.dir,"Res_",idx,".rds"))
	
	cluster.marker <- FindAllMarkers(object=seu.obj, only.pos=T, min.pct=0.25, thresh.use=0.25)
	saveRDS(cluster.marker, file=paste0(work.dir, "Res_Epi_cluster_marker.rds"))
	
	epi.by.sample[[idx]] <- seu.obj
}
