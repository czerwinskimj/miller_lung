library(Seurat)
library(dplyr)
library(phateR)
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

# prepare matrix for SPRING
# filter out EpCAM- clusters
selected.sample <- rownames(epi@meta.data)[!epi@meta.data$res.0.8%in%c(16,18,19)]

# merge clusters
cl.res <- epi@meta.data$res.2
merged.cl <- rep(0, length(cl.res))
merged.cl[which(cl.res%in%c(3,9,8,15))] <- 1
merged.cl[which(cl.res==11)] <- 2
merged.cl[which(cl.res==26)] <- 3
merged.cl[which(cl.res%in%c(2,20,22,27))] <- 4
merged.cl[which(cl.res==19)] <- 5
merged.cl[which(cl.res%in%c(1,6))] <- 6
merged.cl[which(cl.res%in%c(23,25,28))] <- 0
merged.cl[which(cl.res==21)] <- 7
merged.cl[which(cl.res%in%c(4,5,16))] <- 8
merged.cl[which(cl.res%in%c(10,12,13))] <- 9
merged.cl[which(cl.res%in%c(7,14,18))] <- 10
merged.cl[which(cl.res==24)] <- 11
merged.cl[which(cl.res==0)] <- 12
merged.cl[which(cl.res==17)] <- 13
epi@meta.data <- cbind(epi@meta.data, merged.cl)
colnames(epi@meta.data)[ncol(epi@meta.data)] <- "merged.cl"
epi@ident <- as.factor(epi@meta.data$merged.cl)
names(epi@ident) <- epi@cell.names

registerDoParallel(13)
marker.by.cluster <- foreach(k=seq(13), .multicombine=T, .combine='list')%dopar%{
	mat <- FindMarkers(epi, ident.1=k, ident.2=seq(13)[-k], only.pos=T, min.pct=0.25, logfc.threshold=0.25)
	return(mat)
}

combined.cm <- foreach(k=seq(13), .multicombine=T, .combine='rbind')%dopar%{
	mat <- marker.by.cluster[[k]]
	dat <- data.frame(mat, "cluster"=rep(k, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
	rownames(dat) <- paste(rownames(dat), k, sep="_")
	return(dat)
}
saveRDS(combined.cm, file="epi_by_sample/Res_Epi_merged_cluster_marker.rds")

# prepare matrix for SPRING
# use cluster markers
top <- combined.cm %>% group_by(cluster) %>% top_n(50, avg_logFC)
selected.markers <- unique(top$gene_name)
selected.samples <- rownames(epi@meta.data)[which(epi@meta.data$merged.cl!=0)]
#spring.dir <- "SPRING_merged_cm_scaled/"
spring.dir <- "SPRING_merged_cm_scaled_regressed/"
dir.create(file.path(spring.dir))
#hvg.expr <- as.matrix(epi@data[selected.markers, selected.samples])
#hvg.expr <- t(scale(t(hvg.expr)))
hvg.expr <- as.matrix(epi@scale.data[selected.markers, selected.samples])
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(epi@meta.data[selected.samples,])
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
disMat <- 1-cor(hvg.expr)


spring.coor <- read.csv(paste0(spring.dir,"SPRING_merged_cm_scaled_regressed_k50/coordinates.txt"),head=F,row.names=1)
rownames(spring.coor) <- colnames(hvg.info)
# get the kNN network
k <- 15
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))

# plot the marker gene expression along kNN network
#g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g1 <- c("LTF", "TTYH1", "AZGP1", "CALML5", "SCGB3A1", "IL33", "ACTA2", "MUC5A", "MUC5B", "TFF2","SPDEF","CHGA","ASCL1","FOXI1","CFTR","BPIFA1","SCGB1A1","TSLP","IL25","DLCK1","ASCL2","RGS13")
g2 <- intersect(g1, rownames(epi@data))
expressed.idx <- which(apply(epi@data[g2,], 1, max)>0)
g2 <- g2[expressed.idx]
#spring.file <- "epi_by_sample/Plot_SPRING_merged_cm_scaled_regressed_k50-2.png"
#spring.file <- "epi_by_sample/Plot_SPRING_merged_cm_scaled_regressed_k50-2-bo.png"
spring.file <- "epi_by_sample/Plot_SPRING_merged_cm_scaled_regressed_k50-cl11-rg.png"
png(spring.file, width=4000, height=3000)
par(mfrow=c(3,4), mar=c(6,6,6,6))
for(x in g2){
	#plotFeature(coor=spring.coor, values=epi@data[x,selected.samples], nCols=c("navy", "darkorange1"), knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2)
	plotFeature(coor=spring.coor, values=epi@data[x,selected.samples], knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2)
}
dev.off()


spring.file <- "epi_by_sample/Plot_cluster_Info_along_SPRING_merged_cm_scaled_regressed_k50-noClIdx.png"
cluster.res=as.character(as.numeric(as.matrix(hvg.info["merged.cl",])))
png(spring.file, height=1000, width=1000)
par(mar=c(6,6,6,6))
plotFeature2(spring.coor, values=cluster.res, knn.pairs=knn.idx, xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=1.5, bty="n", xaxt="n", yaxt="n")
#for(i in unique(cluster.res)){
#	sample.idx <- colnames(hvg.info)[which(cluster.res==i)]
#	g.coor <- as.numeric(apply(spring.coor[sample.idx,], 2, mean))
#	text(g.coor[1], g.coor[2], labels=i, cex=1.5)
#}
dev.off()

#sample.cols <- c("#49A999","#197636","#A84798","#8BCCEC","#CA6778","#872555")
sample.cols <- c("#49A999","#197636","#332786","#A84798","#CA6778","#872555")
sample.cols <- c("#8BCCEC", "#49A999", "#197636", "#DCCB7C", "#CA6778", "#A84798")
names(sample.cols) <- c("D103dis", "D125dis", "D103air", "D125prox", "D103tra", "D145tra")
spring.file <- "epi_by_sample/Plot_sample_Info_along_SPRING_merged_cm_scaled_regressed_k50-4.png"
cluster.res=as.character(as.matrix(hvg.info["orig.ident",]))
png(spring.file, height=1000, width=1000)
par(mar=c(6,6,6,6))
plotFeature2(spring.coor, values=cluster.res, gCols=sample.cols, knn.pairs=knn.idx, xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=1.5, bty="n", xaxt="n", yaxt="n")
legend("bottomleft", legend=names(sample.cols), text.col=sample.cols, bty="n",cex=2)
dev.off()


spring.file <- "epi_by_sample/Plot_byCluster_Info_along_SPRING_merged_cm_scaled_regressed_k50.png"
png(spring.file, height=4000, width=4000)
par(mfrow=c(4,4), mar=c(6,6,6,6))
for(x in seq(13)){
	plotFeature(spring.coor, values=cluster.res, knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", emphasize=which(hvg.info["merged.cl",]==x))
}
dev.off()

spring.file <- "epi_by_sample/Plot_bupTips2basalCells.png"
png(spring.file, height=1000, width=1000)
par(mar=c(6,6,6,6))
plotFeature(spring.coor, values=cluster.res, knn.pairs=knn.idx, xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", emphasize=which(hvg.info["merged.cl",]==13 & spring.coor[,2]>=350 | hvg.info["merged.cl",]==9 | hvg.info["merged.cl",]==8 | hvg.info["merged.cl",]==10 | hvg.info["merged.cl",]==4))
dev.off()

# get the distance between selected cells and the mass center of the belonging cluster, and filter out those closer to clusters other than the cluster it belongs to
cluster.res=hvg.info["merged.cl",]
cluster.center.coor <- sapply(seq(13), function(i){
	sample.idx <- colnames(hvg.info)[which(cluster.res==i)]
	g.coor <- as.numeric(apply(spring.coor[sample.idx,], 2, mean))
})
selected.cells <- colnames(hvg.info)[which(hvg.info["merged.cl",]==13 & spring.coor[,2]>=350 | hvg.info["merged.cl",]==9 | hvg.info["merged.cl",]==8 | hvg.info["merged.cl",]==10 | hvg.info["merged.cl",]==4)]
cell.coor <- spring.coor[selected.cells,]
cellDis2CC <- sapply(seq(13), function(cl.idx){
	sqrt((cell.coor[,1]-cluster.center.coor[1,cl.idx])^2+(cell.coor[,2]-cluster.center.coor[2,cl.idx])^2)
})
rownames(cellDis2CC) <- selected.cells
devIdx <- rep(T, length(selected.cells))
names(devIdx) <- selected.cells
min.dis.cl <- apply(cellDis2CC[selected.cells,],1,which.min)
c1 <- selected.cells[which(min.dis.cl%in%c(4,8,9,10,13))]
devIdx[c1] <- F

c13 <- intersect(names(which(devIdx)), names(which(hvg.info["merged.cl",]==13)))
devIdx[c13[order(spring.coor[c13,2])[-c(1,2)]]] <- F
c4 <- intersect(names(which(devIdx)), names(which(hvg.info["merged.cl",]==4)))
devIdx[c4[order(spring.coor[c4,2])[-c(1,2,3)]]] <- F

spring.file <- "epi_by_sample/Plot_non_deviated_cells-2.png"
cluster.res=hvg.info["merged.cl",]
png(spring.file, height=1000, width=1000)
par(mar=c(6,6,6,6))
plotFeature2(spring.coor, values=cluster.res, knn.pairs=knn.idx, xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", emphasize=which(colnames(hvg.info)%in%names(which(!devIdx))))
for(i in c(4,8,9,10,13)){
	sample.idx <- colnames(hvg.info)[which(cluster.res==i)]
	g.coor <- as.numeric(apply(spring.coor[sample.idx,], 2, mean))
	text(g.coor[1], g.coor[2], labels=i, cex=1.5)
}
dev.off()

bl.cells <- names(which(!devIdx))
bl.expr <- hvg.expr[,bl.cells]
basal <- SubsetData(epi, cells.use=bl.cells)
basal <- FindVariableGenes(basal, do.plot=F)
bl.expr.hvg <- basal@scale.data[basal@var.genes,]
saveRDS(bl.expr.hvg, file="epi_by_sample/mergedCMScaledRegressed/Dat_basal_lineage_cell_hvg_expr.rds")
# perform diffusion map for extracted basal lineage cells
library(destiny)
dm <- DiffusionMap(t(bl.expr.hvg))
saveRDS(dm, file="Res_diffusion_map_basal_hvg.rds")
dm.coor <- dm@eigenvectors
rownames(dm.coor) <- colnames(bl.expr.hvg)

#color cells by clusters
source("~/Work/commonScript/Script_functions.R")
cluster.res <- basal@meta.data$merged.cl
names(cluster.res) <- basal@cell.names
png("epi_by_sample/mergedCMScaledRegressed/Plot_dm_hvg_cluster_info.png", width=1200, height=400)
par(mfrow=c(1,3))
plotFeature2(dm.coor[,1:2], values=cluster.res, xlab="DC1", ylab="DC2")
for(i in unique(cluster.res)){
	sample.idx <- names(which(cluster.res==i))
	g.coor <- as.numeric(apply(dm.coor[sample.idx,1:2], 2, mean))
	text(g.coor[1], g.coor[2], labels=i, cex=1.5)
}

plotFeature2(dm.coor[,c(1,3)], values=cluster.res, xlab="DC1", ylab="DC3")
for(i in unique(cluster.res)){
	sample.idx <- names(which(cluster.res==i))
	g.coor <- as.numeric(apply(dm.coor[sample.idx,c(1,3)], 2, mean))
	text(g.coor[1], g.coor[2], labels=i, cex=1.5)
}

plotFeature2(dm.coor[,c(2,3)], values=cluster.res, xlab="DC2", ylab="DC3")
for(i in unique(cluster.res)){
	sample.idx <- names(which(cluster.res==i))
	g.coor <- as.numeric(apply(dm.coor[sample.idx,c(2,3)], 2, mean))
	text(g.coor[1], g.coor[2], labels=i, cex=1.5)
}
dev.off()

time <- rank(-dm.coor[,1])
all.cell.time <- rep(0, length(selected.samples))
names(all.cell.time) <- selected.samples
all.cell.time[names(time)] <- time
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
spring.file <- "epi_by_sample/mergedCMScaledRegressed/Plot_SPRING_basalLineage_pseudoTime-2.png"
png(spring.file, height=1000, width=1000)
par(mar=c(6,6,6,6))
#plotFeature2(spring.coor, values=all.cell.time, xlab="", ylab="", main="pseudotime", xaxt="n", yaxt="n", cex.lab=1.3, cex.axis=1.3, cex.main=3, bty="n")
plotFeature2(spring.coor, values=all.cell.time, xlab="", ylab="", main="pseudotime", xaxt="n", yaxt="n", cex.lab=1.3, cex.axis=1.3, cex.main=3, bty="n")
for(i in 1:30) rect(600+3*(i-1),300,600+3*i,310, border=NA, col=colorPal(30)[i])
dev.off()

clusters <- c(9,8,10,13,4)
png("epi_by_sample/mergedCMScaledRegressed/Plot_pseudotime_and_cluster.png")
plot(time,rep(1,length(time)),type="n",ylim=c(1,5))
for(i in seq(length(clusters))){
	cl <- clusters[i]
	idx <- which(basal@meta.data$merged.cl==cl)
	points(time[idx], rep(6-i, length(idx)), pch=16)
}
dev.off()

time.vec <- basal@meta.data$time
stage.vec <- basal@meta.data$stage
#cols <- c("#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177")
#colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
#cols <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(1:5, breaks=30, right=F,include.lowest=T))]
gCols <- setNames(scales::hue_pal()(length(unique(hvg.info["merged.cl",]))), paste0("Cluster", unique(hvg.info["merged.cl",])))
cols <- gCols[paste0("Cluster",c(9,8,10,13,4))]
used.bw=100
used.lwd=5
pdf("epi_by_sample/mergedCMScaledRegressed/Plot_density_pseudotime_per_cluster-1.pdf")
plot(density(time.vec[which(stage.vec==1)],from=1,to=max(time.vec),bw=used.bw), lwd=used.lwd, ylim=c(0,0.002),col=cols[1], main="", xlab="Pseudotime",cex.lab=1.3)
for(i in 2:max(stage.vec)){
	lines(density(time.vec[which(stage.vec==i)],from=1,to=max(time.vec),bw=used.bw), lwd=used.lwd, col=cols[i])
}
legend("topleft", legend=c("Bud tip progenitors","Newly differentiating cells adjacent to bud tip","Hub progenitor cells","Basal cells","Differentiating basal cells"), text.col=cols, bty="n")
dev.off()

# hierarchical clustering for the expression profiles along pseudotime for highly variable genes
hc <- hclust(as.dist(1-cor(t(bl.expr.hvg))), method="ward.D2")
pdf("epi_by_sample/mergedCMScaledRegressed/Plot_hc_wardD2_basal_lineage_highly_variable_genes.pdf")
plot(hc)
abline(h=2.6)
dev.off()
g.num <- 6
g <- cutree(hc, k=g.num)
g.list <- lapply(seq(g.num), function(i){
	names(which(g==i))
})
g.size <- sapply(g.list, length)
g.mean <- matrix(NA, nrow=ncol(bl.expr.hvg), ncol=g.num)
g.sd <- matrix(NA, nrow=ncol(bl.expr.hvg), ncol=g.num)
for(i in seq(g.num)){
	g.mean[,i] <- apply(bl.expr.hvg[g.list[[i]],],2,mean)
	g.sd[,i] <- apply(bl.expr.hvg[g.list[[i]],],2,sd)
}

spline.time <- seq(from=1,to=max(time),length=50)
png("epi_by_sample/mergedCMScaledRegressed/Plot_basal_lineage_hvg_expr_cluster-6.png", height=800, width=1200)
par(mfrow=c(2,3))
for(i in seq(g.num)){
	y.max <- max(g.mean[,i]+g.sd[,i])
	y.min <- min(g.mean[,i]-g.sd[,i])
	up.bound <- g.mean[,i]+g.sd[,i]
	low.bound <- g.mean[,i]-g.sd[,i]
	plot(time, g.mean[,i], ylim=c(y.min, y.max), type="n", main=paste("Cluster",i, "(",g.size[i],")"),xlab="Pseudotime",ylab="Relative expression")
	polygon(c(1:max(time), max(time):1), c(up.bound[order(time)],low.bound[order(time,decreasing=T)]), border="NA", col="#a1d99b")
	lines(spline.time, predict(smooth.spline(time, g.mean[,i], df=8),spline.time)$y, lwd=5, col="#31a354")
	points(time, g.mean[,i], pch=16, col="#31a354")
}
dev.off()

# plot expression pattern of specific genes along pseudotime
pdf("epi_by_sample/mergedCMScaledRegressed/Plot_marker_gene_scaled_expression_along_basal_lineage_pseudotime.pdf", height=10, width=20)
par(mfrow=c(2,4))
for(x in c("SFTPC", "HOPX", "SFTPB", "SCGB3A2", "TP63", "KRT15", "KRT5")){
	plot(time, basal@scale.data[x,], pch=16, main=x, xlab="Pseudotime", ylab="Relative expression", col="#a1d99b50")
	lines(seq(from=1,to=max(time),length=50), predict(smooth.spline(time, basal@scale.data[x,], df=8), seq(from=1,to=max(time),length=50))$y, lwd=5, col="#31a354", cex.main=3, cex.lab=2)
}
dev.off()

pdf("epi_by_sample/mergedCMScaledRegressed/Plot_marker_gene_normed_expression_along_basal_lineage_pseudotime.pdf", height=10, width=20)
par(mfrow=c(2,4))
for(x in c("SFTPC", "HOPX", "SFTPB", "SCGB3A2", "TP63", "KRT15", "KRT5")){
	plot(time, basal@data[x,], pch=16, main=x, xlab="Pseudotime", ylab="Relative expression", col="#a1d99b50")
	lines(seq(from=1,to=max(time),length=50), predict(smooth.spline(time, basal@data[x,], df=8), seq(from=1,to=max(time),length=50))$y, lwd=5, col="#31a354", cex.main=3, cex.lab=2)
}
dev.off()

symbol2ID <- read.table("~/Work/Annotation/Table_10X_gene_symbols_and_ensemblID.txt",stringsAsFactors=F,sep="\t")
expressed.gene.id <- unique(symbol2ID[symbol2ID[,2]%in%rownames(basal@data),1])
write.table(expressed.gene.id, file="epi_by_sample/mergedCMScaledRegressed/List_expressed_gene_ID.txt",row.names=F,col.names=F,quote=F)
for(i in seq(length(g.list))){
	file.name <- paste0("epi_by_sample/mergedCMScaledRegressed/List_cluster_",i,"_gene.txt")
	cl.id <- unique(symbol2ID[symbol2ID[,2]%in%g.list[[i]],1])
	write.table(cl.id, file=file.name, row.names=F,col.names=F,quote=F)
}


cell.ident <- rep(NA, nrow(basal@meta.data))
cell.ident[which(basal@meta.data$merged.cl==9)] <- "Bud tip progenitors"
cell.ident[which(basal@meta.data$merged.cl==8)] <- "Newly differentiating cells adjacent to bud tip"
cell.ident[which(basal@meta.data$merged.cl==10)] <- "Hub progenitor cells"
cell.ident[which(basal@meta.data$merged.cl==13)] <- "Basal cells"
cell.ident[which(basal@meta.data$merged.cl==4)] <- "Differentiating basal cells"
basal@meta.data <- cbind(basal@meta.data, cell.ident)
basal@meta.data <- cbind(basal@meta.data, time)

stage <- rep(NA, nrow(basal@meta.data))
stage[which(basal@meta.data$merged.cl==9)] <- 1
stage[which(basal@meta.data$merged.cl==8)] <- 2
stage[which(basal@meta.data$merged.cl==10)] <- 3
stage[which(basal@meta.data$merged.cl==13)] <- 4
stage[which(basal@meta.data$merged.cl==4)] <- 5
basal@meta.data <- cbind(basal@meta.data, stage)

gCols <- setNames(scales::hue_pal()(length(unique(hvg.info["merged.cl",]))), paste0("Cluster",unique(hvg.info["merged.cl",])))
cols <- gCols[paste0("Cluster", c(9,8,10,13,4))]
pdf("epi_by_sample/mergedCMScaledRegressed/Plot_basal_lineage_pseudotime_and_cluster-2.pdf")
VlnPlot(object = basal, features.plot = "time", group.by="stage", cols.use=cols, point.size.use=0.5)
dev.off()

# enrichment of functional gene sets
gs1 <- read.table("~/Work/Annotation/MSigDB/List_HALLMARK_TGF_BETA_SIGNALING_geneset.txt",stringsAsFactors=F,sep="\t",skip=2)[,1]
gs2 <- read.table("~/Work/Annotation/MSigDB/List_SMAD4_Q6_geneset.txt",stringsAsFactors=F,sep="\t",skip=2)[,1]
function.gene.set <- list("tgfb1"=gs1, "smad4binding"=gs2)
smad.pval <- sapply(seq(length(function.gene.set)), function(j){
	sapply(seq(length(g.list)), function(i){
		d <- nrow(basal@data)
		c <- length(g.list[[i]])
		g0 <- intersect(function.gene.set[[j]], rownames(basal@data))
		b <- length(g0)
		a <- sum(g0%in%g.list[[i]])
		fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")$p.value
	})
})
colnames(smad.pval) <- c("HALLMARK_TGF_BETA_SIGNALING", "SMAD4_binding")
library(qusage)
kegg <- read.gmt("~/Work/Annotation/MSigDB/c2.all.v6.1.symbols.gmt")
d <- nrow(basal@data)
keggPval <- t(sapply(seq(length(kegg)), function(j){
	pathway.genes <- intersect(rownames(basal@data), kegg[[j]])
	c <- length(pathway.genes)
	sapply(seq(length(g.list)), function(i){
		cluster.genes <- g.list[[i]]
		overlap <- intersect(cluster.genes, pathway.genes)
		b <- length(cluster.genes)
		a <- length(overlap)
		if(a<round(b*0.05)){
			return(NA)
		}else{
			pval <- fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")$p.value
			return(pval)
		}
		
	})
}))
keggAdj <- apply(keggPval,2,p.adjust, method="bonferroni")
keggAdj <- matrix(p.adjust(as.vector(keggPval), method="bonferroni"), ncol=ncol(keggPval))
rownames(keggAdj) <- names(kegg)
apply(keggAdj<0.05, 2, sum, na.rm=T)


hallmark <- read.gmt("~/Work/Annotation/MSigDB/h.all.v6.1.symbols.gmt")
hallmarkPval <- t(sapply(seq(length(hallmark)), function(j){
	pathway.genes <- intersect(rownames(basal@data), hallmark[[j]])
	c <- length(pathway.genes)
	sapply(seq(length(g.list)), function(i){
		cluster.genes <- g.list[[i]]
		overlap <- intersect(cluster.genes, pathway.genes)
		b <- length(cluster.genes)
		a <- length(overlap)
		if(a<round(b*0.05)){
			return(NA)
		}else{
			pval <- fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")$p.value
			return(pval)
		}
		
	})
}))
rownames(hallmarkPval) <- names(hallmark)
hallmarkOR <- t(sapply(seq(length(hallmark)), function(j){
	pathway.genes <- intersect(rownames(basal@data), hallmark[[j]])
	c <- length(pathway.genes)
	sapply(seq(length(g.list)), function(i){
		cluster.genes <- g.list[[i]]
		overlap <- intersect(cluster.genes, pathway.genes)
		b <- length(cluster.genes)
		a <- length(overlap)
		if(a<round(b*0.05)){
			return(NA)
		}else{
			pval <- fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")$estimate
			return(pval)
		}
		
	})
}))
rownames(hallmarkOR) <- names(hallmark)
hallmarkAdj <- matrix(p.adjust(as.vector(hallmarkPval), method="bonferroni"), ncol=ncol(hallmarkPval))
rownames(hallmarkAdj) <- names(hallmark)
selected.gene.set <- unique(c(which(apply(hallmarkAdj<0.05,1,sum,na.rm=T)>0), which(hallmarkPval[,5]<0.05)))
usedMat <- rbind(hallmarkPval[selected.gene.set,],t(smad.pval))
vec <- as.vector(usedMat)
vec[which(is.na(vec))] <- 1
usedMat <- matrix(vec, ncol=ncol(usedMat))
input <- sqrt(-log10(usedMat))
rownames(input) <- c(names(hallmark)[selected.gene.set], colnames(smad.pval))

library(gplots)
colors <- colorRampPalette(c("#fde0dd","#49006a"), space="rgb")(50)
pdf("Plot_heatmap_functional_gene_set_enrichment.pdf")
heatmap.2(input,trace="none",density.info="none",scale="none",col=colors)
dev.off()

genes <- c("EPCAM","CDH1","VIM","DCN","PDPN","AGER","ASCL1","CHGA","SFTPC","SFTPB","SFTPD","SCGB1A1","SCGB3A2","TP63","KRT5","FOXJ1")
g.expr <- cbind(vivo.cl.ave.expr[genes, vivo.order], vitro.comb.expr[genes, vitro.order])
normed.expr <- t(apply(g.expr, 1, function(vec){
	(vec-min(vec))/(max(vec)-min(vec))
}))
pdf("Plot_cmOverlap_combinedVitro2Vivo_cellTypeMarker_expression_quantile-2.pdf")
heatmap.2(normed.expr,trace="none",density.info="none",dendrogram="none",scale="none",col=g.colors,Rowv=NA, Colv=NA)
dev.off()




spring.file <- "epi_by_sample/Plot_age_and_location_Info_along_SPRING_merged_cm_scaled_regressed_k50.png"
age.idx <- rep(NA, ncol(hvg.info))
age.idx[grep("D145", hvg.info["orig.ident",])] <- "Day145"
age.idx[grep("D103", hvg.info["orig.ident",])] <- "Day103"
age.idx[grep("D125", hvg.info["orig.ident",])] <- "Day125"
pos.idx <- rep("Proximal", ncol(hvg.info))
pos.idx[grep("dis", hvg.info["orig.ident",])] <- "Distal"
names(pos.idx) <- colnames(hvg.info)
names(age.idx) <- colnames(hvg.info)

png(spring.file, width=3000, height=3000)
par(mfrow=c(3,3), mar=c(6,6,6,6))
plotFeature(spring.coor, values=age.idx, knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Age")
legend("bottomright", text.col=age.cols, legend=unique(age.idx), bty="n", cex = 4)
plotFeature(spring.coor, values=pos.idx, knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Position")
legend("bottomright", text.col=pos.cols, legend=unique(pos.idx), bty="n", cex = 
4)
cplotFeature(spring.coor, values=hvg.info["orig.ident",], knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Sample")
legend("bottomright", text.col=sample.cols, legend=unique(hvg.info["orig.ident",]), bty="n", cex = 4)
for(x in unique(hvg.info["orig.ident",]))
	plotFeature(spring.coor, values=hvg.info["orig.ident",], knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main=x, emphasize=which(hvg.info["orig.ident",]==x))
dev.off()
