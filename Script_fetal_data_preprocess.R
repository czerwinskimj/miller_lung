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


c10.bp <- names(which(dm.coor[names(which(cluster.res=="10")),"DC3"]<0.01))
png("epi_by_sample/mergedCMScaledRegressed/Plot_test.png", width=1200, height=400)
par(mfrow=c(1,3))
plotFeature2(dm.coor[,1:2], values=cluster.res, xlab="DC1", ylab="DC2", emphasize=which(rownames(dm.coor)%in%c10.bp))
plotFeature2(dm.coor[,c(1,3)], values=cluster.res, xlab="DC1", ylab="DC3", emphasize=which(rownames(dm.coor)%in%c10.bp))
plotFeature2(dm.coor[,c(2,3)], values=cluster.res, xlab="DC2", ylab="DC3", emphasize=which(rownames(dm.coor)%in%c10.bp))
dev.off()

cluster.res=hvg.info["merged.cl",]
png("epi_by_sample/mergedCMScaledRegressed/Plot_SPRING_DC_c10_basal_lineage_primed.png", height=1000, width=1000)
par(mar=c(6,6,6,6))
plotFeature2(spring.coor, values=cluster.res, xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", emphasize=which(colnames(hvg.info)%in%c10.bp))
dev.off()

c10.bulb <- setdiff(intersect(names(which(cluster.res=="10")), rownames(spring.coor)[which(spring.coor[,1]>780)]), c10.bp)
cluster.res=hvg.info["merged.cl",]
png("epi_by_sample/mergedCMScaledRegressed/Plot_SPRING_c10_bulb-2.png", height=1000, width=1000)
par(mar=c(6,6,6,6))
plotFeature2(spring.coor, values=cluster.res, xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", emphasize=which(colnames(hvg.info)%in%c10.bulb))
dev.off()

cluster.res <- hvg.info["merged.cl", bl.cells]
png("epi_by_sample/mergedCMScaledRegressed/Plot_DM_SPRING_c10_bulb-2.png", width=1200, height=400)
par(mfrow=c(1,3))
plotFeature2(dm.coor[,1:2], values=cluster.res, xlab="DC1", ylab="DC2", emphasize=which(rownames(dm.coor)%in%c10.bulb))
plotFeature2(dm.coor[,c(1,3)], values=cluster.res, xlab="DC1", ylab="DC3", emphasize=which(rownames(dm.coor)%in%c10.bulb))
plotFeature2(dm.coor[,c(2,3)], values=cluster.res, xlab="DC2", ylab="DC3", emphasize=which(rownames(dm.coor)%in%c10.bulb))
dev.off()

# get the edge distributions between c10 and others
k <- 50
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))
c10 <- intersect(names(which(hvg.info["merged.cl",]==10)), bl.cells)
c10.idx <- which(rownames(disMat)%in%c10)
edge.num <- sapply(seq(13), function(x){
	cx <- intersect(names(which(hvg.info["merged.cl",]==x)), bl.cells)
	cx.idx <- which(rownames(disMat)%in%cx)
	sum((knn.idx[,1]%in%cx.idx & knn.idx[,2]%in%c10.idx) | (knn.idx[,1]%in%c10.idx & knn.idx[,2]%in%cx.idx))
})
names(edge.num) <- paste0("Cluster",seq(13))

#get the cells in cluster 10 linked to cluster 13 or vise vesa
c13 <- intersect(names(which(hvg.info["merged.cl",]==13)), bl.cells)
c13.idx <- which(rownames(disMat)%in%c13)
c10.1 <- rownames(disMat)[unique(knn.idx[which(knn.idx[,1]%in%c13.idx & knn.idx[,2]%in%c10.idx),2])]
c10.2 <- rownames(disMat)[unique(knn.idx[which(knn.idx[,1]%in%c10.idx & knn.idx[,2]%in%c13.idx),1])]
c10 <- union(c10.1, c10.2)

c13.1 <- rownames(disMat)[unique(knn.idx[which(knn.idx[,1]%in%c10.idx & knn.idx[,2]%in%c13.idx),2])]
c13.2 <- rownames(disMat)[unique(knn.idx[which(knn.idx[,1]%in%c13.idx & knn.idx[,2]%in%c10.idx),1])]
c13 <- union(c13.1, c13.2)

png("epi_by_sample/mergedCMScaledRegressed/Plot_SPRING_c10_cells_link2C13.png", height=500, width=1000)
par(mfrow=c(1,2))
plotFeature2(spring.coor, values=hvg.info["merged.cl",], xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", emphasize=which(colnames(hvg.info)%in%c10))
plotFeature2(spring.coor, values=hvg.info["merged.cl",], xlab="", ylab="", main="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", emphasize=which(colnames(hvg.info)%in%c13))
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



# construct cluster-level network
k <- 50
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))

knn.sample.pairs <- cbind(rownames(disMat)[knn.idx[,1]], rownames(disMat)[knn.idx[,2]])
cluster.link <- sapply(1:13, function(c0){
	sapply(1:13, function(c1){
		if(c0==c1){
			return(1)
		}else{
			c0.sample <- colnames(hvg.info)[which(hvg.info["merged.cl",]==c0)]
			c1.sample <- colnames(hvg.info)[which(hvg.info["merged.cl",]==c1)]
			other.sample <- setdiff(rownames(disMat), union(c0.sample, c1.sample))
			c0.to.c1.links <- sum(knn.sample.pairs[,2]%in%c0.sample & knn.sample.pairs[,1]%in%c1.sample)
			c0.to.other.links <- sum(knn.sample.pairs[,2]%in%c0.sample & knn.sample.pairs[,1]%in%other.sample)
			c1.size <- length(c1.sample)
			other.size <- length(other.sample)
			pval <- fisher.test(matrix(c(c0.to.c1.links, c0.to.other.links, c1.size, other.size),c(2,2)), alternative="g")$p.value
			return(pval)
		}
	})
})
rownames(cluster.link) <- paste0("Cluster", 1:13)
colnames(cluster.link) <- paste0("Cluster", 1:13)
pAdj <- cluster.link
diag(pAdj) <- NA
pAdj <- matrix(p.adjust(as.vector(pAdj), method="bonferroni"),ncol=ncol(pAdj))
diag(pAdj) <- 1
rownames(pAdj) <- paste0("Cluster", 1:13)
colnames(pAdj) <- paste0("Cluster", 1:13)

link.freq <- sapply(1:13, function(c0){
	sapply(1:13, function(c1){
		c0.sample <- colnames(hvg.info)[which(hvg.info["merged.cl",]==c0)]
		c1.sample <- colnames(hvg.info)[which(hvg.info["merged.cl",]==c1)]
		c0.to.c1.links <- sum(knn.sample.pairs[,2]%in%c0.sample & knn.sample.pairs[,1]%in%c1.sample)
		return(c0.to.c1.links)	
	})
})
rownames(link.freq) <- paste0("Cluster", 1:13)
colnames(link.freq) <- paste0("Cluster", 1:13)

cutoff <- 0.00001
adj <- "corrected"
if(adj=="corrected"){
	L <- pAdj<cutoff
}else{
	L <- cluster.link<cutoff
}

method <- "either"
if(method=="mutual"){
	usedL <- L & t(L)
}else{
	usedL <- L	
}

usedL <- matrix(as.numeric(usedL), nrow=nrow(usedL))
rownames(usedL) <- paste0("Cluster", 1:13)
colnames(usedL) <- paste0("Cluster", 1:13)
library(igraph)
g <- graph.adjacency(
	usedL,
	mode="undirected",
	weighted=TRUE,
	diag=FALSE
)
edgeweights <- E(g)$weight * 2.0
l1 <- layout_nicely(g)
#Plot the tree object
dir <- "/r1/people/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/epi_by_sample/mergedCMScaledRegressed/"
pdf(paste0(dir,"Plot_epi_vivo_cluster_link_",method,"_",adj,"_",cutoff,"_k",k,".pdf"))
plot(g, layout=l1, edge.width=edgeweights)
dev.off()
saveRDS(l1, file=paste0(dir,"Res_layout_",method,"_",adj,"_",cutoff,"_k",k,".rds"))




# try MNN
# check the average transcriptome similarity between samples
sample.ave.expr <- foreach(s.idx=unique(epi@meta.data$orig.ident), .combine='cbind', .multicombine=T)%dopar%{
	cat(paste(s.idx, "start\n"))
	sample.idx <- grep(s.idx, selected.samples, value=T)
	e <- apply(epi@data[, sample.idx],1,mean)
	return(e)
}
colnames(sample.ave.expr) <- unique(epi@meta.data$orig.ident)
hc.all <- hclust(as.dist(1-cor(sample.ave.expr)), method="ward.D2")
hc.cm <- hclust(as.dist(1-cor(sample.ave.expr[selected.markers,])), method="ward.D2")
pdf("epi_by_sample/Plot_ave_transcriptome_similarity_by_sample.pdf", height=7, width=10)
par(mfrow=c(1,2))
plot(hc.all)
plot(hc.cm)
dev.off()

library(scran)
blood.idx <- rep("Bloodless", length=nrow(epi@meta.data))
blood.idx[epi@meta.data$orig.ident%in%c("D103tra", "D145tra")] <- "Blood"
epi@meta.data <- cbind(epi@meta.data, blood.idx)
colnames(epi@meta.data)[ncol(epi@meta.data)] <- "blood.idx"

selected.samples <- rownames(epi@meta.data)[which(epi@meta.data$merged.cl!=0)]
comb.rbc <- SubsetData(epi, cells.use=intersect(epi@cell.names[which(blood.idx=="Blood")], selected.samples), do.clean=T)
comb.non.rbc <- SubsetData(epi, cells.use=intersect(epi@cell.names[which(blood.idx=="Bloodless")], selected.samples), do.clean=T)

comb.rbc <- ScaleData(object=comb.rbc, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
comb.rbc <- FindVariableGenes(object=comb.rbc, do.plot=FALSE)
comb.non.rbc <- ScaleData(object=comb.non.rbc, vars.to.regress=c("nUMI", "percent.mito"), do.par=T)
comb.non.rbc <- FindVariableGenes(object=comb.non.rbc, do.plot=FALSE)
saveRDS(comb.rbc, file="Dat_epi_comb_rbc.rds")
saveRDS(comb.non.rbc, file="Dat_epi_comb_non_rbc.rds")
hvg.sec <- intersect(comb.rbc@var.genes, comb.non.rbc@var.genes)
hvg.union <- union(comb.rbc@var.genes, comb.non.rbc@var.genes)

# try MNN
library(scran)
mnn.data <- mnnCorrect(comb.rbc@scale.data, comb.non.rbc@scale.data)
corrected.comb.raw <- do.call("cbind", mnn.data$corrected)
colnames(corrected.comb.raw) <- c(colnames(comb.rbc@scale.data), colnames(comb.non.rbc@scale.data))
scale.data <- t(scale(t(corrected.comb.raw)))

mnn <- MergeSeurat(object1=comb.rbc, object2=comb.non.rbc)
mnn@scale.data <- scale.data
mnn <- FindVariableGenes(object=mnn, do.plot=F)
#mnn <- RunPCA(object=mnn, pc.genes=mnn@var.genes, do.print=FALSE)
#mnn <- RunPCA(object=mnn, pc.genes=hvg.sec, do.print=FALSE)
mnn <- RunPCA(object=mnn, pc.genes=hvg.union, do.print=FALSE)
mnn <- FindClusters(object=mnn, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
mnn <- RunTSNE(object=mnn, dims.use=1:20, do.fast=TRUE)
p2 <- TSNEPlot(object=mnn, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=mnn, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE_res0.8_afterMNN_hvg_union.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(mnn, file="Dat_mnn_hvg_union.rds")

# plot the marker gene expression along kNN network
g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g2 <- intersect(g1, rownames(mnn@data))
png("Plot_marker_gene_expression_on_tSNE_afterMNN_hvg_union.png", height=1600, width=1600)
FeaturePlot(object=mnn, features.plot=g2, cols.use=c("gray", "blue"), no.legend=F)
dev.off()

cluster.num <- max(as.numeric(mnn@meta.data$res.0.8))
library(doParallel)
registerDoParallel(cluster.num+1)
combined.cm <- foreach(k=0:cluster.num, .multicombine=T, .combine='rbind')%dopar%{
	mat <- FindMarkers(mnn, ident.1=k, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
	dat <- data.frame(mat, "cluster"=rep(k, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
	rownames(dat) <- paste(rownames(dat), k, sep="_")
	return(dat)
}
stopImplicitCluster()

# prepare matrix for SPRING
# use cluster markers
library(dplyr)
top <- combined.cm %>% group_by(cluster) %>% top_n(50, avg_logFC)
selected.markers <- unique(top$gene_name)
spring.dir <- "../SPRING_mnn_cm/"
dir.create(file.path(spring.dir))
hvg.expr <- as.matrix(mnn@scale.data[selected.markers, ])
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(mnn@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
disMat <- 1-cor(hvg.expr)

# use PCs
spring.dir <- "../SPRING_mnn_pc/"
dir.create(file.path(spring.dir))
hvg.expr <- as.matrix(t(mnn@dr$pca@cell.embeddings))
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(mnn@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
disMat <- as.matrix(dist(t(hvg.expr)))

spring.coor <- read.csv(paste0(spring.dir,"SPRING_mnn_pc_k20/coordinates.txt"),head=F,row.names=1)
rownames(spring.coor) <- colnames(hvg.info)
# get the kNN network
k <- 15
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))

# plot the marker gene expression along kNN network
g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g2 <- intersect(g1, rownames(mnn@data))
spring.file <- "Plot_SPRING_mnn_pc_k20.png"
png(spring.file, width=4000, height=4000)
par(mfrow=c(4,4), mar=c(6,6,6,6))
for(x in g2){
	plotFeature(coor=spring.coor, values=mnn@scale.data[x,], knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2)
}
dev.off()

spring.file <- "Plot_age_and_location_Info_along_SPRING_mnn_pc_k20.png"
age.idx <- rep(NA, ncol(hvg.info))
age.idx[grep("D145", hvg.info["orig.ident",])] <- "Day145"
age.idx[grep("D103", hvg.info["orig.ident",])] <- "Day103"
age.idx[grep("D125", hvg.info["orig.ident",])] <- "Day125"
pos.idx <- rep("Proximal", ncol(hvg.info))
pos.idx[grep("dis", hvg.info["orig.ident",])] <- "Distal"
names(pos.idx) <- colnames(hvg.info)
names(age.idx) <- colnames(hvg.info)
age.cols <- c("#F8766D", "#00BA38", "#619CFF")
names(age.cols) <- c("Day125", "Day103", "Day145")
pos.cols <- c("#F8766D", "#00BFC4")
names(pos.cols) <- c("Distal", "Proximal")
sample.cols <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
names(sample.cols) <- c("D125dis", "D125prox", "D103dis", "D103air", "D103tra", "D145tra")

png(spring.file, width=3000, height=3000)
par(mfrow=c(3,3), mar=c(6,6,6,6))
plotFeature(spring.coor, values=age.idx, knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Age")
legend("bottomright", text.col=age.cols, legend=unique(age.idx), bty="n", cex = 4)
plotFeature(spring.coor, values=pos.idx, knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Position")
legend("bottomright", text.col=pos.cols, legend=unique(pos.idx), bty="n", cex = 4)
plotFeature(spring.coor, values=hvg.info["orig.ident",], knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Sample")
legend("bottomright", text.col=sample.cols, legend=unique(hvg.info["orig.ident",]), bty="n", cex = 4)
for(x in unique(hvg.info["orig.ident",]))
	plotFeature(spring.coor, values=hvg.info["orig.ident",], knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main=x, emphasize=which(hvg.info["orig.ident",]==x))
dev.off()


# try CCA
cca <- RunCCA(object=comb.rbc, object2=comb.non.rbc, num.cc=50, genes.use=hvg.sec)
pdf("epi_by_sample/Plot_CC_num.pdf", width=10, height=8)
MetageneBicorPlot(cca, grouping.var = "blood.idx", dims.eval = 1:50, display.progress = TRUE)
dev.off()
cca <- AlignSubspace(object=cca, reduction.type="cca", grouping.var="blood.idx", dims.align=1:15)
cca <- FindClusters(object=cca, reduction.type="cca.aligned", dims.use=1:15, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
cca <- RunTSNE(object=cca, reduction.use="cca.aligned", dims.use=1:15, do.fast=TRUE)

p2 <- TSNEPlot(object=cca, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=cca, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("epi_by_sample/Plot_tSNE_res0.8_afterCCA.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(cca, file="epi_by_sample/Dat_cca_hvg_sec.rds")

# plot the marker gene expression along kNN network
g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g2 <- intersect(g1, rownames(epi@data))
png("epi_by_sample/Plot_marker_gene_expression_on_tSNE_afterCCA.png", height=1600, width=1600)
FeaturePlot(object=cca, features.plot=g2, cols.use=c("gray", "blue"), no.legend=F)
dev.off()

cluster.num <- max(as.numeric(cca@meta.data$res.0.8))
registerDoParallel(cluster.num+1)
combined.cm <- foreach(k=0:cluster.num, .multicombine=T, .combine='rbind')%dopar%{
	mat <- FindMarkers(epi, ident.1=k, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
	dat <- data.frame(mat, "cluster"=rep(k, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
	rownames(dat) <- paste(rownames(dat), k, sep="_")
	return(dat)
}
stopImplicitCluster()

# prepare matrix for SPRING
# use cluster markers
top <- combined.cm %>% group_by(cluster) %>% top_n(50, avg_logFC)
selected.markers <- unique(top$gene_name)
spring.dir <- "SPRING_cca_cm/"
dir.create(file.path(spring.dir))
hvg.expr <- as.matrix(cca@scale.data[selected.markers, ])
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(cca@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
disMat <- 1-cor(hvg.expr)

# use PCs
spring.dir <- "SPRING_cca_pc/"
dir.create(file.path(spring.dir))
hvg.expr <- as.matrix(t(cca@dr$cca.aligned@cell.embeddings))
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(cca@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
disMat <- as.matrix(dist(t(hvg.expr)))

spring.coor <- read.csv(paste0(spring.dir,"SPRING_cca_k20/coordinates.txt"),head=F,row.names=1)
rownames(spring.coor) <- colnames(hvg.info)
# get the kNN network
k <- 15
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))

# plot the marker gene expression along kNN network
g1 <- c("PDPN", "AGER", "HOPX", "SFTPC", "SOX9", "SFTPB", "CFTR", "SFTPD", "FOXJ1", "SCGB1A1", "SCGB3A1", "SCGB3A2", "TP63", "KRT5", "KRT15", "ASCL1")
g2 <- intersect(g1, rownames(cca@data))
spring.file <- "epi_by_sample/Plot_SPRING_cca_pc_k20.png"
png(spring.file, width=4000, height=4000)
par(mfrow=c(4,4), mar=c(6,6,6,6))
for(x in g2){
	plotFeature(coor=spring.coor, values=cca@data[x,], knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2)
}
dev.off()

spring.file <- "epi_by_sample/Plot_age_and_location_Info_along_SPRING_cca_pc_k20.png"
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
legend("bottomright", text.col=pos.cols, legend=unique(pos.idx), bty="n", cex = 4)
plotFeature(spring.coor, values=hvg.info["orig.ident",], knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main="Sample")
legend("bottomright", text.col=sample.cols, legend=unique(hvg.info["orig.ident",]), bty="n", cex = 4)
for(x in unique(hvg.info["orig.ident",]))
	plotFeature(spring.coor, values=hvg.info["orig.ident",], knn.pairs=knn.idx, xlab="", ylab="", cex.lab=1.3, cex.axis=1.3, cex.main=5, bty="n", xaxt="n", yaxt="n", cex=1.2, main=x, emphasize=which(hvg.info["orig.ident",]==x))
dev.off()
