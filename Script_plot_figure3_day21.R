library(Seurat)
source("~/Work/commonScript/Script_functions.R")
library(doParallel)
library(gplots)

d21 <- readRDS("Res_Day21_B1.rds")

cl.ave.expr <- getAveExpr(d21@meta.data, feature.to.calc="res.0.8", expr=d21@data, core.num=20)

cell.type.markers <- c("SOX9", "SFTPC", "SCGB3A2", "SFTPB", "CFTR", "TP63", "KRT15", "LTF", "FOXJ1", "SPDEF", "MUC5AC", "CHGA","CDK1","MKI67")
cl.marker.expr <- cl.ave.expr[cell.type.markers,]
png("Plot_cannonical_cell_type_markers_across_clusters.png", height=2800, width=800)
par(mfrow=c(14,1))
for(i in cell.type.markers){
	barplot(cl.marker.expr[i,], main=i)
}
dev.off()

d21.cm <- findAllMarkers(seu.obj=d21, selected.column="res.0.8")
saveRDS(d21.cm, file="Res_Day21_cluster_marker_res0.8_20PC.rds")
library(dplyr)
top <- d21.cm %>% group_by(cluster) %>% top_n(3, avg_logFC)
top.cm <- top$gene_name
png("Plot_cluster_markers_across_clusters.png", height=6000, width=800)
par(mfrow=c(30,1))
for(i in top.cm){
	barplot(cl.ave.expr[i,], main=i)
}
dev.off()

vivo.cl.expr <- readRDS("/home/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vitro/only_D0_B1_1000_remove_D0_B2/Dat_in_vivo_merged_cluster_average_expr.rds")
in.vivo.cm <- readRDS("~/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/epi_by_sample/Res_Epi_merged_cluster_marker.rds")
top.vivo.cm <- read.table("/home/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/SPRING_merged_cm_scaled_regressed/List_top50_cluster_markers.txt",stringsAsFactors=F)[,1]
 
# calculate the transcriptome similarity between 
## use interset of in vivo cluster markers and in vitro highly variable genes
#g1 <- intersect(in.vivo.cm$gene_name, d21@var.genes)
## use interset of in vivo cluster markers and in vitro cluster markers
#g1 <- intersect(in.vivo.cm$gene_name, d21.cm$gene_name)
## use all in vivo cluster markers
#g1 <- unique(in.vivo.cm$gene_name)
## use top in vivo cluster markers
g1 <- top.vivo.cm

# 

#de.gene <- readRDS("~/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vitro/only_D0_B1_1000_remove_D0_B2/decompose/Res_SMAD_activation_related_gene.rds")
#g1 <- setdiff(top.vivo.cm, de.gene)

#cl.cols <- setNames(scales::hue_pal()(13), paste0("Cluster", seq(13)))
tsne.coor <- d21@dr$tsne@cell.embeddings

registerDoParallel(20)
#vitroCell2vivoCluster.cm.hvg <- foreach(k=seq(ncol(d21@data)), .multicombine=T, .combine='cbind')%dopar%{
#vitroCell2vivoCluster.cm.cm <- foreach(k=seq(ncol(d21@data)), .multicombine=T, .combine='cbind')%dopar%{
#vitroCell2vivoCluster.cm <- foreach(k=seq(ncol(d21@data)), .multicombine=T, .combine='cbind')%dopar%{
vitroCell2vivoCluster.top.cm <- foreach(k=seq(ncol(d21@data)), .multicombine=T, .combine='cbind')%dopar%{
	if(k%%500==0){cat(paste0(k,"start\n"))}
	cor.vec <- sapply(seq(ncol(vivo.cl.expr)), function(cl){
		cor(d21@data[g1,k], vivo.cl.expr[g1,cl], method="spearman")
	})
	return(cor.vec)
}
stopImplicitCluster()
saveRDS(vitroCell2vivoCluster.top.cm, file="Res_vitroCell2vivoCluster_topVivoClusterMarkers_SCC.rds")
saveRDS(vitroCell2vivoCluster.cm, file="Res_vitroCell2vivoCluster_VivoClusterMarkers.rds")
saveRDS(vitroCell2vivoCluster.cm.cm, file="Res_vitroCell2vivoCluster_vivo_and_vitro_ClusterMarkers.rds")
saveRDS(vitroCell2vivoCluster.cm.hvg, file="Res_vitroCell2vivoCluster_VivoClusterMarkers_vitroHvg.rds")

#bestVivoCluster2VItroCell.cm.hvg <- apply(vitroCell2vivoCluster.cm.hvg, 2, which.max)
#values=paste0("Cluster", bestVivoCluster2VItroCell.cm.hvg)

#bestVivoCluster2VItroCell.cm.cm <- apply(vitroCell2vivoCluster.cm.cm, 2, which.max)
#values=paste0("Cluster", bestVivoCluster2VItroCell.cm.cm)

#bestVivoCluster2VItroCell.cm <- apply(vitroCell2vivoCluster.cm, 2, which.max)
#values=paste0("Cluster", bestVivoCluster2VItroCell.cm)

#bestVivoCluster2VItroCell.top.cm <- apply(vitroCell2vivoCluster.top.cm, 2, which.max)
bestVivoCluster2VItroCell.top.cm <- apply(vitroCell2vivoCluster.top.cm, 2, function(vec){order(vec, decreasing=T)[2]})
values=paste0("Cluster", bestVivoCluster2VItroCell.top.cm)


# plot the coorelation distribution
vitro.cluster.res <- as.numeric(d21@meta.data$res.0.8)
#selected.cor.mat <- vitroCell2vivoCluster.top.cm
selected.cor.mat <- vitroCell2vivoCluster.cm.hvg
selected.cor.dist <- lapply(0:max(vitro.cluster.res), function(i){
	cl.sample.idx <- which(vitro.cluster.res==i)
	cl.cor.dist <- selected.cor.mat[,cl.sample.idx]
	return(cl.cor.dist)
})
pdf("Plot_PCC_distribution_across_vivo_clusters_for_each_vitro_cluster_cm_hvg.pdf", height=30, width=8)
par(mfrow=c(10,1))
for(i in seq(length(selected.cor.dist))){
	boxplot(t(selected.cor.dist[[i]]), main=paste("In Vitro cluster",i-1), xlab="In vivo cluster", ylab="PCC")
}
dev.off()

# plot the heatmap showing distribution of correlation to different in vivo clusters for each day in vitro cells
selected.cor.mat <- t(vitroCell2vivoCluster.top.cm)
#scaled.cor.mat <- t(scale(t(selected.cor.mat)))
#max.cor.p <- apply(scaled.cor.mat,1,function(vec){
#	pnorm(max(vec))
#})
#cor.diff <- apply(selected.cor.mat, 1, function(vec){
#	max(vec)-min(vec)
#})
#pcc.cutoff <- 0
#selected.cell.idx <- which(cor.diff>pcc.cutoff)
#used.mat <- selected.cor.mat[selected.cell.idx,]
normed.mat <- t(apply(selected.cor.mat, 1, function(vec){
	(vec-min(vec))/(max(vec)-min(vec))
}))
#hc.col <- hclust(as.dist(1-cor(normed.mat)), method="ward.D2")
hc.col <- hclust(as.dist(1-cor(vivo.cl.expr[top.vivo.cm,])), method="ward.D2")
hc.col$order <- c(1,2,7,8,9,10,11,13)
hc.row <- hclust(as.dist(1-cor(t(normed.mat), method="spearman")), method="ward.D2")
vitro.cell.idx <- hc.row$order
values=paste0("Cluster", bestVivoCluster2VItroCell.top.cm)
vitro.cell.idx.2 <- unlist(lapply(paste0("Cluster",1:13), function(x){
	intersect(vitro.cell.idx, which(values==x))
}))
colors <- colorRampPalette(c("#2166ac","#67a9cf","#d1e5f0","#f7f7f7","#fddbc7","#ef8a62","#b2182b"), space="rgb")(50)

#cl.cols <- setNames(scales::hue_pal()(13), paste0("Cluster", seq(13)))
#cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6","#4DA7B0")
cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6","#FDB164")
names(cl.cols) <- paste0("Cluster", seq(13))
cell.col <- cl.cols[values[vitro.cell.idx.2]]
cell.types <- c("Ciliated", "Proliferating ciliated", "Undefined", "Differentiating basal cell", "Secretory progenitor", "Mixture", "Neuroendocrine", "Bud tip adjacent", "Bud tip progenitor", "Hub progenitor", "Undefined proliferating cell", "Submucosal gland",  "Basal cell")
names(cell.types) <- paste0("Cluster", seq(13))

file.name <- "Plot_heatmap_vivoCell2VitroCluster_topCMexcludeDE-1.pdf"
pdf(file.name)
par(mar=c(6,6,5,3))
heatmap.2(normed.mat[vitro.cell.idx.2,hc.col$order],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE,labRow=NA,labCol=cell.types[hc.col$order],ColSideColors=cl.cols[hc.col$order], RowSideColors=cell.col)
dev.off()

# plot the best correlated in vivo cluster distribution on the in vitro tSNE coordinate
tsne.coor <- d21@dr$tsne@cell.embeddings
values=paste0("Cluster", bestVivoCluster2VItroCell.top.cm)
pdf("Plot_tSNE_vitroCell_vs_vivoCluster_topCM_SCC_2nd.pdf", height=20, width=20)
plot(tsne.coor, pch=16, col=cl.cols[values], bty="n", xaxt="n",yaxt="n", xlab="", ylab="",cex=2)
legend("bottomleft", legend=cell.types[intersect(names(cl.cols), values)], text.col=cl.cols[names(cl.cols)%in%values], bty="n")
dev.off()

png("Plot_tSNE_vitroCell_vs_vivoCluster_top_cm-byCluster.png", height=1200, width=2000)
par(mfrow=c(3,5))
for(i in seq(13)){
	x <- paste0("Cluster",i)
	plot(tsne.coor, pch=16, col=rep("gray", length(values)), bty="n", xaxt="n",yaxt="n", xlab="", ylab="",cex.main=2,main=x)
	if(sum(values==x)>1){
		points(tsne.coor[which(values==x),], pch=16, col=cl.cols[values][which(values==x)])
	}else{
		x.coor <- tsne.coor[which(values==x),1]
		y.coor <- tsne.coor[which(values==x),2]
		points(x.coor, y.coor, pch=16, col=cl.cols[values][which(values==x)])
		
	}
}
dev.off() 

# feature plots
markers <- c("FOXJ1","TMEM190","CDK1","MKI67","ASCL1","SFTPC", "SCGB3A2", "SFTPB", "TP63", "KRT15")
markers <- c("MUC5AC", "MUC5B", "SCGB1A1", "CHGA", "SPDEF", "CFTR", "EPCAM", "SOX9", "KRT5", "EGFR", "F3", "PDPN")

colorPal <- grDevices::colorRampPalette(c("navy", "darkorange1"))
png("Plot_tSNE_marker_gene_expression_bo-additional.png", height=4500, width=6000)
par(mfrow=c(3,4))
for(x in markers){
	plotFeature2(tsne.coor, values=d21@data[x, ], xaxt="n", yaxt="n", bty="n", main=x, xlab="", ylab="",cex.main=5,nCols=c("navy", "darkorange1"),cex=4)
}
#for(i in 1:30) rect(20+0.5*(i-1),30,20+0.5*i,33, border=NA, col=colorPal(30)[i])
dev.off()
