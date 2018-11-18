library(Seurat)
library(gplots)
library(doParallel)
source("~/Work/commonScript/Script_functions.R")
d0 <- readRDS("Res_Day0_B1.rds")
d21 <- readRDS("../Day21_B1/Res_Day21_B1.rds")

cell.type.markers <- c("SOX9", "SFTPC", "SCGB3A2", "SFTPB", "CFTR", "TP63", "KRT15", "LTF", "FOXJ1", "SPDEF", "MUC5AC", "CHGA","CDK1","MKI67")

png("Plot_cell_type_marker_expression.png", width=1200, height=1200)
FeaturePlot(d0, features.plot=cell.type.markers, cols.use=c("gray", "blue"), no.legend=T)
dev.off()

vivo.cl.expr <- readRDS("/home/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vitro/only_D0_B1_1000_remove_D0_B2/Dat_in_vivo_merged_cluster_average_expr.rds")
in.vivo.cm <- readRDS("~/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/epi_by_sample/Res_Epi_merged_cluster_marker.rds")
top.vivo.cm <- read.table("/home/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/SPRING_merged_cm_scaled_regressed/List_top50_cluster_markers.txt",stringsAsFactors=F)[,1]
 
# calculate the transcriptome similarity between 
## use interset of in vivo cluster markers and in vitro highly variable genes
#g1 <- intersect(in.vivo.cm$gene_name, d0@var.genes)
## use interset of in vivo cluster markers and in vitro cluster markers
#g1 <- intersect(in.vivo.cm$gene_name, d0.cm$gene_name)
## use all in vivo cluster markers
#g1 <- unique(in.vivo.cm$gene_name)
## use top in vivo cluster markers
g1 <- top.vivo.cm

cl.cols <- setNames(scales::hue_pal()(13), paste0("Cluster", seq(13)))
tsne.coor <- d0@dr$tsne@cell.embeddings


registerDoParallel(20)
#vitroCell2vivoCluster.cm.hvg <- foreach(k=seq(ncol(d0@data)), .multicombine=T, .combine='cbind')%dopar%{
#vitroCell2vivoCluster.cm.cm <- foreach(k=seq(ncol(d0@data)), .multicombine=T, .combine='cbind')%dopar%{
#vitroCell2vivoCluster.cm <- foreach(k=seq(ncol(d0@data)), .multicombine=T, .combine='cbind')%dopar%{
vitroCell2vivoCluster.top.cm <- foreach(k=seq(ncol(d0@data)), .multicombine=T, .combine='cbind')%dopar%{
	if(k%%500==0){cat(paste0(k,"start\n"))}
	cor.vec <- sapply(seq(ncol(vivo.cl.expr)), function(cl){
		cor(d0@data[g1,k], vivo.cl.expr[g1,cl])
	})
	return(cor.vec)
}
stopImplicitCluster()
saveRDS(vitroCell2vivoCluster.top.cm, file="Res_vitroCell2vivoCluster_topVivoClusterMarkers.rds")

#bestVivoCluster2VItroCell.cm.hvg <- apply(vitroCell2vivoCluster.cm.hvg, 2, which.max)
#values=paste0("Cluster", bestVivoCluster2VItroCell.cm.hvg)

#bestVivoCluster2VItroCell.cm.cm <- apply(vitroCell2vivoCluster.cm.cm, 2, which.max)
#values=paste0("Cluster", bestVivoCluster2VItroCell.cm.cm)

#bestVivoCluster2VItroCell.cm <- apply(vitroCell2vivoCluster.cm, 2, which.max)
#values=paste0("Cluster", bestVivoCluster2VItroCell.cm)

bestVivoCluster2VItroCell.top.cm <- apply(vitroCell2vivoCluster.top.cm, 2, which.max)
values=paste0("Cluster", bestVivoCluster2VItroCell.top.cm)

# plot the heatmap showing distribution of correlation to different in vivo clusters for each day in vitro cells
selected.cor.mat <- t(vitroCell2vivoCluster.top.cm)

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
cell.types <- c("Ciliated", "Proliferating ciliated", "Undefined", "Differentiating basal cell", "Secretory progenitor", "Mixture", "Neuroendocrine", "Bud tip adjacent", "Bud tip progenitor", "Hub progenitor", "Undefined proliferating cell", "Submucosal gland",  "Basal cell")
names(cell.types) <- paste0("Cluster", seq(13))

cell.col <- cl.cols[values[vitro.cell.idx.2]]
file.name <- "Plot_heatmap_vivoCell2VitroCluster-10.pdf"
pdf(file.name)
par(mar=c(10,6,5,3))
heatmap.2(normed.mat[vitro.cell.idx.2,hc.col$order],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE,labRow=NA,labCol=cell.types[hc.col$order],ColSideColors=cl.cols[hc.col$order], RowSideColors=cell.col)
dev.off()

# plot the best correlated in vivo cluster distribution on the in vitro tSNE coordinate
tsne.coor <- d0@dr$tsne@cell.embeddings
values=paste0("Cluster", bestVivoCluster2VItroCell.top.cm)
pdf("Plot_tSNE_vitroCell_vs_vivoCluster_top_cm-9.pdf", height=20, width=20)
plot(tsne.coor, pch=16, col=cl.cols[values], bty="n", xaxt="n",yaxt="n", xlab="", ylab="",cex=3)
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
markers <- c("SFTPC","SOX9", "TP63")
colorPal <- grDevices::colorRampPalette(c("navy", "darkorange1"))
sample.vec <- rep(c("Day0","Day21"), c(ncol(d0@data), ncol(d21@data)))
png("Plot_tSNE_marker_gene_expression_day0.png", height=3000, width=3000)
par(mfrow=c(2,2))
for(x in markers){
	values <- c(d0@data[x,], d21@data[x,])
	cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F,include.lowest=T))]
	cellColor[values == 0] <- "#bdbdbd30"
	used.cellColor <- cellColor[which(sample.vec=="Day0")]
	plot(tsne.coor, col=used.cellColor, pch=16, bty="n", xaxt="n", yaxt="n", main=x, xlab="", ylab="",cex.main=4,cex=7)
}
#for(i in 1:30) rect(7+0.2*(i-1),12,7+0.2*i,14, border=NA, col=colorPal(30)[i])
dev.off()

png("Plot_scale_bar.png")
plot(seq(from=5, to=15, length=10), seq(from=11, to=15, length=10), type="n")
for(i in 1:30) rect(7+0.2*(i-1),12,7+0.2*i,14, border=NA, col=colorPal(30)[i])
dev.off()

markers <- c("SFTPC", "SOX9", "FOXJ1", "TP63")
png("Plot_tSNE_day0_SFTPC_SOX9_FOXJ1_TP63.png", height=1200, width=1200)
par(mfrow=c(2,2))
for(x in markers){
	plotFeature2(tsne.coor, values=d0@data[x,], nCol=c("navy", "darkorange1"), bty="n", xaxt="n", yaxt="n", main=x, xlab="", ylab="",cex.main=4,cex=5)
}
#for(i in 1:30) rect(7+0.2*(i-1),12,7+0.2*i,14, border=NA, col=colorPal(30)[i])
dev.off()
