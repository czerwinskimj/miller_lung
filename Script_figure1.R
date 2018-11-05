library(Seurat)
library(dplyr)
library(doParallel)
source("~/Work/commonScript/Script_functions.R")

epi <- readRDS("../Res_epi.rds")
fetal <- readRDS("/home/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/Res_Master_1500.rds")
png("Plot_test_all.png")
FeaturePlot(fetal, features.plot=c("EPCAM", "DCN"), cols.use=c("gray", "blue"))
dev.off()

g1 <- c("AIF1", "CD79B", "MS4A1", "CD3E", "CD3D", "DCN", "ACTA2", "COL1A2","ASCL1","CHGA", "CDH1","EPCAM")
cl.expr <- getAveExpr(fetal@meta.data, feature.to.calc="res.0.8", expr=fetal@data)

png("Plot_tSNE_fetal_DCN_VIM_EPCAM_PTPRC.png", height=1500, width=1000)
FeaturePlot(fetal, features.plot=c("DCN", "VIM", "EPCAM", "PTPRC","HBA1"), no.legend=F, cols.use=c("gray", "blue"))
dev.off()

current.cluster.ids <- c(0, seq(26))
new.cluster.ids <- c(12,13,1,2,21,14,26,15,3,16,17,18,19,4,20,5,6,24,27,11,25,22,7,8,9,23,10)
cell.idx <- rep(NA, nrow(fetal@meta.data))
for(new.idx in seq(27)){
	old.idx=current.cluster.ids[which(new.cluster.ids==new.idx)]
		cell.idx[which(fetal@meta.data$res.0.8==old.idx)] <- new.idx
	
}
fetal@meta.data[,"new_cl_idx"] <- cell.idx
lineage <- rep("unknown", nrow(fetal@meta.data))
epi.cluster <- seq(11)
meso.cluster <- c(12:17, 19, 20)
immu.cluster <- c(21,23,24)
rbc.cluster <- 27
unknown <- c(18,22,25,26)
lineage[which(fetal@meta.data[,"new_cl_idx"]%in%epi.cluster)] <- "epithelial"
lineage[which(fetal@meta.data[,"new_cl_idx"]%in%meso.cluster)] <- "mesenchymal"
lineage[which(fetal@meta.data[,"new_cl_idx"]%in%immu.cluster)] <- "immune"
lineage[which(fetal@meta.data[,"new_cl_idx"]%in%rbc.cluster)] <- "rbc"
 
 
expr <- as.matrix(fetal@data)
registerDoParallel(30)
res.epi <- foreach(i=seq(nrow(expr)), .multicombine=F, .combine='rbind')%dopar%{
	if(i%%500==0){cat(paste(i,"done\n"))}
	e0 <- exp(expr[i,])-1
	e1 <- e0[which(lineage=="epithelial")]
	e2 <- e0[which(lineage!="epithelial")]
	p1 <- sum(e1>0)/length(e1)
	p2 <- sum(e2>0)/length(e2)
	logfc <- log(mean(e1)/mean(e0))
	pval <- wilcox.test(e1, e0, alternative="g")$p.value
	return(c(p1, p2, p1-p2, logfc, pval))
}
rownames(res.epi) <- rownames(expr)
colnames(res.epi) <- c("p1", "p2", "p_dif", "avg_logFC", "pval")

pan.epi.markers <- rownames(res.epi)[which(res.epi[,"p_dif"]>0.7 & res.epi[,"avg_logFC"]>0.5 & res.epi[,"pval"]<0.05)]
pan.epi.markers <- rev(c("EPCAM", "KRT8", "KRT18", "FXYD3", "PERP", "CDH1"))

get.pan.markers <- function(group.idx.vec=lineage, selected.group="epithelial", core.num=40){
	registerDoParallel(core.num)
	res.mat <- foreach(i=seq(nrow(expr)), .multicombine=F, .combine='rbind')%dopar%{
		if(i%%500==0){cat(paste(i,"done\n"))}
		e0 <- exp(expr[i,])-1
		e1 <- e0[which(group.idx.vec==selected.group)]
		e2 <- e0[which(group.idx.vec!=selected.group)]
		p1 <- sum(e1>0)/length(e1)
		p2 <- sum(e2>0)/length(e2)
		logfc <- log(mean(e1)/mean(e0))
		pval <- wilcox.test(e1, e0, alternative="g")$p.value
		return(c(p1, p2, p1-p2, logfc, pval))
	}
	stopImplicitCluster()
	rownames(res.mat) <- rownames(expr)
	colnames(res.mat) <- c("p1", "p2", "p_dif", "avg_logFC", "pval")
	return(res.mat)
}

res.imm <- get.pan.markers(group.idx.vec=lineage, selected.group="immune", core.num=40)
res <- res.imm
pan.immune.markers <- setdiff(rownames(res)[which(res[,"p_dif"]>0.7 & res[,"avg_logFC"]>0.5 & res[,"pval"]<0.05)],c("LSP1","PTPRCAP","TBC1D10C","TRAF3IP3","LCK","IL2RG","RAC2","GMFG","LMD2","ARHGDIB","LTB","CD52","LMD2","GPSM3","LIMD2"))
pan.meso.markers <- c("COL6A2","MFAP4","DCN","FHL1","COL1A2","COL3A1")
pan.imm.markers <- c("CD37","CORO1A","LCP1","PTPRC","CD53","LAPTM5")
endothelial.markers <- c("FLT1", "CDH5", "CLDN5", "EGFL7", "ESAM","IFI27")

pdf("Plot_fetal_dotplot-all.pdf", width=10, height=7)
DotPlot(fetal, genes.plot=c(pan.imm.markers, pan.meso.markers, pan.epi.markers), x.lab.rot=TRUE, plot.legend=TRUE, group.by="new_cl_idx")
dev.off()

# rearrange the order again #
current.cluster.ids <- seq(27)
new.cluster.ids <- c(seq(20),23,22,24,25,26,21,27)
cell.idx <- rep(NA, nrow(fetal@meta.data))
for(new.idx in new.cluster.ids){
	old.idx=current.cluster.ids[which(new.cluster.ids==new.idx)]
		cell.idx[which(fetal@meta.data$new_cl_idx==old.idx)] <- new.idx
	
}
fetal@meta.data[,"cl_idx_3"] <- cell.idx
pdf("Plot_fetal_dotplot_all_final.pdf", width=10, height=7)
DotPlot(fetal, genes.plot=c(pan.imm.markers, pan.meso.markers, pan.epi.markers), x.lab.rot=TRUE, plot.legend=TRUE, group.by="cl_idx_3")
dev.off()

# rearrange the order again #
new.cluster.ids <- seq(27)
current.cluster.ids <- c(11,seq(5),7:10,6,12:27)
cell.idx <- rep(NA, nrow(fetal@meta.data))
for(new.idx in new.cluster.ids){
	old.idx=current.cluster.ids[which(new.cluster.ids==new.idx)]
		cell.idx[which(fetal@meta.data$cl_idx_3==old.idx)] <- new.idx
	
}
fetal@meta.data[,"cl_idx_4"] <- cell.idx
neuron.mix.markers <- c("CHGA","ASCL1","NNAT","STMN2","GRP","MPZ")


pdf("Plot_fetal_dotplot_all_add_neuron_add_endothelial.pdf", width=10, height=7)
DotPlot(fetal, genes.plot=c(pan.imm.markers, endothelial.markers, pan.meso.markers, pan.epi.markers, neuron.mix.markers), x.lab.rot=TRUE, plot.legend=TRUE, group.by="cl_idx_4", cols.use=c("#d9d9d9", "#252525"))
dev.off()

fetal@ident <- as.factor(fetal@meta.data[,"cl_idx_4"])
names(fetal@ident) <- fetal@cell.names
png("Plot_tSNE_fetal_cluster_final.png")
TSNEPlot(fetal, do.label=TRUE)
dev.off()

lineage <- rep("RBC", nrow(fetal@meta.data))
lineage[which(fetal@meta.data$cl_idx_4%in%seq(11))] <- "Epithelial"
lineage[which(fetal@meta.data$cl_idx_4%in%c(12:20))] <- "Mesenchymal"
lineage[which(fetal@meta.data$cl_idx_4==21)] <- "Endothelial"
lineage[which(fetal@meta.data$cl_idx_4%in%c(22:26))] <- "Immune"
fetal@meta.data[,"lineage"] <- lineage
lineage.cols <- c("#f0f0f0","#e3e3e3", "#c9c9c9", "#a9a9a9", "#ef6548")
names(lineage.cols) <- c("RBC","Immune", "Endothelial", "Mesenchymal", "Epithelial")
cell.cols <- lineage.cols[lineage]
png("Plot_tSNE_lineage_assignment.png", height=2000, width=2000)
plot(tsne.coor, col=cell.cols, pch=19, main="", bty="n", cex.main=5,  cex=3.5, xaxt="n", yaxt="n", xlab="", ylab="")
legend("topleft", pch=19, col=rev(lineage.cols), legend=rev(names(lineage.cols)),cex=3.5,bty="n")
dev.off()

pan.lineage.markers <- list("epithelial"=pan.epi.markers, "mesenchymal"=pan.meso.markers, "immune"=pan.imm.markers)
saveRDS(pan.lineage.markers, file="Res_pan_lineage_markers.rds")

tsne.coor <- fetal@dr$tsne@cell.embeddings
colors <- c("#005a32", "#41ab5d", "#78c679", "#ffffb2", "#fc4e2a", "#e31a1c","#b10026")
g1 <- pan.lineage.markers[["epithelial"]]
e1 <- rank(apply(fetal@scale.data[g1,], 2, mean))
normed.e <- apply(fetal@data[g1,], 1, function(vec){
	(vec-min(vec))/(max(vec)-min(vec))
})
e1 <- apply(normed.e,1,mean)
png("Plot_tSNE_selected_pan_epithelial_lineage_markers_expression-quantile.png", height=2000, width=2000)
plotFeature2(tsne.coor, values=e1, main="Epithelial", nCols=colors, bty="n", cex.main=5,  cex=3, xaxt="n", yaxt="n", xlab="", ylab="", zeroAsGray=FALSE)
for(i in seq(27)){
	idx <- which(fetal@meta.data[,"cl_idx_4"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=4.5)
}
dev.off()

g1 <- pan.lineage.markers[["mesenchymal"]]
e1 <- rank(apply(fetal@scale.data[g1,], 2, mean))
png("Plot_tSNE_selected_pan_mesenchymal_lineage_markers_expression-quantile.png", height=2000, width=2000)
plotFeature2(tsne.coor, values=e1, main="Mesenchymal", bty="n", cex.main=5, nCols=colors, cex=3, xaxt="n", yaxt="n", xlab="", ylab="")
for(i in seq(27)){
	idx <- which(fetal@meta.data[,"cl_idx_4"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=4.5)
}
dev.off()

g1 <- pan.lineage.markers[["immune"]]
e1 <- rank(apply(fetal@scale.data[g1,], 2, mean))
png("Plot_tSNE_selected_pan_immune_lineage_markers_expression-quantile.png", height=2000, width=2000)
plotFeature2(tsne.coor, values=e1, main="Immune", bty="n", cex.main=5, nCols=colors, cex=3, xaxt="n", yaxt="n", xlab="", ylab="")
for(i in seq(27)){
	idx <- which(fetal@meta.data[,"cl_idx_4"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=4.5)
}
dev.off()

# identify potential doublet using doubletFinder
library(doubletFinder)
fetal.doub <- doubletFinder(fetal, expected.doublets=1000, proportion.NN = 0.01)
fetal.doub.2 <- doubletFinder(fetal, expected.doublets=1000, proportion.NN = 0.005)
fetal@meta.data[,"pANN"] <- fetal.doub.2@meta.data$pANN
saveRDS(fetal@meta.data[,"pANN"], file="Res_fetal_doubletFinder_pANN_under_pnn0.005.rds")

 
seu <- fetal.doub.2
tsne.coor <- seu@dr$tsne@cell.embeddings
proportion.mat <- sapply(seq(from=100, to=2000, by=200), function(top.num){
	doublet.idx <- order(seu@meta.data$pANN,decreasing=T)[1:top.num]
	sapply(seq(27), function(i){
		cl.idx <- which(seu@meta.data[,"cl_idx_3"]==i)
		sum(cl.idx%in%doublet.idx)/length(cl.idx)
	})
})
png("Plot_tSNE_doubletFinder_pnn0.005.png", height=1000, width=2000)
par(mfrow=c(1,2), mar=c(5,5,5,5))
plot(tsne.coor, col="#bdbdbd30", pch=19, bty="n", xaxt="n", yaxt="n", xlab="", ylab="", cex=2)
for(top.num in seq(from=100, to=2000, by=200)){
	doublet.idx <- order(seu@meta.data$pANN,decreasing=T)[1:top.num]
	points(tsne.coor[doublet.idx,], col="#3182bd10", pch=19, cex=2)
}
for(i in seq(27)){
	idx <- which(seu@meta.data[,"cl_idx_3"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=2)
}
legend("topleft", legend="doublet", pch=19, bty="n", col="#3182bd", cex=2)
plot(seq(27), proportion.mat[,1], ylim=c(0,1), type="l", xaxt="n", xlab="Cluster", ylab="Doublet proportion under different cutoffs", cex.axis=1.5, cex.lab=2, bty="n")
axis(1, at=seq(27), las=1, cex.axis=1.5)
for(j in seq(ncol(proportion.mat))){
	lines(seq(27), proportion.mat[,j])
}
points(seq(27), proportion.mat[,5], pch=19)
abline(h=900/nrow(seu@meta.data), lty=2, col="gray", lwd=3)
legend("topleft", pch=19, legend="Top900 cells ranked by proportions of artificial doublet \n as nearest neighbor as putative doublet", bty="n", cex=2)
dev.off()

fetal.cm.6 <- FindMarkers(fetal, ident.1=6, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
g1 <- rownames(fetal.cm.6)[order(fetal.cm.6$avg_logFC, decreasing=T)[1:10]]
pdf("Plot_fetal_dotplot_cl6.pdf", width=10, height=7)
DotPlot(fetal, genes.plot=g1, x.lab.rot=TRUE, plot.legend=TRUE, group.by="cl_idx_3")
dev.off()

fetal.cm <- findAllMarkers(seu.obj=fetal, selected.column="cl_idx_3", cluster.to.test=c(2, 7, 8, 22))
saveRDS(fetal.cm, file="Res_all_fetal_cells_specific_cluster_markers_cl_idx_3_cluster2_7_8_22.rds")
fetal.top10 <- fetal.cm %>% group_by(cluster) %>% top_n(10, avg_logFC)
g2 <- unique(fetal.top10$gene_name)
pdf("Plot_fetal_dotplot_cl6_2_7_8_22.pdf", width=15, height=7)
DotPlot(fetal, genes.plot=union(g1,g2), x.lab.rot=TRUE, plot.legend=TRUE, group.by="cl_idx_3")
dev.off()

lineage <- rep("RBC", nrow(fetal@meta.data))
lineage[which(fetal@meta.data$cl_idx_4%in%seq(11))] <- "epithelial"
lineage[which(fetal@meta.data$cl_idx_4%in%seq(11))]


fetal.cm <- findAllMarkers(seu.obj=fetal, selected.column="res.0.8", cluster.to.test=c(26, 21, 18, 19, 6))
saveRDS(fetal.cm, file="Res_all_fetal_cells_specific_cluster_markers.rds")

fetal.top10 <- fetal.cm %>% group_by(cluster) %>% top_n(10, avg_logFC)
g1 <- unique(fetal.top10$gene_name)
png("Plot_dot_plot_denovo.png", width=1000)
DotPlot(fetal, genes.plot=g1, x.lab.rot=TRUE, plot.legend=TRUE)
dev.off()

# plot fetal tSNE showing sample information
g.cols <- c("#49A999", "#332786", "#CA6778", "#197636", "#A84798", "#872555")
names(g.cols) <- c("D103dis","D103air","D103tra","D125dis","D125prox","D145tra")
sample.names <- c("15 week distal", "15 week airway", "15 week trachea", "18 week distal", "18 week proximal", "21 week trachea")
group.res <- as.character(fetal@meta.data$orig.ident)
pdf("Plot_tSNE_fetal_all_no_sample_identity.pdf", height=15, width=30)
par(mfrow=c(1,2))
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=g.cols, cex=1.8)

plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=g.cols, type="n")
legend("topleft", legend=sample.names, text.col=g.cols, bty="n")
dev.off()

# plot fetal tSNE showing cluster information
head(fetal@meta.data)
tsne.coor <- fetal@dr$tsne@cell.embeddings
group.res <- as.character(fetal@meta.data$cl_idx_4)
pdf("Plot_tSNE_fetal_all_cluster_default_cols.pdf", height=15, width=15)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1.8)
for(i in seq(27)){
	idx <- which(fetal@meta.data[,"cl_idx_4"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=2)
}
dev.off()

unique(fetal@meta.data$lineage)
epi.cl <- unique(fetal@meta.data$cl_idx_4[which(fetal@meta.data$lineage=="Epithelial")])
epi.cols <- colorRampPalette(c("#fee0d2", "#fb6a4a","#a50f15"))(length(epi.cl))
names(epi.cols) <- paste0("Cluster", epi.cl)

mes.cl <- unique(fetal@meta.data$cl_idx_4[which(fetal@meta.data$lineage=="Mesenchymal")])
mes.cols <- colorRampPalette(c("#e5f5e0", "#74c476","#006d2c"))(length(mes.cl))
names(mes.cols) <- paste0("Cluster", mes.cl)

imm.cl <- unique(fetal@meta.data$cl_idx_4[which(fetal@meta.data$lineage=="Immune")])
imm.cols <- colorRampPalette(c("#deebf7", "#6baed6","#08519c"))(length(imm.cl))
names(imm.cols) <- paste0("Cluster", imm.cl)

endo.cols <- "#6a51a3"
names(endo.cols) <- paste0("Cluster", unique(fetal@meta.data$cl_idx_4[which(fetal@meta.data$lineage=="Endothelial")]))
rbc.cols <- "#bdbdbd"
names(rbc.cols) <- paste0("Cluster", unique(fetal@meta.data$cl_idx_4[which(fetal@meta.data$lineage=="RBC")]))
cell.type.cols <- c(epi.cols, mes.cols, imm.cols, endo.cols, rbc.cols)
group.res <- paste0("Cluster", as.character(fetal@meta.data$cl_idx_4))
pdf("Plot_tSNE_fetal_all_cluster_sequential_cols.pdf", height=15, width=15)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1.8,gCols=cell.type.cols)
for(i in seq(27)){
	idx <- which(fetal@meta.data[,"cl_idx_4"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=2)
}
dev.off()


## part 2. only focus on epithelial cells
new.cluster.ids <- 0:13
current.cluster.ids <- c(0,9,8,11,10,13,4,6,5,7,12,2,3,1)
cell.idx <- rep(NA, nrow(epi@meta.data))
for(new.idx in new.cluster.ids){
	old.idx=current.cluster.ids[which(new.cluster.ids==new.idx)]
		cell.idx[which(epi@meta.data$merged.cl==old.idx)] <- new.idx
	
}
epi@meta.data[,"cl_idx_2"] <- cell.idx

cluster.idx <- 0:13
group.idx <- c(0,1,1,2,3,4,4,5,6,8,7,9,9,9)
pan.cluster.idx <- rep(NA, nrow(epi@meta.data))
for(idx1 in group.idx){
	idx2=cluster.idx[which(group.idx==idx1)]
		pan.cluster.idx[which(epi@meta.data$cl_idx_2%in%idx2)] <- idx1
	
}
epi@meta.data[,"pan_cluster_idx"] <- pan.cluster.idx

pan.cluster.names <- c("0:Non-epithelial", "1:Bud tip", "2:Proliferative intermediate", "3:Hub progenitor", "4:Basal cells", "5:Early secretory progenitor", "6:Goblet/club", "7:SMG", "8:NE", "9:Ciliated")
names(pan.cluster.names) <- paste0("Cluster",0:9)
pan.cluster.identity <- pan.cluster.names[paste0("Cluster", epi@meta.data[,"pan_cluster_idx"])]
epi@meta.data[,"pan_cluster_identity"] <- pan.cluster.identity
epi@ident <- as.factor(pan.cluster.identity)
names(epi@ident) <- epi@cell.names
png("Plot_tSNE_epi_pan_cluster_idx.png")
TSNEPlot(epi, do.label=TRUE)
dev.off()

epi.pan.cm <- findAllMarkers(epi, selected.column="pan_cluster_idx", core.num=10, p.value.cutoff=0.05, bg=1:9, cluster.to.test=1:9)
saveRDS(epi.pan.cm, file="Res_epi_pan_cluster_markers.rds")
cl.ave.expr <- getAveExpr(meta.data=epi@meta.data, feature.to.calc="pan_cluster_idx", expr=epi@data, core.num=20)
idx <- epi.pan.cm$cluster
idx[which(epi.pan.cm$cluster==7)] <- 8
idx[which(epi.pan.cm$cluster==8)] <- 7
epi.pan.cm.2 <- epi.pan.cm 
epi.pan.cm.2$cluster <- idx
saveRDS(epi.pan.cm.2, file="Res_epi_pan_cluster_markers.rds")
saveRDS(epi, file="Res_epi_final.rds")
saveRDS(fetal, file="Res_fetal_final.rds")



#top.cm <- epi.pan.cm %>% group_by(cluster) %>% top_n(20, avg_logFC)
#g2 <- unique(top.cm$gene_name)
g1 <- unique(unlist(lapply(seq(9), function(i){
	cl.idx <- which(epi.pan.cm.2$cluster==i)
	mat <- epi.pan.cm.2[cl.idx,]
	p.dif <- mat$pct.1-mat$pct.2
	fc <- mat$avg_logFC
	g0.p <- mat$gene_name[order(p.dif, decreasing=T)[1:20]]
	#g0.fc <- mat$gene_name[order(fc, decreasing=T)[1:20]]
	#g0 <- union(g0.p, g0.fc)
})))
g1 <- c(g1, "FOXJ1")
g1 <- setdiff(g1, c("NREP", "SPARC", "CPE", "VIM", "CXCL17", "HES6", "SERPINF1", "IFITM2", "TGM2", "MYL9","CAV1"))
g1 <- c(g1[1:60], "TP63", g1[61:length(g1)])
res <- epi.pan.cm.2[epi.pan.cm.2$gene_name%in%g1,]
write.table(res, file="Table_figure1c_heatmap_genes.txt", sep="\t", row.names=F,col.names=T,quote=F)
png("Plot_epi_denovo_cluster_markers-4.png", width=2000, height=1500)
par(mar=c(6,2,2,2))
DoHeatmap(epi, data.use=epi@scale.data[, which(epi@meta.data$pan_cluster_idx!=0)], genes.use=g1, slim.col.label=TRUE, group.by="pan_cluster_identity", group.label.rot=TRUE)
dev.off()

cl.ave.expr.2 <- cl.ave.expr[,-1]
cl.ave.expr.2[,7] <- cl.ave.expr[,9]
cl.ave.expr.2[,8] <- cl.ave.expr[,8]
colors <- colorRampPalette(c("#FF00FF","#000000","#FFFF00"), space="rgb")(50)
colors <- colorRampPalette(c("#8073ac","#f7f7f7","#e08214"), space="rgb")(50)
colors <- colorRampPalette(c("#f7f7f7", "#d9d9d9", "#252525"), space="rgb")(50)


cl.cols <- c("#9E0142","#A84798", "#5E4FA2", "#EED951","#FDD884","#DCCB7C","#3C7AB6","#4DA7B0","#197636")
names(cl.cols) <- c("1:Bud tip", "2:Proliferative intermediate", "3:Hub progenitor", "4:Basal cells", "5:Early secretory progenitor", "6:Goblet/club", "7:SMG", "8:NE", "9:Ciliated")


file.name <- "Plot_epi_denovo_cluster_markers_by_cluster-5.pdf"
pdf(file.name, height=8)
heatmap.2(cl.ave.expr.2[g1,],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE, labRow=FALSE, labCol=names(cl.cols), ColSideColors=cl.cols, srtCol=30, margins = c(8, 5))
dev.off()

res <- read.table("Table_top50_cluster_marker_genes_for_vivo_epi_SPRING.txt",sep="\t",stringsAsFactors=F,head=T)
cluster.order <- c(9,8,11,10,13,4,6,3,1,2,5,12,7)
g1 <- unique(res$gene_name)
cl.ave.expr <- getAveExpr(meta.data=epi@meta.data, feature.to.calc="merged.cl", expr=epi@data, core.num=20,genes=g1,specified.order=seq(13))
scale.expr <- t(apply(cl.ave.expr,1,scale))
colnames(scale.expr) <- colnames(cl.ave.expr)
g.ordered <- t(sapply(g1, function(g){
	cl.idx <- res$cluster[which(res$gene_name==g)]
	vec <- scale.expr[g, cl.idx]
	res <- c(cl.idx[which.max(vec)], max(vec))
}))
colnames(g.ordered) <- c("max.cluster.idx", "max.cluster.scale.expr")
saveRDS(g.ordered, file="Res_springMarkers_max_cluster_idx_and_expr.rds")

res.ave <- lapply(cluster.order, function(i){
	mat <- g.ordered[which(g.ordered[,"max.cluster.idx"]==i),]
	mat <- mat[order(mat[,"max.cluster.scale.expr"], decreasing=T),]
	return(mat)
})
res.ave <- do.call('rbind', res.ave)

cluster.names <- c("0:Non-epithelial", "1:Multiciliated","2:Proliferating ciliated","3:Intermediate ciliated","4:Differentiating basal cells", "5:Goblet-like secretory","6:Club-like secretory", "7:NE", "8:Bud tip adjacent", "9:Bud tip progenitor","10:Hub progenitor", "11:Proliferative intermediate","12:SMG",   "13:Basal cells")
names(cluster.names) <- paste0("Cluster",c(0,seq(13)))
cluster.identity <- cluster.names[paste0("Cluster", epi@meta.data[,"merged.cl"])]
epi@meta.data[,"cluster_identity"] <- cluster.identity
epi@ident <- as.factor(cluster.identity)

cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6","#FDB164")
names(cl.cols) <- paste0("Cluster", seq(13))
cl.cols <- cl.cols[paste0("Cluster",cluster.order)]
colors <- colorRampPalette(c("#f0f0f0", "#d9d9d9", "#252525"), space="rgb")(50)
out <- scale.expr[rownames(res.ave), paste0("Cluster", cluster.order)]
file.name <- "Plot_epi_denovo_cluster_markers_spring_markers.pdf"
pdf(file.name, height=8)
heatmap.2(out,trace="none",density.info="none",scale="none",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE, labRow=FALSE, labCol=cluster.names[paste0("Cluster",cluster.order)], ColSideColors=cl.cols, srtCol=30, margins = c(8, 5))
dev.off()
  
known.marker <- c("SFTPC", "AGER", "MKI67", "SFTPB", "SCGB3A2","TP63", "KRT5", "SCGB1A1", "FOXJ1","CDK1", "MUC5B", "LTF", "CHGA")
e1 <- lapply(known.marker, function(x){
	lapply(cluster.order, function(i){
		epi@data[x,which(epi@meta.data$merged.cl==i)]
	})	
})


pdf("Plot_boxplot_known_markers_13clusters-1.pdf", height=50)
par(mfrow=c(length(known.marker),1))
for(i in seq(length(e1))){
	boxplot(e1[[i]], main=known.marker[i], outline=F, bty="n", col=cl.cols, border=cl.cols, lwd=1, names=cluster.order)
}
dev.off()


colorPal <- grDevices::colorRampPalette(c("navy", "darkorange1"))
png("Plot_gene_expression_scale.png")
plot(seq(from=18,to=50,length=10),seq(from=28,to=35,length=10),type="n",main="",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
for(i in 1:3) rect(20+0.5*(i-1),30,20+0.5*i,33, border=NA, col="#bdbdbd")
for(i in 4:33) rect(20+0.5*(i-1),30,20+0.5*i,33, border=NA, col=colorPal(30)[i-3])
dev.off()


pdf("Plot_boxplot_EGFR_F3.pdf", height=20)
par(mfrow=c(3,1))
boxplot(e1, main="F3+EGFR", outline=F, bty="n", col=cl.cols, border=cl.cols, lwd=4)
boxplot(e2, main="EGFR", outline=F, bty="n", col=cl.cols, border=cl.cols, lwd=4)
boxplot(e3, main="F3", outline=F, bty="n", col=cl.cols, border=cl.cols, lwd=4)
dev.off()


png("Plot_epi_dot_plot.png")
DotPlot(epi, genes.plot=c("EPCAM", "CDH1", "VIM", "DCN"))
dev.off()

png("Plot_epi_tSNE_TM4SF1_SFTPC_SOX2_KRT5.png", width=800, height=800)
FeaturePlot(epi, features.plot=c("TM4SF1", "SFTPC", "SOX2", "KRT5"), cols.use=c("gray", "blue"), no.legend=F)
dev.off()

png("Plot_epi_tSNE_neuronal_lineage.png", width=1600, height=1200)
FeaturePlot(epi, features.plot=fetal.top10$gene_name[fetal.top10$cluster==19], cols.use=c("gray", "blue"), no.legend=F)
dev.off()

# bud tip progenitor/basal cell markers
#epi@ident <- as.factor(epi@meta.data$merged.cl)
#names(epi@ident) <- epi@cell.names
#cm <- findAllMarkers(seu.obj=epi, pos.only=F, selected.column="merged.cl", cluster.to.test=c(9,13), bg=1:13)
#write.table(cm, file="/r1/people/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/epi_by_sample/Table_cluster9_13_marker_compared2Epi.txt",sep="\t",quote=F)

res <- read.table("Table_top50_cluster_marker_genes_for_vivo_epi_SPRING.txt",sep="\t",stringsAsFactors=F,head=T)
cl.idx <- which(res$cluster==9)
mat <- res[cl.idx,]
up.9 <- mat$gene_name[order(mat$avg_logFC, decreasing=T)[1:20]]

cl.idx <- which(res$cluster==13)
mat <- res[cl.idx,]
up.13 <- union(mat$gene_name[order(mat$avg_logFC, decreasing=T)[1:20]], c("EGFR", "F3"))

cl.ave.expr <- getAveExpr(meta.data=epi@meta.data, feature.to.calc="merged.cl", expr=epi@data, genes=union(up.9, up.13))
cl.ave.expr <- cl.ave.expr[,-1]
colors <- colorRampPalette(c("midnightblue","white","darkorange1"))(n=299)
plot.file <- "Plot_heatmap_cluster_marker_basalCell_features_rankByFC.pdf"
pdf(plot.file)
heatmap.2(cl.ave.expr[up.13,],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE)
dev.off() 

plot.file <- "Plot_heatmap_cluster_marker_budTipProgenitor_features_rankByFC.pdf"
pdf(plot.file)
heatmap.2(cl.ave.expr[up.9,],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE)
dev.off() 


tsne.coor <- epi@dr$tsne@cell.embeddings
idx <- which(epi@meta.data$merged.cl!=0)
png("Plot_tSNE_featureplot_fig1d.png", height=4000, width=2000)
par(mfrow=c(2,1))
for(x in c("SFTPC", "CA2")){
	plotFeature2(tsne.coor[idx,], values=epi@data[x,idx], xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=6, cex.main=5, nCols=c("navy", "darkorange1"))
}
dev.off()

png("Plot_tSNE_featureplot_fig1d-3.png", height=1500, width=1500)
par(mfrow=c(5,5))
for(x in c(up.9,"SOX9")){
	plotFeature2(tsne.coor[idx,], values=epi@data[x,idx], xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", nCols=c("navy", "darkorange1"))
}
dev.off()

png("Plot_tSNE_featureplot_fig1h.png", height=4000, width=4000)
par(mfrow=c(2,2))
for(x in c("TP63", "IL33", "KRT5", "KRT15")){
	plotFeature2(tsne.coor[idx,], values=epi@data[x,idx], xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=6, cex.main=5, nCols=c("navy", "darkorange1"))
}
dev.off()

png("Plot_tSNE_featureplot_supFig2-noName.png", height=3000, width=7500)
par(mfrow=c(2,5))
for(x in c("MUC5B", "MUC5AC", "SPDEF", "SCGB1A1", "SCGB3A2", "LTF", "FOXJ1", "MKI67","CDK1")){
	plotFeature2(tsne.coor[idx,], values=epi@data[x,idx], xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=6, cex.main=5, nCols=c("navy", "darkorange1"))
}
dev.off()
