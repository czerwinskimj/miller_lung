epi <- readRDS("../Res_epi.rds")
source("~/Work/commonScript/Script_functions.R")

cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6","#FDB164","#bdbdbd")
names(cl.cols) <- paste0("Cluster", c(seq(13),0))
cell.types <- c("Multiciliated cells", "Proliferating ciliated cells", "Undefined", "Differentiating basal cells", "Early secretory progenitors", "Late secretory/ differentiating club/ goblet cells", "Neuroendocrine", "Bud tip adjacent/newly differentiating bud tip progenitors", "Bud tip progenitors", "Hub progenitors", "Proliferative intermediates", "Submucosal gland",  "Basal cells","Mesenchymal cells")
names(cell.types) <- paste0("Cluster", c(seq(13),0))

# tSNE
tsne.coor <- epi@dr$tsne@cell.embeddings
cluster.res <- paste0("Cluster", epi@meta.data$merged.cl)
pdf("Plot_tSNE_no_cluster_identity.pdf", height=15, width=30)
par(mfrow=c(1,2))
plotFeature2(tsne.coor, values=cluster.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=cl.cols, cex=1.8)

plotFeature2(tsne.coor, values=cluster.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=cl.cols, type="n")
legend("topleft", legend=cell.types, text.col=cl.cols, bty="n")
dev.off()

# SPRING
spring.coor <- read.csv("/r1/people/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/SPRING_merged_cm_scaled_regressed/SPRING_merged_cm_scaled_regressed_k50/coordinates.txt",head=F,row.names=1)
hvg.info <- readRDS("/r1/people/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/epi_by_sample/Res_mergedCluster_top50_cluster_marker_hvg_info.rds")
cluster.res <- as.numeric(as.matrix(hvg.info["merged.cl",]))
hvg.expr <- readRDS("/r1/people/qianhui_yu/Work/Lung/1500gene_Proximal_distal_separate_analysis/in_vivo/not_to_regress_out_rbc_gene_expression/epi_by_sample/Res_mergedCluster_top50_cluster_marker_hvg_scaled_expr.rds")
disMat <- 1-cor(hvg.expr)
k <- 15
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))

cluster.res <- as.numeric(as.matrix(hvg.info["merged.cl",]))
cluster.res <- paste0("Cluster", cluster.res)
pdf("Plot_SPRING_no_cluster_identity-2.pdf", height=25, width=50)
#png("Plot_SPRING_no_cluster_identity-2.png", height=1500, width=3000)
par(mfrow=c(1,2))
plotFeature2(spring.coor, values=cluster.res, knn.pairs=knn.idx, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=cl.cols, cex=4)

plotFeature2(spring.coor, values=cluster.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=cl.cols, type="n")
legend("topleft", legend=cell.types[-length(cell.types)], text.col=cl.cols[-length(cell.types)], bty="n", cex=3)
dev.off()

g1 <- c("TP63", "SFTPC", "PDPN", "SCGB3A2", "SFTPB", "KRT5",  "FOXJ1",  "SCGB1A1", "SOX9", "CHGA")
pdf("Plot_SPRING_cell_type_marker_gene_expr_main.pdf", height=30, width=75)
par(mfrow=c(2,5))
for(x in g1){
	expr <- epi@data[x,colnames(hvg.info)]
	plotFeature2(spring.coor, values=expr, knn.pairs=knn.idx, xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", nCols=c("navy", "darkorange1"), cex=4, cex.main=5)
}
dev.off()

g2 <- c("CFTR", "MUC5AC",  "LTF", "SPDEF")
pdf("Plot_SPRING_cell_type_marker_gene_expr_supp.pdf", height=15, width=60)
par(mfrow=c(1,4))
for(x in g2){
	expr <- epi@data[x,colnames(hvg.info)]
	plotFeature2(spring.coor, values=expr, knn.pairs=knn.idx, xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", nCols=c("navy", "darkorange1"), cex=4, cex.main=5)
}
dev.off()
