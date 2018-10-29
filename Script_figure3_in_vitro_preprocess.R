samples <- unique(samples)
seu.by.sample <- list()
for(idx in samples){
	cat(paste(idx, "start\n"))
	work.dir <- paste0(idx, "/")
	dir.create(file.path(work.dir))
	
	seu.obj <- SubsetData(InVitro, cells.use = which(InVitro@meta.data$samples==idx), do.clean=T)
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
	
	seu.by.sample[[idx]] <- seu.obj
}
