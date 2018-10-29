#' Plot feature's gradients/classes across all samples (cells) given the plotting coordinates
#'
#' This function plot the samples based on the given coordinates, coloring
#' each dot based on the value of its corresponding feature value/class.
#'
#' @param coord 	The plotting coordinates, expected to be of two columns. Each row represents one dot for one sample (cell).
#' @param value 	The values to be shown, with each value representing one sample (cell) with the same order as coord.
#' @param emphasize 	The Indices of samples (cells) to be emphasized. When it is set, the colors of all the other cells are set to be #bdbdbd30.
#' @param col 	The customized colors of dots. Either col or value must be given.
#' @param ... 	Other arguments passing to the plot function.
#' @export
plotFeature <- function(coord, values = NULL, emphasize = NULL, colorScheme="gyr", col = NULL, ...){
	if (is.null(col)){
		if (is.numeric(values)){
			if(colorScheme=="gyr"){
				colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
			}else if (colorScheme=="beach"){
				colorPal <- grDevices::colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
			}else if (colorScheme=="pink"){
				colorPal <- grDevices::colorRampPalette(c("#cdcdcd", "#cdcdcd", "#fcc5c0", "#dd3497","#49006a"))
			}else if (colorScheme=="bo"){
				colorPal <- grDevices::colorRampPalette(c("navy", "darkorange1"))
			}else if (colorScheme=="gy"){
				colorPal <- grDevices::colorRampPalette(c("#ffffe5","#f7fcb9","#d9f0a3", "#addd8e", "#78c679","#41ab5d","#238443"))
			}else if(colorScheme=="rb1"){
				colorPal <- grDevices::colorRampPalette(c("#a8c1c8", "#e9d9e5", "#e7c5dc", "#de97b1", "#ee538d"))
			}else if (colorScheme=="rb2"){
				colorPal <- grDevices::colorRampPalette(c("#49616d", "#ffccb0", "#ffd6de", "#e56a6b", "#a80500"))
			}else if (colorScheme=="pink2"){
				colorPal <- grDevices::colorRampPalette(c("#2d0000", "#cec2cd", "#de91c5", "#e854a5", "#8b0043"))
			}else if (colorScheme=="gy2"){
				colorPal <- grDevices::colorRampPalette(c("#043007", "#005c13", "#e2e2eb", "#ffd500", "#e59800"))
			}else if (colorScheme=="rg"){
				colorPal <- grDevices::colorRampPalette(c("#4a9590", "#bcdbdb", "#bd6578", "#b72b3d","#540101"))
			}else if (colorScheme=="purple"){
				colorPal <- grDevices::colorRampPalette(c("#e6e1dd", "#c9abc3", "#864d7f", "#4e1e4a"))
			}else if (colorScheme=="by"){
				colorPal <- grDevices::colorRampPalette(c("#000066", "#006699", "#006666", "#FFDE00","#CC9900"))
			}
			cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F,include.lowest=T))]
			if(min(values, na.rm=T) == 0) cellColor[values == 0] <- "#bdbdbd30"
		}else{
			if (colorScheme=="ggplot"){
				cols <- setNames(scales::hue_pal()(length(unique(values))), unique(values))
			}else if (colorScheme=="jet"){
				colorPal <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
				cols <- setNames(colorPal(length(unique(values))), unique(values))
			}
			
			cellColor <- cols[values]
		}
	} else{
		if (length(col) == 1) col <- rep(col, nrow(coord))
		cellColor <- col
	}
	
	if (is.null(emphasize)){
		plot(coord, col = cellColor, pch=16, ...)
	} else{
		plot(coord, col = "#bdbdbd30", pch=16, ...)
		points(coord[emphasize, ], col = cellColor[emphasize], pch=16, ...)
	}
}

plotFeature2 <- function(coor, values, knn.pairs=NULL, emphasize=NULL, nCols=NULL, gCols=NULL, zeroAsGray=T, ...){
	if(is.numeric(values)){
		if(is.null(nCols)){
			colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
		}else{
			colorPal <- grDevices::colorRampPalette(nCols)
		}
		cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F,include.lowest=T))]
		if(zeroAsGray){
			if(min(values, na.rm=T) == 0) cellColor[values == 0] <- "#bdbdbd30"
		}
	}else{
		if(is.null(gCols)){
			gCols <- setNames(scales::hue_pal()(length(unique(values))), unique(values))
		}
		cellColor <- gCols[values]
	}
	if(is.null(emphasize)){
		plot(coor, col=cellColor, pch=16, ...)
	}else{
		plot(coor, col="#bdbdbd30", pch=16, ...)
		points(coor[emphasize,], col=cellColor[emphasize], pch=16, ...)
	}
	
	if(!is.null(knn.pairs)){
		for(i in seq(nrow(knn.pairs))){
			x.coor <- coor[knn.pairs[i,],1]
			y.coor <- coor[knn.pairs[i,],2]
			lines(x.coor, y.coor, col="#bdbdbd30", lwd=0.5, ...)
		}
		if(is.null(emphasize)){
			points(coor, col=cellColor, pch=16, ...)
		}else{
			points(coor[emphasize,], col=cellColor[emphasize], pch=16, ...)
		}
		
	}

}

getAveExpr <- function(meta.data, feature.to.calc, specified.order=NULL, feature.type="numeric", expr, genes=NULL, core.num=5){
	if(sum(rownames(meta.data)!=colnames(expr))>0){
		return(cat("Inconsistent meta.data and expression sample names\n"))
	}
	if(is.null(genes)){
		genes=rownames(expr)
	}
	library(doParallel)
	registerDoParallel(core.num)
	if(is.null(specified.order)){
		if(feature.type=="numeric"){
			sorted.features <- sort(as.numeric(unique(meta.data[,feature.to.calc])))
		}else{
			sorted.features <- sort(unique(meta.data[,feature.to.calc]))
		}
	}else{
		sorted.features <- specified.order
	}
	ave.expr <- foreach(k=sorted.features, .multicombine=T, .combine='cbind')%dopar%{
		sample.idx <- which(meta.data[,feature.to.calc]==k)
		apply(expr[genes, sample.idx], 1, mean)
	}
	stopImplicitCluster()
	rownames(ave.expr) <- genes
	if(feature.type=="numeric"){
		sorted.features <- paste0("Cluster", sorted.features)
	}
	colnames(ave.expr) <- sorted.features
	return(ave.expr)
} 

findAllMarkers <- function(seu.obj, selected.column, core.num=5, p.value.cutoff=1, bg=NULL, cluster.to.test=NULL, pos.only=T){
	meta.data <- seu.obj@meta.data
	cluster.num <- max(as.numeric(meta.data[,selected.column]))
	if(is.null(cluster.to.test)){
		cluster.to.test <- 0:cluster.num
	}
	if(is.null(bg)){
		bg=0:cluster.num
	}
	library(doParallel)
	registerDoParallel(core.num)
	combined.cm <- foreach(k=cluster.to.test, .multicombine=T, .combine='rbind')%dopar%{
		bg.2 = setdiff(bg, k)
		mat <- FindMarkers(seu.obj, ident.1=k, ident.2=bg.2, only.pos=pos.only, min.pct=0.25, logfc.threshold=0.25)
		dat <- data.frame(mat, "cluster"=rep(k, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
		rownames(dat) <- paste(rownames(dat), k, sep="_")
		return(dat)
	}
	stopImplicitCluster()
	combined.cm <- combined.cm[which(combined.cm$p_val<p.value.cutoff),]
	return(combined.cm)
}
