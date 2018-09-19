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

plotFeature2 <- function(coor, values, knn.pairs=NULL, emphasize=NULL, nCols=NULL, gCols=NULL, zeroAsGray=T, ...){
  if(is.numeric(values)){
    if(is.null(nCols)){
      colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
    }else{
      colorPal <- grDevices::colorRampPalette(nCols)
    }
    cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F,include.lowest=T))]
    if(zeroAsGray){
#      if(min(values, na.rm=T) == 0) cellColor[values == 0] <- "#bdbdbd"
      if(min(values, na.rm=T) == 0) cellColor[values == 0] <- "#bdbdbd"
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
    points(coor[emphasize,], col=cellColor[emphasize], pch=16)
  }
  
  if(!is.null(knn.pairs)){
    for(i in seq(nrow(knn.pairs))){
      x.coor <- coor[knn.pairs[i,],1]
      y.coor <- coor[knn.pairs[i,],2]
      lines(x.coor, y.coor, col="#bdbdbd30", lwd=0.5)
    }
    if(is.null(emphasize)){
      points(coor, col=cellColor, pch=16)
    }else{
      points(coor[emphasize,], col=cellColor[emphasize], pch=16)
    }
    
  }
  
}
