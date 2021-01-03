
#' Interactive visualization of the simulated tumor
#' 
#' @description Interactive visualization of the simulated tumor using the \code{rgl} package (if available). 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}()}.
#' @param plot.type Which type of plot to draw. "Normal" assigns a random rgb value to each genotype while
#' "heat" colors cells with more mutations red and cells with fewer mutations blue. 
#' @param background If rgl is installed, this will set the color of the background
#' @param axes Will include axes (rgl only). 
#' 
#' @details If \pkg{rgl} is installed, then the plots will be interactive. If \pkg{rgl} is unavailable, static plots will
#' be created with \pkg{scatterplot3d}. Since plotting performance with \pkg{scatterplot3d} is reduced, it is strongly
#' recommended that \pkg{rgl} is installed for optimal use of this function. 
#' 
#' @return None. 
#' 
#' @author Phillip B. Nicol 
#' 
#' 
visualizeTumor <- function(tumor, plot.type = "normal", background = "black", axes = FALSE) {
  if(!requireNamespace("rgl", quietly = TRUE)) {
    warning("Installing package 'rgl' is recommended for interactive visualization. Feautures and performance limited.")
    vistum_scatter(tumor, plot.type = plot.type)
  } else{
    min_x <- min(tumor$cell_ids[,1])
    max_x <- max(tumor$cell_ids[,1])
    min_y <- min(tumor$cell_ids[,2])
    max_y <- max(tumor$cell_ids[,2])
    min_z <- min(tumor$cell_ids[,3])
    max_z <- max(tumor$cell_ids[,3])
  
    if(plot.type == "heat") {
    col.pal <- colorRampPalette(c("blue", "red"))
    hotcold <- col.pal(max(tumor$cell_ids$nmuts) + 1)
    rgl::open3d()
    rgl::bg3d(background)
    rgl::plot3d(tumor$cell_ids[,1], tumor$cell_ids[,2], tumor$cell_ids[,3], 
            col = hotcold[tumor$cell_ids$nmuts+1], box = axes, axes = axes,
            xlim = c(min_x - 10, max_x + 10), ylim = c(min_y - 10, max_y + 10), size = 5,
            zlim = c(min_z - 10, max_z + 10), xlab = ifelse(axes,"X",""), ylab = ifelse(axes,"Y",""), 
            zlab = ifelse(axes,"Z",""))
    } else if(plot.type == "normal") {
      rgl::open3d()
      rgl::bg3d(background)
      rgl::plot3d(tumor$cell_ids[,1], tumor$cell_ids[,2], tumor$cell_ids[,3], 
          col = tumor$color_scheme[tumor$cell_ids[,4]+1],  box = axes, axes = axes,
          xlim = c(min_x - 10, max_x + 10), ylim = c(min_y - 10, max_y + 10), size = 5,
          zlim = c(min_z - 10, max_z + 10), xlab = ifelse(axes,"X",""), ylab = ifelse(axes,"Y",""), 
          zlab = ifelse(axes,"Z",""))
    }
  }
}

#' 2D cross section of the simulated tumor
#' 
#' @description 2D cross section of the simulated tumor.
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}()}.
#' @param slice.dim One of "x", "y" or "z", which denotes the dimension which will be fixed to obtain a 2D cross section.
#' @param level Which value will the dimension given in \code{slice.dim} be fixed at? 
#' @param plot.type Which type of plot to draw. "Normal" assigns a random rgb value to each genotype while
#' "heat" colors cells with more mutations red and cells with fewer mutations blue. This is exactly the same as \code{plot.type}
#' in \code{visualizeTumor}. 
#' 
#' @return None. 
#' 
#' @author Phillip B. Nicol 
#' 
#' 
plotSlice <- function(tumor, slice.dim = "x", level = 0, plot.type = "normal") {
  slice <- switch(slice.dim, "x" = 1, "y" = 2, "z" = 3)
  
  df_slice <- tumor$cell_ids[tumor$cell_ids[,slice] == level,-slice]
  if(plot.type == "heat") {
    col.pal <- colorRampPalette(c("blue", "red"))
    hotcold <- col.pal(max(tumor$cell_ids$nmuts) + 1)
    plot(df_slice[,1], df_slice[,2], col = hotcold[df_slice$nmuts+1], pch = 16,
         xlab = NA, ylab = NA)
    
  }
  else if(plot.type == "normal") {
    plot(df_slice[,1], df_slice[,2], col = tumor$color_scheme[df_slice$genotype+1], pch = 16,
         xlab = NA, ylab = NA)
  }    
}

vistum_scatter <- function(tumor, plot.type = "normal") {
  min_x <- min(tumor$cell_ids[,1])
  max_x <- max(tumor$cell_ids[,1])
  min_y <- min(tumor$cell_ids[,2])
  max_y <- max(tumor$cell_ids[,2])
  min_z <- min(tumor$cell_ids[,3])
  max_z <- max(tumor$cell_ids[,3])
  
  if(plot.type == "heat") {
    col.pal <- colorRampPalette(c("blue", "red"))
    hotcold <- col.pal(max(tumor$cell_ids$nmuts) + 1)
    scatterplot3d::scatterplot3d(tumor$cell_ids[,1], tumor$cell_ids[,2], tumor$cell_ids[,3], 
                color = hotcold[tumor$cell_ids$nmuts+1], 
                xlim = c(min_x - 10, max_x + 10), ylim = c(min_y - 10, max_y + 10), pch = 16,
                zlim = c(min_z - 10, max_z + 10), xlab = "x", ylab = "y", zlab = "z")
  }
  else if(plot.type == "normal") {
    scatterplot3d::scatterplot3d(tumor$cell_ids[,1], tumor$cell_ids[,2], tumor$cell_ids[,3], 
                                 color = tumor$color_scheme[tumor$cell_ids[,4]+1], 
                                 xlim = c(min_x - 10, max_x + 10), ylim = c(min_y - 10, max_y + 10), pch = 16,
                                 zlim = c(min_z - 10, max_z + 10), xlab = "x", ylab = "y", zlab = "z")
  }  
}