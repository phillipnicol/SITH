

visualizeTumor <- function(tumor_instance, plot.type = "normal", background = "black") {
  if(!requireNamespace("rgl", quietly = TRUE)) {
    stop("rgl package required for interactive tumor visualization. Please install this package to use this feature.")
  }
  min_x <- min(tumor_instance$cell_ids[,1])
  max_x <- max(tumor_instance$cell_ids[,1])
  min_y <- min(tumor_instance$cell_ids[,2])
  max_y <- max(tumor_instance$cell_ids[,2])
  min_z <- min(tumor_instance$cell_ids[,3])
  max_z <- max(tumor_instance$cell_ids[,3])
  
  if(plot.type == "heat") {
    col.pal <- colorRampPalette(c("blue", "red"))
    hotcold <- col.pal(max(tumor_instance$cell_ids$nmuts) + 1)
    rgl::open3d()
    rgl::bg3d(background)
    rgl::plot3d(tumor_instance$cell_ids[,1], tumor_instance$cell_ids[,2], tumor_instance$cell_ids[,3], 
           col = hotcold[tumor_instance$cell_ids$nmuts+1], box = F, axes = F,
           xlim = c(min_x - 10, max_x + 10), ylim = c(min_y - 10, max_y + 10), size = 5,
           zlim = c(min_z - 10, max_z + 10), pch = 16, xlab = "", ylab = "", zlab = "")
  }
  else if(plot.type == "normal") {
  rgl::open3d()
  rgl::bg3d(background)
  rgl::plot3d(tumor_instance$cell_ids[,1], tumor_instance$cell_ids[,2], tumor_instance$cell_ids[,3], 
         col = tumor_instance$color_scheme[tumor_instance$cell_ids[,4]+1], box = F, axes = F,
         xlim = c(min_x - 10, max_x + 10), ylim = c(min_y - 10, max_y + 10), size = 5,
         zlim = c(min_z - 10, max_z + 10), pch = 16, xlab = "", ylab = "", zlab = "")
  }
}