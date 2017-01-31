#' Interactively extract points from ggplot2 objects
#'
#' Adds "locator" functionallity to plots made with mmplot (ggplot2). 
#'
#' @usage mmplot_locator(p, sig)
#'
#' @param p (required) A ggplot2 object.
#' @param sig The number of significant digits (default: 3).
#' 
#' @return X and Y coordinates marked in the plot.
#' 
#' @export
#' @import grid
#' @import ggplot2
#' @author Soren M. Karst \email{sorenkarst@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' 
#' p <- mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "phylum", minlength = 3000)
#' p
#' sel <- mmplot_locator(p)
#' sel <- data.frame(C13.12.03  =  c(1.39, 2.07, 16.8, 19.4, 7.72, 1.76),
#'                   C14.01.09  =  c(29.4, 67.6, 85.9, 43.6, 16.7, 14.9))
#' mmplot_selection(p, sel)
#' }

mmplot_locator <- function(p, trans_x="", trans_y = "", sig=3, network = F){
  # Arguments
  point_max <- 100
  
  # Build ggplot object
  suppressWarnings(ggobj <- ggplot_build(p))
  
  # Extract coordinates
  xr <- ggobj$layout$panel_ranges[[1]]$x.range
  yr <-ggobj$layout$panel_ranges[[1]]$y.range
  
  # Variable for selected points
  sp <- data.frame(x = as.numeric(), y = as.numeric())
  colnames(sp) <- c(ggobj$plot$mapping$x, ggobj$plot$mapping$y)
  
  # Plot and move to plot area viewport
  suppressWarnings(print(ggobj$plot))
  panels <- unlist(current.vpTree()) %>%
    grep("panel", ., fixed = TRUE, value = TRUE)
  p_n <- length(panels)
  if (p_n == 0){
    stop("No plot detected")
  }
  if (p_n > 1){
    stop("Multiple-panel plot is not supported")  
  }
  seekViewport(panels, recording=TRUE)
  pushViewport(viewport(width=1, height=1))
  
  # Select point, plot, store and repeat
  for (i in 1:point_max){
    tmp <- grid.locator('native')
    if (is.null(tmp))break
    grid.points(tmp$x,tmp$y, pch = 16, gp=gpar(cex=0.5, col="darkred"))
    sp[i, ] <- as.numeric(tmp)
  }
  grid.polygon(x= unit(sp[,1], "native"), y= unit(sp[,2], "native"), gp=gpar(fill=NA))
  
  # Scale data values
  sp[,1] <- (max(xr)- min(xr)) * sp[,1] + min(xr)
  sp[,2] <- (max(yr)- min(yr)) * sp[,2] + min(yr)
  if (trans_x == "log10" & network == F) sp[,1] <- 10^(sp[,1])
  if (trans_y == "log10" & network == F) sp[,2] <- 10^(sp[,2])
  if (trans_x == "sqrt" & network == F) sp[,1] <- (sp[1, ])^(sp[,1])
  if (trans_y == "sqrt" & network == F) sp[,2] <- (sp[2, ])^(sp[,2])
  
  if (network == F){
    # Write values to console
    show(paste(colnames(sp)[1], " = ", list(signif(sp[,1], sig))))
    show(paste(colnames(sp)[2], " = ", list(signif(sp[,2], sig))))
    # Return selected points dataframe
    return(sp)
  }
  if (network == T){
    # Extract numbers
    in.selection <- point.in.polygon(ggobj$plot$data$x, ggobj$plot$data$y, sp[,1], sp[,2], mode.checked=T) > 0
    scaffolds <- as.character(ggobj$plot$data$scaffold[in.selection])
    dput(scaffolds)
    return(scaffolds)
  }
}
