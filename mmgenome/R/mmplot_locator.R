#' Interactively define a selection subspace on mmplot or mmplot_network
#'
#' Interactively define a selection subspace on \code{mmplot} and \code{mmplot_network}. The selection can be used to extract subsets of mmgenome data objects using \code{mmextract}. The function works on generic ggplot2 objects,  
#'
#' @usage mmplot_locator(p, sig)
#'
#' @param p (required) A mmplot or mmplot_network object.
#' @param sig The number of significant digits (default: 3).
#' @param trans Specify x and y-axis transformation used in mmplot. The format is trans=c(x-axis, y-axis) i.e. trans=c("", "log10").  The following transformations are supported: No transformation (default) = "". Log10 = "log10". Square root transformation = "sqrt".  
#' 
#' @return X and Y coordinates marked in the plot or scaffold names.
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
#' p <- mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "pps_phylum", minlength = 5000)
#' sel <- mmplot_locator(p, trans=c("log10", "log10"))
#' sel <- data.frame(C13.12.03  =  c(1.39, 2.07, 16.8, 19.4, 7.72, 1.76),
#'                   C14.01.09  =  c(29.4, 67.6, 85.9, 43.6, 16.7, 14.9))
#' mmplot_selection(p, sel)
#' 
#' dA <- mmextract(d, sel)
#' 
#' pn <- mmplot_network(data = dA, network = pe, color = "essential", nconnections = 10)
#' sel <- mmplot_locator(pn)
#' sel <- c("3951", "6944")
#' mmplot_network(data = dA, network = pe, color = "essential", nconnections = 10, highlight= sel)
#' }

mmplot_locator <- function(p, sig=3, trans=c("","")){
  ### Arguments
  point_max <- 100
  plot_type <- attr(p, "comment")
  if (is.null(plot_type)) {plot_type <- "other"}
  trans_x = trans[1]
  trans_y = trans[2]

  ### Build ggplot object
  suppressWarnings(ggobj <- ggplot_build(p))
  
  ### Extract coordinates
  xr <- ggobj$layout$panel_ranges[[1]]$x.range
  yr <-ggobj$layout$panel_ranges[[1]]$y.range
  
  ### Variable for selected points
  sp <- data.frame(x = as.numeric(), y = as.numeric())
  colnames(sp) <- c(ggobj$plot$mapping$x, ggobj$plot$mapping$y)
  
  ### Detect and move to plot area viewport
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
  
  ### Select point, plot, store and repeat
  for (i in 1:point_max){
    tmp <- grid.locator('native')
    if (is.null(tmp))break
    grid.points(tmp$x,tmp$y, pch = 16, gp=gpar(cex=0.5, col="darkred"))
    sp[i, ] <- as.numeric(tmp)
  }
  grid.polygon(x= unit(sp[,1], "native"), y= unit(sp[,2], "native"), gp=gpar(fill=NA))
  
  ### Scale data values
  sp[,1] <- (max(xr)- min(xr)) * sp[,1] + min(xr)
  sp[,2] <- (max(yr)- min(yr)) * sp[,2] + min(yr)
  if (trans_x == "log10" & plot_type != "mmplot_network") sp[,1] <- 10^(sp[,1])
  if (trans_y == "log10" & plot_type != "mmplot_network") sp[,2] <- 10^(sp[,2])
  if (trans_x == "sqrt" & plot_type != "mmplot_network") sp[,1] <- (sp[1, ])^(sp[,1])
  if (trans_y == "sqrt" & plot_type != "mmplot_network") sp[,2] <- (sp[2, ])^(sp[,2])
  if (!(trans_x %in% c("", "log10", "sqrt")) & plot_type != "mmplot_network") warning("x-axis transformation [", trans_x, "] not recognized. Defaulting to no transformation.")
  if (!(trans_y %in% c("", "log10", "sqrt")) & plot_type != "mmplot_network") warning("y-axis transformation [", trans_y, "] not recognized. Defaulting to no transformation.")
  if (!(trans_y == "" | trans_x == "") & plot_type == "mmplot_network") warning("Axis transformation of mmplot_network not supported. Defaulting to no transformation")

  
  ### Return output based on plottype
  if(plot_type == "mmplot_network"){
    # Extract scaffold names
    in.selection <- point.in.polygon(ggobj$plot$data$x, ggobj$plot$data$y, sp[,1], sp[,2], mode.checked=T) > 0
    scaffolds <- as.character(ggobj$plot$data$scaffold[in.selection])
    dput(scaffolds)
    return(scaffolds)
  } else if (plot_type == "mmplot"){
    # Write values to console
    show(paste(colnames(sp)[1], " = ", list(signif(sp[,1], sig))))
    show(paste(colnames(sp)[2], " = ", list(signif(sp[,2], sig))))
    # Return selected points dataframe
    return(sp)
  } else {
    warning("Plot not recognized as mmplot or mmplot_network. mmplot_locator might not work as intended.")
    # Write values to console
    show(paste(colnames(sp)[1], " = ", list(signif(sp[,1], sig))))
    show(paste(colnames(sp)[2], " = ", list(signif(sp[,2], sig))))
    # Return selected points dataframe
    return(sp)
  }
}
