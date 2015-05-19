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
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
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

mmplot_locator <- function(p, sig=3){
  pi <- ggplot_build(p)
  x <- unlist(lapply(pi$data, function(l){l$x}))
  y <- unlist(lapply(pi$data, function(l){l$y})) 
  x <- x[!is.infinite(x)]
  y <- y[!is.infinite(y)]
  n = 100
  d <- data.frame(matrix(as.numeric(), 0, 2))
  colnames(d) <- c(pi$plot$mapping$x, pi$plot$mapping$y)
  
  for (i in 1:n){    
    seekViewport('panel.3-4-3-4', recording=TRUE) 
    pushViewport(dataViewport(x,y))
    tmp <- grid.locator('native')
    if (is.null(tmp))break
    grid.points(tmp$x,tmp$y, pch = 16, gp=gpar(cex=0.5, col="darkred"))
    d[i, ] <- as.numeric(tmp)
  }
  grid.polygon(x= unit(d[,1], "native"), y= unit(d[,2], "native"), gp=gpar(fill=NA))
  
  if (pi$panel$x_scales[[1]]$trans$name == "log-10") d[,1] <- 10^(d[,1])
  if (pi$panel$y_scales[[1]]$trans$name == "log-10") d[,2] <- 10^(d[,2])
  cat(paste0(colnames(d)[1]," = ",list(signif(d[,1],sig)),',\n'))
  cat(paste0(colnames(d)[2]," = ",list(signif(d[,2],sig)),'\n\n'))
  return(d)
}
