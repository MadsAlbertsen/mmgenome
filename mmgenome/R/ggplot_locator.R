#' Interactively extract points ggplot2 figures
#'
#' Adds "locator" functionallity to plots made with ggplot2
#'
#' @usage ggplot_locator(p)
#'
#' @param p (required). A ggplot2 object.
#' 
#' @return X and Y coordinates marked in the plot.
#' 
#' @export
#' @import grid
#' @import ggplot2
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

ggplot_locator <- function(p){
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
  show(paste(colnames(d)[1]," = ",list(round(d[,1],4))))
  show(paste(colnames(d)[2]," = ",list(round(d[,2],4))))
  return(d)
}