#' Wraps ggplot2 for pretty default plots
#'
#' A nice long description.
#'
#' @usage mmplot(data, x, y, logx, logy, color)
#'
#' @param data A dataframe with scaffold information.
#' @param x x-axis variable.
#' @param y y-axis variable.
#' @param logx log10 scale the x-axis (default True)
#' @param logy log10 scale the y-axis (default True)
#' @param color Color by "phylum", "gc" or a coverage dataset (default: phylum)
#' 
#' @return a ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmplot <- function(data, x, y, logx=F, logy=F, color = "phylum"){
  if (color == "phylum"){
    p <- ggplot(data=data, aes_string(x = x, y = y, size = "length", color = "phylum")) + 
      geom_point(alpha=0.1, color = 'black') +
      geom_point(data=subset(data, phylum != "NA"), shape = 1, alpha = 0.7) +
      scale_size_area(name= "Scaffold length", max_size=20) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
  }
  if (color == "gc"){
    p <- ggplot(data=data, aes_string(x = x, y = y, size = "length", color = "gc")) +       
      geom_point(alpha = 0.3) +
      scale_size_area(name = "Scaffold length", max_size = 20) +
      scale_colour_gradientn(colours = c("red", "green", "blue"))
  }
  if(color != "gc" & color != "phylum"){
    options(digits=2)
    p <- ggplot(data=data, aes_string(x = x, y = y, size = "length", color = color)) +       
      geom_point(alpha = 0.3) +
      scale_size_area(name = "Scaffold length", max_size = 20) +
      scale_colour_gradientn(colours = c("red", "green", "blue"), trans = "log")
  }
  
  if (logx == T){p <- p + scale_x_log10()}
  if (logy == T){p <- p + scale_y_log10()}
  return(p)
}
