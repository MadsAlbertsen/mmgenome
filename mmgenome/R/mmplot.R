#' Wraps ggplot2 for pretty default plots
#'
#' A nice long description.
#'
#' @usage mmplot(data, x, y)
#'
#' @param data (required) A dataframe with scaffold information.
#' @param x (required) x-axis variable.
#' @param y (required) y-axis variable.
#' @param logx log10 scale the x-axis (default True)
#' @param logy log10 scale the y-axis (default True)
#' @param color Color by "phylum", "gc" or a coverage dataset (default: phylum)
#' @param minlength Minimum length of plotted scaffolds.
#' @param network Network used to plot connections between scaffolds.
#' @param nconnections Minimum number of connections to plot (default: 0).
#' @param duplicates Mark scaffolds with duplicated essential genes (default: F)
#' @param labels If scaffold names are to be plotted (default: F)
#' 
#' @return a ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmplot <- function(data, x, y, logx=F, logy=F, color = "phylum", minlength = NULL, network = NULL, nconnections = 0, duplicates = F, labels = F){
  
  ## Subset based on length constrain
  
  if (!is.null(minlength)){
    data$scaffolds <- subset(data$scaffolds, length >= minlength)
    data$essential <- subset(data$essetinal, data$essential$scaffold %in% data$scaffolds$scaffold)    
  }
  
  ## Extract connections between scaffolds
  if (!is.null(network)) {
    snetwork <- subset(network, network$scaffold1 %in% data$scaffolds$scaffold & network$scaffold2 %in% data$scaffolds$scaffold & connections >= nconnections)
    
    links <- merge(snetwork, data$scaffolds[,c("scaffold",x,y)], by.x = "scaffold1", by.y = "scaffold") 
    links <- merge(links, data$scaffolds[,c("scaffold",x,y)], by.x = "scaffold2", by.y = "scaffold") 
    colnames(links)[4:7] <- c("x","y","xend","yend") 
  }
  
  if (duplicates == T){
    dup<-data$essential[which(duplicated(data$essential$hmm.id) | duplicated(data$essential$hmm.id, fromLast=TRUE)),] 
    dup <- dup[order(dup$hmm.id),c(1,3)]
    
    t2 <- data.frame(scaffold1 = numeric(0), scaffold2 = numeric(0))
    for (i in 1:nrow(dup)){
      for (j in 1:nrow(dup)){
        if (dup[i,1] != dup[j,1] & dup[i,2] == dup[j,2]){
          t2[nrow(t2)+1, ] <- c(dup[i,1], dup[j,1])
        }
      }
    }
    
    dlinks <- merge(t2, data$scaffolds[,c("scaffold",x,y)], by.x = "scaffold1", by.y = "scaffold") 
    dlinks <- merge(dlinks, data$scaffolds[,c("scaffold",x,y)], by.x = "scaffold2", by.y = "scaffold") 
    colnames(dlinks)[3:6] <- c("x","y","xend","yend") 
  }
  
  if (color == "phylum"){
    p <- ggplot(data=data$scaffolds, aes_string(x = x, y = y, size = "length", color = "phylum")) + 
      geom_point(alpha=0.1, color = 'black') +
      geom_point(data=subset(data$scaffolds, phylum != "NA"), shape = 1, alpha = 0.7) +
      scale_size_area(name= "Scaffold length", max_size=20) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
  }
  if (color == "gc"){
    p <- ggplot(data=data$scaffolds, aes_string(x = x, y = y, size = "length", color = "gc")) +       
      geom_point(alpha = 0.3) +
      scale_size_area(name = "Scaffold length", max_size = 20) +
      scale_colour_gradientn(colours = c("red", "green", "blue"))
  }
  if(color != "gc" & color != "phylum"){
    options(digits=2)
    p <- ggplot(data=data$scaffolds, aes_string(x = x, y = y, size = "length", color = color)) +       
      geom_point(alpha = 0.3) +
      scale_size_area(name = "Scaffold length", max_size = 20) +
      scale_colour_gradientn(colours = c("red", "green", "blue"), trans = "log")
  }
  
  if (logx == T){p <- p + scale_x_log10()}
  if (logy == T){p <- p + scale_y_log10()}
  
  if (!is.null(network)){
    p <- p +  
      geom_segment(data = links, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgrey", size = 1) +
      geom_point(data = links, aes(x = x, y = y), size = 2, color = "darkgrey") +
      geom_point(data = links, aes(x = xend, y = yend), size = 2, color = "darkgrey")
  }
  
  if (duplicates == T){
    p <- p +
      geom_segment(data = dlinks, aes(x = x, y = y, xend = xend, yend = yend), color = "darkred", size = 1) +
      geom_point(data = dlinks, aes(x = x, y = y), size = 2, color = "darkred") +
      geom_point(data = dlinks, aes(x = xend, y = yend), size = 2, color = "darkred")
  }
  
  if (labels == T){
    p <- p + geom_text(label = data$scaffolds$scaffold, size = 4, color = "black")
  }
  
  
  return(p)
}
