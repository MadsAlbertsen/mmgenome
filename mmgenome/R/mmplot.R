#' Function to visualise metagenomes
#' 
#' Versitile plotting function to visualise metagnome assemblies loaded through \code{mmload}. Wraps ggplot2 for pretty default plots.
#'
#' @usage mmplot(data, x, y)
#'
#' @param data (required) A dataframe with scaffold information.
#' @param x (required) x-axis variable.
#' @param y (required) y-axis variable.
#' @param color Color by a specific variable or "none" (default: "essential").
#' @param log.x log10 scale the x-axis (default: F).
#' @param log.y log10 scale the y-axis (default: F).
#' @param log.color log10 scale the colors (default: F).
#' @param minlength Minimum length of plotted scaffolds.
#' @param network Network used to plot connections between scaffolds.
#' @param nconnections Minimum number of connections to plot (default: 0).
#' @param duplicates Mark scaffolds with duplicated essential genes (default: F).
#' @param labels If scaffold names are to be plotted (default: F).
#' @param resize Constant to rescale the size of the scaffolds (default: 1).
#' @param point.size Use a fixed size for points instead of scaffold length.
#' @param highlight Mark selected scaffolds on the plot. Either as a vector of scaffold names or as a full subset of data.
#' @param hightlight.color Color of the highlighted scaffolds (default: "darkred").
#' @param alpha Transparency of the plotted points.
#' @param factor.shape Plot factor colors as "outline" or "solid" shapes (default: "outline").
#' @param esom.map Contour map of the esom analysis.
#' @param esom.cut Cutoff value for the esom map data.
#' @param esom.color Color of the esom contour boundaries.
#' 
#' @return a ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "essential", minlength = 10000)
#' }

mmplot <- function(data, x, y, log.x=F, log.y=F, color = "essential", minlength = NULL, network = NULL, nconnections = 0, duplicates = F, labels = F, log.color = F,  resize = 1, point.size = NULL, highlight = NULL, highlight.color = "darkred", alpha = NULL, esom.map = NULL, esom.cut = 0.15, esom.color = "darkred", factor.shape = "outline"){
  
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
  
  ## Extract duplicates
  
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
  
  ## Plot the data
  
  ### Colors: none
  if (color == "none"){
    if (is.null(alpha)){alpha <- 0.1}
    p <- ggplot(data=data$scaffolds, aes_string(x = x, y = y, size = "length"))              
    if (is.null(point.size)){
      p <- p + geom_point(alpha = alpha, color = "black") +
              scale_size_area(name = "Scaffold length", max_size = 20*resize)
    } else {
      p <- p + geom_point(alpha = alpha, color = "black", size = point.size)
    }
  } else {
  ### Colors: factors    
    if (class(data$scaffolds[,color]) == "factor"){
      if(factor.shape == "solid"){fs <- 16} else {fs <- 1}
      
      if (is.null(alpha)){alpha <- 0.1}
      p <- ggplot(data=data$scaffolds, aes_string(x = x, y = y, size = "length", color = color))
      if (is.null(point.size)){
        p <- p + geom_point(alpha=alpha, color = 'black') +
                 geom_point(data=subset(data$scaffolds, data$scaffolds[, color] != "NA"), shape = fs, alpha = 0.7) +
                 scale_size_area(name = "Scaffold length", max_size = 20*resize)
      } else{
        p <- p + geom_point(alpha=alpha, color = 'black', size = point.size) +
                 geom_point(data=subset(data$scaffolds, data$scaffolds[, color] != "NA"), shape = fs, alpha = 0.7, size = point.size)
      }
      p <- p + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
    }
  
  ### Colors: numeric
    if (class(data$scaffolds[,color]) != "factor"){
      if (is.null(alpha)){alpha <- 0.3}
      options(digits=2)
      p <- ggplot(data=data$scaffolds, aes_string(x = x, y = y, size = "length", color = color))
      if (is.null(point.size)){
        p <- p + geom_point(alpha = alpha) +
                 scale_size_area(name = "Scaffold length", max_size = 20*resize)
      } else{
        p <- p + geom_point(alpha = alpha, size = point.size)
      }
      if (log.color == F){p <- p + scale_colour_gradientn(colours = c("red", "green", "blue"))}
      if (log.color == T){p <- p + scale_colour_gradientn(colours = c("red", "green", "blue"), trans = "log10")} 
    }
  }
    
  ### Scale axis
  if (log.x == T){p <- p + scale_x_log10()}
  if (log.y == T){p <- p + scale_y_log10()}
  
  ### Plot connections between scaffolds
  if (!is.null(network)){
    p <- p +  
      geom_segment(data = links, aes(x = x, y = y, xend = xend, yend = yend), color = "darkgrey", size = 1, alpha = 0.5) +
      geom_point(data = links, aes(x = x, y = y), size = 2, color = "darkgrey") +
      geom_point(data = links, aes(x = xend, y = yend), size = 2, color = "darkgrey")
  }
  
  ### Plot duplicates
  if (duplicates == T){
    p <- p +
      geom_segment(data = dlinks, aes(x = x, y = y, xend = xend, yend = yend), color = "darkred", size = 1) +
      geom_point(data = dlinks, aes(x = x, y = y), size = 2, color = "darkred") +
      geom_point(data = dlinks, aes(x = xend, y = yend), size = 2, color = "darkred")
  }
  
  ### Plot labels
  if (labels == T){
    p <- p + geom_text(label = data$scaffolds$scaffold, size = 4, color = "black")
  }
  
  ### Highlight selected scaffolds
  if (!is.null(highlight)){
    if (class(highlight) != "list"){
      highlight <- as.character(highlight)
      sdata <- subset(data$scaffolds, scaffold %in% highlight)
      p <- p + geom_point(data = sdata, color = highlight.color, shape = 1)  
    } else{
      p <- p + geom_point(data = highlight$scaffold, color = highlight.color, shape = 1)  
    }
  }
  
  p <- p + theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey95"),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            legend.key = element_blank()
            )
  
  if(!is.null(esom.map)){
    es <- subset(esom.map, density > esom.cut)
    p <- p + annotate(geom = "point", x = es$esomx, y = es$esomy, size = 2, color = esom.color) +
      theme(panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  }
  
  ### Add identifier-tag for utility functions
  attr(p, "comment") <- "mmplot"
  
  ### Output
  return(p)
}
