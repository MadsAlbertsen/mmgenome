#' Plots connected scaffolds in a network
#'
#' Plots connected scaffolds in a network using igraph and ggplot2.
#'
#' @usage mmplot_network(data, network)
#'
#' @param data (required) The dataframe with the network to be plotted.
#' @param network (required) A network file with connections between scaffolds.
#' @param nconnections The minimum number of connections (default: 2).
#' @param color Color by gc, phylum or a specific coverage dataset (default: "gc").
#' @param log.color log10 scale the colors (default: F)
#' @param labels If scaffold names are to be plotted (default: F).
#' @param highlight Mark selected scaffolds on the plot. Either as a vector of scaffold names or as a full subset of data.
#' @param hightlight.color Color of the highlighted scaffolds (default: "darkred").
#' @param print.nolinks Print scaffolds with no links to the console (default: F).
#' @param links.scale Scale the width of the links between scaffolds by a constant (default: 1).
#' @param seed Set the seed to obtain reproducible graph layout. Usefull when using mmplot_locator.
#' 
#' @return A ggplot2 object.
#' 
#' @import ggplot2
#' @import igraph
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' 
#' p <- mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "essential", minlength = 3000)
#' sel <- data.frame(C13.12.03  =  c(4.94, 8.99, 10.8, 7.22, 5.79),
#'                   C14.01.09  =  c(43.9, 57.4, 35.7, 25.7, 32.6))
#' mmplot_selection(p, sel) 
#' 
#' dA <- mmextract(d, sel)
#' mmplot(data = dA, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "essential", minlength = 3000)
#' mmplot_network(data = dA, network = pe, color = "essential", nconnections = 10)
#' 
#' dB <- mmextract_network(subset = dA, original = d, network = pe, nconnections = 10, type = "direct")
#' mmplot(data = dB, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "essential", minlength = 3000)
#' mmplot_network(data = dB, network = pe, color = "essential", nconnections = 10) 
#' 
#' }


mmplot_network <- function(data, network, nconnections = 2, labels = F, color = "gc", log.color = F, highlight = NULL, highlight.color = "darkred", print.nolinks = F, scale.links = 1, seed = NULL){

  ## Subset the network
  snetwork <- subset(network, network$scaffold1 %in% data$scaffolds$scaffold & network$scaffold2 %in% data$scaffolds$scaffold & network$connections >= nconnections)
  
  ## Convert to graph 
  g <- graph.data.frame(snetwork, directed = F)
  
  ## Calculate a layout
  if (!is.null(seed)) set.seed(seed)
  t <- layout_with_fr(g)  

  ## Extract layout coordinates
  gpoints <- data.frame( "scaffold" = V(g)$name, "x" = t[,1], "y" = t[,2])
  gpoints1 <- merge(gpoints, data$scaffolds)
  
  ## Extract link coordinates 
  links <- merge(snetwork, gpoints1[,1:3], by.x = "scaffold1", by.y = "scaffold")
  links1 <- merge(links, gpoints1[,1:3], by.x = "scaffold2", by.y = "scaffold")
  colnames(links1)[4:7] <- c("x", "y", "xend", "yend")
  
  ## Plot the data
  if (class(data$scaffolds[,color]) == "factor"){
    p <- ggplot(data = gpoints1, aes_string(x = "x", y = "y", size = "length", color = color)) +
      geom_segment(data=links1, aes(x=x, y=y, xend=xend, yend=yend), color = "darkgrey", size = log10(links1$connections)*scale.links) +
      geom_point(alpha=0.1, color = 'black') +
      geom_point(data=subset(gpoints1, gpoints1[,color] != "NA"), shape = 1) +
      scale_size_area(name= "Scaffold length", max_size=20) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
  }
  
  
  if (class(data$scaffolds[,color]) != "factor"){
    p <- ggplot(data = gpoints1, aes_string(x = "x", y = "y", size = "length", color = color)) +
      geom_segment(data=links1, aes(x=x, y=y, xend=xend, yend=yend), color = "darkgrey", size = log10(links1$connections)*scale.links) +
      geom_point(alpha=0.7) +
      scale_size_area(name= "Scaffold length", max_size=20)
      if (log.color == F){p <- p + scale_colour_gradientn(colours = c("red", "green", "blue"))}
    if (log.color == T){p <- p + scale_colour_gradientn(colours = c("red", "green", "blue"), trans = "log10")}
  }
  
  if (labels == T){
    p <- p + geom_text(label = gpoints1$scaffold, size = 4, color = "black")
  }
  
  p <- p + 
    #theme_bw() + 
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank(),
          panel.background = element_blank(),
          legend.key = element_blank())
  
  ### Highlight selected scaffolds
  if (!is.null(highlight)){
    if (class(highlight) != "list"){
      highlight <- as.character(highlight)
      sdata <- subset(gpoints1, scaffold %in% highlight)
      p <- p + geom_text(data = sdata, color = highlight.color, size = 4, label = sdata$scaffold)  
    } else{
      sdata <- subset(gpoints1, scaffold %in% highlight$scaffolds$scaffold)
      p <- p + geom_text(data = sdata, color = highlight.color, size = 4, label = sdata$scaffold)  
    }
  }
  
  if(print.nolinks == T){
    nolinks <- data$scaffolds$scaffold[!(data$scaffolds$scaffold %in% gpoints$scaffold)]
    print("The following scaffolds have no links to other scaffolds:")
    print(nolinks)
  }
  
  ### Add identifier-tag for utility functions
  attr(p, "comment") <- "mmplot_network"
  
  ### Output
  return(p)
}
