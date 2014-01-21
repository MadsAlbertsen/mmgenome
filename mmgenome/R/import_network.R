#' Imports a connection network and adds attributes
#'
#' A nice long description
#'
#' @usage import_network(network, data, nconnections, nsize, ewidth)
#'
#' @param network Flat network file.
#' @param data Data frame of scaffold information.
#' @param nconnections Minumum number of connections.
#' @param nsize Scale for node size.
#' @param ewidth Scale for edge width.
#' 
#' @return A igraph object with scaffold attributes attached.
#' 
#' @export
#' @import igraph
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

import_network <- function (network, data , nconnections=0, nsize=.1, ewidth=.3) {
  snetwork <- subset(network, connections >= nconnections)
  g<-graph.data.frame(snetwork, directed = F)  
  g.stats <- data[V(g)$name,]
  V(g)$label <- ""
  V(g)$size <- sqrt(g.stats$length)*nsize
  V(g)$color <- g.stats$gc-min(d$gc)
  V(g)$length <- g.stats$length
  V(g)$frame.color <- "black"
  E(g)$width <- log(E(g)$connections)*ewidth 
  return(g)
}