#' Identifies clusters of scaffolds from a graph
#'
#' A nice long description
#'
#' @usage merge_network(graph, data)
#'
#' @param graph A igraph object.
#' @param data Data frame of scaffold information.
#' 
#' @return A data frame with scaffold cluster information attached.
#' 
#' @export
#' @import igraph
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


merge_network <- function (graph, data) {  
  g.clust <- clusters(graph)
  V(graph)$membership <- g.clust$membership  
  tm <- cbind(1:length(V(graph)),g.clust$membership)
  
  colnames(tm) <- c("id","cluster")
  ts <- cbind(1:g.clust$no, g.clust$csize)
  colnames(ts) <- c("cluster", "size")
  ta <- merge(x=tm, y=ts, by="cluster", all = T)
  ta <- ta[with(ta, order(id)), ]
  V(graph)$csize <- ta$size
  tc <- cbind.data.frame(1:length(V(graph)),as.numeric(V(graph)$name),V(graph)$membership,V(graph)$csize)
  colnames(tc) <- c("vertex.id","scaffold","cluster","cluster.size")
  d<- merge(x = data, y = tc , by = "scaffold", all = T)
  return(d)
}
