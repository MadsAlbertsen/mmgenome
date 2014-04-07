#' Extract scaffolds connected by paired-end reads
#'
#' A nice long description.
#'
#' @usage mmextract_network(subset, network, original, nconnections, type)
#'
#' @param subset The subset dataframe to use for extraction.
#' @param network The network of connected scaffold.
#' @param original The original dataframe (default: d).
#' @param nconnections The minumum number of connections needed (default: 2).
#' @param type Extract direct or complete (all) linked scaffolds (default: direct). 
#' 
#' @return A dataframe with scaffold information.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 

mmextract_network <- function(subset, network, original = d, nconnections = 2, type = "direct"){
  
  if (type == "direct"){
    ns <- subset(network, (network$scaffold1 %in% subset$scaffolds$scaffold | network$scaffold2 %in% subset$scaffolds$scaffold) & network$connections >= nconnections)
    out <- subset(original$scaffolds, original$scaffolds$scaffold %in% ns$scaffold1 | original$scaffolds$scaffold %in% ns$scaffold2 | original$scaffolds$scaffold %in% subset$scaffolds$scaffold)
  }
  
  if (type == "complete"){
    snetwork <- subset(network, connections >= nconnections)
    g <- graph.data.frame(snetwork, directed = F)
    g.clust <- clusters(g)
    clusters <- cbind.data.frame(V(g)$name,g.clust$membership)
    colnames(clusters) <- c("scaffold", "cluster")  
    ext.clusters <- subset(clusters, clusters$scaffold %in% subset$scaffold)
    ext.scaffolds <- subset(clusters, clusters$cluster %in% ext.clusters$cluster)
    out <- subset(original$scaffolds, original$scaffolds$scaffold %in% subset$scaffold | original$scaffolds$scaffold %in% ext.scaffolds$scaffold)
  }

  es <- subset(original$essential, original$essential$scaffold %in% out$scaffold)
  
  outlist <- list(scaffolds = out, essential = es)
  
  return(outlist)  
}