#' Extract scaffolds connected by paired reads
#'
#' A nice long description
#'
#' @usage extract_nodes(original, subset)
#'
#' @param original The original dataframe with all scaffolds (d).
#' @param subset The subset of extracted scaffolds.
#' 
#' @return A list with nodes and scaffolds
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
extract_nodes <- function (original, subset) {
  clusters<-as.numeric(as.character(unique(na.omit(subset$cluster))))
  original.clusters<-original[original$cluster %in% clusters,]
  nodes<-as.numeric(as.character(na.omit(original.clusters$vertex.id)))
  scaffolds <- unique(c(subset$scaffold, original.clusters$scaffold))
  out <- list(nodes, scaffolds)
  names(out) <- c("nodes", "scaffolds")
  return(out)
} 

