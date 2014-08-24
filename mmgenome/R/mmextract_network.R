#' Extract scaffolds connected by paired-end reads
#'
#' Extracts scaffolds that are connected by paired-end or mate-pair reads.
#'
#' @usage mmextract_network(subset, network, original, nconnections, type)
#'
#' @param subset The subset dataframe to use for network extraction or a vector of scaffold names.
#' @param network The network of connected scaffold.
#' @param original The original dataframe (default: d).
#' @param nconnections The minumum number of connections needed (default: 2).
#' @param type Extract "direct" or "complete" (all) linked scaffolds (default: direct). 
#' 
#' @return A dataframe with scaffold information.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' 
#' p <- mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "phylum", minlength = 3000)
#' sel <- data.frame(C13.12.03  =  c(4.94, 8.99, 10.8, 7.22, 5.79),
#'                   C14.01.09  =  c(43.9, 57.4, 35.7, 25.7, 32.6))
#' mmplot_selection(p, sel) 
#' 
#' dA <- mmextract(d, sel)
#' mmplot(data = dA, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "phylum", minlength = 3000)
#' mmplot_network(data = dA, network = pe, color = "phylum", nconnections = 10)
#' 
#' dB <- mmextract_network(subset = dA, original = d, network = pe, nconnections = 10, type = "direct")
#' mmplot(data = dB, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "phylum", minlength = 3000)
#' mmplot_network(data = dB, network = pe, color = "phylum", nconnections = 10) 
#' 
#' }

mmextract_network <- function(subset, network, original = d, nconnections = 2, type = "direct"){
  
  if (class(subset) != "list"){
    subset <- list(scaffolds = data.frame(scaffold = subset))
  }
  
  
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
    ext.clusters <- subset(clusters, clusters$scaffold %in% subset$scaffolds$scaffold)
    ext.scaffolds <- subset(clusters, clusters$cluster %in% ext.clusters$cluster)
    out <- subset(original$scaffolds, original$scaffolds$scaffold %in% subset$scaffold | original$scaffolds$scaffold %in% ext.scaffolds$scaffold)
  }

  if (length(subset) == 2){es <- subset(original$essential, original$essential$scaffold %in% out$scaffold)}
  
  if (length(subset) == 2){
    outlist <- list(scaffolds = out, essential = es)
  } else {
    outlist <- list(scaffolds = out)
  }  
  
  return(outlist)  
}