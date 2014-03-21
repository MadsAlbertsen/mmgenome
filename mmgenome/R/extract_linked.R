#' Extract scaffolds connected by paired-end reads
#'
#' A nice long description.
#'
#' @usage extract_linked(subset, original, network, minimum)
#'
#' @param subset The subset of extracted scaffolds
#' @param original The original dataframe with all scaffolds (default: d)
#' @param network The network of connected scaffold (default: n)
#' @param minumum The minumum number of connections needed (default: 5)
#' 
#' 
#' @return A dataframe with scaffold information.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 


extract_linked <- function(subset, original = d, network = n, minimum = 5){
  ns <- subset(network, (network$scaffold1 %in% subset$scaffold | network$scaffold2 %in% subset$scaffold) & network$connections >= minimum)
  extract <- unique(c(ns$scaffold1, ns$scaffold2, subset$scaffold))
  out <- subset(original, original$scaffold %in% extract)
  return(out)
}