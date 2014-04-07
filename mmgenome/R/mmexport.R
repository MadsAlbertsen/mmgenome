#' Exports scaffolds in fasta format
#'
#' A nice long description
#'
#' @usage mmexport(data, assembly, file)
#'
#' @param data The dataframe containing all data.
#' @param assembly The assembly
#' @param file Name of output file.
#' 
#' @return An overview of key statistics of the scaffolds.
#' 
#' @export
#' @import Biostrings
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmexport <- function(data, assembly, file){
  writeXStringSet(assembly[names(assembly) %in% as.character(data$scaffolds$scaffold)], file = file)
}
