#' Exports scaffolds in fasta format
#'
#' Exports scaffolds in fasta format using the Biostings package.
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
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' 
#' p <- mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "phylum", minlength = 3000)
#' 
#' sel <- data.frame(C13.12.03  =  c(1.39, 2.07, 16.8, 19.4, 7.72, 1.76),
#'                   C14.01.09  =  c(29.4, 67.6, 85.9, 43.6, 16.7, 14.9))
#' mmplot_selection(p, sel)
#' 
#' dA <- mmextract(d, sel)
#' 
#' mmexport(data=dA, assembly=assembly, file = "my_subset.fa")
#' }


mmexport <- function(data, assembly, file){
  writeXStringSet(assembly[names(assembly) %in% as.character(data$scaffolds$scaffold)], file = file)
}
