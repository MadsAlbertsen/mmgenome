#' Imports and formats silva rRNA results
#'
#' A nice long description
#'
#' @usage import_rrna(data, type)
#'
#' @param data The dataframe containing all data.
#' @param type Either 16S or 23S. Just adds the header.
#' 
#' @return Nicely formatted dataframe.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


import_rrna <- function(data, type = "16S") {
  temp <- sapply(data$sequence_identifier, function (x) strsplit(as.character(x), "\\.", perl = T)[[1]][1])
  temp <-cbind(temp,as.character(data$lca_tax_slv))  
  colnames(temp) <- c("scaffold",type)
  return(temp)
}
