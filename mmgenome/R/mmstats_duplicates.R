#' Identifies dublicates essential genes
#'
#' A nice long description
#'
#' @usage mmstats_duplicates(essential)
#'
#' @param essential The dataframe with information on essential genes.
#' 
#' @return A overview of duplicated essential genes and their parrent scaffold.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmstats_duplicates <- function (essential) {
  dub<-essential[which(duplicated(essential$hmm.id) | duplicated(essential$hmm.id, fromLast=TRUE)),] 
  dub <- dub[order(dub$hmm.id),c(1,3)]
  return(dub)
}