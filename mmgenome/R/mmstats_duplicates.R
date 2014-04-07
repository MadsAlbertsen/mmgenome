#' Identifies duplicates essential genes
#'
#' A nice long description
#'
#' @usage mmstats_duplicates(data)
#'
#' @param data The dataframe.
#' 
#' @return A overview of duplicated essential genes and their parrent scaffold.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmstats_duplicates <- function (data) {
  dup<-essential[which(duplicated(data$essential$hmm.id) | duplicated(data$essential$hmm.id, fromLast=TRUE)),] 
  dup <- dup[order(dup$hmm.id),c(1,3)]
  return(dup)
}