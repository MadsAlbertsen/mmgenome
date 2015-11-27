#' Add data to any mmgenome object
#'
#' Add new dataframes to an excisting mmgenome object.
#'
#' @usage mmadd(data, new)
#'
#' @param data (required) The existing mmgenome object.
#' @param newdata (required) The new dataframe that is to be added.
#' 
#' @return A list with 2 dataframes: scaffolds and essential (if provided).
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmadd <- function(data, newdata){
  colnames(newdata)[1] <- "scaffold"
  newdata$scaffold <- as.character(newdata$scaffold)
  
  data$scaffolds <- merge(data$scaffolds, newdata, by = "scaffold", all.y = F, all.x = T)
  
  return(data)
}