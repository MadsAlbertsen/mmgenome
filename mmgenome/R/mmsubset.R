#' Subset scaffolds
#'
#' Subset scaffolds using the default subset syntax.
#'
#' @usage mmsubset(data, ...)
#'
#' @param data (required) The dataframe containing all data.
#' @param ... (required) Further arguments are passed on to the subset function.
#' 
#' @return The subset of scaffolds.
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
#' subset <- mmsubset(data = d, gc > 60)
#' }

mmsubset <- function(data, ...){
  d_new <- subset(data$scaffolds, ...)
  e_new <- subset(data$essential, data$essential$scaffold %in% d_new$scaffold) 
  return(list(scaffolds = d_new, 
              essential = e_new))
}