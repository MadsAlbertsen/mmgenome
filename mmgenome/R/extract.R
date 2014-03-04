#' Extracts a subset of scaffolds
#'
#' A nice long description
#'
#' @usage extract(data, selection)
#'
#' @param data The dataframe containing all data.
#' @param selection The subspace to extract.
#' 
#' @return The subset of scaffolds in the original dataframe within the defined selection.
#' 
#' @export
#' @import sp
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

extract <-
function(data, selection){
  xname <- names(selection[1])
  yname <- names(selection[2])
  in.selection <- point.in.polygon(unlist(data[xname]), unlist(data[yname]), unlist(selection[1]), unlist(selection[2]), mode.checked=T)
  in.selection <- in.selection > 0
  return(data[in.selection, ])
}