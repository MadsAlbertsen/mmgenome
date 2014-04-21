#' Extracts a subset of scaffolds
#'
#' A nice long description
#'
#' @usage mmextract(data, selection)
#'
#' @param data The dataframe containing all data.
#' @param selection The subspace to extract.
#' @param minlength Minimum scaffold length.
#' 
#' @return The subset of scaffolds in the original dataframe within the defined selection.
#' 
#' @export
#' @import sp
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmextract <-  function(data, selection, minlength = NULL){
  if (!is.null(minlength)){
    data$scaffolds <- subset(data$scaffolds, length >= minlength)
    data$essential <- subset(data$essential, data$essential$scaffold %in% data$scaffolds$scaffold)    
  } 
  xname <- names(selection[1])
  yname <- names(selection[2])
  in.selection <- point.in.polygon(unlist(data$scaffolds[xname]), unlist(data$scaffolds[yname]), unlist(selection[1]), unlist(selection[2]), mode.checked=T)
  in.selection <- in.selection > 0
  
  out <- data$scaffolds[in.selection, ]
  es <- subset(data$essential, data$essential$scaffold %in% out$scaffold)
  
  outlist <- list(scaffolds = out, essential = es)  
  
  return(outlist)
}