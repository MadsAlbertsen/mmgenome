#' Extracts a subset of scaffolds
#'
#' Extracts a subset of scaffolds defined on a \code{mmplot} and selected using \code{mmplotlocator}.
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
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' 
#' p <- mmplot(data = d, x = "C13.12.03", y = "C14.01.09", log.x = T, log.y = T, color = "phylum", minlength = 3000)
#' p
#' sel <- mmplot_locator(p)
#' sel <- data.frame(C13.12.03  =  c(1.39, 2.07, 16.8, 19.4, 7.72, 1.76),
#'                   C14.01.09  =  c(29.4, 67.6, 85.9, 43.6, 16.7, 14.9))
#' mmplot_selection(p, sel)
#' 
#' dA <- mmextract(d, sel)
#' }

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