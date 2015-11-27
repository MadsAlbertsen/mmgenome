#' Extracts a subset of scaffolds
#'
#' Extracts a subset of scaffolds defined on a \code{mmplot} and selected using \code{mmplot_locator}.
#'
#' @usage mmextract(data, selection)
#'
#' @param data (required) The dataframe containing all data.
#' @param selection (required) Extract the scaffolds within the subspace.
#' @param minlength Minimum scaffold length.
#' @param exclude Vector of scaffold names to exclude.
#' @param include Vector of scaffold names to include.
#' @param original The original dataset, required when using include.
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

mmextract <-  function(data, selection = NULL, minlength = NULL, exclude = NULL, include = NULL, original = NULL){
### To extract only using scaffold names
  if (is.null(selection) & !is.null(include)){
    out <- subset(original$scaffold, scaffold %in% c(include, data$scaffolds$scaffold))
    if (length(data) == 2){es <- subset(data$essential, data$essential$scaffold %in% out$scaffold)}
    if (length(data) == 2){
      outlist <- list(scaffolds = out, essential = es)
    } else {
      outlist <- list(scaffolds = out)
    }  
### To extract using selection
  } else{
    if (!is.null(minlength)){
      data$scaffolds <- subset(data$scaffolds, length >= minlength)
      if (length(data) == 2){data$essential <- subset(data$essential, data$essential$scaffold %in% data$scaffolds$scaffold)}
    } 
    xname <- names(selection[1])
    yname <- names(selection[2])
    in.selection <- point.in.polygon(unlist(data$scaffolds[xname]), unlist(data$scaffolds[yname]), unlist(selection[1]), unlist(selection[2]), mode.checked=T)
    in.selection <- in.selection > 0
  
    out <- data$scaffolds[in.selection, ]  
  
    if (!is.null(exclude)){
      out <- subset(out, !(scaffold %in% exclude))
    }
    if (length(data) == 2){es <- subset(data$essential, data$essential$scaffold %in% out$scaffold)}
  
    if (!is.null(include)){
      include_unique <- include[!(include %in% out$scaffold)]
      out_in <- subset(original$scaffold, scaffold %in% include_unique)
      out <- rbind.data.frame(out, out_in)
      if (length(data) == 2){es <- subset(original$essential, original$essential$scaffold %in% out$scaffold)}
    }
  
    if (length(data) == 2){
      outlist <- list(scaffolds = out, essential = es)
    } else {
      outlist <- list(scaffolds = out)
    }  
  }
  return(outlist)
}