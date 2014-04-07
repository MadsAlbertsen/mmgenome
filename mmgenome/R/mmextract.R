#' Extracts a subset of scaffolds
#'
#' A nice long description
#'
#' @usage mmextract(data, seletion)
#'
#' @param data The dataframe containing all data.
#' @param selection The subspace to extract.
#' 
#' @return The subset of scaffolds in the original dataframe within the defined selection.
#' 
#' @export
#' @import alphahull
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmextract <- function(data, selection){
  xname <- names(selection[1])
  yname <- names(selection[2])
  xr <- range(selection[xname])
  yr <- range(selection[yname])
  sel.hull <- ahull(selection, alpha=100000)
  ds <- subset(data$scaffolds, data$scaffolds[xname] > min(xr) & data$scaffolds[xname] < max(xr) & data$scaffolds[yname] > min(yr) & data$scaffolds[yname] < max(yr))
  in.selection <- apply(ds[c(xname, yname)],1,function(x){inahull(ahull.obj=sel.hull, p=x)})
  
  out <- ds[in.selection, ]
  es <- subset(data$essential, data$essential$scaffold %in% out$scaffold)
  
  outlist <- list(scaffolds = out, essential = es)  
  return(outlist)  
  
}