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
  ds <- subset(data, data[xname] > min(xr) & data[xname] < max(xr) & data[yname] > min(yr) & data[yname] < max(yr))
  in.selection <- apply(ds[c(xname, yname)],1,function(x){inahull(ahull.obj=sel.hull, p=x)})
  return(ds[in.selection, ])
}