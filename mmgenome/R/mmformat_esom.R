#' Format raw esom data for use in mmgenome
#'
#' Format raw esom data generated using databionics ESOM software.
#'
#' @usage mmformat_esom(names, bm, umx)
#'
#' @param names (required) The .names esom file.
#' @param bm (required) The .bm esom file.
#' @param umx (required) The .umx esom file.
#' @param xname Name of the esom x-coordinate (default: "esomx")
#' @param yname Name of the esom y-coordinate (default: "esomy")
#' 
#' @return A list with two dataframes. The esom scaffold coordinates and the esom map boundaries.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmformat_esom <- function(names, bm, umx, xname = "esomx", yname = "esomy"){
  
  colnames(names) <- c("ename", "splitname", "scaffold")
  colnames(bm) <- c("ename", "esomx", "esomy")
  
  esom_coord <- merge(names, bm) %>% select(scaffold, esomx, esomy)
  esom_coord$esomx <- as.numeric(esom_coord$esomx)
  esom_coord$esomy <- as.numeric(esom_coord$esomy)
  colnames(esom_coord) <- c("scaffold", yname, xname)
  
  umx$R <- 1:nrow(umx)
  
  esom_map <- melt(umx, id.vars = "R")
  esom_map$variable <- sub(pattern = "V", replacement = "", x = esom_map$variable)
  colnames(esom_map) <- c("esomy", "esomx", "density")
  esom_map$esomx <- as.numeric(esom_map$esomx)
  
  
  return(list(esom_coord = esom_coord, 
              esom_map = esom_map))
}