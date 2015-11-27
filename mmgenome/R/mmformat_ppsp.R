#' Format ppsp data for use in mmgenome
#'
#' Format PhyloPythiaS+ data for use in mmgenome
#'
#' @usage mmformat_ppsp(data)
#'
#' @param data (required) The PP.pOUT file from PPSP.
#' @param ranks A vector of which phylogenetic ranks to keep (default: c("phylum", "class))
#' @param ranknames A vector with names of the kept ranks.
#' 
#' @return A dataframe with cleaned PPSP data.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmformat_ppsp <- function(data, ranks = c("phylum", "class"), ranknames = NULL){
  
  colnames(data)[1] <- c("scaffold")
  data[data == ""] <- NA
  data$scaffold <- as.character(data$scaffold)
  data <- data[,c("scaffold",ranks)]
  if(!is.null(ranknames)){colnames(data) <- c("scaffold",ranknames)}
  
  return(data)
}