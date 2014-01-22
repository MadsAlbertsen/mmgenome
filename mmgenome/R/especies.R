#' Extracts overview of essential genes at Phylum level
#'
#' A nice long description
#'
#' @usage especies(eg)
#'
#' @param eg Pre-formatted file of genomes and essential genes: data(eg)
#' 
#' @return A list with two objects. $plot a ggplot object, $table a table with summery data.
#' 
#' @export
#' @import ggplot2
#' @import reshape
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' @examples
#' 
#' data(eg)
#' eg_alpha <- subset(eg, Class == "Alphaproteobacteria")
#' species_data <- especies(eg_alpha)
#' species_data$plot
#' speces_data$table[1,]

especies <- function(eg){
  em <- melt(eg,id.vars=c(colnames(eg)[1:9]))
  colnames(em)[10:11] <- c("HMM", "Count")
  cHMM <- aggregate(Count ~ HMM, sum, data = em)
  cHMM <- cHMM[ order(cHMM$Count, decreasing=T) ,]
  em$HMM <- factor(em$HMM, levels = cHMM$HMM)
  cAcc <- aggregate(Count ~ Accession, sum, data = em)
  cAcc <- cAcc[ order(cAcc$Count, decreasing=T) ,]
  em$Accession <- factor(em$Accession, levels = cAcc$Accession)
  
  p <- ggplot(em, aes(x= Accession, y = HMM, fill = Count)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3))
  
  tHMM <- aggregate(Count ~ Accession, sum, data = em)
  uHMM <- aggregate(HMM ~ Accession, length, data = subset(em, Count >= 1))
  mHMM <- merge(tHMM, uHMM)
  colnames(mHMM) <- c("Accession","Total.HMM", "Unique.HMM")
  summary <- merge(eg, mHMM, by ="Accession")
  out <- list(p, summary, em) 
  names(out) <- c("plot", "table", "data")
  return(out)
}