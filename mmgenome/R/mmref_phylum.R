#' Extracts overview of essential genes at Phylum level
#'
#' A nice long description
#'
#' @usage mmref_phylum(eg)
#'
#' @param eg Pre-formatted file of genomes and essential genes: data(eg)
#' 
#' @return A list with three objects. $plot a ggplot object, $table a table with summery data and $data with the data used for plotting.
#' 
#' @export
#' @import ggplot2
#' @import reshape2
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' @examples
#' 
#' data(eg)
#' phylum_data <- mmref_phylum(eg)
#' phylum_data$plot
#' phylum_data$table

mmref_phylum <- function(eg){
  em <- melt(eg,id.vars=c(colnames(eg)[1:9]))
  colnames(em)[10:11] <- c("HMM", "Count")
  cHMM <- aggregate(Count ~ HMM, sum, data = em)
  cHMM <- cHMM[ order(cHMM$Count, decreasing=T) ,]
  em$HMM <- factor(em$HMM, levels = cHMM$HMM)
  cAcc <- aggregate(Count ~ Accession, sum, data = em)
  cAcc <- cAcc[ order(cAcc$Count, decreasing=T) ,]
  em$Accession <- factor(em$Accession, levels = cAcc$Accession)
  emp <- aggregate(Count ~ Kingdom + Phylum + HMM, median, data = em)
  cPhylum <- aggregate(Accession ~ Phylum, length, data = eg)
  cPHMM <- aggregate(Count ~ Phylum, sum, data = emp)
  cPHMM$Count <- round(cPHMM$Count,1)
  uPHMM <- aggregate(HMM ~ Phylum, length, data = subset(emp, Count > .8))
  tPHMM <- merge(uPHMM, cPHMM)
  tPHMM <- merge(tPHMM, cPhylum)
  colnames(tPHMM) <- c("Phylum", "Unique.HMM", "Total.HMM", "n.Genomes")
  tPHMM <- tPHMM[ order(tPHMM$Unique, decreasing=T) ,]
  emp$Phylum <- factor(emp$Phylum, levels = tPHMM$Phylum)
  phylum.label <- paste(tPHMM$Phylum, " (", tPHMM$n.Genomes , ",", tPHMM$Unique.HMM, ", ", round(tPHMM$Total.HMM,0), ")", sep = "")
  
 p <-  ggplot(emp, aes(x= Phylum, y = HMM, fill = Count)) + 
    geom_tile() +
    scale_x_discrete(labels=phylum.label) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) +
    coord_flip()  
 
 out <- list(p, tPHMM, emp) 
 names(out) <- c("plot", "table", "data")
 return(out)
}