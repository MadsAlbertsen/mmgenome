#' Exports scaffolds in fasta format
#'
#' A nice long description
#'
#' @usage mmref_compare(data, compare, level, output, display.all)
#'
#' @param data The dataframe containing all data.
#' @param compare Which phyla to compare to.
#' @param level What level to compare on. Currently only supports "phyla".
#' @param output Either a table or plot (default: plot)
#' @param display.all Show all HMMs even if they are not contained in the selected references (default: F). 
#' 
#' @return A table or figure comparing the essential genes in the bin to selected reference genomes.
#' 
#' @export
#' @import Biostrings
#' @import reshape2
#' @import ggplot2
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmref_compare <- function(data, compare, level="phylum", output = "plot", display.all = F){
  data(eg)
  if (level == "phylum"){
    phylum.data <- mmref_phylum(eg)$data
    t <- subset(phylum.data, Phylum %in% compare)[,2:4]
    if(display.all == F){ t <- subset(t, Count != 0) }
    t2 <- cast(t, HMM~Phylum, value = "Count")
    ess <- cbind(data$essential[,c(1,3)],rep(1,nrow(data$essential)))
    colnames(ess) <- c("scaffold", "HMM", "mmbin")
    cHMM <- aggregate(mmbin ~ HMM, sum, data = ess)
    cHMM$HMM <- sapply(cHMM$HMM, function (x) strsplit(as.character(x), "\\.", perl = T)[[1]][1])
    out <- merge(t2, cHMM, by = "HMM", all.x = T, all.y = T)
  }
  
  if(output != "table"){
    out[is.na(out)] <- 0 
    temp <- melt(out, id.var = "HMM")
    colnames(temp) <- c("HMM","Phylum", "Count")
    p <-  ggplot(temp, aes(x= Phylum, y = HMM, fill = Count)) + 
      geom_tile() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) #+
    out <- p
  }  
  return(out)
}



