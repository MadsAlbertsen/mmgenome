#' Calculates statistics on a set of scaffolds
#'
#' A nice long description
#'
#' @usage sstats(scaffolds, essential, ncov)
#'
#' @param scaffolds The dataframe containing all data on scaffolds.
#' @param essential The dataframe containing all data on essential genes.
#' @param ncov The number of coverage datasets (default = 2).
#' 
#' @return An overview of key statistics of the scaffolds.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

sstats <- function(scaffolds, essential, ncov=2) {
  coverage <- list()
  for (i in 1:ncov){
    coverage[i] <- round(sum((scaffolds[,i+3]*scaffolds$length))/sum(scaffolds$length),1)
    names(coverage)[i] <- paste("Coverage", i, sep="")
  }
  coverage<-t(data.frame(coverage))
  total.length <- sum(scaffolds$length)
  mean.length <- round(mean(scaffolds$length),1)
  n.scaffolds <- nrow(scaffolds)
  max.length <- max(scaffolds$length)
  mean.gc <- round(sum((scaffolds$gc*scaffolds$length))/sum(scaffolds$length),1)
  total.ess <- nrow(essential)
  unique.ess <- length(unique(essential$hmm.id))
  
  out<-rbind(total.length,n.scaffolds, mean.length, max.length, mean.gc, coverage, total.ess, unique.ess)
  colnames(out) <- "General Stats"
  return(out)
}

