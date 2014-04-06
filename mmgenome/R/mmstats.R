#' Calculates statistics on a set of scaffolds
#'
#' A nice long description
#'
#' @usage mmstats(scaffolds, essential, ncov)
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

mmstats <- function(scaffolds, essential, ncov=2) {
  coverage <- list()
  for (i in 1:ncov){
    coverage[i] <- round(sum((scaffolds[,i+3]*scaffolds$length))/sum(scaffolds$length),1)
    names(coverage)[i] <- paste("Coverage",names(scaffolds)[i+3])
  }
  coverage<-t(data.frame(coverage))
  Length.total <- sum(scaffolds$length)
  Length.mean <- round(mean(scaffolds$length),1)
  n.scaffolds <- nrow(scaffolds)
  Length.max <- max(scaffolds$length)
  GC.mean <- round(sum((scaffolds$gc*scaffolds$length))/sum(scaffolds$length),1)
  Ess.total <- nrow(essential)
  Ess.unique <- length(unique(essential$hmm.id))
  
  sorted_length <- sort(scaffolds$length)
  cum_sum <- 0
  for (i in 1:length(sorted_length)){
    cum_sum <- cum_sum + sorted_length[i]
    if (cum_sum >= Length.total/2){
      N50 <- sorted_length[i]
      break
    }
  }
  
  out<-rbind(n.scaffolds, GC.mean, N50, Length.total, Length.max, Length.mean, coverage, Ess.total, Ess.unique)
  colnames(out) <- "General Stats"
  return(out)
}