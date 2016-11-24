#' Calculates summary statistics on a set of scaffolds
#'
#' Calculates summary statistics on a set of scaffolds. Can be used on any object loaded using \code{mmload}.
#'
#' @usage mmstats(data, ncov)
#'
#' @param data The dataframe containing the bin data (e.g. dA).
#' @param ncov The number of coverage datasets (default = 2).
#' @param data_compare The dataframe containing the all data (e.g. d).
#' 
#' @return An overview of key statistics of the scaffolds.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' data(rocco)
#' mmstats(data = d, ncov = 4)
#' }

mmstats <- function(data, ncov=2,data_compare=NULL) {
  options(scipen=8)
  data$scaffolds$length <- as.numeric(data$scaffolds$length)
  coverage <- list()
  for (i in 1:ncov){
    coverage[i] <- round(sum((data$scaffolds[,i+3]*data$scaffolds$length))/sum(data$scaffolds$length),2)
    names(coverage)[i] <- paste("Coverage",names(data$scaffolds)[i+3])
  }
  coverage<-t(data.frame(coverage))
  Length.total <- sum(data$scaffolds$length)
  Length.mean <- round(mean(data$scaffolds$length),1)
  n.scaffolds <- nrow(data$scaffolds)
  Length.max <- max(data$scaffolds$length)
  GC.mean <- round(sum((data$scaffolds$gc*data$scaffolds$length))/sum(data$scaffolds$length),1)
  if (length(data) == 2){Ess.total <- nrow(data$essential)}
  if (length(data) == 2){Ess.unique <- length(unique(data$essential$hmm.id))}
  
  sorted_length <- sort(data$scaffolds$length)
  cum_sum <- 0
  for (i in 1:length(sorted_length)){
    cum_sum <- cum_sum + sorted_length[i]
    if (cum_sum >= Length.total/2){
      N50 <- sorted_length[i]
      break
    }
  }
  # If there is  a data object to compare it to:
  if (!is.null(data_compare)) {
    data_compare$scaffolds$length <- as.numeric(data_compare$scaffolds$length)
    coverage_compare <- list()
    for (i in 1:ncov){
      coverage_compare[i] <- round(sum((data_compare$scaffolds[,i+3]*data_compare$scaffolds$length))/sum(data_compare$scaffolds$length),2)
      names(coverage_compare)[i] <- paste("Coverage",names(data$scaffolds)[i+3])
    }
    coverage_compare<-t(data.frame(coverage_compare))
    Length.total_compare <- sum(data_compare$scaffolds$length)
    fractions<-round(Length.total/Length.total_compare*(coverage/coverage_compare)*100,2)
    row.names(fractions)<-gsub(pattern = "Coverage",replacement = "Fraction(%)",x = row.names(fractions))
    coverage<-fractions
  }
  
  
  if (length(data) == 2){
    out<-rbind(n.scaffolds, GC.mean, N50, Length.total, Length.max, Length.mean, coverage, Ess.total, Ess.unique)
  } else{
    out<-rbind(n.scaffolds, GC.mean, N50, Length.total, Length.max, Length.mean, coverage)
  }
  colnames(out) <- "General Stats"
  return(out)
}
