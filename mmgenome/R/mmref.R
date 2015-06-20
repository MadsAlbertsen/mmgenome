#' Compare essential genes in a bin to a group of references
#'
#' Compare essential genes in an extracted genome bin to a group of references
#'
#' @usage mmref(data, tax.level, tax.compare)
#'
#' @param data The dataframe containing the genome bin of interest (optional).
#' @param tax.level What level to compare on. Supports Kingdom, Phylum, Class, Genus, Species (default: Phylum)
#' @param tax.compare A vector with names on what should be used for comparison (e.g. c("Proteobacteria", "Actinobacteria")) or "all".
#' @param tax.aggregate Summarise the data to a fixed taxonomic rank. Supports Kingdom, Phylum, Class, Genus, Species (default: tax.level).
#' @param summarise Report "mean" or "median" HMM counts (default: mean).  
#' @param output Either complete or plot (default: plot)
#' @param plot.display Show all HMMs even if they are not contained in the selected references (default: T).
#' @param plot.sort Sort the plot (default: T)
#' 
#' @return An overview of essential genes in the bin and selected reference genomes.
#' 
#' @export
#' @import Biostrings
#' @import reshape2
#' @import ggplot2
#' @import plyr
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' mmref(tax.level = "Phylum", tax.compare = c("Proteobacteria", "Actinobacteria")) 
#' 
#' res <- mmref(tax.level = "Phylum", tax.compare = c("Proteobacteria", "Actinobacteria"), output = "complete") 
#' res$table
#' 
#' mmref(tax.level = "Phylum", tax.compare = "all") 
#' mmref(tax.level = "Phylum", tax.compare = "Proteobacteria", tax.aggregate = "Class")
#' mmref(tax.level = "Phylum", tax.compare = "Proteobacteria", tax.aggregate = "Class", summarise = "median") 
#' mmref(tax.level = "Class", tax.compare = "Betaproteobacteria", tax.aggregate = "Strain", summarise = "median")
#' 
#' mmref(tax.level = "Phylum", tax.compare = "Proteobacteria", tax.aggregate = "Class", summarise = "median", data = dA)}

mmref <- function(data = NULL, tax.level="Phylum", tax.compare, tax.aggregate = tax.level, summarise = "mean", output = "plot", plot.display = T, plot.sort = T){
  
  ### Load the reference data
  data(eg)
  
  ### Subset to the relevant dataset
  suppressWarnings(if(tax.compare != "all"){ es <- subset(eg, eg[,tax.level] %in% tax.compare) } else {es <- eg})
  
  ### Convert the dataframe to long format
  em <- melt(es,id.vars=c(colnames(eg)[1:9]), value.name = "Count", variable.name = "HMM")
  
  ### Summarise the HMM counts as either mean or medians
  hcount <- ddply(em, c(tax.aggregate,"HMM"), summarize, median = median(Count), mean = mean(Count))  
  
  ### Load essential gene informatio from the extracted genome bin
  if (!is.null(data)){
    if (length(data) != 2){
      stop(paste("There is no data on essential genes in the supplied data!"))
    }
    data$essential$hmm.id <- as.factor(gsub("\\..*","", as.character(data$essential$hmm.id)))
    mmbin <- cbind.data.frame(HMM = data$essential$hmm.id, name = "mmgenomebin")
    colnames(mmbin)[ncol(mmbin)] <- tax.aggregate
    mmcount <- ddply(mmbin, c(tax.aggregate,"HMM"), summarize, median = length(HMM), mean = length(HMM))      
    mmcount_missing <- data.frame(name = "mmgenomebin", HMM = levels(hcount$HMM)[!(levels(hcount$HMM) %in% levels(mmcount$HMM))],median = 0, mean = 0)
    colnames(mmcount_missing)[1] <- tax.aggregate
    hcount <- rbind.data.frame(hcount, mmcount, mmcount_missing)
  }
  
  ### Remove any HMMs not seen in the datasets
  if (plot.display == F){
    if (summarise == "mean"){hcount <- subset(hcount, mean != 0)}
    if (summarise == "median"){hcount <- subset(hcount, median != 0)}
  }
  
  ### Sort the data for prettyness
  if (plot.sort == T){
    cHMM <- ddply(hcount, ~HMM, summarize, mean = mean(mean))  
    cHMM <- cHMM[order(cHMM$mean),]
    hcount$HMM <- factor(hcount$HMM, levels = rev(cHMM$HMM))
    
    cTax <- ddply(hcount, tax.aggregate, summarize, Total = sum(mean))  
    cTax <- cTax[order(cTax$Total),]
    hcount[,tax.aggregate] <- factor(hcount[,tax.aggregate], levels = rev(cTax[,tax.aggregate]))
  }
  
  ### Plot the data    
  p <- ggplot(data = hcount, aes_string(x = tax.aggregate, y = "HMM", fill = summarise)) +
    geom_tile() +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3))
  
  ### Count the number of genomes in each group
  gcount <- ddply(es, tax.aggregate, summarize, Genomes = length(Accession))
  ems <- subset(em, Count != 0)
  ems1 <- ddply(ems, c(tax.aggregate,"Accession"), summarize, Total = sum(Count), Unique = length(HMM))
  ems_table <- ddply(ems1, tax.aggregate, summarize, Genomes = length(Accession), Total.HMM = median(Total), Unique.HMM = median(Unique))
  
  if (!is.null(data)){
    mmbin_table <- data.frame(name = "mmgenomebin", Genomes = 1, Total.HMM = sum(mmcount$median), Unique.HMM = length(mmcount$HMM))
    colnames(mmbin_table)[1] <- tax.aggregate
    ems_table <- rbind.data.frame(ems_table, mmbin_table)
  }
  
  if(output == "complete") {return(list(plot = p, data = hcount, table = ems_table))}
  if(output == "plot") {return(p)}
}