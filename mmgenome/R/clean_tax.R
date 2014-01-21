#' Cleans the taxonomy output
#'
#' A nice long description
#'
#' @usage clean_tax(tax, occurence, expand, clean)
#'
#' @param tax The taxonomy dataframe (required).
#' @param occurence Remove phyla with less than X occurences.
#' @param expand Given a specific phylum, use the class level assignments instead.
#' @param clean Simplify some of the phylum names.
#' 
#' @return A cleaned taxonomy file.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

clean_tax <- function (tax, occurence=NULL, expand=NULL, clean=FALSE){
  tax$phylum <- as.character(tax$phylum)
  tax$class <- as.character(tax$class)
  
  # Given a specific phyla - uses class assignment instead of phyla (often usefull with Proteobacteria)
  if (!is.null(expand)){
    for (i in 1:nrow(tax)){
      if(tax$phylum[i] == expand){
        tax$phylum[i] <- tax$class[i]   
      }
    }        
  }

  # Make some of the phyla names more pretty
  if (clean == TRUE){
    tax$phylum<-gsub(" <phylum>", "", tax$phylum)
    tax$phylum<-gsub("unclassified Bacteria", "Unclassified Bacteria", tax$phylum)
    tax$phylum<-gsub("Fibrobacteres/Acidobacteria group", "Acidobacteria", tax$phylum)
    tax$phylum<-gsub("Bacteroidetes/Chlorobi group", "Bacteroidetes", tax$phylum)
    tax$phylum<-gsub("delta/epsilon subdivisions", "Deltaproteobacteria", tax$phylum)
    tax$phylum<-gsub("Chlamydiae/Verrucomicrobia group", "Verrucomicrobia", tax$phylum)  
  }
  
  tax<-subset(tax, phylum != "NA")
  tax$phylum <- as.factor(tax$phylum)
  tax$class <- as.factor(tax$class)  

  # Remove phyla seen less than X times
  if (!is.null(occurence)){
    uniquetaxa<-as.character(unique(tax$phylum)) 
    for (i in 1:length(uniquetaxa)) {
      tax.match <- which(tax$phylum==uniquetaxa[i]) 
      no.occurences <- sum(tax$count[tax.match])
      if (no.occurences < occurence) { 
        tax$phylum[tax.match] <- NA
      }
    }     
    tax<-subset(tax, phylum != "NA") 
  }  

  tax<-droplevels(tax)
  return(tax)
}