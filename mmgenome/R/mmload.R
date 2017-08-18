#' Data loading function
#'
#' This function is used to load information on scaffolds into a single object that is used by all subsequent mmgenome functions.
#'
#' @usage mmload(assembly, coverage, essential)
#'
#' @param assembly (required) The assembly.
#' @param coverage (required) A vector with the names of the coverage datasets.
#' @param essential A dataframe with information on essential genes. 
#' @param tax A dataframe with taxonomic classification.
#' @param tax.freq Remove phyla seen less than X times (default: 20).
#' @param tax.expand Use class level assignments for specific phyla.
#' @param tax.clean Clean up some of the phyla names (default: T)
#' @param pca Add PCA analyis of tetranucleotide frequencies (default: T)
#' @param pca.ncomp Number of principal components to return (default: 3).
#' @param other A vector with the names of additional datasets.
#' 
#' @return A list with 2 dataframes: scaffolds and essential (if provided).
#' 
#' @export
#' @import Biostrings
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' d <- mmload(assembly = assembly, 
#'             coverage = c("A", "B", "C", "D"), 
#'             essential = ess,           
#'             tax = tax,
#'             tax.expand = "Proteobacteria",
#'             tax.freq = 50,
#'             other = "rRNA16S"
#'             )
#' }

mmload <- function(assembly, coverage, essential = NULL, pca = T, pca.ncomp = 3, tax = NULL, tax.freq = 20, tax.clean = T, tax.expand = NULL, other = NULL, BH_SNE = F, threads = 2){
  
  ##### Load coverage, gc and length #####
  
  print("Loading scaffold length, coverage and gc.")  
  out <- cbind.data.frame(names(assembly), width(assembly), letterFrequency(assembly, letters = c("CG"), as.prob=T)*100)  
  colnames(out) <- c("scaffold","length","gc")
  out$scaffold <- as.character(out$scaffold)
  
  for (i in 1:length(coverage)){
    data <- get(coverage[i])
    if (ncol(data) != 2){
      stop(paste("The coverage file", coverage[i], "does not have 2 columns, the first column needs to be scaffold name and the second scaffold coverage."))
    }
    if (nrow(data) > nrow(out)) {
      stop(paste("The coverage file", coverage[i], "contains more scaffolds than the assembly."))
    }
    if (nrow(data) < nrow(out)) {
      warning(paste("The coverage file", coverage[i], "contains less scaffolds than the assembly. Setting the coverage of the missing scaffolds to 0."))
    }
    colnames(data) <- c("scaffold", "coverage")
    data$scaffold <- as.character(data$scaffold)
    cout <- nrow(out)
    out <- merge(x = out, y = data, by = "scaffold", all = T)
    if(nrow(out) > cout){
      stop(paste("The coverage file", coverage[i], "contains scaffolds with names not found in the assembly. Make sure that the names of the scaffolds in the assembly and coverage files are identical."))
    }
    colnames(out)[length(out)] <- coverage[i]
  }
  out[is.na(out)] <- 0
  
  ##### Calculate PCA on tetranucleotide frequencies #####
  
  if (pca == T){
    print("Calculating PCA.")
    kmer_for<-oligonucleotideFrequency(assembly, 4, as.prob=T)
    kmer_revC <- oligonucleotideFrequency(reverseComplement(assembly), 4, as.prob=T)
    kmer <- (kmer_for + kmer_revC)/2*100
    rda <- rda(kmer[,2:ncol(kmer)])
    res <- cbind.data.frame(scaffold = as.character(names(assembly)), scores(rda,choices=1:pca.ncomp)$sites)
    out <- merge(x = out, y = res, by = "scaffold")
  }
  
  if (BH_SNE == T){
    set.seed(42) # Sets seed for reproducibility
    if (pca == F) {
      print("Be patient calculating kmer sitributions")
      kmer_for<-oligonucleotideFrequency(assembly, 4, as.prob=T)
      kmer_revC <- oligonucleotideFrequency(reverseComplement(assembly), 4, as.prob=T)
      kmer <- (kmer_for + kmer_revC)/2*100
    }
    print("Be patient calculating BH-SNE")
    tsne_out <- Rtsne.multicore(kmer,num_threads = threads,check_duplicates = F)
    tsne_plot <- data.frame(x_tsne = tsne_out$Y[,1], y_tsne = tsne_out$Y[,2])
    res <- cbind.data.frame(scaffold = as.character(names(assembly)), tsne_plot)
    out <- merge(x = out, y = res, by = "scaffold")
  }
  
  ##### Load and clean the taxonomic classification #####
  
  if (!is.null(tax)){
    print("Loading taxonomy.")
    tax$phylum <- as.character(tax$phylum)
    tax$class <- as.character(tax$class)
    if (!is.null(tax.expand)){
      for (i in 1:nrow(tax)){
        if(tax$phylum[i] %in% tax.expand){
          tax$phylum[i] <- tax$class[i]   
        }
      }        
    }
    
    # Make some of the NCBI phyla names more pretty
    if (tax.clean == TRUE){
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
    if (!is.null(tax.freq)){
      uniquetaxa<-as.character(unique(tax$phylum)) 
      for (i in 1:length(uniquetaxa)) {
        tax.match <- which(tax$phylum==uniquetaxa[i]) 
        no.occurences <- sum(tax$count[tax.match])
        if (no.occurences < tax.freq) { 
          tax$phylum[tax.match] <- NA
        }
      }     
      tax<-subset(tax, phylum != "NA") 
    }  
    tax<-droplevels(tax)
    tax$scaffold <- as.character(tax$scaffold)
    colnames(tax)[2] <- "essential"
    out <- merge(x = out, y = tax[,c(1,2)], by = "scaffold", all = T)
  }
  
  ##### Load other datasets #####
  
  if (!is.null(other)){
    print("Loading additional datasets.")
    length_check <- nrow(out)
    for (i in 1:length(other)){
      data <- get(other[i])
      colnames(data)[1] <- c("scaffold")
      if (ncol(data) < 2){
        stop(paste("The file", other[i], "does not have 2 columns. At least a column with scaffold names and a coloumn with additional information is required."))
      }
      if (nrow(data) > nrow(out)) {
        stop(paste("The file", other[i], "contains more scaffolds than the assembly."))
      }
      if (nrow(data) < nrow(out)) {
        warning(paste("The file", other[i], "contains less scaffolds than the assembly. Missing values are treated as NA."))
      }      
      data$scaffold <- as.character(data$scaffold)
      out <- merge(x = out, y = data, by = "scaffold", all.x = T)
      #colnames(out)[length(out)] <- other[i]
      if(nrow(out) > length_check){
        stop(paste("The file", other[i], "contains scaffolds with names not found in the assembly. Make sure that the names of the scaffolds in the assembly and other files are identical."))
      }
    }
  } 
  
  ##### Output the data
  if(!is.null(essential)){
    outlist <- list(scaffolds = out, essential = ess)  
  } else {
    outlist <- list(scaffolds = out)  
  }
  
  
  return(outlist)
}
