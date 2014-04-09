#' Data loading function
#'
#' A nice long description.
#'
#' @usage mmload(assembly, coverage, essential, tax)
#'
#' @param assembly (required) The assembly.
#' @param coverage (required) A list of coverage datasets.
#' @param essential (required) A dataframe with information on essential genes. 
#' @param tax (required) A dataframe with taxonomic classification.
#' @param coverge.type The type of coverage files either clc or simple (default: simple).
#' @param tax.freq Remove phyla seen less than X times (default: 20).
#' @param tax.expand Use class level assignments for specific phyla.
#' @param tax.clean Clean up some of the phyla names (default: T)
#' @param pca Add pca analyis of tetranucleotide frequencies (default: T)
#' @param pca.ncomp Number of principal components to return (default: 3).
#' @param rRNA16S 16S taxonomic information.
#' @param rRNA23S 23S taxonomic information.
#' @param rRNA.type Currently only supports "silva" .csv files.
#' 
#' @return A list with 2 dataframes: scaffolds and essential
#' 
#' @export
#' @import Biostrings
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmload <- function(assembly, coverage, coverage.type = "simple", essential, pca = T, pca.ncomp = 3, tax = NULL, tax.freq = 20, tax.clean = T, tax.expand = NULL, rRNA16S = NULL, rRNA23S = NULL, rRNA.type = "silva"){

##### Load coverage, gc and length #####
  if (coverage.type == "clc"){
    for (i in 1:length(coverage)){  
      temp <- cbind.data.frame(coverage[[i]]$Reference.sequence,coverage[[i]]$Average.coverage)
      colnames(temp) <- c("scaffold", names(coverage)[i])
      coverage[[i]] <- temp
    }
  }
  
  b <- cbind.data.frame(names(assembly), width(assembly), letterFrequency(assembly, letters = c("CG"), as.prob=T)*100)
  colnames(b) <- c("scaffold","length","gc")
  l <- c(list(b), coverage)
  out <- Reduce(function(...) merge(..., all=T), l)
  out[is.na(out)] <- 0
  out$scaffold <- as.integer(as.character(out$scaffold))
  out <- out[with(out, order(scaffold)), ]

##### Calculate PCA on tetranucleotide frequencies #####

  if (pca == T){
    kmer_for<-oligonucleotideFrequency(assembly, 4, as.prob=T)
    kmer_revC <- oligonucleotideFrequency(reverseComplement(assembly), 4, as.prob=T)
    kmer <- (kmer_for + kmer_revC)/2*100
    rda <- rda(kmer[,2:ncol(kmer)])
    out<-cbind(out,scores(rda,choices=1:pca.ncomp)$sites)
  }

##### Load and clean the taxonomic classification #####

  if (!is.null(tax)){
      tax$phylum <- as.character(tax$phylum)
      tax$class <- as.character(tax$class)
      if (!is.null(tax.expand)){
        for (i in 1:nrow(tax)){
          if(tax$phylum[i] == tax.expand){
            tax$phylum[i] <- tax$class[i]   
          }
        }        
      }
      
      # Make some of the phyla names more pretty
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
    out <- merge(out,tax[,c(1,2)], by = "scaffold", all = T)
  }

##### Load rRNA taxonomic classifications #####

   if (!is.null(rRNA16S)){
     if (rRNA.type == "silva"){
       temp <- sapply(rRNA16S$sequence_identifier, function (x) strsplit(as.character(x), "\\.", perl = T)[[1]][1])
       temp <-cbind(temp,as.character(rRNA16S$lca_tax_slv))  
       colnames(temp) <- c("scaffold","16S")
       out <- merge(out, temp , by = "scaffold", all = T)
     }
   }
    
  if (!is.null(rRNA23S)){
    if (rRNA.type == "silva"){
      temp <- sapply(rRNA23S$sequence_identifier, function (x) strsplit(as.character(x), "\\.", perl = T)[[1]][1])
      temp <-cbind(temp,as.character(rRNA23S$lca_tax_slv))  
      colnames(temp) <- c("scaffold","23S")
      out <- merge(out, temp , by = "scaffold", all = T)
    }
  }     
  colnames(ess) = c("scaffold","orf","hmm.id")
  outlist <- list(scaffolds = out, essential = ess)
  return(outlist)
}