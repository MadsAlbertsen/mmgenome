\name{mmload}
\alias{mmload}
\title{Data loading function}
\usage{
mmload(assembly, coverage, essential, tax)
}
\arguments{
  \item{assembly}{(required) The assembly.}

  \item{coverage}{(required) A list of coverage datasets.}

  \item{essential}{(required) A dataframe with information
  on essential genes.}

  \item{tax}{(required) A dataframe with taxonomic
  classification.}

  \item{coverge.type}{The type of coverage files either clc
  or simple (default: simple).}

  \item{tax.freq}{Remove phyla seen less than X times
  (default: 20).}

  \item{tax.expand}{Use class level assignments for
  specific phyla.}

  \item{tax.clean}{Clean up some of the phyla names
  (default: T)}

  \item{pca}{Add pca analyis of tetranucleotide frequencies
  (default: T)}

  \item{pca.ncomp}{Number of principal components to return
  (default: 3).}

  \item{other}{A list of additional datasets.}
}
\value{
A list with 2 dataframes: scaffolds and essential
}
\description{
A nice long description.
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

