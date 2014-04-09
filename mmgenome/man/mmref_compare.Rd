\name{mmref_compare}
\alias{mmref_compare}
\title{Exports scaffolds in fasta format}
\usage{
mmref_compare(data, compare, level, output, display.all)
}
\arguments{
  \item{data}{The dataframe containing all data.}

  \item{compare}{Which phyla to compare to.}

  \item{level}{What level to compare on. Currently only
  supports "phyla".}

  \item{output}{Either a table or plot (default: plot)}

  \item{display.all}{Show all HMMs even if they are not
  contained in the selected references (default: F).}
}
\value{
A table or figure comparing the essential genes in the bin
to selected reference genomes.
}
\description{
A nice long description
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

