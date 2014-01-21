\name{clean_tax}
\alias{clean_tax}
\title{Cleans the taxonomy output}
\usage{
clean_tax(tax, occurence, expand, clean)
}
\arguments{
  \item{tax}{The taxonomy dataframe (required).}

  \item{occurence}{Remove phyla with less than X
  occurences.}

  \item{expand}{Given a specific phylum, use the class
  level assignments instead.}

  \item{clean}{Simplify some of the phylum names.}
}
\value{
A cleaned taxonomy file.
}
\description{
A nice long description
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

