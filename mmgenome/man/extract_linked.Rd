\name{extract_linked}
\alias{extract_linked}
\title{Extract scaffolds connected by paired-end reads}
\usage{
extract_linked(subset, original, network, minimum)
}
\arguments{
  \item{subset}{The subset of extracted scaffolds}

  \item{original}{The original dataframe with all scaffolds
  (default: d)}

  \item{network}{The network of connected scaffold
  (default: n)}

  \item{minumum}{The minumum number of connections needed
  (default: 5)}
}
\value{
A dataframe with scaffold information.
}
\description{
A nice long description.
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

