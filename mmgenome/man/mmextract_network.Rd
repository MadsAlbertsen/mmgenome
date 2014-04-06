\name{mmextract_network}
\alias{mmextract_network}
\title{Extract scaffolds connected by paired-end reads}
\usage{
mmextract_network(subset, network, original, nconnections, type)
}
\arguments{
  \item{subset}{The subset dataframe to use for
  extraction.}

  \item{network}{The network of connected scaffold.}

  \item{original}{The original dataframe (default: d).}

  \item{nconnections}{The minumum number of connections
  needed (default: 2).}

  \item{type}{Extract direct or complete (all) linked
  scaffolds (default: direct).}
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

