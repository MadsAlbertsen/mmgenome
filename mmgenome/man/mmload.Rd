\name{mmload}
\alias{mmload}
\title{Import wrapper for coverage files}
\usage{
mmimport_coverage(assembly, coverage, type)
}
\arguments{
  \item{assembly}{The assembly.}

  \item{coverage}{A list of coverage datasets.}

  \item{type}{The type of coverage files either clc or
  simple (default: simple).}

  \item{pca}{Add pca analyis of tetranucleotide frequencies
  (default: T)}
}
\value{
a dataframe with scaffold coverage, gc and lenght
information.
}
\description{
A nice long description.
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

