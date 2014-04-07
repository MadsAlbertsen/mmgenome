\name{mmload}
\alias{mmload}
\title{Data loading function}
\usage{
mmload(assembly, coverage, essential, tax)
}
\arguments{
  \item{assembly}{The assembly.}

  \item{coverage}{A list of coverage datasets.}

  \item{coverge.type}{The type of coverage files either clc
  or simple (default: simple).}

  \item{essential}{}

  \item{tax}{}

  \item{tax.freq}{}

  \item{tax.expand}{}

  \item{tax.clean}{}

  \item{pca}{Add pca analyis of tetranucleotide frequencies
  (default: T)}

  \item{pca.ncomp}{}

  \item{rRNA16S}{}

  \item{rRNA23S}{}

  \item{rRNA.type}{}
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

