\name{mmplot_network}
\alias{mmplot_network}
\title{Plots connected scaffolds in a network graph}
\usage{
mmplot_network(data, network)
}
\arguments{
  \item{data}{(required) The dataframe with the network to
  be plotted.}

  \item{network}{(required) A network file with connections
  between scaffolds.}

  \item{nconnections}{The minimum number of connections
  (default: 2).}

  \item{color}{Color by gc, phylum or a specific coverage
  dataset (default: "gc").}

  \item{log.color}{log10 scale the colors (default: F)}

  \item{labels}{If scaffold names are to be plotted
  (default: F).}
}
\value{
An igraph network plot.
}
\description{
A nice long description
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

