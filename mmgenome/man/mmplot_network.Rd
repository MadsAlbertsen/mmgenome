\name{mmplot_network}
\alias{mmplot_network}
\title{Plots connected scaffolds in a network graph}
\usage{
mmplot_network(data, network, nconnections, ewidth, nsize, labels, color)
}
\arguments{
  \item{data}{The dataframe with the network to be
  plotted.}

  \item{network}{A network file with connections between
  scaffolds.}

  \item{nconnections}{The minimum number of connections.}

  \item{ewidth}{Scale the width of edges by connections.}

  \item{nsize}{Scale the size of nodes by length.}

  \item{labels}{If scaffold names are to be plotted.}

  \item{color}{Color by gc, phylum or coverage.}
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

