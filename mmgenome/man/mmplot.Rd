\name{mmplot}
\alias{mmplot}
\title{Wraps ggplot2 for pretty default plots}
\usage{
mmplot(data, x, y)
}
\arguments{
  \item{data}{(required) A dataframe with scaffold
  information.}

  \item{x}{(required) x-axis variable.}

  \item{y}{(required) y-axis variable.}

  \item{color}{Color by a specific variable (default:
  phylum)}

  \item{log.x}{log10 scale the x-axis (default: T)}

  \item{log.y}{log10 scale the y-axis (default: T)}

  \item{log.color}{log10 scale the colors (default: F)}

  \item{minlength}{Minimum length of plotted scaffolds.}

  \item{network}{Network used to plot connections between
  scaffolds.}

  \item{nconnections}{Minimum number of connections to plot
  (default: 0).}

  \item{duplicates}{Mark scaffolds with duplicated
  essential genes (default: F)}

  \item{labels}{If scaffold names are to be plotted
  (default: F)}

  \item{resize}{Rescale the size of the scaffolds (default:
  1)}
}
\value{
a ggplot2 object
}
\description{
A nice long description.
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

