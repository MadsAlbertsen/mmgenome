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

  \item{logx}{log10 scale the x-axis (default True)}

  \item{logy}{log10 scale the y-axis (default True)}

  \item{color}{Color by "phylum", "gc" or a coverage
  dataset (default: phylum)}

  \item{minlength}{Minimum length of plotted scaffolds.}

  \item{network}{Network used to plot connections between
  scaffolds.}

  \item{nconnections}{Minimum number of connections to plot
  (default: 0).}

  \item{duplicates}{Mark scaffolds with duplicated
  essential genes (default: F)}

  \item{labels}{If scaffold names are to be plotted
  (default: F)}
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

