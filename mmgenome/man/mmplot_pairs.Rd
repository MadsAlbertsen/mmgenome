\name{mmplot_pairs}
\alias{mmplot_pairs}
\title{Generates a ggplot2 style pairs plot}
\usage{
mmplot_pairs(data, variables, color, nsize)
}
\arguments{
  \item{data}{(required) The dataframe to plot.}

  \item{variables}{(required) The vaiables to plot as a
  list.}

  \item{color}{Coloring of the plot (default: gc)}

  \item{minlength}{Minimum length of the plotted
  scaffolds.}

  \item{log}{A vector of variables to log scale.}

  \item{log.color}{Log scale the color (default: F)}

  \item{resize}{Rescale the size of the scaffolds (default:
  1)}
}
\value{
A pairs plot object.
}
\description{
A nice long description
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

