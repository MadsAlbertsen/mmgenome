\name{mmplot}
\alias{mmplot}
\title{Wraps ggplot2 for pretty default plots}
\usage{
mmplot(data, x, y, logx, logy, color)
}
\arguments{
  \item{data}{A dataframe with scaffold information.}

  \item{x}{x-axis variable.}

  \item{y}{y-axis variable.}

  \item{logx}{log10 scale the x-axis (default True)}

  \item{logy}{log10 scale the y-axis (default True)}

  \item{color}{Color by "phylum", "gc" or a coverage
  dataset (default: phylum)}

  \item{minlength}{Minimum length of plotted scaffolds.}
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

