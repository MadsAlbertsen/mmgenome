\docType{data}
\name{data}
\alias{data}
\title{Wraps ggplot2 for pretty default plots}
\format{\preformatted{List of 2
 $ scaffolds:'data.frame':	91 obs. of  10 variables:
  ..$ scaffold: int [1:91] 192 209 216 355 442 596 636 677 707 709 ...
  ..$ length  : int [1:91] 37855 99827 132251 23810 39343 61203 102931 27270 21924 56336 ...
  ..$ gc      : num [1:91] 63 64.8 63.2 64.2 61.9 ...
  ..$ S26     : num [1:91] 2.81 3.16 3.34 2.86 8.38 ...
  ..$ S27     : num [1:91] 12.8 14.1 14.1 11.9 25.7 ...
  ..$ PC1     : num [1:91] 0.136 0.16 0.14 0.148 0.132 ...
  ..$ PC2     : num [1:91] 0.07809 0.02156 0.06879 0.06958 0.00621 ...
  ..$ PC3     : num [1:91] -0.227 -0.239 -0.243 -0.293 -0.017 ...
  ..$ phylum  : Factor w/ 9 levels "Actinobacteria",..: NA NA NA NA NA NA NA NA NA NA ...
  ..$ 16S     :List of 91
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
  .. ..$ : logi NA
 $ essential:'data.frame':	187 obs. of  3 variables:
  ..$ scaffold: int [1:187] 1341 36731 713 26473 26473 713 713 26473 26473 713 ...
  ..$ orf     : int [1:187] 43 125 59 147 137 49 61 149 134 46 ...
  ..$ hmm.id  : Factor w/ 109 levels "PF00162.14","PF00276.15",..: 1 1 2 2 3 3 4 4 5 5 ...
}}
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
\keyword{datasets}

