\name{especies}
\alias{especies}
\title{Extracts overview of essential genes at Phylum level}
\usage{
especies(eg)
}
\arguments{
  \item{eg}{Pre-formatted file of genomes and essential
  genes: data(eg)}
}
\value{
A list with two objects. $plot a ggplot object, $table a
table with summery data.
}
\description{
A nice long description
}
\examples{
data(eg)
eg_alpha <- subset(eg, Class == "Alphaproteobacteria")
species_data <- especies(eg_alpha)
species_data$plot
speces_data$table[1,]
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

