\name{mmref_phylum}
\alias{mmref_phylum}
\title{Extracts overview of essential genes at Phylum level}
\usage{
mmref_phylum(eg)
}
\arguments{
  \item{eg}{Pre-formatted file of genomes and essential
  genes: data(eg)}
}
\value{
A list with three objects. $plot a ggplot object, $table a
table with summery data and $data with the data used for
plotting.
}
\description{
A nice long description
}
\examples{
data(eg)
phylum_data <- mmref_phylum(eg)
phylum_data$plot
phylum_data$table
}
\author{
Soren M. Karst \email{smk@bio.aau.dk}

Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

