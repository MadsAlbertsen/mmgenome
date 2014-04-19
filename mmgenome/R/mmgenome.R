#' mmgenome: multi-metagenome
#'
#' Tools for extracting individual genomes from metagenomes
#'
#' @section Individual functions:
#'
#' \itemize{
#'   \item \code{\link{mmload}}: Merge data function.
#'   \item \code{\link{mmplot}}: Default plot function.
#'   \item \code{\link{mmplot_selction}}: Highlights a selection made with mmplot_locater on a mmplot.
#'   \item \code{\link{mmplot_locator}}: Adds "locator" functionallity to ggplot2 objects.
#'   \item \code{\link{mmplot_pairs}}: Pairs plot.
#'   \item \code{\link{mmplot_network}}: Network plot.
#'   \item \code{\link{mmextract}}: Extracts a subset of scaffolds based on input from \code{\link{mmplot_locator}}.
#'   \item \code{\link{mmextract_network}}: Extracts PE/MP connected scaffolds.
#'   \item \code{\link{mmstats}}: An overview of key statistics of the scaffolds.
#'   \item \code{\link{mmstats_duplicates}}: Identifies duplicated essential single copy genes.
#'   \item \code{\link{mmref_especies}}: Reference data on essential genes.
#'   \item \code{\link{mmref_ephylum}}: Reference data on essential genes on phylum level.
#'   \item \code{\link{mmref_compare}}: Compare essential gene content to selected references.
#'   \item \code{\link{mmexport}}: Export selected scaffolds to fasta file.
#' }
#'
#' @docType package
#' @name mmgenome
NULL