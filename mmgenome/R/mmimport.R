#' Import data
#'
#' Imports data through a Rmarkdown template.
#'
#' @usage mmimport(file)
#'
#' @param file The Rmarkdown file that contains the load definitions.
#' @param save Save an image of the importated data for fast loading (default: T).
#' @param force Forces re-import of data if an previous image is present (default: F).
#' @param save.filename Filename of the image to save or load (default: "data.RData")
#' 
#' @return Loaded data as detailed in the Rmarkdown file.
#' 
#' @export
#' @import knitr
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' mmimport(file = "Load_data.Rmd")
#' }


mmimport <- function(file, save = T, force = F, save.filename = "data.RData"){
  
  if (force == T) {
    knit(input = file, tangle = TRUE)
    Rfile <- gsub(".Rmd", ".R", file)
    source(file = Rfile)
    if (save == T){save.image(file = save.filename)}
  } else{
    if (file.exists(save.filename)){
      load(file=save.filename, envir = globalenv())
      print("Data already generated, loading it directly!")  
    }  else {
      knit(input = file, tangle = TRUE)
      Rfile <- gsub(".Rmd", ".R", file)
      source(file = Rfile)
      if (save == T){save.image(file = save.filename)}
    }
  }  
}
