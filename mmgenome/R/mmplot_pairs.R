#' Generates a pairs plot
#'
#' A nice long description
#'
#' @usage mmplot_pairs(data, variables, color, nsize)
#'
#' @param data The dataframe to plot.
#' @param variables The vaiables to plot as a list.
#' @param color Coloring of the plot (default: gc)
#' @param nsize Resize the points (default: 1)
#' 
#' @return A pairs plot object.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmplot_pairs <- function(data, variables, color = "gc", nsize = 1){
  
  if(color == "gc"){
    rgb.c <- colorRampPalette(c("red", "green", "blue"))
    rgb.a <- adjustcolor(rgb.c(max(data$gc)-min(data$gc)+1), alpha.f = 0.3)
    palette(rgb.a)
    pairs(data[,variables], upper.panel=NULL, col = data$gc-min(data$gc), cex = sqrt(data$length)/75*nsize, pch=20)
  }
  
  if(color == "phylum"){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length=n+1)
      hcl(h=hues, l=65, c=100)[1:n]
    }
    data$phylum <- factor(data$phylum, levels=c(levels(data$phylum),"None"))
    data$phylum[is.na(data$phylum)] <- "None"
    pal <- c(gg_color_hue(length(levels(data$phylum))-1),"#0000000D")
    palette(pal)
    data <- data[with(data, order(rev(phylum))), ]
    pairs(data[,variables], upper.panel=NULL, col = data$phylum, cex = sqrt(data$length)/75*nsize, pch=20)
  }
  
  if(color != "gc" & color != "phylum"){
    rgb.c <- colorRampPalette(c("red", "green", "blue"))
    rgb.a <- adjustcolor(rgb.c(100), alpha.f = 0.5)
    palette(rgb.a)
    diff <- log10(max(data[,color])+1) - log10(min(data[,color])+1)
    tcol <- (log10(data[,color]+1) -log10(min(data[,color])+1)) * (99 / diff) +1
    pairs(data[,variables], upper.panel=NULL, col = tcol, cex = sqrt(data$length)/75*nsize, pch=20)
  }  
}
