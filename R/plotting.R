#' Function that plots two types of data points against it's chromosomal location.
#' Note: This is a plot PER chromosome.
#' @noRd
#'
create.haplotype.plot = function(chrom.position, points_blue, points_red, x_min, x_max, title, xlab, ylab) {
  par(pch=".", cex=1, cex.main=0.8, cex.axis = 0.6, cex.lab=0.7,yaxp=c(-0.05,1.05,6))
  plot(c(x_min,x_max), c(0,1), type="n",, main=title, xlab=xlab, ylab=ylab)
  points(chrom.position, points_blue, col="blue")
  points(chrom.position, points_red, col="red")
}