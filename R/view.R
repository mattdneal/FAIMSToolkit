#' Plot FAIMS data matrices
#'
#' @param dataMatrix a matrix of FAIMS data
#' @param rowsToPlot the rows to plot
#' @param plotLayout a matrix to specify the layout for multiple plots per page (default is one plot per page)
#' @param absolute whether to plot the absolute value of the FAIMS data.
#' @param ... additional parameters to pass to \code{image}
#'
#' @return NULL
#' @export
plotFAIMSdata <- function(FAIMSObject, rowsToPlot, plotLayout=1, absolute=T, ...) {
  if (!requireNamespace("viridis", quietly = TRUE)) {
    warning('Base-R colour palettes are not colour-blind friendly and transfer
poorly to greyscale. Install viridis [install.packages("viridis")]
to use a better palette.')
    color.palette <- topo.colors(512)
  } else {
    color.palette <- viridis::plasma(512)
  }

  if (!inherits(FAIMSObject, "FAIMS")) stop("FAIMSObject must inherit class FAIMS")
  dataMatrix <- FAIMSObject$data
  faimsDim <- FAIMSObject$faimsDim
  if (absolute) {
    dataMatrix <- abs(dataMatrix)
  }
  layout(plotLayout)

  mar.original <- par("mar")
  par(mar=c(0.2, 0.5, 1.1, 0.5))

  for(i in rowsToPlot) {
    image(matrix(dataMatrix[i, ], nrow=faimsDim[1]),
          zlim=range(dataMatrix),
          main=rownames(dataMatrix)[i],
          axes=FALSE,
          useRaster=TRUE,
          col=color.palette,
          ...)
    box()
  }
  layout(1)
  par(mar=mar.original)
  invisible(NULL)
}

