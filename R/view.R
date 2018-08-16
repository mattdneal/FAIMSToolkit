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


#' Plot FAIMS data matrices with axes and ion current scale
#'
#' @param dataMatrix a matrix of FAIMS data
#' @param rowToPlot the row to plot
#' @param runToPlot the run to plot
#'
#' @return NULL
#' @import lattice
#' @importFrom grid grid.text
#' @export
prettyFAIMSPlot <- function(FAIMSObject, rowToPlot, runToPlot=1, title="") {
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
  plotMatrix <- abs(matrix(dataMatrix[rowToPlot, ], nrow=faimsDim[1])[, 1:faimsDim[2] + (runToPlot - 1) * faimsDim[2]])

  zlim <- range(plotMatrix)
  col.at <- seq(zlim[1], zlim[2], length.out = 256)


  colnames(plotMatrix) <- seq(0, 100, length.out = faimsDim[2])
  rownames(plotMatrix) <- round(seq(-6, 6, length.out=faimsDim[1]), 2)
  rownames(plotMatrix)[round(seq(1, faimsDim[1], length.out=5))] <- c(-6, -3, 0, 3, 6)

  print(levelplot(plotMatrix,
            aspect=1,
            col.regions=color.palette,
            at=col.at,
            contour=F,
            cuts=256,
            pretty=T,
            xlab="Compensation Voltage (V)",
            ylab="Dispersion Field (%)",
            scales=list(x=list(at=round(seq(1, faimsDim[1], length.out=5))),
                        y=list(at=round(seq(1, faimsDim[2] ,length.out=6)))),
            main=title,
            useRaster=T))
  trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
  grid::grid.text("Ion count\n(arb. units)", 0.2, -0.02, hjust=0.35, vjust=1.1)
  trellis.unfocus()
}
