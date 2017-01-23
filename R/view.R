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
plotFAIMSdata <- function(dataMatrix, rowsToPlot, plotLayout=1, absolute=T, numRows=512, ...) {
  if (absolute) {
    dataMatrix <- abs(dataMatrix)
  }
  layout(plotLayout)

  for(i in rowsToPlot) {
    image(matrix(dataMatrix[i, ], nrow=numRows),
          zlim=range(dataMatrix),
          main=rownames(dataMatrix)[i],
          axes=FALSE,
          ...)
  }
  layout(1)
  return(NULL)
}

