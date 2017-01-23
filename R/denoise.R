#' Remove background noise from FAIMS data
#'
#' Denoises FAIMS data. First identifies pixels likely to be solely noise (mean
#' value across all samples close to zero). Then identifies a value to subtract
#' from all pixels by taking the fractionToRemove quantile of all the values in
#' the noise pixels, and subtracts this from the absolute value of each pixels,
#' zeroing any values which are negative after subtracting background noise.
#'
#' @param dataMatrix FAIMS data matrix
#' @param fractionNoise The percentage of pixels to identify as "noise"
#' @param fractionToRemove The quantile value of the noise pixels to subtract from
#' all pixels
#'
#' @return a denoised FAIMS data matrix
#' @export
denoiseFaimsData <- function(dataMatrix, fractionNoise=0.05, fractionToRemove=0.99) {
  dataColSums <- abs(colSums(dataMatrix))
  threshold <- quantile(dataColSums, fractionNoise)
  noise <- which(dataColSums < threshold)

  dataMatrix <- abs(dataMatrix)
  noiseThreshold <- quantile(dataMatrix[, noise], fractionToRemove)

  dataMatrix <- dataMatrix - noiseThreshold
  dataMatrix[dataMatrix < 0] <- 0

  return(dataMatrix)
}

getNeighbourIndices <- function(index, width, height, faimsDim) {
  rowNeighbours <- (index - width):(index + width)
  keepRowNeighbours <- ceiling(rowNeighbours / 512) == ceiling(index / 512)
  rowNeighbours <- rowNeighbours[keepRowNeighbours]

  neighbours <- c()
  for (i in 0:height) {
    neighbours <- c(neighbours, rowNeighbours + faimsDim[1] * i)
    if (i != 0) {
      neighbours <- c(neighbours, rowNeighbours - faimsDim[1] * i)
    }
  }

  numPixels <- prod(faimsDim)
  matrixNum <- ceiling(index / numPixels)

  keep <- which(neighbours > numPixels * (matrixNum - 1)
                & neighbours <= numPixels * matrixNum
                & neighbours != index)

  neighbours <- neighbours[keep]
  return(neighbours)
}

#' Remove background noise from FAIMS data using local correlations
#'
#' The signal in FAIMS data exhibits a high degree of local correlation. This
#' function zeroes out a fraction of pixels with the lowest degree of local
#' correlation.
#'
#' @param dataMatrix FAIMS data matrix
#' @param neighbourhoodSize the size of the neighbourhood in which to look for local correlations
#' @param alpha p-value required for a pixel to be considered correlated to its neighbours
#'
#' @return A FAIMS data matrix
#' @export
denoiseFaimsData.localCorr <- function(dataMatrix,
                                       neighbourhoodSizeCol=1, neighbourhoodSizeRow=1,
                                       alpha=0.05, plot=TRUE) {

  n <- nrow(dataMatrix)
  scores <- numeric(ncol(dataMatrix))
  faimsDim <- attr(dataMatrix, faimsDimName)
  dataMatrix <- scale(dataMatrix, center=T, scale=F)
  for (i in 1:ncol(dataMatrix)) {
    var.cor <- cor(rowMeans(dataMatrix[,getNeighbourIndices(i,
                                                            neighbourhoodSizeCol,
                                                            neighbourhoodSizeRow,
                                                            faimsDim)]),
                   dataMatrix[,i])
    t.value <- var.cor * sqrt(n - 2) / sqrt(1 - var.cor^2)
    scores[i] <- (1 - pt(abs(t.value), n - 2)) * 2
  }

  signal <- (1:ncol(dataMatrix) %in% SGoF(scores, alpha))

  dataMatrix <- abs(dataMatrix)
  if (length(signal) > 0) {
    half.index <- logical(ncol(dataMatrix))
    half.index[1:(ncol(dataMatrix) / 2)] <- TRUE
    max.1 <- as.numeric(quantile(dataMatrix[, half.index & !signal], 0.99))
    max.2 <- as.numeric(quantile(dataMatrix[, !half.index & !signal], 0.99))
    dataMatrix[, half.index] <- dataMatrix[, half.index] - max.1
    dataMatrix[, !half.index] <- dataMatrix[, !half.index] - max.2
    dataMatrix[as.logical(dataMatrix < 0)] <- 0
  } else {
    stop("No non-noise pixels found.")
  }

  attr(dataMatrix, faimsDimName) <- faimsDim
  return(dataMatrix)
}
