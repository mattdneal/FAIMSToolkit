#' Find mean k-nearest neighbour distance for each sample
#'
#' Compute mean k-nn distances for outlier detection
#'
#' @param dataMatrix a data matrix
#' @param max.k max k for mean k-nn distance
#'
#' @return a \code{nrow(dataMatrix)} by \code{k} matrix with the mean k-nn
#'   distance for each sample, with k along the columns and samples along the
#'   rows
#' @export
knnMeanDistance <- function(dataMatrix, max.k=10) {
  numSamples <- nrow(dataMatrix)
  distMatrix <- as.matrix(dist(dataMatrix))
  scores <- matrix(0, nrow=numSamples, ncol=max.k)
  rownames(scores) <- rownames(dataMatrix)
  for (i in 1:numSamples) {
    dists <- distMatrix[i, -i]
    for (k in 1:max.k) {
      scores[i, k] <- mean(dists[order(dists, decreasing = F)[1:k]])
    }
  }
  return(scores)
}
