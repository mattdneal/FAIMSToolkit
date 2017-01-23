#' Find mean k-nearest neighbour distance for each sample
#'
#' @param dataMatrix
#' @param max.k
#'
#' @return
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
