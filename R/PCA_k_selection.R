evidence_for_k <- function(prcompObj, k) {
  log.evidence <- numeric(length(k))
  # Using formula from Minka 2000, "Automatic Choice of dimensionality for PCA"
  eigenvalues <- (prcompObj$sdev)^2
  N <- ncol(prcompObj$rotation)
  d <- length(eigenvalues)
  for (i in 1:length(k)) {
    m <- d * k[i] - k[i] * (k[i] + 1) / 2
    if (k[i]==d) {
      stop("Cannot assess evidence when k==d using Minka's approach")
    } else if (k[i] > d) {
      stop("k > d specified. You can't choose k PCs from d PCs if k > d")
    }
    noise <- mean(eigenvalues[(k[i]+1):length(eigenvalues)])
    log.evidence[i] <- (  (-N / 2) * sum(log(eigenvalues[1:k[i]]))
                        - (N * (d - k[i]) / 2) * log(noise)
                        - ((m + k[i]) / 2) * log(N))
  }
  return(log.evidence)
}

select_k <- function(prcompObj) {
  max.k <- length(prcompObj$sdev) - 1
  if (max.k == 1) return(2)
  evidence <- evidence_for_k(prcompObj, 1:max.k)
  continue <- TRUE
  for (i in 1:(max.k - 1)) {
    if (evidence[i] > evidence[i + 1]) {
      return(i)
    }
  }
  return(length(prcompObj$sdev))
}
