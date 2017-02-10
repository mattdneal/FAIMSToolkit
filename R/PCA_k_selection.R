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
  best.k <- 1
  best.k.evidence <- evidence_for_k(prcompObj, best.k)
  continue <- TRUE
  for (k in 2:(length(prcompObj$sdev) - 1)) {
    new.k <- k
    new.k.evidence <- evidence_for_k(prcompObj, new.k)
    if (new.k.evidence > best.k.evidence) {
      best.k <- new.k
      best.k.evidence <- new.k.evidence
    } else {
      return(best.k)
    }
  }
  return(length(prcompObj$sdev))
}
