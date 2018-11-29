#' Check the minimum flow rates for all the runs in a dataset
#'
#' @param dir Directory pat to check
#' @param filePattern Regex pattern of exported FAIMS files
#'
#' @return data frame listing minimum flow rate per folder
#' @export
checkFlowRate <- function(dir, threshold=1.9, filePattern='.*[.](txt|asc)', moveTo=NULL) {
  out <- data.frame(sample.id=character(),
                    min.flow.rate=numeric(),
                    which.min=numeric(),
                    num.below.threshold=numeric(),
                    num.runs=numeric(),
                    stringsAsFactors = F)
  count <- 0

  if (!is.null(moveTo)) {
    dir.create(moveTo, recursive=T)
  }

  for (id in list.dirs(dir, full.names=F, recursive = F)) {
    count <- count + 1
    fileList  = list.files(file.path(dir, id), pattern=filePattern, full.names=TRUE, recursive=F)
    min.flow.rate <- Inf
    out[count, 1] <- id
    out[count, 4] <- 0
    out[count, 5] <- length(fileList)
    for (i in seq_along(fileList)) {
      monitoringData <- matrix(scan(fileList[i], skip=63, nlines=51, quiet=T, what=character()), 51, byrow=T)
      min.flow.rate.file <- min(monitoringData[, 15])
      if (min.flow.rate.file < min.flow.rate) {
        out[count, 3] <- i
        min.flow.rate <- min.flow.rate.file
      }
      if (min.flow.rate.file < threshold) {
        out[count, 4] <- out[count, 4] + 1
        if (!is.null(moveTo)) {
          dir.create(file.path(moveTo, id), recursive=T)
          file.rename(fileList[i], file.path(moveTo, id, basename(fileList[i])))
        }
      }
    }
    out[count, 2] <- min.flow.rate
  }
  return(out)
}

#' Check for FAIMS drift
#'
#' @param faims a faims object
#' @param index the samples to include. NULL (the default) includes all samples
#' @param ... other parameters to pass to plot
#'
#' @return a list with the correlation, a 95% confidence interval, and the matched rankings.
#' @export
checkDrift <- function(faims, index=NULL, ...) {
  if (is.null(index)) {
    index <- 1:nrow(faims$data)
  }
  faims_dist <- as.numeric(dist(faims$data[index,]))
  time_dist <- as.numeric(dist(faims$timestamps[index, 1]))
  include <- time_dist != 0
  faims_dist <- faims_dist[include]
  time_dist <- time_dist[include]
  data <- data.frame(time_rank=rank(time_dist, ties.method="random"), dist_rank=rank(faims_dist, ties.method = "random"))
  plot(data[,1], data[,2], xlab="Time distance rank", ylab="Data distance rank", ...)
  boot.out <- boot::boot(data, function(data, index) cor(data[index,1], data[index,2]), R=10000)
  cor <- cor(data[,1], data[,2])
  cor.ci <- boot::boot.ci(boot.out, conf = 0.95, type = "bca")
  return(list(cor=cor, cor.ci=cor.ci, data=data))
}
