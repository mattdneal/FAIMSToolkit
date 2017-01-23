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
