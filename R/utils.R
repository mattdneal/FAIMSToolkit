progBarInit <- function(range) {
  cat("0%")
  cat(rep(" ", 17), sep="")
  cat("25%")
  cat(rep(" ", 16), sep="")
  cat("50%")
  cat(rep(" ", 16), sep="")
  cat("75%")
  cat(rep(" ", 16), sep="")
  cat("100%\n")
  cat("|")
  out <- list(range=range, current=range[1], done=FALSE)
}

progBarUpdate <- function(bar.obj, new.value) {
  if (!bar.obj$done) {
    while(bar.obj$current < new.value) {
      if (bar.obj$current >= bar.obj$range[2]
                             - (bar.obj$range[2] - bar.obj$range[1]) / 80) {
        cat("|\n")
        bar.obj$done <- TRUE
      } else {
        cat("-")
      }
      bar.obj$current <- bar.obj$current + (bar.obj$range[2] - bar.obj$range[1]) / 80
    }
  }
  return(bar.obj)
}

#' Delete samples from a FAIMS object
#'
#' @param FAIMSObject a FAIMS object
#' @param deleteIndices indices to delete
#'
#' @return a FAIMS object
#' @export
deleteFAIMSSample <- function(FAIMSObject, deleteIndices) {
  if (length(deleteIndices) != 0) {
    data <- FAIMSObject$data[-deleteIndices, ]
    minFlowRate <- FAIMSObject$minFlowRate[-deleteIndices]
    out <- FAIMSObjectFactory(data=data,
                              faimsDim=FAIMSObject$faimsDim,
                              minFlowRate=minFlowRate)
    if (!is.null(FAIMSObject$arrayData)) {
      out$arrayData <- FAIMSObject$arrayData[-deleteIndices,,,,]
      class(out) <- c(class(out), "FAIMSArray")
    }
  }
}
