#' perform a feature selection using a Wilcoxon rank-sum test
#'
#' @param dataMatrix
#' @param targetValues
#' @param threshold
#' @param nKeep
#'
#' @return
#' @export
FeatureSelection <- function(dataMatrix, targetValues){

  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  nDataItems = nrow(dataMatrix)
  nFeatures  = ncol(dataMatrix)
  ##----------------------------------------------------------------------
  ## FEATURE SELECTION ---------------------------------------------------
  ##----------------------------------------------------------------------
  pValues       = numeric(nFeatures)
  index.disease = which(targetValues==TRUE)
  index.control = which(targetValues==FALSE)
  data.disease  = dataMatrix[index.disease,]
  data.control  = dataMatrix[index.control,]
  cat("Wilcoxon testing", nFeatures, "features", fill=TRUE)
  bar <- progBarInit(c(1, nFeatures))
  for (i in 1:nFeatures){
    bar <- progBarUpdate(bar, i)
    current.disease = data.disease[, i]
    current.control = data.control[, i]
    outputObject    = wilcox.test(current.disease, current.control, exact=FALSE)
    pValues[i]      = outputObject$p.value
  }

  return(pValues)
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------
