#' perform a feature selection using a Wilcoxon rank-sum test
#'
#' @param dataMatrix a numeric matrix of training data
#' @param targetValues class labels (must be a logical vector)
#'
#' @return a list of scores, one for each sample. Lower is better.
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
