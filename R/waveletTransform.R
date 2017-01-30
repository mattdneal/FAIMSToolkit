#' Perform a 1D wavelet transform on a FAIMS data matrix
#'
#' @param FAIMSObject a FAIMS object
#' @param discardTopLevels number of levels to discard as noise (0 keeps all
#'   data)
#'
#' @return a matrix of wavelet-transformed FAIMS data
#' @export
#'
#' @importFrom wavethresh wd accessD
WaveletTransform <- function(FAIMSObject, discardTopLevels=2){
  if (!inherits(FAIMSObject, "FAIMS")) stop("FAIMSObject must inherit class FAIMS")

  dataMatrix <- FAIMSObject$data
  faimsDim <- FAIMSObject$faimsDim
  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  nDataItems         = nrow(dataMatrix)
  nFeatures          = ncol(dataMatrix)
  numRows            = faimsDim[1]
  numCols            = nFeatures / faimsDim[1]
  nWavelets          = 2^(ceiling(log2(numRows)) - discardTopLevels)
  data.wd            = matrix(0, nDataItems, nWavelets * numCols)
  row.names(data.wd) = row.names(dataMatrix)
  ##----------------------------------------------------------------------
  ## GENERATE THE WAVELET TRANSFORMED DATA -------------------------------
  ##----------------------------------------------------------------------

  pb <- progBarInit(c(0, nDataItems * numCols))
  count <- 0

  for (i in 1:nDataItems){
    waveColIndex <- 1
    for (j in 1:numCols) {
      count <- count + 1
      pb <- progBarUpdate(pb, count)
      colIndex             = 1:numRows + (j - 1) * numRows
      currentData          = dataMatrix[i, colIndex]
      working              = numeric(2^(ceiling(log2(numRows))))
      working[1:numRows]   = currentData
      current.wd           = wd(working, filter.number=1, bc="symmetric", family = "DaubExPhase")
      for (k in 0:(current.wd$nlevels - 1 - discardTopLevels)) {
        data.wd[i, waveColIndex:(waveColIndex + 2^k - 1)] <- accessD(current.wd, k)
        waveColIndex <- waveColIndex + 2^k
      }
    }
  }
  ##----------------------------------------------------------------------
  ## REMOVE ANY ZERO-VARIANCE FEATURES -----------------------------------
  ##----------------------------------------------------------------------
  sigmaValues = apply(data.wd, 2, sd)
  keep        = which(sigmaValues>0)
  data.wd     = data.wd[, keep]

  attr(data.wd, "keptColumns") <- keep

  return(data.wd)
}

#' Perform a 2D wavelet transform on a FAIMS data matrix
#'
#' This function expects an input matrix of dimension (nDataItems * nFeatures)
#' each row is therefore the 1D representation of the 2D data for a single item
#' This means that this function needs to know the underlying dimensionality of
#' the 2D FAIMS run
#'
#' @param FAIMSObject A FAIMS object
#' @param cropSize size to crop image to
#' @param discardHighestNLevels number of levels to discard as noise
#'
#' @return A matrix of wavelet-transformed FAIMS data
#' @export
#'
#' @importFrom wavethresh imwd lt.to.name
WaveletTransform_2D <- function(FAIMSObject, cropSize=NULL, discardHighestNLevels=2){
  if (!inherits(FAIMSObject, "FAIMS")) stop("FAIMSObject must inherit class FAIMS")

  dataMatrix <- FAIMSObject$data
  faimsDim <- FAIMSObject$faimsDim
  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  dataMatrix <- as.matrix(dataMatrix)
  nDataItems <- nrow(dataMatrix)
  nFeatures <- ncol(dataMatrix)

  pixelsPerRun <- prod(faimsDim)

  nFaimsRuns <- nFeatures / pixelsPerRun

  ##----------------------------------------------------------------------
  ## GENERATE THE WAVELET TRANSFORMED DATA -------------------------------
  ##----------------------------------------------------------------------

  pb <- progBarInit(c(0, nDataItems * nFaimsRuns))
  count <- 0

  for (i in 1:nDataItems){
    working.data <- c()
    for (runNum in 1:nFaimsRuns) {
      count <- count + 1
      pb <- progBarUpdate(pb, count)
      currentData <- dataMatrix[i, 1:pixelsPerRun + (runNum - 1) * pixelsPerRun]
      dim(currentData) <- faimsDim
      ##PAD IMAGE WITH ZEROS, TO MAKE A SQUARE MATRIX OF SIZE 2^N
      nSize <- 2^ceiling(log2(max(faimsDim)))
      data.padded <- matrix(0, nSize, nSize)
      data.padded[1:faimsDim[1], 1:faimsDim[2]] <- currentData

      ##OPTION TO CROP THE SQUARE MATRIX
      ##this is pretty specific to the FAIMS format....
      if (is.null(cropSize)==FALSE){
        edgeSize    = (nSize - cropSize) / 2
        lower       = 1     + edgeSize
        upper       = lower + cropSize - 1
        data.padded = data.padded[lower:upper, 1:cropSize]
      }
      ##WAVELET TRANSFORM
      current.wd            = imwd(data.padded,
                                   filter.number=1,
                                   bc="symmetric",
                                   family = "DaubExPhase")
      nLevels               = current.wd$nlevels-1-discardHighestNLevels
      if (nLevels < 0) stop("Can't discard more levels than are in the transform")
      working.data          = c(working.data, current.wd$w0Lconstant)
      for (j in 0:nLevels) {
        for (k in 1:4){
          working.data      = c(working.data,current.wd[[ lt.to.name(j, k) ]])
        }
      }
    }
    if (i==1){
      data.wd            = matrix(0, nDataItems, length(working.data))
      row.names(data.wd) = row.names(dataMatrix)
    }
    data.wd[i, ]         = working.data
  }

  return(data.wd)
}
