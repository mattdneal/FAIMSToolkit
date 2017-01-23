#' Perform a 1D wavelet transform on a FAIMS data matrix
#'
#' @param dataMatrix FAIMS data
#'
#' @return a matrix of wavelet-transformed FAIMS data
#' @export
#'
#' @importFrom wavethresh wd accessD
WaveletTransform <- function(dataMatrix, discardTopLevels=1){
  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  nDataItems         = nrow(dataMatrix)
  nFeatures          = ncol(dataMatrix)
  faimsDim           = attr(dataMatrix, faimsDimName)
  numRows            = faimsDim[1]
  numCols            = nFeatures / faimsDim[1]
  nWavelets          = 2^(ceiling(log2(numRows)) - discardTopLevels)
  data.wd            = matrix(0, nDataItems, nWavelets * numCols)
  row.names(data.wd) = row.names(dataMatrix)
  ##----------------------------------------------------------------------
  ## GENERATE THE WAVELET TRANSFORMED DATA -------------------------------
  ##----------------------------------------------------------------------
  for (i in 1:nDataItems){
    waveColIndex <- 1
    for (j in 1:numCols) {
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
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------

##Function to perform a 2D wavelet transform on an input FAIMS data matrix
##
##This function expects an input matrix of dimension (nDataItems * nFeatures)
##each row is therefore the 1D representation of the 2D data for a single item
##This means that this function needs to know the underlying dimensionality of the 2D FAIMS run
##
#' Perform a 2D wavelet transform on a FAIMS data matrix
#'
#' This function expects an input matrix of dimension (nDataItems * nFeatures)
#' each row is therefore the 1D representation of the 2D data for a single item
#' This means that this function needs to know the underlying dimensionality of the 2D FAIMS run
#'
#' @param dataMatrix
#' @param cropSize
#'
#' @return
#' @export
#'
#' @importFrom wavethresh imwd
WaveletTransform_2D <- function(dataMatrix, cropSize=NULL, discardHighestNLevels=0){

  library(wavethresh)
  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  dataMatrix <- as.matrix(dataMatrix)
  nDataItems <- nrow(dataMatrix)
  nFeatures <- ncol(dataMatrix)

  faimsDim <- attr(dataMatrix, faimsDimName)

  pixelsPerRun <- prod(faimsDim)

  nFaimsRuns <- nFeatures / pixelsPerRun

  ##----------------------------------------------------------------------
  ## GENERATE THE WAVELET TRANSFORMED DATA -------------------------------
  ##----------------------------------------------------------------------
  for (i in 1:nDataItems){
    working.data <- c()
    for (runNum in 1:nFaimsRuns) {
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
                                   filter.number=10,
                                   family="DaubLeAsymm",
                                   type="wavelet",
                                   bc="symmetric")
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

  # length(current.wd2$w4L4)   ## 8 levels, level 7 has length 16384 (4 components)
  ## level 6 has length 4096
  ## level 5 has length 1024 ....
  # Number of columns for wavelet transform according to matrix size:
  #   512x512--> 349525
  #   256x256--> 87381
  #   128x128--> 21845 ## additional benefit of dimensionality reduction

  return(data.wd)
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------
