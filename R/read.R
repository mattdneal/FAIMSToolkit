faimsDimName <- "faimsDim"

#' Read in FAIMS data from a list of files
#'
#' @param fileList A list of files containing FAIMS data (in ASCII format)
#'
#' @return a data matrix containing the FAIMS data (with one row per file)
#' @export
ReadInFaimsData <- function(fileList, dec=".", minFlowRate=1.9) {
  ##----------------------------------------------------------------------
  ## FIND ALL THE DATA FILES TO READ IN ----------------------------------
  ##----------------------------------------------------------------------
  nFiles     = length(fileList)

  nlines <- 51

  numCVMeasurements <- 512

  numPixels <- numCVMeasurements * nlines * 2

  dataMatrix = matrix(0, nFiles, numPixels)

  colnames(dataMatrix)

  ##----------------------------------------------------------------------
  ## READ IN ALL THE DATA FILES ------------------------------------------
  ##----------------------------------------------------------------------
  for (i in 1:nFiles){
    monitoringData <- matrix(scan(fileList[i], skip=63, nlines=nlines, quiet=T, what=character(), dec=dec), nlines, byrow=T)
    min.flow = min(monitoringData[, 15])
    paste(min.flow)
    if (min.flow < minFlowRate) stop(paste(fileList[i], "has flow rate below minimum of", minFlowRate))
    data.positiveIon = scan(fileList[i], skip=116, nlines=nlines, quiet=TRUE, dec=dec)
    data.negativeIon = scan(fileList[i], skip=168, nlines=nlines, quiet=TRUE, dec=dec)
    dataVector       = c(data.positiveIon, data.negativeIon)
    dataMatrix[i,]   = dataVector
  }
  ##----------------------------------------------------------------------
  ## ITEM, FEATURE NAMES -------------------------------------------------
  ##----------------------------------------------------------------------

  attr(dataMatrix, faimsDimName) <- c(numCVMeasurements, nlines)

  return(dataMatrix)
}

#' Read in FAIMS files from a set of directories
#'
#' Search for directories containing FAIMS data below a specified directory, and return a named matrix with one row per directory.
#' This function will check that there is at least one ingested file for each directory in /code{dataPath}, and return a warning
#' if not. If multiple matching files are found an error will be thrown.
#'
#' @param dataPath the directory containing the data
#' @param filePattern a regular expression. All files matching this regular expression will be included.
#'
#' @return A named matrix with one row per data file read
#' @export
ReadInFaimsDirectories <- function(dataPath, filePattern, dec=".", minFlowRate=1.9){

  ##----------------------------------------------------------------------
  ## FIND DATA FILES, READ IN THE DATA -----------------------------------
  ##----------------------------------------------------------------------
  dirList    = list.dirs(dataPath, recursive=FALSE)
  dataFiles  = list.files(dirList, pattern=filePattern, full.names=TRUE)

  folderNames <- sapply(strsplit(dirList, "/"), function(x) x[length(x)])
  fileFolderNames <- sapply(strsplit(dataFiles, "/"), function(x) x[length(x) - 1])

  for (folder in folderNames) {
    fileCount <- sum(fileFolderNames==folder)
    if (fileCount > 1) {
      stop(paste("Too many matching files in", folder, "for filePattern", filePattern))
    }
    if (fileCount == 0) {
      warning(paste("No matching file found in", folder, "for filePattern", filePattern))
    }
  }

  dataMatrix = ReadInFaimsData(dataFiles, dec=dec, minFlowRate=minFlowRate)
  row.names(dataMatrix) = fileFolderNames

  return(dataMatrix)
}

#' Read and combine multiple FAIMS data files per sample
#'
#' Search for directories containing FAIMS data below a specified directory, and return a named matrix with one row per directory.
#' Multiple file patterns may be passed as a character vector to combine multiple files per directory.
#'
#'
#' @param dataPath the directory containing the data
#' @param filePatterns a character vector of regular expressions. All files matching this regular expression will be included.
#'
#' @return a named matrix with one row per directory
#' @export
ReadInFaimsDirectoriesMultiFile <- function(dataPath, filePatterns, dec=".", minFlowRate=1.9) {
  dataMatrix <- NULL
  for (filePattern in filePatterns) {
    newData <- ReadInFaimsDirectories(dataPath, filePattern, dec=dec, minFlowRate=minFlowRate)
    if (is.null(dataMatrix)) {
      dataMatrix <- newData
      faimsDim <- attr(dataMatrix, faimsDimName)
    } else {
      # Check that the rownames match
      currentNames <- rownames(dataMatrix)
      newNames <- rownames(newData)
      if (all(currentNames %in% newNames)
          & all(newNames %in% currentNames)) {
        dataMatrix <- cbind(dataMatrix[currentNames, ], newData[currentNames, ])
        rownames(dataMatrix) <- currentNames
        attr(dataMatrix, faimsDimName) <- faimsDim
      } else {
        print("Folder names so far:")
        print(currentNames)
        print("New folder names:")
        print(newNames)
        stop(paste("Change in sample names for", filePattern))
      }
    }
  }

  return(dataMatrix)
}

#' Read FAIMS data generated using an autosampler
#'
#' @param dataPath The directory to scan for data
#' @param numRuns The number of runs to take
#' @param filePattern a regex which should match the names of the files you want to ingest
#' @param runNumPattern a regex which should match the run number of your data.
#'        Any characters found will be stripped and the result converted using as.numeric
#'
#' @return
#' @export
#'
#' @examples
ReadInFaimsDirectoriesAutosampler <- function(dataPath, numRuns,
                                              filePattern='.*[.](txt|asc)',
                                              runNumPattern='[0-9]+[.](txt|asc)',
                                              dec=".", minFlowRate=1.9) {

  dirList    = list.dirs(dataPath, recursive=FALSE)
  folderNames <- sapply(strsplit(dirList, "/"), function(x) x[length(x)])

  numSamples <- length(dirList)
  dataMatrix <- NULL

  rowNum <- 1
  for (dir in dirList) {
    cat("Processing folder", folderNames[rowNum], fill=TRUE)
    dataFiles  = list.files(dir, pattern=filePattern, full.names=TRUE)
    runNumMatches <- regexpr(runNumPattern, dataFiles)
    runNums <- as.numeric(gsub('[^0-9]', '', regmatches(dataFiles, runNumMatches)))
    # Reorder the files so we read them in the order they were measured
    dataFiles <- dataFiles[order(runNums)]
    if (length(dataFiles) < numRuns) {
      stop(paste("Insufficient runs for sample", dir))
    }
    if (length(dataFiles) > numRuns) {
      warning(paste("Not using all", length(dataFiles), "runs for sample", dir))
      cat(paste("Not using all", length(dataFiles), "runs for sample", dir), fill=TRUE)
      dataFiles <- dataFiles[1:numRuns]
    }
    sampleRow <- as.numeric(t(ReadInFaimsData(dataFiles, dec=dec, minFlowRate=minFlowRate)))
    if (is.null(dataMatrix)) {
      dataMatrix <- matrix(0, nrow=numSamples, ncol=length(sampleRow))
    }
    dataMatrix[rowNum, ] <- sampleRow
    rowNum <- rowNum + 1
  }

  rownames(dataMatrix) <- folderNames

  nlines <- 51
  numCVMeasurements <- 512
  attr(dataMatrix, faimsDimName) <- c(numCVMeasurements, nlines)

  return(dataMatrix)
}

#' Convert FAIMS data to an array
#'
#' @param data Faims data matrix
#'
#' @return
#' @export
convertToArray <- function(data) {
  nrow <- nrow(data)
  faimsDim <- attr(data, faimsDimName)
  nruns <- prod(dim(data)) / (nrow * prod(faimsDim) * 2)
  dim <- c(nrow, faimsDim, 2, nruns)
  dimNamesPrefix <- c("sampleID", "CV", "dispersionStrength", "polarity", "runNumber")
  dimNames <- list()
  dimNames[[1]] <- rownames(data)
  for (i in 2:5) {
    dimNames[[i]] <- paste(dimNamesPrefix[i], 1:dim[i], sep=".")
  }

  out <- array(data, dim=dim, dimnames = dimNames)
  return(out)
}
