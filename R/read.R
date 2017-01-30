
FAIMSObjectFactory <- function(data, faimsDim, minFlowRate) {
  out <- list()
  class(out) <- "FAIMS"
  out$data <- data
  out$faimsDim <- faimsDim
  out$minFlowRate <- minFlowRate
  return(out)
}

#' Read in FAIMS data from a list of files
#'
#' @param fileList A list of files containing FAIMS data (in ASCII format)
#' @param dec Decimal symbol
#' @param minFlowRate minimum acceptable flow rate. If min flow rate falls below
#'   this threshold an error is thrown
#'
#' @return A FAIMS object
#' @export
ReadInFaimsData <- function(fileList,
                            dec=".",
                            minFlowRate=1.9) {
  ##----------------------------------------------------------------------
  ## FIND ALL THE DATA FILES TO READ IN ----------------------------------
  ##----------------------------------------------------------------------
  nFiles     = length(fileList)

  #Read the number of lines and CV measurements from the first file
  sample.file <- scan(file=fileList[1], what="character", sep = "\n", quiet=TRUE)

  count <- 0
  line=53
  notAComment <- TRUE
  while (notAComment) {
    if (substr(sample.file[line], 1, 2) != "//") {
      count <- count + 1
      line <- line + 1
    } else {
      notAComment <- FALSE
    }
  }
  nlines <- count

  numCVMeasurements <- as.numeric(strsplit(sample.file[21], "\t")[[1]][3])

  faimsDim <- c(numCVMeasurements, nlines)
  names(faimsDim) <- c("CompensationVoltage", "DispersionFieldStrength")

  numPixels <- numCVMeasurements * nlines * 2

  data = matrix(0, nFiles, numPixels)

  minFlowRateVec <- numeric(nFiles)
  ##----------------------------------------------------------------------
  ## READ IN ALL THE DATA FILES ------------------------------------------
  ##----------------------------------------------------------------------
  for (i in 1:nFiles){
    monitoringData <- matrix(scan(fileList[i],
                                  skip=63,
                                  nlines=nlines,
                                  quiet=T,
                                  what=character(),
                                  dec=dec),
                             nlines, byrow=T)
    minFlowRateVec[i] = min(monitoringData[, 15])
    if (minFlowRateVec[i] < minFlowRate) {
      stop(paste(fileList[i], "has flow rate below minimum of", minFlowRate))
    }
    data.positiveIon = scan(fileList[i], skip=116, nlines=nlines, quiet=TRUE, dec=dec)
    data.negativeIon = scan(fileList[i], skip=168, nlines=nlines, quiet=TRUE, dec=dec)
    dataVector = c(data.positiveIon, data.negativeIon)
    data[i,] = dataVector
  }
  out <- FAIMSObjectFactory(data=data, faimsDim=faimsDim, minFlowRate=minFlowRateVec)
  out$ReadInFaimsData <- list()
  out$ReadInFaimsData$argList <- as.list(match.call())
  return(out)
}

#' Read in FAIMS files from a set of directories
#'
#' Search for directories containing FAIMS data below a specified directory, and
#' return a named matrix with one row per directory. This function will check
#' that there is at least one ingested file for each directory in
#' /code{dataPath}, and return a warning if not. If multiple matching files are
#' found an error will be thrown.
#'
#' @param dataPath the directory containing the data
#' @param filePattern a regular expression. All files matching this regular
#'   expression will be included.
#' @param dec Decimal symbol
#' @param minFlowRate minimum acceptable flow rate. If min flow rate falls below
#'   this threshold an error is thrown
#'
#' @return A named matrix with one row per data file read
#' @export
ReadInFaimsDirectories <- function(dataPath,
                                   filePattern='.*[.](txt|asc)',
                                   dec=".",
                                   minFlowRate=1.9){

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

  out = ReadInFaimsData(dataFiles, dec=dec, minFlowRate=minFlowRate)
  out$ReadInFaimsDirectories <- list()
  out$ReadInFaimsDirectories$argList <- as.list(match.call())
  out$ReadInFaimsDirectories$dirList <- dirList
  out$ReadInFaimsDirectories$dataFiles <- dataFiles
  row.names(out$data) = fileFolderNames
  names(out$minFlowRate) = fileFolderNames

  return(out)
}

#' Read and combine multiple FAIMS data files per sample
#'
#' Search for directories containing FAIMS data below a specified directory, and
#' return a named matrix with one row per directory. Multiple file patterns may
#' be passed as a character vector to combine multiple files per directory.
#'
#' @param dataPath the directory containing the data
#' @param filePatterns a character vector of regular expressions. All files
#'   matching this regular expression will be included.
#' @param dec Decimal symbol
#' @param minFlowRate minimum acceptable flow rate. If min flow rate falls below
#'   this threshold an error is thrown
#'
#' @return a named matrix with one row per directory
#' @export
ReadInFaimsDirectoriesMultiFile <- function(dataPath,
                                            filePatterns='.*[.](txt|asc)',
                                            dec=".",
                                            minFlowRate=1.9) {
  ReadInFaimsData <- list()
  ReadInFaimsDirectories <- list()
  dataMatrix <- NULL

  for (filePattern in filePatterns) {

    newData <- ReadInFaimsDirectories(dataPath, filePattern, dec=dec, minFlowRate=minFlowRate)

    ReadInFaimsData[[filePattern]] <- newData$ReadInFaimsData
    ReadInFaimsDirectories[[filePattern]] <- newData$ReadInFaimsDirectories

    if (is.null(dataMatrix)) {

      dataMatrix <- newData$data
      faimsDim <- newData$faimsDim
      minFlowRateMatrix <- as.matrix(newData$minFlowRate)

    } else {

      # Check that the rownames match
      currentNames <- rownames(dataMatrix)
      newNames <- rownames(newData$data)

      if (all(currentNames %in% newNames)
          & all(newNames %in% currentNames)) {

        dataMatrix <- cbind(dataMatrix[currentNames, ],
                            newData$data[currentNames, ])
        rownames(dataMatrix) <- currentNames
        minFlowRateMatrix <- cbind(minFlowRateMatrix[currentNames, ],
                                   newData$minFlowRate[currentNames])
        rownames(minFlowRateMatrix) <- currentNames

      } else {

        message("Folder names so far:")
        message(currentNames)
        message("New folder names:")
        message(newNames)
        stop(paste("Change in sample names for", filePattern))

      }
    }
  }

  colnames(minFlowRateMatrix) <- filePatterns

  out <- FAIMSObjectFactory(data=dataMatrix, faimsDim=faimsDim, minFlowRate=minFlowRateMatrix)
  out$ReadInFaimsData <- ReadInFaimsData
  out$ReadInFaimsDirectories <- ReadInFaimsDirectories

  out$ReadInFaimsDirectoriesMultiFile <- list()
  out$ReadInFaimsDirectoriesMultiFile$argList <- as.list(match.call())

  return(out)

}

#' Read FAIMS data generated using an autosampler
#'
#' @param dataPath The directory to scan for data
#' @param numRuns The number of runs to take
#' @param filePattern a regex which should match the names of the files you want to ingest
#' @param runNumPattern a regex which should match the run number of your data.
#'        Any characters found will be stripped and the result converted using as.numeric
#' @param dec Decimal symbol
#' @param minFlowRate minimum acceptable flow rate. If min flow rate falls below
#'   this threshold an error is thrown
#'
#' @return A FAIMS object
#' @export
ReadInFaimsDirectoriesAutosampler <- function(dataPath, numRuns,
                                              filePattern='.*[.](txt|asc)',
                                              runNumPattern='[0-9]+[.](txt|asc)',
                                              dec=".", minFlowRate=1.9) {

  dirList    = list.dirs(dataPath, recursive=FALSE)
  folderNames <- sapply(strsplit(dirList, "/"), function(x) x[length(x)])

  numSamples <- length(dirList)
  dataMatrix <- NULL
  numFiles <- numeric(numSamples)
  names(numFiles) <- dirList
  minFlowRateMatrix <- matrix(0, numSamples, numRuns)

  rowNum <- 1
  for (dir in dirList) {
    cat("Processing folder", folderNames[rowNum], fill=TRUE)
    dataFiles  = list.files(dir, pattern=filePattern, full.names=TRUE)
    runNumMatches <- regexpr(runNumPattern, dataFiles)
    runNums <- as.numeric(gsub('[^0-9]', '', regmatches(dataFiles, runNumMatches)))
    # Reorder the files so we read them in the order they were measured
    dataFiles <- dataFiles[order(runNums)]
    numFiles[dir] <- length(dataFiles)
    if (length(dataFiles) < numRuns) {
      stop(paste("Insufficient runs for sample", dir))
    }
    if (length(dataFiles) > numRuns) {
      warning(paste("Not using all", length(dataFiles), "runs for sample", dir))
      cat(paste("Not using all", length(dataFiles), "runs for sample", dir), fill=TRUE)
      dataFiles <- dataFiles[1:numRuns]
    }
    sampleData <- ReadInFaimsData(dataFiles, dec=dec, minFlowRate=minFlowRate)
    sampleRow <- as.numeric(t(sampleData$data))
    minFlowRateMatrix[rowNum, ] <- sampleData$minFlowRate
    if (is.null(dataMatrix)) {
      dataMatrix <- matrix(0, nrow=numSamples, ncol=length(sampleRow))
      faimsDim <- sampleData$faimsDim
    }
    dataMatrix[rowNum, ] <- sampleRow
    rowNum <- rowNum + 1
  }

  rownames(dataMatrix) <- folderNames
  rownames(minFlowRateMatrix) <- folderNames

  out <- FAIMSObjectFactory(data=dataMatrix, faimsDim=faimsDim, minFlowRate=minFlowRateMatrix)
  out$ReadInFaimsDirectoriesAutosampler$dirList <- dirList
  out$ReadInFaimsDirectoriesAutosampler$argList <- as.list(match.call())

  return(out)
}

#' Convert FAIMS data to an array
#'
#' Converts a FAIMS object to a FAIMSArray object, with the data stored in a
#' 5-dimensional array. The dimensions correspond to sample ID, CV, dispersion
#' field strength, polarity, and runNumber.
#'
#' If keepMatrixData=TRUE, the output inherits both FAIMS and FAIMSArray.
#' Otherwise, it only inherits FAIMSArray and will not work with functions which
#' require FAIMS objects.
#'
#' @param FAIMSObject A FAIMS object
#' @param keepMatrixData Whether to keep the matrix-shaped FAIMS data
#'
#' @return A FAIMSArray object (also inheriting FAIMS if
#'   \code{keepMatrixData=TRUE})
#' @export
convertToArray <- function(FAIMSObject, keepMatrixData=TRUE) {
  if (!inherits(FAIMSObject, "FAIMS")) stop("FAIMSObject must inherit class FAIMS")
  data <- FAIMSObject$data
  nrow <- nrow(data)
  faimsDim <- FAIMSObject$faimsDim
  nruns <- prod(dim(data)) / (nrow * prod(faimsDim) * 2)
  dim <- c(nrow, faimsDim, 2, nruns)
  dimNamesPrefix <- c("SampleID", "CompensationVoltage", "DispersionFieldStrength", "Polarity", "RunNumber")
  dimNames <- list()
  dimNames[[1]] <- rownames(data)
  for (i in 2:5) {
    dimNames[[i]] <- paste(dimNamesPrefix[i], 1:dim[i], sep=".")
  }

  arrayData <- array(data, dim=dim, dimnames = dimNames)
  names(dim(arrayData)) <- dimNamesPrefix
  out <- FAIMSObject

  out$arrayData <- arrayData

  currentClass <- class(out)
  if (!keepMatrixData) {
    out$data <- NULL
    class(out) <- unique(c(currentClass[currentClass!= "FAIMS"], "FAIMSArray"))
  } else {
    class(out) <- unique(c(currentClass, "FAIMSArray"))
  }

  return(out)
}
