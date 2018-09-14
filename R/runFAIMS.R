setParams <- function(bestModel) {
  if (bestModel$selectionMethod == "SGoF") {
    nKeep = NULL
    SGoF = as.numeric(bestModel$threshold)
    threshold=NULL
  } else if (bestModel$selectionMethod == "nKeep") {
    nKeep = as.numeric(bestModel$threshold)
    SGoF = NULL
    threshold=NULL
  } else if (bestModel$selectionMethod == "threshold"){
    nKeep = NULL
    SGoF = NULL
    threshold=as.numeric(bestModel$threshold)
  } else {
    stop(paste("Unknown selection method", bestModel$selectionMethod))
  }
  model <- bestModel$classifierMethod
 return(list(nKeep=nKeep, SGoF=SGoF, threshold=threshold, model=model))
}

#' Run a standard analysis given a FAIMS object and class labels
#'
#' @param FAIMSObject a FAIMS object
#' @param targetValues class labels
#' @param models a list of \link{caret::train} models
#' @param modelSelectFolds pre-generated folds for model selection
#' @param modelSelectScores pre-generated scores for model selection
#' @param bestModelFolds pre-generated folds for best model assessment
#' @param bestModelScores pre-generated scores for best model assessment
#' @param waveletData pre-computed wavelet data
#' @param SGoF Select variables using sequential goodness of fit? (only for PCA
#'   analysis)
#' @param nKeep Select variables using keep top N?
#' @param extraData Additional data to feed to the classifier
#'
#' @return A list of results (see out$bestModelSummary and
#'   out$modelSelectSummary for a summary of results)
#' @export
runFAIMS <- function(FAIMSObject, targetValues,
                     models=c("rf", "glmnet", "svmRadial", "svmLinear", "gbm", "nnet", "glm"),
                     modelSelectFolds=NULL, modelSelectScores=NULL,
                     bestModelFolds=NULL, bestModelScores=NULL,
                     waveletData=NULL, SGoF=TRUE, nKeep=TRUE,
                     extraData=NULL) {
  out <- list()
  if (is.null(waveletData)) {
    if (is.null(FAIMSObject)) stop("FAIMSObject not supplied")
    out$waveletData <- as.data.frame(WaveletTransform(FAIMSObject))
  } else {
    out$waveletData <- waveletData
  }
  out$targetValues.scrambled <- sample(targetValues)

  out$modelSelect <- list()
  out$bestModel <- list()

  # Select a classifier method
  if (is.null(modelSelectFolds)) {
    out$modelSelect$folds <- generateFolds(targetValues, 10, T)
  } else {
    out$modelSelect$folds <- modelSelectFolds
  }

  if (nKeep) {
    nKeep.pca <- 2^(0:10)
    nKeep.nopca <- 2^(0:8)
  } else {
    nKeep.pca <- NULL
    nKeep.nopca <- NULL
  }

  if (SGoF) {
    SGoF.pca <- c(0.2, 0.1, 0.05, 0.01, 0.001)
  } else {
    SGoF.pca <- NULL
  }

  out$modelSelect$pca.cv <- CrossValidation(out$waveletData,
                                            targetValues,
                                            models=models,
                                            PCA=T,
                                            nKeep=nKeep.pca,
                                            SGoF=SGoF.pca,
                                            folds=out$modelSelect$folds,
                                            nFolds=10,
                                            precomputedScores=modelSelectScores,
                                            extraData=extraData)

  out$modelSelect$scores <- out$modelSelect$pca.cv$scores
  out$modelSelect$nopca.cv <- CrossValidation(out$waveletData,
                                              targetValues,
                                              models=models,
                                              PCA=F,
                                              nKeep=nKeep.nopca,
                                              folds=out$modelSelect$folds,
                                              nFolds=10,
                                              precomputedScores=out$modelSelect$scores,
                                              extraData=extraData
                                              )

  out$modelSelect$pca.results <- CrossValRocCurves(out$modelSelect$pca.cv)
  out$modelSelect$nopca.results <- CrossValRocCurves(out$modelSelect$nopca.cv)


  out$modelSelect$pca.bestModel <-
    out$modelSelect$pca.results$summary[which.max(out$modelSelect$pca.results$summary$auc), ]
  out$modelSelect$nopca.bestModel <-
    out$modelSelect$nopca.results$summary[which.max(out$modelSelect$nopca.results$summary$auc), ]

  # Take the best classifier for PCA and no PCA and retrain with new fold selection
  if (is.null(bestModelFolds)) {
    out$bestModel$folds <- generateFolds(targetValues, 10, T)
  } else {
    out$bestModel$folds <- bestModelFolds
  }

  out$bestModel$scores <- bestModelScores
  out$bestModel$scrambled.scores <- NULL

  out$modelSelectSummary <-
    rbind(cbind(PCA=TRUE, out$modelSelect$pca.results$summary),
          cbind(PCA=FALSE, out$modelSelect$nopca.results$summary))

  out$bestModelSummary <- NULL
  out$bestModelSummary.scrambled <- NULL



  params <- tryCatch(setParams(out$modelSelect$pca.bestModel),
                     error = function(e) {
                       warning(e)
                       return(NA)
                     }
  )
  if (!identical(NA, params)) {
    out$bestModel$pca.cv <-
      CrossValidation(out$waveletData,
                      targetValues,
                      models=unique(c("glm", params$model)),
                      PCA=T,
                      nKeep=params$nKeep,
                      SGoF=params$SGoF,
                      threshold=params$threshold,
                      folds=out$bestModel$folds,
                      nFolds=10,
                      precomputedScores=bestModelScores,
                      extraData=extraData)
    out$bestModel$scores <- out$bestModel$pca.cv$scores

    out$bestModel$pca.results <- CrossValRocCurves(out$bestModel$pca.cv)

    #rerun with scrambled classes
    out$bestModel$scrambled.pca.cv <-
      CrossValidation(out$waveletData,
                      out$targetValues.scrambled,
                      models=unique(c("glm", params$model)),
                      PCA=T,
                      nKeep=params$nKeep,
                      SGoF=params$SGoF,
                      threshold=params$threshold,
                      folds=out$bestModel$folds,
                      nFolds=10,
                      extraData=extraData)
    out$bestModel$scrambled.scores <- out$bestModel$scrambled.pca.cv$scores
    out$bestModel$scrambled.pca.results <- CrossValRocCurves(out$bestModel$scrambled.pca.cv)

    out$bestModelSummary <- rbind(out$bestModelSummary,
                                  cbind(PCA=TRUE, out$bestModel$pca.results$summary)
    )
    out$bestModelSummary.scrambled <- rbind(out$bestModelSummary.scrambled,
                                            cbind(PCA=TRUE, out$bestModel$scrambled.pca.results$summary)
    )
  }
  params <- tryCatch(setParams(out$modelSelect$nopca.bestModel),
                     error = function(e) {
                       warning(e)
                       return(NA)
                     }
  )
  if (!identical(NA, params)) {
    out$bestModel$nopca.cv <-
      CrossValidation(out$waveletData,
                      targetValues,
                      models=unique(c("glm", params$model)),
                      PCA=F,
                      nKeep=params$nKeep,
                      SGoF=params$SGoF,
                      threshold=params$threshold,
                      folds=out$bestModel$folds,
                      nFolds=10,
                      precomputedScores=out$bestModel$scores,
                      extraData=extraData)

    out$bestModel$nopca.results <- CrossValRocCurves(out$bestModel$nopca.cv)

    out$bestModel$scrambled.nopca.cv <-
      CrossValidation(out$waveletData,
                      out$targetValues.scrambled,
                      models=unique(c("glm", params$model)),
                      PCA=F,
                      nKeep=params$nKeep,
                      SGoF=params$SGoF,
                      threshold=params$threshold,
                      folds=out$bestModel$folds,
                      nFolds=10,
                      precomputedScores=out$bestModel$scrambled.scores,
                      extraData=extraData)
    out$bestModel$scrambled.nopca.results <- CrossValRocCurves(out$bestModel$scrambled.nopca.cv)


    out$bestModelSummary <- rbind(out$bestModelSummary,
                                  cbind(PCA=FALSE, out$bestModel$nopca.results$summary)
    )
    out$bestModelSummary.scrambled <- rbind(out$bestModelSummary.scrambled,
                                            cbind(PCA=FALSE, out$bestModel$scrambled.nopca.results$summary)
    )
  }

  return(out)
}
