SGoF <- function(scores, alpha) {
  observed <- sum(scores < alpha)
  expected <- qbinom(prob=alpha, size=length(scores), p=1 - alpha)
  diff <- max(observed - expected + 1, 2)
  return(order(scores, decreasing=F)[1:diff])
}

threshold <- function(scores, threshold) {
  selected <- which(scores < threshold)
  selected <- selected[order(scores[selected], decreasing=F)]
  if (length(selected) == 0) selected <- order(scores, decreasing=F)[1:2]
  return(selected)
}

nKeep <- function(scores, nKeep) {
  return(order(scores, decreasing=F)[1:nKeep])
}

#' Title Generate a division of a data set into folds
#'
#' @param targetValues the class labels
#' @param nFolds number of folds
#' @param stratified whether folds should be stratified by \code{targetValues}
#'
#' @return a vector giving the fold for each sample
#' @export
#'
#' @examples
#' classes <- sample(c(rep(TRUE, 20), rep(FALSE, 50)), 70)
#' folds <- generateFolds(classes, 10, T)
generateFolds <- function(targetValues, nFolds, stratified=T) {
  nDataItems <- length(targetValues)
  if (stratified) {
    folds <- numeric(nDataItems)
    nPos <- sum(targetValues)
    nNeg <- sum(!targetValues)
    folds[targetValues] <- sample(rep(1:nFolds, length.out=nPos), nPos)
    folds[!targetValues] <- sample(rep(1:nFolds, length.out=nNeg), nNeg)
  } else {
    folds <- sample(rep(1:nFolds, length.out=nDataItems), nDataItems)
  }
  return(folds)
}

#' Run a set of classification models on input training, test data sets
#'
#' @param data.train a data frame of training data
#' @param targetValues a logical vector
#' @param data.test a data frame of test data (columns must match \code{data.train})
#' @param models a list of \link{caret::train} models to run
#' @param kFolds number of folds for model selection within each fold
#' @param repeats number of repeats for model selection within each fold
#' @param tuneLength number of parameters to tune
#' @param verbose verbose output if TRUE
#'
#' @return a data frame containing prediction probabilities for each classification algorithm.
#'These are the predicted probabililty of \code{targetValues==TRUE}
#'
#' @importFrom caret train trainControl
ClassifierModels <- function(data.train, targetValues, data.test,
                             models, kFolds=2, repeats=5, tuneLength=5,
                             verbose=F){

  stopifnot(is.logical(targetValues))
  stopifnot(is.data.frame(data.train))
  stopifnot(is.data.frame(data.test))
  ##convert target values to a more convenient form
  targetValues = as.factor(targetValues)
  levels(targetValues) <- c("Neg", "Pos")

  tc <- trainControl(method="repeatedcv",
                     number=kFolds, repeats = repeats,
                     verboseIter = verbose,
                     savePredictions = "final",
                     search="grid",
                     classProbs=T)

  predictions <- matrix(0, nrow=nrow(data.test), ncol=length(models))

  colnames(predictions) <- models
  rownames(predictions) <- rownames(data.test)

  for (model in models) {
    trainedModel <- tryCatch(train(data.train, targetValues,
                                   method=model,
                                   trControl = tc,
                                   tuneLength=tuneLength),
                             error=function(e) {
                               # convert error to a warning and return NA so we can keep going
                               warning(e)
                               return(NA)
                               }
                             )
    if (!is.na(trainedModel)) {
      predictions[, model] <- predict(trainedModel, data.test, type="prob")[, 1]
    } else {
      predictions[, model] <- NA
    }
  }

  return(predictions)
}

#' Utility function to create predictions for a given fold and set of variables
#'
#' @param trainingData A data frame of training data
#' @param trainingTargets A logical vector of training class labels
#' @param testData A data frame of test data
#' @param models A vector of \link{caret::train} models
#' @param tuneKFolds number of folds for parameter tuning
#' @param tuneRepeats number of repeats for parameter tuning
#' @param PCA apply PCA?
#' @param verbose verbose output?
#' @param tuneLength number of parameters to try when tuning
#' @param keep numeric vector of column indices to keep in training and test data
#' @param extraTrainingData additional training data
#' @param extraTestData additional test data
#' @param heatmap plot a heatmap before training models?
#'
#' @return predictions for this fold
runFold <- function(trainingData,
                    trainingTargets,
                    testData,
                    models,
                    tuneKFolds,
                    tuneRepeats,
                    tuneLength,
                    keep,
                    PCA,
                    verbose,
                    heatmap,
                    extraTrainingData,
                    extraTestData
                    ) {

  if (is.null(keep)==FALSE){
    trainingData = trainingData[, keep, drop=F]
    testData     = testData[, keep, drop=F]
    if (heatmap){
      working           = as.matrix(trainingData)
      rownames(working) = trainingTargets
      heatmap(t(working), trace="none")
    }
  }
  ##OPTION TO ADD IN EXTRA FEATURES, AFTER THE SELECTION
  if (!is.null(extraTrainingData) & !is.null(extraTestData)){
    colnames(extraTrainingData) <- colnames(extraTestData) <- paste("ED", seq(ncol(extraTrainingData)), sep="")
    trainingData = cbind(trainingData, extraTrainingData)
    testData     = cbind(testData,     extraTestData)
  }
  ##OPTION TO RUN PCA ON THE FEATURES
  if (PCA){
    model.pca = prcomp(trainingData)
    cum.var   = 0
    k         = 1
    while (cum.var<0.95) {  ## Choose the principal components that explain 95% of the variance
      k       = k+1        # ALWAYS KEEP AT LEAST TWO PCs
      cum.var = sum(model.pca$sdev[1:k]^2)/sum(model.pca$sdev^2)
    }
    trainingData = as.data.frame(model.pca$x[,1:k])
    cat('No. Principal Components= ', k, fill=TRUE)
    ## TRANSFORM TEST DATA USING PCA MODEL
    object.testData = predict(model.pca, testData)
    testData        = as.data.frame(object.testData[,1:k])
  }
  ##TRAIN MODELS; MAKE PREDICTIONS
  working = ClassifierModels(trainingData,
                             trainingTargets,
                             testData,
                             models=models,
                             kFolds=tuneKFolds,
                             repeats=tuneRepeats,
                             tuneLength=tuneLength,
                             verbose=verbose)
  return(working)
}


#' Cross-validation for classification models
#'
#' @param data.train Training data to be divided into folds (must be a data
#'   frame)
#' @param targetValues Target responses
#' @param models list of \link{caret::train} models to train
#' @param nFolds number of folds
#' @param stratified TRUE for stratified folds
#' @param threshold list of threshold p-values for selecting features to keep
#' @param nKeep list of number of features to keep
#' @param SGoF alpha for Sequential Goodness of Fit selection of features
#' @param verbose TRUE for verbose output
#' @param heatmap TRUE to plot a heatmap of selected features for each fold
#' @param PCA TRUE to apply PCA to selected features in each fold
#' @param extraData any additional data to add to training data after feature
#'   selection
#' @param tuneKFolds number of folds for tuning models within each fold
#' @param tuneRepeats number of repeats for tuning models within each fold
#' @param tuneLength number of parameters to test when tuning models
#' @param folds optionally pre-specify which samples go in which fold. Should be
#'   NULL to select folds randomly, or a vector of length
#'   \code{nrow(data.train)} containing values in \code{seq(nFolds)}
#'
#' @return A list of predictions for the given model combinations, ready to be
#'   passed to \link{CrossValRocCurves}
#' @export
CrossValidation <- function(data.train,
                            targetValues,
                            models=c("rf", "glmnet", "svmRadial", "svmLinear", "gbm", "nnet", "glm"),
                            nFolds=10,
                            stratified=TRUE,
                            threshold=NULL,
                            nKeep=NULL,
                            SGoF=NULL,
                            verbose=FALSE,
                            heatmap=FALSE,
                            PCA=FALSE,
                            extraData=NULL,
                            tuneKFolds=2, tuneRepeats=5, tuneLength=5,
                            folds=NULL){
  ##----------------------------------------------------------------------
  ## ASSERTIONS ABOUT THE INPUT ------------------------------------------
  ##----------------------------------------------------------------------
  stopifnot(is.logical(targetValues))
  stopifnot(is.data.frame(data.train))
  if (is.null(threshold) & is.null(nKeep) & is.null(SGoF)) {
    stop("You need to set at least one of threshold, nKeep and SGoF.")
  }
  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  nDataItems   = length(targetValues)
  featureNames = colnames(data.train)
  ##----------------------------------------------------------------------
  ## RUN THE CROSS-VALIDATION --------------------------------------------
  ##----------------------------------------------------------------------
  nModels <- length(models)

  if (is.null(folds)) {
    folds <- generateFolds(targetValues=targetValues,
                           nFolds=nFolds,
                           stratified=stratified)
  } else {
    # do some checks
    if(nFolds != max(folds)) stop("folds does not match nfolds")
    if(!all(seq(max(folds)) %in% folds)) stop("Some folds are empty")
    if(length(folds) != length(targetValues)) stop("length(folds) != length(targetValues")
  }


  out <- list()

  out$targetValues <- targetValues

  out$folds <- folds

  out$scores <- list()


  predictions <- matrix(0, nDataItems, nModels)
  colnames(predictions) <- models
  predictions <- as.data.frame(predictions)
  row.names(predictions) <- row.names(data.train)

  if (!is.null(SGoF)) {
    for (i in SGoF) {
      out$predictions$SGoF[[as.character(i)]] <- predictions
      out$selected$SGoF[[as.character(i)]] <- list()
    }
  }

  if (!is.null(nKeep)) {
    for (i in nKeep) {
      out$predictions$nKeep[[as.character(i)]] <- predictions
      out$selected$nKeep[[as.character(i)]] <- list()
    }
  }

  if (!is.null(threshold)) {
    for (i in threshold) {
      out$predictions$threshold[[as.character(i)]] <- predictions
      out$selected$threshold[[as.character(i)]] <- list()
    }
  }

  foldProgressBar <- tcltk::tkProgressBar(title="Fold progress",
                                          label=paste("0/", nFolds, sep=""),
                                          min = 0,
                                          max = nFolds,
                                          initial = 0,
                                          width = 300)
  for (fold in 1:nFolds){
    cat('Starting fold ', fold, " of ", nFolds, fill=TRUE)

    tcltk::setTkProgressBar(pb=foldProgressBar,
                            label=paste(fold, "/", nFolds, sep=""),
                            value = fold)
    ##SEPARATE DATA INTO TRAINING, TEST SETS
    index           = which(folds == fold)
    testData        = data.train[index,]
    trainingData    = data.train[-index,]
    trainingTargets = targetValues[-index]
    extraTrainingData = extraData[-index, , drop=FALSE]
    extraTestData = extraData[index, , drop=FALSE]
    ##OPTION TO RUN SIMPLE FEATURE PRE-SELECTION USING THE CURRENT TRAINING SET
    scores <- FeatureSelection(trainingData, trainingTargets)

    out$scores[[fold]] <- scores

    for (alpha in SGoF) {
      selected <- SGoF(scores, alpha)
      out$selected$SGoF[[as.character(alpha)]][[fold]] <- selected

      cat("SGoF ", alpha, ": ", length(selected), "variables selected", fill=TRUE)

      working <- runFold(trainingData=trainingData,
                         trainingTargets=trainingTargets,
                         testData=testData,
                         models=models,
                         tuneKFolds=tuneKFolds,
                         tuneRepeats=tuneRepeats,
                         tuneLength=tuneLength,
                         keep=selected,
                         PCA=PCA,
                         verbose=verbose,
                         heatmap=heatmap,
                         extraTrainingData=extraTrainingData,
                         extraTestData=extraTestData
      )

      out$predictions$SGoF[[as.character(alpha)]][index, ] = working
    }

    for (keep in nKeep) {
      selected <- nKeep(scores, keep)
      out$selected$nKeep[[as.character(keep)]][[fold]] <- selected

      cat("Keeping ", keep, " top variables", fill=TRUE)

      working <- runFold(trainingData,
                         trainingTargets,
                         testData,
                         models,
                         tuneKFolds,
                         tuneRepeats,
                         tuneLength,
                         selected,
                         PCA,
                         verbose,
                         heatmap,
                         extraTrainingData=extraTrainingData,
                         extraTestData=extraTestData
      )

      out$predictions$nKeep[[as.character(keep)]][index, ] = working
    }

    for (i in threshold) {
      selected <- threshold(scores, i)
      out$selected$threshold[[as.character(i)]][[fold]] <- selected

      cat("Threshold ", i, ": ", length(selected), "variables selected", fill=TRUE)

      working <- runFold(trainingData,
                         trainingTargets,
                         testData,
                         models,
                         tuneKFolds,
                         tuneRepeats,
                         tuneLength,
                         selected,
                         PCA,
                         verbose,
                         heatmap,
                         extraTrainingData=extraTrainingData,
                         extraTestData=extraTestData
      )

      out$predictions$threshold[[as.character(i)]][index, ] = working
    }


  }
  close(foldProgressBar)
  return(out)
}

#EXPERIMENTAL
GPLVMCrossValidation <- function(data.train,
                                 targetValues,
                                 models=c("rf", "glmnet", "svmRadial", "svmLinear", "gbm", "nnet"),
                                 nFolds=10,
                                 stratified=TRUE,
                                 sigma_d=10^(c(-5,-3,-1,0)),
                                 q=3,
                                 verbose=FALSE,
                                 heatmap=FALSE,
                                 PCA=TRUE,
                                 extraData=NULL,
                                 tuneKFolds=2, tuneRepeats=5, tuneLength=5,
                                 ...) {
  if (!requireNamespace("GPLVM", quietly = TRUE)) {
    stop("Package GPLVM needed for this function to work. Please install it.",
         call. = FALSE)
  }

  stopifnot(is.logical(targetValues))
  stopifnot(is.array(data.train))

  ##----------------------------------------------------------------------
  ## FIND USEFUL VALUES --------------------------------------------------
  ##----------------------------------------------------------------------
  nDataItems   = length(targetValues)
  featureNames = colnames(data.train)
  classes <- factor(targetValues)
  data.train.unstructured <- t(apply(data.train, 1, as.numeric))
  ##----------------------------------------------------------------------
  ## RUN THE CROSS-VALIDATION --------------------------------------------
  ##----------------------------------------------------------------------
  nModels <- length(models)

  folds <- generateFolds(targetValues=targetValues,
                         nFolds=nFolds,
                         stratified=stratified)


  out <- list()

  out$targetValues <- targetValues

  out$folds <- folds

  out$scores <- list()


  predictions <- matrix(0, nDataItems, nModels)
  colnames(predictions) <- models
  predictions <- as.data.frame(predictions)
  row.names(predictions) <- row.names(data.train)
  for (i in sigma_d) {
    out$predictions$GPLVM[[as.character(i)]] <- predictions
    out$model$GPLVM[[as.character(i)]] <- list()
  }

  for (fold in 1:nFolds){
    cat('Starting fold ', fold, " of ", nFolds, fill=TRUE)
    ##SEPARATE DATA INTO TRAINING, TEST SETS
    index           = which(folds == fold)
    testData.unstructured        = data.train.unstructured[index,]
    attr(testData.unstructured, faimsDimName) =c(512, 51)
    trainingData.unstructured    = data.train.unstructured[-index,]
    attr(trainingData.unstructured, faimsDimName) =c(512, 51)
    trainingTargets = targetValues[-index]
    trainingClasses = classes[-index]
    testData <- convertToArray(testData.unstructured)
    trainingData <- convertToArray(trainingData.unstructured)
    extraTrainingData = extraData[-index, ]
    extraTestData = extraData[index, ]

    for (sigma in sigma_d) {

      out$model$GPLVM[[as.character(sigma)]]$train[[fold]] <-
        GPLVM::fit.lsa_bcsgplvm(X=trainingData,
                                q=q,
                                classes=as.numeric(trainingClasses),
                                Z.prior="discriminative",
                                Z.prior.params=list(classes=trainingClasses,
                                                    sigma_d=sigma),
                                verbose=verbose,
                                save.X=FALSE,
                                ...)

      out$model$GPLVM[[as.character(sigma)]]$test[[fold]] <-
        GPLVM:::predict.LSA_BCSGPLVM(out$model$GPLVM[[as.character(sigma)]]$train[[fold]],
                                    testData,
                                    training.data=trainingData)

      pairs(rbind(out$model$GPLVM[[as.character(sigma)]]$train[[fold]]$final.Z,
                  out$model$GPLVM[[as.character(sigma)]]$test[[fold]]$predictions),
            col=c(as.numeric(trainingClasses), as.numeric(classes[index]) + 2))

      working <- runFold(trainingData=out$model$GPLVM[[as.character(sigma)]]$train[[fold]]$final.Z,
                         trainingTargets=trainingTargets,
                         testData=out$model$GPLVM[[as.character(sigma)]]$test[[fold]]$predictions,
                         models=models,
                         tuneKFolds=tuneKFolds,
                         tuneRepeats=tuneRepeats,
                         tuneLength=tuneLength,
                         keep=1:q,
                         PCA=PCA,
                         verbose=verbose,
                         heatmap=heatmap,
                         extraTrainingData=extraTrainingData,
                         extraTestData=extraTestData
      )

      out$predictions$GPLVM[[as.character(sigma)]][index, ] = working
    }
  }
  return(out)
}
