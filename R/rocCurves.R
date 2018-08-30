#' Plot a ROC curve with CI
#'
#' @param rocCurve A ROC object
#' @param titleString Title for plot
#' @param ci Plot confidence intervals if TRUE
#'
#' @return NULL
#'
#' @export
#' @importFrom pROC plot.roc
plotRocCurve <- function(rocCurve, titleString="", ci=TRUE) {
  if (!ci) {
    plot(rocCurve, lwd=4, cex=4, mar=c(5, 4, 4, 2)+.1)
  } else {
    rocobj = plot.roc(rocCurve, lwd=4, cex=4, mar=c(5, 4, 4, 2)+.1,
                      ci=TRUE, print.auc=FALSE)
    ciobj  = ci.se(rocobj, specificities=seq(0, 1, length=100))
    plot(ciobj, type="shape", col="gray")
  }

  auc               = as.numeric(rocCurve$auc)
  confidence        = as.numeric(rocCurve$ci)

  title(paste(titleString, " (auc=", format(auc, digits=2), ") (95% CI: ",
              format(confidence[1], digits=2), ", ", format(confidence[3], digits=2), ")", sep=""))
  return(NULL)
}

#' Return a ROC curve object, given input class probabilities
#'
#' @param predictions vector of predictions
#' @param targetValues vector of class labels
#' @param titleString title string if \code{plotCurve==TRUE}
#' @param ci logical - compute confidence interval?
#' @param lr logical - compute likelihood ratios?
#' @param plotCurve logical - plot the resulting ROC curve?
#'
#' @return a \link{pROC::roc} object
#' @export
#'
#'@importFrom pROC roc ci.se
getRocCurve <- function(predictions, targetValues, titleString="", ci=FALSE, lr=FALSE, plotCurve=FALSE) {

  ##----------------------------------------------------------------------
  ## GENERATE A ROC CURVE ------------------------------------------------
  ##----------------------------------------------------------------------
  rocCurve          = roc(as.numeric(targetValues), predictions, ci=TRUE, direction="auto")
  auc               = as.numeric(rocCurve$auc)
  confidence <- as.numeric(rocCurve$ci)
  sensitivityValues = rocCurve$sensitivities
  specificityValues = rocCurve$specificities
  ##----------------------------------------------------------------------
  ## PLOT THE ROC CURVE --------------------------------------------------
  ##----------------------------------------------------------------------
  if (plotCurve) {
    plotRocCurve(rocCurve, titleString=titleString, ci=ci)
  }
  ##----------------------------------------------------------------------
  ## FIND A THRESHOLD TO BALANCE SENSITIVITY, SPECIFICITY ----------------
  ##----------------------------------------------------------------------
  ##use distance between null line and actual curve
  ##Obviously, this simplifies... :-)
  #  metricScore = (sensitivityValues - 1 + specificityValues)^2 +
  #                (specificityValues - 1 + sensitivityValues)^2
  #  index       = which.max(metricScore)

  ##this finds sens/spec values that are approximately equal
  metricScore = abs(sensitivityValues - specificityValues)
  index       = which.min(metricScore)

  sensitivity = sensitivityValues[index]
  specificity = specificityValues[index]
  threshold   = rocCurve$thresholds[index]
  ##----------------------------------------------------------------------
  ## FIND THE POSITIVE, NEGATIVE LIKELIHOOD RATIOS -----------------------
  ##----------------------------------------------------------------------
  lr.pos = sensitivity     / (1-specificity)
  lr.neg = (1-sensitivity) / specificity
  ##----------------------------------------------------------------------
  ## COMPUTE CONFIDENCE INTERVALS FOR SENSITIVITY, SPECIFICITY -----------
  ##----------------------------------------------------------------------
  if (rocCurve$direction=="<")
    predictedClass = predictions>threshold
  else
    predictedClass = predictions<threshold

  truePositive   = sum((predictedClass==TRUE)  & (targetValues==TRUE))
  trueNegative   = sum((predictedClass==FALSE) & (targetValues==FALSE))
  falsePositive  = sum((predictedClass==TRUE)  & (targetValues==FALSE))
  falseNegative  = sum((predictedClass==FALSE) & (targetValues==TRUE))
  sensitivity.ci = as.numeric(binom.test(truePositive, truePositive + falseNegative)$conf.int)
  specificity.ci = as.numeric(binom.test(trueNegative, trueNegative + falsePositive)$conf.int)
  precision      = truePositive / (truePositive + falsePositive)
  recall         = sensitivity
  f.measure      = 2 * precision*recall / (precision + recall)
  ##also recompute sens/spec, as a sanity check
  sens.2         = truePositive / (truePositive + falseNegative)
  spec.2         = trueNegative / (trueNegative + falsePositive)
  ##----------------------------------------------------------------------
  ## PRINT USEFUL INFO TO SCREEN -----------------------------------------
  ##----------------------------------------------------------------------
  cat(paste(titleString, "(ROC)"), fill=TRUE)
  cat("------------------", fill=TRUE)
  cat(paste("AUC         = ", format(auc,         digits=2)), "   (", format(confidence[1], digits=2), " - ",
      format(confidence[3], digits=2), ")", fill=TRUE, sep="")
  cat(paste("sensitivity = ", format(sensitivity, digits=2)), "   (", format(sensitivity.ci[1], digits=2), " - ",
      format(sensitivity.ci[2], digits=2), ")", fill=TRUE, sep="")
  cat(paste("specificity = ", format(specificity, digits=2)), "   (", format(specificity.ci[1], digits=2), " - ",
      format(specificity.ci[2], digits=2), ")", fill=TRUE, sep="")
  cat("\n")

  if (lr==TRUE){
    cat(paste("LR+ = ", format(lr.pos, digits=2)), fill=TRUE, sep="")
    cat(paste("LR- = ", format(lr.neg, digits=2)), fill=TRUE, sep="")
    cat("\n")
  }

  return(rocCurve)
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------

#' Compute ROC objects for the results of a cross validation
#'
#' @param crossVal.obj the output from \link{CrossValidation}
#'
#' @return A list of \link{pROC::roc} objects and a summary of the results
#' @export
CrossValRocCurves <- function(crossVal.obj) {
  out <- list()

  out$ROC.list <- list()

  out$summary <- data.frame(selectionMethod=character(0),
                            threshold=character(0),
                            classifierMethod=character(0),
                            auc=numeric(0),
                            auc.ci.lower=numeric(0),
                            auc.ci.higher=numeric(0),
                            sens=numeric(0),
                            sens.ci.lower=numeric(0),
                            sens.ci.higher=numeric(0),
                            spec=numeric(0),
                            spec.ci.lower=numeric(0),
                            spec.ci.higher=numeric(0),
                            ci.perc=numeric(0),
                            classifier.threshold=numeric(0),
                            stringsAsFactors = FALSE)

  data <- crossVal.obj$predictions

  targetValues <- crossVal.obj$targetValues

  for (selectionMethod in names(data)) {
    out$ROC.list[[selectionMethod]] <- list()
    for (threshold in names(data[[selectionMethod]])) {
      out$ROC.list[[selectionMethod]][[threshold]] <- list()
      preds <- data[[selectionMethod]][[threshold]]
      for (classifierMethod in colnames(preds)) {
        if (!any(is.na(preds[, classifierMethod]))) {
          rocCurve <- getRocCurve(predictions = preds[, classifierMethod],
                                  targetValues = targetValues,
                                  ci=TRUE,
                                  plotCurve=FALSE)

          out$ROC.list[[selectionMethod]][[threshold]][[classifierMethod]] <- rocCurve

          auc <- as.numeric(rocCurve$auc)
          confidence <- as.numeric(rocCurve$ci)

          sensitivityValues = rocCurve$sensitivities
          specificityValues = rocCurve$specificities

          ##this finds sens/spec values that are approximately equal
          metricScore = abs(sensitivityValues - specificityValues)
          index       = which.min(metricScore)

          sensitivity = sensitivityValues[index]
          specificity = specificityValues[index]
          classifier.threshold   = rocCurve$thresholds[index]

          if (rocCurve$direction=="<")
            predictedClass = preds[, classifierMethod]>classifier.threshold
          else
            predictedClass = preds[, classifierMethod]<classifier.threshold

          truePositive   = sum((predictedClass==TRUE)  & (targetValues==TRUE))
          trueNegative   = sum((predictedClass==FALSE) & (targetValues==FALSE))
          falsePositive  = sum((predictedClass==TRUE)  & (targetValues==FALSE))
          falseNegative  = sum((predictedClass==FALSE) & (targetValues==TRUE))
          sensitivity.ci = as.numeric(binom.test(truePositive, truePositive + falseNegative)$conf.int)
          specificity.ci = as.numeric(binom.test(trueNegative, trueNegative + falsePositive)$conf.int)

          line <- c(selectionMethod=selectionMethod,
                    threshold=threshold,
                    classifierMethod=classifierMethod,
                    auc=auc,
                    auc.ci.lower=confidence[1],
                    auc.ci.higher=confidence[3],
                    sens=sensitivity,
                    sens.ci.lower=sensitivity.ci[1],
                    sens.ci.higher=sensitivity.ci[2],
                    spec=specificity,
                    spec.ci.lower=specificity.ci[1],
                    spec.ci.higher=specificity.ci[2],
                    classifier.threshold=classifier.threshold,
                    ci.perc=0.95)
        } else {
          line <- c(selectionMethod=selectionMethod,
                    threshold=threshold,
                    classifierMethod=classifierMethod,
                    auc=NA,
                    auc.ci.lower=NA,
                    auc.ci.higher=NA,
                    sens=NA,
                    sens.ci.lower=NA,
                    sens.ci.higher=NA,
                    spec=NA,
                    spec.ci.lower=NA,
                    spec.ci.higher=NA,
                    classifier.threshold=NA,
                    ci.perc=NA)
        }
        out$summary <- rbind.data.frame(out$summary, line, stringsAsFactors=FALSE, make.row.names=FALSE)
      }
    }
  }
  colnames(out$summary) <- c("selectionMethod",
                             "threshold",
                             "classifierMethod",
                             "auc",
                             "auc.ci.lower",
                             "auc.ci.higher",
                             "sens",
                             "sens.ci.lower",
                             "sens.ci.higher",
                             "spec",
                             "spec.ci.lower",
                             "spec.ci.higher",
                             "classifier.threshold",
                             "ci.perc")
  return(out)
}

#' Calculate the mean ROC curve from a list of predictors and targets
#'
#' @param predictions_list list of predictors
#' @param targets_list list of targets
#' @param points points to calculate the ROC curves at
#' @param fixed_axis whether \code{points} fixes sensitivity or specificity
#' @param boot_index list of bootstrap indices for each fold. Generated automatically if NULL
#' @param return_folds if TRUE, return the ROC curve and AUC for each fold as well as the mean
#'
#' @return list containing mean ROC curve and AUC
#' @import pROC
#'
calc_mean_roc <- function(predictions_list,
                          targets_list,
                          points,
                          fixed_axis=c("specificity", "sensitivity"),
                          boot_index=NULL,
                          return_folds=FALSE) {
  num_folds <- length(predictions_list)
  fixed_axis_index <- which(c("specificity", "sensitivity") == fixed_axis) + 1

  if (is.null(boot_index)) {
    boot_index <- list()
    for (fold in 1:num_folds) {
      boot_index[[fold]] <- seq_along(predictions_list[[fold]])
    }
  }

  fold_roc <- matrix(0, nrow=num_folds, ncol=length(points))
  fold_auc <- numeric(num_folds)

  for (fold in 1:num_folds) {
    temp_roc <- roc(targets_list[[fold]][boot_index[[fold]]], predictions_list[[fold]][boot_index[[fold]]])

    fold_roc[fold, ] <- coords(temp_roc, points, fixed_axis)[-c(1, fixed_axis_index), ]

    fold_auc[fold] <- auc(temp_roc)
  }

  out_roc <- apply(fold_roc, 2, mean)

  out_auc <- mean(fold_auc)

  out <- list(roc=out_roc, auc=out_auc)

  if (return_folds) {
    out$folds_roc <- fold_roc
    out$folds_auc <- fold_auc
  }

  return(out)
}

#' Compute bootstraps of the mean ROC
#'
#' Compute the ROC curve separately for each fold, and compute the mean
#' ROC curve, instead of pooling predictions.
#'
#' @param predictions
#' @param folds
#' @param targets
#' @param num_samples
#' @param alt_predictions
#' @param fixed_axis
#' @param points
#'
#' @return
boot_auc_folds <- function(predictions,
                           folds,
                           targets,
                           num_samples=1000,
                           alt_predictions=NULL,
                           fixed_axis=c("sensitivity", "specificity"),
                           points=seq(0, 1, by=0.01),
                           stratified=T) {
  fixed_axis <- match.arg(fixed_axis)
  fixed_axis_index <- which(c("specificity", "sensitivity") == fixed_axis) + 1
  predictions_list <- list()
  targets_list <- list()
  if (!is.null(alt_predictions)) alt_predictions_list <- list()

  num_folds <- max(folds)
  if (any(seq(num_folds) != sort(unique(folds)))) {
    stop("folds is missing a fold number")
  }
  for (i in 1:num_folds) {
    predictions_list[[i]] <- predictions[folds==i]
    targets_list[[i]] <- targets[folds==i]
    if (!is.null(alt_predictions)) alt_predictions_list[[i]] <- alt_predictions[folds==i]
  }

  booted_auc <- numeric(num_samples)
  fold_auc <- numeric(num_folds)

  if (!is.null(alt_predictions)) {
    alt_booted_auc <- numeric(num_samples)
    alt_fold_auc <- numeric(num_folds)
  }

  booted_roc <- matrix(nrow=num_samples, ncol=length(points))
  fold_roc <- matrix(0, nrow=num_folds, ncol=length(points))

  if (!is.null(alt_predictions)) {
    alt_booted_roc <- matrix(nrow=num_samples, ncol=length(points))
    alt_fold_roc <- matrix(0, nrow=num_folds, ncol=length(points))
  }

  boot_index <- list()
  for (iteration in 1:num_samples) {
    for (fold in 1:num_folds) {
      if (stratified) {
        boot_index[[fold]] <- c(sample(which(targets_list[[fold]]), replace=T),
                                sample(which(!targets_list[[fold]]), replace=T))
      } else {
        valid <- FALSE
        while (!valid) {
          boot_index[[fold]] <- sample(seq_along(predictions_list[[fold]]), replace = T)
          if (any(targets_list[[fold]][boot_index[[fold]]]) & !all(targets_list[[fold]][boot_index[[fold]]])) {
            valid <- TRUE
          }
        }
      }
    }

    it_roc <- calc_mean_roc(predictions_list=predictions_list,
                            targets_list=targets_list,
                            points=points,
                            fixed_axis=fixed_axis,
                            boot_index=boot_index)

    booted_roc[iteration, ] <- it_roc$roc
    booted_auc[iteration] <- it_roc$auc

    if (!is.null(alt_predictions)) {
      alt_it_roc <- calc_mean_roc(predictions_list=alt_predictions_list,
                                  targets_list=targets_list,
                                  points=points,
                                  fixed_axis=fixed_axis,
                                  boot_index=boot_index)

      alt_booted_roc[iteration, ] <- alt_it_roc$roc
      alt_booted_auc[iteration] <- alt_it_roc$auc
    }


  }

  ROC <- calc_mean_roc(predictions_list=predictions_list,
                       targets_list=targets_list,
                       points=points,
                       fixed_axis=fixed_axis,
                       boot_index=NULL,
                       return_folds=TRUE)

  out <- list(AUC=ROC$auc, ROC=ROC$roc, AUC_boot=booted_auc, ROC_boot=booted_roc,
              fixed_axis=fixed_axis, points=points,
              folds_ROC=ROC$folds_roc, folds_AUC=ROC$folds_auc)

  if (!is.null(alt_predictions)) {
    alt_ROC <- calc_mean_roc(predictions_list=alt_predictions_list,
                             targets_list=targets_list,
                             points=points,
                             fixed_axis=fixed_axis,
                             boot_index=NULL)
    out$alt_AUC <- alt_ROC$auc
    out$alt_ROC <- alt_ROC$roc
    out$alt_AUC_boot <- alt_booted_auc
    out$alt_ROC_boot <- alt_booted_roc
  }

  return(out)
}

plot_booted_roc_curve <- function(roc_list, ci=0.95, ...) {
  if (roc_list$fixed_axis == "sensitivity") {
    sens <- c(0, roc_list$points, 1)
    spec <- c(1, roc_list$ROC, 0)
  } else {
    spec <- c(0, roc_list$points, 1)
    sens <- c(1, roc_list$ROC, 0)
  }
  upper <- apply(roc_list$ROC_boot, 2, function(x) quantile(x, (1 - ci) / 2))
  lower <- apply(roc_list$ROC_boot, 2, function(x) quantile(x, (1 + ci) / 2))
  plot(spec, sens, type="l",
       xlab="Specificity", ylab="Sensitivity",
       xlim=c(1,0), ylim=c(0,1),
       ...)
  if (roc_list$fixed_axis == "sensitivity") {
    upper.sens <- c(0, roc_list$points, 1)
    upper.spec <- c(1, upper, 0)
    lower.sens <- c(0, roc_list$points, 1)
    lower.spec <- c(1, lower, 0)
  } else {
    upper.spec <- c(0, roc_list$points, 1)
    upper.sens <- c(1, upper, 0)
    lower.sens <- c(0, roc_list$points, 1)
    lower.spec <- c(1, lower, 0)
  }
  polygon(c(upper.spec, rev(lower.spec)), c(upper.sens, rev(lower.sens)), col="grey")
  lines(spec, sens)

}

#' Compute a bootstrap of the difference between two ROC curve AUCs
#'
#'
#'
#' @param roc_1 first ROC curve (from pROC)
#' @param roc_2 second ROC curve (from pROC)
#' @param rep number of bootstrap samples
#'
#' @return bootstrap of the difference in AUC
#' @export
#'
#' @import pROC
compare_auc_bootstrap <- function(roc_1, roc_2, rep=10000) {
  auc_diff_bootstrap <- numeric(rep)
  pred_1 <- roc_1$predictor
  pred_2 <- roc_2$predictor
  resp <- roc_1$response
  if (any(roc_1$response != roc_2$response) | length(roc_1$response) != length(roc_2$response)) {
    stop("Responses for the two ROC curves differ")
  }
  for (i in 1:rep) {
    index <- sample(length(resp), replace=T)
    temp_auc_1 <- auc(roc(resp[index], pred_1[index]))
    temp_auc_2 <- auc(roc(resp[index], pred_2[index]))
    auc_diff_bootstrap[i] <- temp_auc_1 - temp_auc_2
  }
  return(auc_diff_bootstrap)
}
