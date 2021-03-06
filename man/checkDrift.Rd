% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{checkDrift}
\alias{checkDrift}
\title{Check for FAIMS drift}
\usage{
checkDrift(faims, index = NULL, ...)
}
\arguments{
\item{faims}{a faims object}

\item{index}{the samples to include. NULL (the default) includes all samples}

\item{...}{other parameters to pass to plot}
}
\value{
a list with the correlation, a 95% confidence interval, and the matched rankings.
}
\description{
Check for FAIMS drift
}
