% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA_k_selection.R
\name{select_k}
\alias{select_k}
\title{Select PCA PCs to retain}
\usage{
select_k(prcompObj)
}
\arguments{
\item{prcompObj}{a prcomp object}
}
\value{
k
}
\description{
Select the number of PCs k to retain based on the evidence for each model. As per Minka 2000, "Automatic Choice of dimensionality for PCA".
}
