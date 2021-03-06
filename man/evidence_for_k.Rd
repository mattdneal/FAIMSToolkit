% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA_k_selection.R
\name{evidence_for_k}
\alias{evidence_for_k}
\title{Calculate the Evidence for k}
\usage{
evidence_for_k(prcompObj, k)
}
\arguments{
\item{prcompObj}{a prcomp object}

\item{k}{k to assess}
}
\value{
log evidence for retaining k PCs
}
\description{
Calculate the evidence for retaining k PCs in a PCA analysis, as in Minka 2000, "Automatic Choice of dimensionality for PCA"
}
