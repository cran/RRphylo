% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geiger_functions.R
\name{sig2BM}
\alias{sig2BM}
\title{Brownian Motion rate computation}
\usage{
sig2BM(tree,y)
}
\arguments{
\item{tree}{a phylogenetic tree. The tree needs not to be ultrametric and
fully dichotomous.}

\item{y}{either a single vector variable or a multivariate dataset. In any
case, \code{y} must be named.}
}
\value{
The Brownian Motion rate of phenotypic evolution for each variable in \code{y}.
}
\description{
The function computes rate of phenotypic evolution along a phylogeny assuming Brownian Motion model of evolution.
}
\examples{

### Univariate data ###
data(DataCetaceans)
DataCetaceans$treecet->treecet
DataCetaceans$masscet->masscet
sig2BM(tree=treecet,y=masscet)

### Multivariate data ###
data(DataUng)
DataUng$treeung->treeung
DataUng$PCscoresung->PCscores
sig2BM(tree=treeung,y=PCscores)
}
\author{
Pasquale Raia, Silvia Castiglione
}
