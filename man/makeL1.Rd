% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeL1.R
\name{makeL1}
\alias{makeL1}
\title{Matrix of branch lengths along a root-to-node path}
\usage{
makeL1(tree)
}
\arguments{
\item{tree}{a phylogenetic tree. The tree needs not to be ultrametric and
fully dichotomous.}
}
\value{
The function returns a \eqn{n * n} matrix of branch lengths for all
  root-to-node paths (one per each node of the tree).
}
\description{
This function produces a \eqn{n * n} matrix, where n=number of
  internal branches. Each row represents the branch lengths aligned along a
  root-to-node path.
}
\examples{
data("DataApes")
DataApes$Tstage->Tstage

makeL1(tree=Tstage)->ml1
}
\author{
Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
  Carotenuto
}
