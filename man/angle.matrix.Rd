% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angle.matrix.R
\name{angle.matrix}
\alias{angle.matrix}
\title{Ontogenetic shape vectors analysis}
\usage{
angle.matrix(RR,node,Y=NULL,select.axes=c("no","yes"),
  type=c("phenotypes","rates"),cova=NULL,clus=0.5)
}
\arguments{
\item{RR}{an object produced by \code{\link{RRphylo}}.}

\item{node}{the number identifying the most recent common ancestor to all the
species the user wants ontogenetic vectors be computed.}

\item{Y}{multivariate trait values at tips.}

\item{select.axes}{if \code{"yes"}, \code{Y} variables are individually
regressed against developmental stages and only significant variables are
retained to compute ontogenetic vectors. All variables are retained
otherwise.}

\item{type}{specifies weather to perform the analysis on phenotypic
(\code{"phenotypes"}) or rate (\code{"rates"}) vectors.}

\item{cova}{the covariate to be indicated if its effect on rate values must
be accounted for. Contrary to \code{\link{RRphylo}}, \code{cova} needs to be as
long as the number of tips in the tree. As the covariate only affects rates
computation, there is no covariate to provide when \code{type =
"phenotypes"}.}

\item{clus}{the proportion of clusters to be used in parallel computing. To
run the single-threaded version of \code{angle.matrix} set \code{clus} = 0.}
}
\value{
A list containing 4 objects:

\enumerate{\item \strong{$regression.matrix} a 'list' including
  'angles between species' and 'angles between species to MRCA' matrices for
  all possible combinations of species pairs from the two sides descending
  from the MRCA. For each matrix, corresponding biogenetic and paedomorphosis
  tests are reported. \item \strong{$angles.2.MRCA.and.vector.size} a
  'data.frame' including angles between the resultant vector of species and
  the MRCA and the size of the resultant vector computed from species to
  MRCA, per stage per species. \item \strong{$ontogenetic.vectors2MRCA} a
  'data.frame' including angle, size, and corresponding x and y components,
  of ontogenetic vectors computed between each species and the MRCA. For both
  angle and size, the p-value for the difference between species pairs is
  reported. \item \strong{$ontogenetic.vectors.to.1st.stage} a 'list'
  containing: \itemize{\item$matrices: for all possible combinations of
  species pairs from the two sides descending form the MRCA, the upper
  triangle of the matrix contains the angles between different ontogenetic
  stages for the first species. The same applies to the lower triangle, but
  for the second species. \item$vectors: for all possible combinations of
  species pairs from the two sides descending form the MRCA, angles and sizes
  of ontogenetic vectors computed to the first stage of each species. For
  both, the p-value for the difference between the species pair is reported.
  } }
}
\description{
This function computes and compares ontogenetic vectors among
  species in a tree.
}
\details{
The \code{angle.matrix} function takes as objects a phylogenetic
  tree (retrieved directly from an \code{\link{RRphylo}} object), including
  the different ontogenetic stages of each species as polytomies. Names at
  tips must be written as species ID and stage number separated by the
  underscore. The \code{RR} object \code{angle.matrix} is fed with is
  just used to extract the dichotomized version of the phylogeny. This is
  necessary because node numbers change randomly at dichotomizing non-binary
  trees. However, when performing \code{angle.matrix} with the covariate the
  \code{RR} object must be produced without accounting for the
  covariate. Furthermore, as the covariate only affects the rates
  computation, it makes no sense to use it when computing vectors for
  phenotypic variables. Once angles and vectors are computed,
  \code{angle.matrix} performs two tests by means of standard major axis
  (SMA) regression. For each species pair, the "biogenetic test" verifies
  whether the angle between species grows during development, meaning that
  the two species becomes less similar to each other during growth. The
  "paedomorphosis test" tells whether there is heterochronic shape change in
  the data. Under paedomorphosis, the adult stages of one (paedomorphic)
  species will resemble the juvenile stages of the other (peramorphic)
  species. The test regresses the angles formed by the shapes at different
  ontogenetic stages of a species to the shape at the youngest stage of the
  other in the pair, against age. Then, it tests whether the two regression
  lines (one per species) have different slopes, and whether they have
  different signs. If the regression lines point to different directions, it
  means that one of the two species in the pair resembles, with age, the
  juveniles of the other, indicating paedomorphosis. Ontogenetic vectors of
  individual species are further computed, in reference to the MRCA of the
  pair, and to the first stage of each species (i.e. intraspecifically).
  Importantly, the size of the ontogenetic vectors of rates tell whether the
  two species differ in terms of developmental rate, which is crucial to
  understand which process is behind paedomorphosis, where it applies.While
  performing the analysis, the function prints messages on-screen informing
  about tests results. If \code{select.axes = "yes"}, informs the user about
  which phenotypic variables are used. Secondly, it specifies whether
  ontogenetic vectors to MRCA, and intraspecific ontogenetic vectors
  significantly differ in angle or size between species pairs. Then, for each
  species pair, it indicates if the biogenetic law and paedomorphosis apply.
}
\examples{
  \dontrun{
  data("DataApes")
  DataApes$PCstage->PCstage
  DataApes$Tstage->Tstage
  DataApes$CentroidSize->CS
  cc<- 2/parallel::detectCores()

  RRphylo(tree=Tstage,y=PCstage,clus=cc)->RRstage
# Case 1. without accounting for the effect of a covariate

 # Case 1.1 selecting shape variables that show significant relationship with age
  # on phenotypic vectors
    angle.matrix(RRstage,node=72,Y=PCstage,select.axes="yes",type="phenotypes",clus=cc)->am1
  # on rates vectors
    angle.matrix(RRstage,node=72,Y=PCstage,select.axes="yes",type="rates",clus=cc)->am2

 # Case 1.2 using all shape variables
  # on phenotypic vectors
    angle.matrix(RRstage,node=72,Y=PCstage,select.axes="no",type="phenotypes",clus=cc)->am3
  # on rates vectors
    angle.matrix(RRstage,node=72,Y=PCstage,select.axes="no",type="rates",clus=cc)->am4


# Case 2. accounting for the effect of a covariate (on rates vectors only)

 # Case 2.1 selecting shape variables that show significant relationship with age
   angle.matrix(RRstage,node=72,Y=PCstage,select.axes="yes",type="rates", cova=CS,clus=cc)->am5


 # Case 2.2 using all shape variables
   angle.matrix(RRstage,node=72,Y=PCstage,select.axes="no",type="rates",cova=CS,clus=cc)->am6
  }
}
\author{
Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
  Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
  Carotenuto
}
