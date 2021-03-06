% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.merger.R
\name{tree.merger}
\alias{tree.merger}
\title{Fast addition of tips and clades on an existing tree}
\usage{
tree.merger(backbone,data,source.tree=NULL,tip.ages = NULL, node.ages
  = NULL,min.branch=NULL,plot=TRUE,filename=NULL)
}
\arguments{
\item{backbone}{the backbone tree to attach tips/clades on.}

\item{data}{a dataset including as columns:\enumerate{\item bind = the
tips/clades to be attached; \item reference = the reference tip or clade
where 'bind' must be attached; \item poly = logical specifying if 'bind'
and 'reference' should form a polytomous clade.} See details for further
explanations.}

\item{source.tree}{the tree where 'bind' clades are to be extracted from. If
no clade has to be attached, it can be left unspecified.}

\item{tip.ages}{as in \code{\link{scaleTree}}, a named vector including the
ages (i.e. the time distance from the youngest tip within the tree) of the
tips. If unspecified, the function assumes all the tips on the backbone
tree are correctly placed and places all the new tips at the maximum
distance from the tree root (i.e. the present if the tips are extant).}

\item{node.ages}{as in \code{scaleTree}, a named vector including the ages
(i.e. the time distance from the youngest tip within the tree) of the
nodes. The nodes must be defined by collating the names of the two
phylogenetically furthest tips it subtends to, separated by a "-" (see
examples). If no calibration date for nodes is supplied, the function may
shift the node position back in time as to place new tips/clades and to fit
tip ages.}

\item{min.branch}{as in \code{scaleTree}, the minimum branch length that will
be imposed for shifted nodes.}

\item{plot}{if \code{TRUE}, the function produces an interactive plotting
device to check the placing of each \code{bind}.}

\item{filename}{if \code{plot=TRUE} and provided a \code{filename} (with or
without the path), the function stores a pdf file showing the plot of the
entire phylogeny.}
}
\value{
Merged phylogenetic tree.
}
\description{
The function attaches new tips and/or clades derived from a
  source phylogeny to a pre-existing backbone tree.
}
\details{
The function attaches tips and/or clades from the \code{source} tree
  to the \code{backbone} tree according to the \code{data} object. Within the
  latter, a clade, either to be binded or to be the reference, must be
  indicated by collating the names of the two phylogenetically furthest tips
  belonging to it, separated by a "-". Duplicated 'bind' produce error.
  Tips/clades set to be attached to the same 'reference' are considered to
  represent a polytomy. Tips set as 'bind' which are already on the backbone
  tree are removed from the latter and placed according to the 'reference'.
  See examples and \href{../doc/Tree-Manipulation.html#tree.merger.html}{vignette} for
  clarifications.
}
\examples{
 \dontrun{
 require(ape)
 require(geiger)
 DataCetaceans$treecet->tree
 data.frame(bind=c("Balaena_mysticetus-Caperea_marginata",
                   "Aetiocetus_weltoni",
                   "Saghacetus_osiris",
                   "Zygorhiza_kochii",
                   "Ambulocetus_natans",
                   "Kentriodon_pernix",
                   "Kentriodon_schneideri",
                   "Kentriodon_obscurus",
                   "Tursiops_truncatus-Delphinus_delphis",
                   "Kogia_sima",
                   "Grampus_griseus"),
            reference=c("Fucaia_buelli-Aetiocetus_weltoni",
                        "Aetiocetus_cotylalveus",
                        "Fucaia_buelli-Tursiops_truncatus",
                        "Saghacetus_osiris-Fucaia_buelli",
                        "Dalanistes_ahmedi-Fucaia_buelli",
                        "Kentriodon_schneideri",
                        "Phocoena_phocoena-Delphinus_delphis",
                        "Kentriodon_schneideri",
                        "Stenella_attenuata-Stenella_longirostris",
                        "Kogia_breviceps",
                        "Globicephala_melas-Pseudorca_crassidens"),
            poly=c(FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE))->dato

 c(Aetiocetus_weltoni=28.0,
   Saghacetus_osiris=33.9,
   Zygorhiza_kochii=34.0,
   Ambulocetus_natans=40.4,
   Kentriodon_pernix=15.9,
   Kentriodon_schneideri=11.61,
   Kentriodon_obscurus=13.65)->tip.ages
 c("Ambulocetus_natans-Fucaia_buelli"=52.6,
   "Balaena_mysticetus-Caperea_marginata"=21.5)->node.ages

 # remove some tips from the original tree and create a source tree
 drop.tip(tree,c(names(tip.ages),
                 tips(tree,131)[-which(tips(tree,131)\%in\%
                               c("Caperea_marginata","Eubalaena_australis"))],
                 tips(tree,195)[-which(tips(tree,195)=="Tursiops_aduncus")]))->backtree
 drop.tip(tree,which(!tree$tip.label\%in\%c(names(tip.ages),
                                          tips(tree,131),
                                          tips(tree,195))))->sourcetree

 plot(backtree,cex=.6)
 plot(sourcetree,cex=.6)

 tree.merger(backbone=backtree,data=dato,source.tree=sourcetree,
             tip.ages=tip.ages,node.ages = node.ages, plot=TRUE)->treeM
   }
}
\references{
aaa
}
\seealso{
\href{../doc/Tree-Manipulation.html#tree.merger.html}{\code{tree.merger} vignette};
  \href{../doc/Tree-Manipulation.html#scaleTree}{\code{scaleTree} vignette};
}
\author{
Silvia Castiglione, Carmela Serio, Pasquale Raia
}
