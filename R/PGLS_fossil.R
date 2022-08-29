#' @title Phylogenetic Generalized Least Square with phylogenies including
#'   fossils
#' @description The function performs pgls for non-ultrametric trees using a
#'   variety of evolutionary models or \code{\link{RRphylo}} rates to change the
#'   tree correlation structure.
#' @usage PGLS_fossil(modform,data,tree=NULL,RR=NULL,...)
#' @param modform the formula for the regression model.
#' @param data a list of named vectors including response and predictor
#'   variables as named in \code{modform}.
#' @param tree a phylogenetic tree to be indicated for any model other than
#'   \code{RRphylo} is used to rescale tree branches. The tree needs not to be
#'   ultrametric and fully dichotomous.
#' @param RR the result of \code{RRphylo} performed on the response variable. If
#'   \code{RR} is specified, tree branches are rescaled to the absolute
#'   branch-wise rate values calculated by \code{RRphylo} to transform the
#'   variance-covariance matrix.
#' @param ... further argument passed to the function
#'   \code{\link[phylolm]{phylolm}}.
#' @importFrom ape corPagel
#' @importFrom stats terms
#' @export
#' @seealso \href{../doc/RRphylo.html}{\code{RRphylo} vignette}
#' @return Fitted pgls parameters and significance.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' \dontrun{
#' library(ape)
#' library(phytools)
#'
#' rtree(100)->tree
#' fastBM(tree)->resp
#' fastBM(tree,nsim=3)->resp.multi
#' fastBM(tree)->pred1
#' fastBM(tree)->pred2
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),tree=tree)->pgls_noRR
#'
#' RRphylo::RRphylo(tree,resp)->RR
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp,x2=pred1,x1=pred2),tree=tree,RR=RR)->pgls_RR
#'
#' # To derive log-likelihood and AIC for PGLS_fossil outputs the function AIC can be applied
#' AIC(pgls_noRR)
#' AIC(pgls_RR)
#'
#'
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree)->pgls2_noRR
#' cc<- 2/parallel::detectCores()
#' RRphylo::RRphylo(tree,resp.multi,clus=cc)->RR
#' PGLS_fossil(modform=y1~x1+x2,data=list(y1=resp.multi,x2=pred1,x1=pred2),tree=tree,RR=RR)->pgls2_RR
#'
#' AIC(pgls2_noRR)
#' AIC(pgls2_RR)
#' }

PGLS_fossil<-function(modform,data,tree=NULL,RR=NULL,...){
  # require(nlme)
  # require(ape)
  # require(phylolm)
  if (!requireNamespace("phylolm", quietly = TRUE)) {
    stop("Package \"phylolm\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(!is.null(RR)) rescaleRR(RR$tree,RR=RR)->tree
  lapply(data,function(k) treedataMatch(tree, k)[[1]])->data

  argum<-list(...)
  if(!is.null(RR)&length(argum)>0){
    warning("The tree is rescaled by using RR rates, no further arguments are allowed",immediate. = TRUE)
    phylolm::phylolm(modform,phy=tree,data=data,model="BM")->res
  } else phylolm::phylolm(modform,phy=tree,data=data,...)->res

  return(res)
}
