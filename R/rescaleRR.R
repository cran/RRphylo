#' @title Rescaling phylogenetic trees according to RRphylo rates
#' @description The function rescales all branches and leaves of the
#'   phylogenetic tree according to branch-wise phenotypic evolutionary rates
#'   fitted by \code{\link{RRphylo}}.
#' @usage rescaleRR(tree,RR)
#' @param tree the phylogenetic tree to be rescaled.
#' @param RR the result of \code{RRphylo} performed on \code{tree}. When a
#'   multivariate phenotype is used, rescaling is operated on the norm-2 vector
#'   of rates
#' @export
#' @return Rescaled phylogenetic tree.
#' @author Silvia Castiglione, Pasquale Raia, Carmela Serio
#' @references Castiglione, S., Serio, C., Piccolo, M., Mondanaro, A.,
#'   Melchionna, M., Di Febbraro, M., Sansalone, G., Wroe, S., & Raia, P.
#'   (2020). The influence of domestication, insularity and sociality on the
#'   tempo and mode of brain size evolution in mammals. \emph{Biological Journal
#'   of the Linnean Society},in press. doi:10.1093/biolinnean/blaa186
#' @examples
#' \dontrun{
#' library(ape)
#'
#' rtree(100)->tree
#' fastBM(tree)->y
#'
#' RRphylo::RRphylo(tree,y)->RR
#' rescaleRR(tree,RR)->treeRes
#' }


rescaleRR<-function(tree,RR){
  # require(ape)

  tree->tree1
  abs(RR$rates[,1])->rts
  sum(tree1$edge.length)->t1ele
  rts[-1]->rts
  names(rts)[Nnode(tree1):length(rts)]<-seq(1,Ntip(tree1))
  rts[match(tree1$edge[,2],names(rts))]->rts
  tree1$edge.length*rts->tree1$edge.length
  t1ele/sum(tree1$edge.length)*tree1$edge.length->tree1$edge.length
  return(tree1)
}
