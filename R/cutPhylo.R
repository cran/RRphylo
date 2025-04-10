#' @title Cut the phylogeny at a given age or node
#' @description The function cuts all the branches of the phylogeny which are
#'   younger than a specific age or node (i.e. the age of the node).
#' @usage cutPhylo(tree,age=NULL,node=NULL,keep.lineage=TRUE)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @param age the age (in terms of time distance from the recent) at which the
#'   tree must be cut
#' @param node the node whose age must be used as cutting limit.
#' @param keep.lineage logical specifying whether lineages with no descendant
#'   tip must be retained (see example below). Default is \code{TRUE}.
#' @export
#' @seealso \href{../doc/Tree-Manipulation.html#cutPhylo}{\code{cutPhylo}
#'   vignette}
#' @importFrom phytools drop.clade
#' @importFrom ape axisPhylo
#' @details When an entire lineage is cut (i.e. one or more nodes along a path)
#'   and \code{keep.lineages = TRUE}, the leaves left are labeled as "l"
#'   followed by a number.
#' @return The function returns the cut phylogeny and plots it into the graphic
#'   device. The time axis keeps the root age of the original tree. Note, tip
#'   labels are ordered according to their position in the tree.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' \dontrun{
#' library(ape)
#'
#' set.seed(22)
#' rtree(100)->tree
#' 3->age
#'
#' cutPhylo(tree,age=age)->t1
#' cutPhylo(tree,age=age,keep.lineage=FALSE)->t1a
#' cutPhylo(tree,node=151)->t2
#' cutPhylo(tree,node=151,keep.lineage=FALSE)->t2a
#' }

cutPhylo<-function(tree,age=NULL,node=NULL,keep.lineage=TRUE){
  # require(ape)
  # require(phytools)

  if(!identical(tree$edge[tree$edge[,2]<=Ntip(tree),2],seq(1,Ntip(tree)))){
    tree$tip.label<-tree$tip.label[tree$edge[tree$edge[,2]<=Ntip(tree),2]]
    tree$edge[tree$edge[,2]<=Ntip(tree),2]<-seq(1,Ntip(tree))
  }


  dist.nodes(tree)[(Ntip(tree)+1),]->dd
  names(dd)[which(as.numeric(names(dd))<=Ntip(tree))]<-tree$tip.label[as.numeric(names(dd)[which(as.numeric(names(dd))<=Ntip(tree))])]
  if(is.null(node)) max(nodeHeights(tree))-age->cutT else dd[match(node,names(dd))]->cutT

  dd[which(dd>=cutT)]->ddcut
  names(ddcut)->cutter

  ### Tips only ###
  tree->tt
  if(all(cutter%in%tree$tip.label)){
    tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]<-
      tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]-(ddcut-cutT)
  }else{
    ### Tips and nodes ###
    #if(any(suppressWarnings(as.numeric(cutter))>Ntip(tree))){
    as.numeric(cutter[which(suppressWarnings(as.numeric(cutter))>Ntip(tree))])->cutn
    i=1
    while(i<=length(cutn)){
      if(any(cutn%in%getDescendants(tree,cutn[i]))) cutn[-which(cutn%in%getDescendants(tree,cutn[i]))]->cutn
      i=i+1
    }

    tree->tt
    i=1
    while(i<=length(cutn)){
      getMRCA(tt,tips(tree,cutn[i]))->nn
      drop.clade(tt,tips(tt,nn))->tt
      tt$tip.label[which(tt$tip.label=="NA")]<-paste("l",i,sep="")
      i=i+1
    }

    diag(vcv(tt))->times
    if(any(times>cutT)){
      times[which(times>cutT)]->times
      tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]<-
        tt$edge.length[match(match(names(times),tt$tip.label),tt$edge[,2])]-(times-cutT)
    }

  }

  if(isFALSE(keep.lineage)){
    drop.tip(tt,which(!tt$tip.label%in%tree$tip.label))->tt
    if(Ntip(tt)<=100) plot(tt,cex=.8) else plot(tt,cex=.5 )
  }else{
    tip.col<-rep("black",Ntip(tt))
    tip.col[which(!tt$tip.label%in%tree$tip.label)]<-"red"
    if(Ntip(tt)<=100) plot(tt,cex=.8,tip.color=tip.col) else plot(tt,cex=.5,tip.color=tip.col)
  }

  axisPhylo(root.time = max(nodeHeights(tree)))

  return(tt)
}
