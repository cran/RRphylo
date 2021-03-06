#' @title Matrix of branch lengths along a root-to-node path
#' @description This function produces a \eqn{n * n} matrix, where n=number of
#'   internal branches. Each row represents the branch lengths aligned along a
#'   root-to-node path.
#' @usage makeL1(tree)
#' @param tree a phylogenetic tree. The tree needs not to be ultrametric and
#'   fully dichotomous.
#' @export
#' @return The function returns a \eqn{n * n} matrix of branch lengths for all
#'   root-to-node paths (one per each node of the tree).
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @examples
#' data("DataApes")
#' DataApes$Tstage->Tstage
#'
#' makeL1(tree=Tstage)


makeL1<-function (tree)
{

  if(!identical(tree$tip.label,tips(tree,(Ntip(tree)+1)))){
    data.frame(tree$tip.label,N=seq(1,Ntip(tree)))->dftips
    tree$tip.label<-tips(tree,(Ntip(tree)+1))
    data.frame(dftips,Nor=match(dftips[,1],tree$tip.label))->dftips
    tree$edge[match(dftips[,2],tree$edge[,2]),2]<-dftips[,3]
  }
  t <- tree
  internals <- unique(c(t$edge[, 1], t$edge[, 2][which(t$edge[,
                                                              2] > Ntip(t))]))
  sort(internals)->internals
  edgedX <- data.frame(t$edge, t$edge.length)
  edged.1 <- edgedX[edgedX$X2 > Ntip(t), ]
  L1 <- matrix(ncol = length(internals), nrow = length(internals))
  colnames(L1) <- internals
  rownames(L1) <- internals
  node.path <- list()
  for (i in 1:length(internals)) {
    a <- getDescendants(t, internals[i])
    a <- a[a > Ntip(t)]
    node.path[[i]] <- data.frame(rep(internals[i], length(a)),
                                 a)
  }
  node.path <- do.call(rbind, node.path)
  pathN <- list()
  for (i in 1:length(edged.1[, 2])) {
    pathN[[i]] <- c(node.path[which(node.path[, 2] == edged.1[i,
                                                              2]), 1], edged.1[i, 2])
  }
  for (j in 1:length(pathN)) {
    a <- list()
    for (i in 2:length(pathN[[j]]) - 1) {
      a[[i]] <- pathN[[j]][c(i, i + 1)]
    }
    b <- do.call(rbind, a)
    L1.match <- b[, 2]
    br.len <- array()
    for (k in 1:dim(b)[1]) {
      br.len[k] <- edged.1[b[k, 1] == edged.1[, 1] & b[k,
                                                       2] == edged.1[, 2], ][, 3]
    }
    d <- data.frame(L1.match, br.len)
    L1[j+1, match(d[, 1], colnames(L1))] <- d[, 2]
  }
  if (is.null(t$root.edge) || t$root.edge == 0)
    L1[, 1] <- 1
  else L1[, 1] <- t$root.edge
  L1[which(is.na(L1))] <- 0
  return(L1)
}
