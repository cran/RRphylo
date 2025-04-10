#' @title Locating shifts in phenotypic evolutionary rates
#' @usage search.shift(RR, status.type = c("clade", "sparse"),node = NULL, state
#'   = NULL, cov = NULL, nrep = 1000, f = NULL)
#' @description The function \code{search.shift} (\cite{Castiglione et al.
#'   2018}) tests whether individual clades or group of tips dispersed through
#'   the phylogeny evolve at different \code{\link{RRphylo}} rates as compared
#'   to the rest of the tree.
#' @param RR an object fitted by the function \code{\link{RRphylo}}.
#' @param status.type whether the \code{"clade"} or \code{"sparse"} condition
#'   must be tested.
#' @param node under the \code{"clade"} condition, the node/s (clades) to be
#'   tested for the rate shift. If \code{node} is left unspecified, the function
#'   performs under the 'auto-recognize' feature, meaning it will automatically
#'   test individual clades for deviation of their rates from the background
#'   rate of the rest of the tree (see details).
#' @param state the named vector of states for each tip, to be provided under
#'   the \code{"sparse"} condition.
#' @param cov the covariate vector to be indicated if its effect on rate values
#'   must be accounted for. Contrary to \code{\link{RRphylo}}, \code{cov} needs to be
#'   as long as the number of tips of the tree.
#' @param nrep the number of simulations to be performed for the rate shift
#'   test, by default \code{nrep} is set at 1000.
#' @param f the size of the smallest clade to be tested. By default, nodes
#'   subtending to one tenth of the tree tips are tested.
#' @importFrom graphics symbols mtext
#' @importFrom stats sd
#' @importFrom utils globalVariables
#' @importFrom grDevices	pdf	dev.off
#' @export
#' @seealso \href{../doc/search.shift.html}{\code{search.shift} vignette}
#' @seealso \code{\link{overfitSS}}; \href{../doc/overfit.html#overfitSS}{\code{overfitSS} vignette}
#' @seealso \code{\link{plotShift}}; \href{../doc/Plotting-tools.html#plotShift}{\code{plotShift} vignette}
#' @details Under the 'auto-recognize' mode, \code{search.shift} automatically
#' tests individual clades (ranging in size from one half of the tree down to
#' \code{f} tips) for deviation of their rates from the background rate of the
#' rest of the tree. An inclusive clade with significantly high rates is likely
#' to include descending clades with similarly significantly high rates. Hence,
#' under 'auto-recognize' \code{search.shift} scans clades individually and
#' selects only the node subtending to the highest difference in mean absolute
#' rates as compared to the rest of the tree. If the argument \code{node}
#' (\code{"clade"} condition) is provided, the function computes the difference
#' between mean rate values of each clade and the rest of the tree, and compares
#' it to a random distribution of differences generated by shuffling rates
#' across tree branches. Additionally, if more than one \code{node} is
#' indicated, the rate difference for one clade is additionally computed by
#' excluding the rate values of the others from the rate vector of the rest of
#' the tree. Also, all the clades are considered as to be under a common rate
#' regime and compared as a single group to the rest of the tree.
#' @return Under \code{"clade"} case without specifying nodes (i.e.
#'   'auto-recognize') a list including:
#' @return \strong{$all.clades} for each detected node, the data-frame includes
#'   the average rate difference (computed as the mean rate over all branches
#'   subtended by the node minus the average rate for the rest of the tree) and
#'   the probability that it do represent a real shift. Probabilities are
#'   contrasted to simulations shuffling the rates across the tree branches for
#'   a number of replicates specified by the argument \code{nrep}. Note that the
#'   p-values refer to the number of times the real average rates are larger (or
#'   smaller) than the rates averaged over the rest of the tree, divided by the
#'   number of simulations. Hence, large rates are significantly larger than the
#'   rest of the tree (at alpha = 0.05), when the probability is > 0.975; and
#'   small rates are significantly small for p < 0.025.
#' @return \strong{$single.clades} the same as with 'all.clades' but restricted
#'   to the largest/smallest rate values along a single lineage (i.e. nested
#'   clades with smaller rate shifts are excluded).
#' @return Under \code{"clade"} condition by specifying the \code{node}
#'   argument:
#' @return \strong{$all.clades.together} if more than one node is tested, this
#'   specifies the average rate difference and the significance of the rate
#'   shift, by considering all the specified nodes as evolving under a single
#'   rate. As with the 'auto-recognize' feature, large rates are significantly
#'   larger than the rest of the tree (at alpha = 0.05), when the probability is
#'   > 0.975; and small rates are significantly small for p < 0.025.
#' @return \strong{$single.clades} gives the significance for individual clades
#'   tested individually against the rest of the tree (\strong{$singles}) and by
#'   excluding the rate values of other shifting clades from the rate vector of the rest of
#'   the tree (\strong{$no.others})
#' @return Under the \code{"sparse"} condition:
#' @return   \strong{$state.results} for each state, the data-frame includes the
#'   average rate difference (computed as the mean rate over all leaves evolving
#'   under a given state, minus the average rate for each other state or the
#'   rest of the tree) and the probability that the shift is real. Large rates
#'   are significantly larger (at alpha = 0.05), when the probability is >
#'   0.975; and small rates are significantly small for p < 0.025. States are
#'   compared pairwise.
#' @return Under all circumstances, if \code{'cov'} values are provided to the
#'   function, \code{search.shift} returns as \strong{$rates} object the vector of
#'   residuals of \code{\link{RRphylo}} rates versus \code{cov} regression.
#' @return The output always has an attribute "Call" which returns an unevaluated call to the function.
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro
#'   Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco
#'   Carotenuto
#' @references Castiglione, S., Tesone, G., Piccolo, M., Melchionna, M.,
#'   Mondanaro, A., Serio, C., Di Febbraro, M., & Raia, P.(2018). A new method
#'   for testing evolutionary rate variation and shifts in phenotypic evolution.
#'   \emph{Methods in Ecology and Evolution}, 9:
#'   974-983.doi:10.1111/2041-210X.12954
#' @examples
#' \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' DataOrnithodirans$statedino->statedino
#' cc<- 2/parallel::detectCores()
#'
#' RRphylo(tree=treedino,y=massdino,clus=cc)->dinoRates
#'
#' # Case 1. Without accounting for the effect of a covariate
#'
#' # Case 1.1 "clade" condition
#' # with auto-recognize
#' search.shift(RR=dinoRates,status.type="clade")->SSauto
#' # testing two hypothetical clades
#' search.shift(RR=dinoRates,status.type="clade",node=c(696,746))->SSnode
#'
#' # Case 1.2 "sparse" condition
#' # testing the sparse condition.
#' search.shift(RR=dinoRates,status.type= "sparse",state=statedino)->SSstate
#'
#'
#' # Case 2. Accounting for the effect of a covariate
#'
#' # Case 2.1 "clade" condition
#' search.shift(RR=dinoRates,status.type= "clade",cov=massdino)->SSauto.cov
#'
#' # Case 2.2 "sparse" condition
#' search.shift(RR=dinoRates,status.type="sparse",state=statedino,cov=massdino)->SSstate.cov
#'     }



search.shift<-function(RR,
                       status.type=c("clade","sparse"),
                       node=NULL,
                       state=NULL,
                       cov=NULL,
                       nrep=1000,
                       f=NULL)
{
  # require(phytools)

  SScore<-function(leaves,rates){
    leaf.rates <- rates[match(leaves, rownames(rates),nomatch = 0),]
    NCrates <- rates[which(!rownames(rates)%in%names(leaf.rates))]
    leaf2NC.diff <- mean(abs(leaf.rates))-mean(abs(NCrates))
    C <- length(leaf.rates)
    NC <- length(rates) - C
    ran.diffR <- replicate(nrep,mean(sample(abs(rates), C)) -
                             mean(sample(abs(rates),NC)))
    p.shift <- rank(c(leaf2NC.diff, ran.diffR[-nrep]))[1]/nrep
    return(cbind(diff=leaf2NC.diff,p=p.shift))
  }

  funcall <- match.call()
  tree <- RR$tree
  rates <- RR$rates[,,drop=FALSE]
  betas<-RR$multiple.rates[,,drop=FALSE]


  if(is.null(f)) f<-round(Ntip(tree)/10)

  if(!is.null(cov)){
    RRphylo(tree,cov,clus=0)->RRcova
    abs(c(RRcova$aces,cov))->Y
    c(rownames(RRcova$aces),names(cov))->names(Y)
    covRates(Y,betas)->betas

    if(ncol(betas)>1) rates <- as.matrix(apply(betas, 1, function(x) sqrt(sum(x^2)))) else rates<-betas
  }

  if (status.type == "clade") {
    if (is.null(node)) {
      st <- subtrees(tree)
      len<-sapply(st,Ntip)
      st <- st[which(len < (Ntip(tree)/2) & len >= round(f))]
      nns <- sapply(st, function(x) getMRCA(tree, x$tip.label))
    } else node->nns


    allres<-do.call(rbind,lapply(nns,function(j) SScore(c(getDescendants(tree, j), tips(tree, j)),rates)))
    rownames(allres)<-nns
    colnames(allres)<-c("rate.difference","p.value")

    l2N.init<-allres[,1]
    p.init<-allres[,2]
    names(l2N.init)<-names(p.init)<-nns
    # data.frame(rate.difference=l2N.init[match(names(p.init),names(l2N.init))],p.value=p.init)->allres

    if(is.null(node)){
      if (length(p.init[p.init>=0.975|p.init<=0.025])==0) p.single <-leaf2N.diff <-NULL else{
        p.single <- p.init[p.init>=0.975|p.init<=0.025]
        leaf2N.diff <- l2N.init[match(names(p.single),names(l2N.init))]
      }
      if (length(p.single)>= 2){
        ups <- p.single[p.single >= 0.975]
        dws <- p.single[p.single <= 0.025]
        ups.sel<-sapply(node.paths(tree,names(ups)),function(x){
          x[which.max(abs(leaf2N.diff[match(x,names(leaf2N.diff))]))]
        })
        dws.sel<-sapply(node.paths(tree,names(dws)),function(x){
          x[which.max(abs(leaf2N.diff[match(x,names(leaf2N.diff))]))]
        })
        ups<-ups[which(names(ups)%in%ups.sel)]
        dws<-dws[which(names(dws)%in%dws.sel)]

        p.single <- p.single[which(names(p.single)%in%names(c(ups, dws)))]
        leaf2N.diff <- leaf2N.diff[match(names(p.single),names(leaf2N.diff))]

        p.single[order(p.single)]->p.single
      }

      if(is.null(p.single)) single<-NULL else
        data.frame(rate.difference=leaf2N.diff[match(names(p.single),names(leaf2N.diff))],
                   p.value=p.single)->single
      res<-list(all.clades=allres,single.clades=single)

    }else{
      if(length(nns)==1)
        res<-list(single.clades=allres) else{
          totshtips<-length(unlist(lapply(nns,function(x) tips(tree,x))))
          if(totshtips>Ntip(tree)*0.5){
            warning("The clades under testing include more than one half of the tree species")
            res<-list(single.clades=allres)
          }else{
            node.paths(tree,nns)->np
            if(any(sapply(np,length)>1)){
              nsp<-np[which(sapply(np,length)>1)]
              warning(paste("Nodes",paste(sapply(nsp,function(x) paste(x,collapse="-")),collapse=" and "),
                            "are on the same path,only one per pair will be tested"),immediate.=TRUE)
              nns<-c(unlist(np[which(sapply(np,length)==1)]),
                     sapply(nsp,function(x) x[which.max(abs(allres[match(x,rownames(allres)),1]))]))

            }

            Cbranch<-unlist(lapply(nns,function(k) getDescendants(tree,k)))
            Cbranch <- unique(Cbranch[which(Cbranch>=Ntip(tree))])
            Ctips<-unique(unlist(lapply(nns,function(k) tips(tree,k))))
            Cleaf <- c(Cbranch, Ctips)
            allclatog<-SScore(Cleaf,rates)
            rownames(allclatog)<-"all"
            colnames(allclatog)<-c("rate.difference","p.value")

            ssc.list<-list()
            for (i in 1:length(nns)) {
              NOD <- nns[-i]
              others <- unlist(lapply(NOD,function(k) tips(tree, k)))
              des<-unlist(lapply(NOD,function(k) getDescendants(tree, k)[which(getDescendants(tree, k)>Ntip(tree))]))
              otdes<- c(des, others)

              Cleaf <- c(getDescendants(tree, nns[i]),unlist(tips(tree, nns[i])))
              NCrates<-rates[which(!rownames(rates)%in%otdes),,drop=FALSE]
              ssc.list[[i]]<-SScore(Cleaf,NCrates)
            }
            noothers<-do.call(rbind,ssc.list)
            rownames(noothers)<-nns
            colnames(noothers)<-c("rate.difference","p.value")

            res<-list(all.clades.together=allclatog,single.clades=list(singles=allres,no.others=noothers))
          }
        }
    }
  } else {
    state<-as.matrix(state)
    state <- treedataMatch(tree, state)[[1]][,1]
    frame <- data.frame(status = as.factor(state),
                        rate = rates[match(names(state),rownames(rates))])
    sta <- tapply(abs(frame$rate), frame$status, mean)
    sta <- sta[match(unique(state), names(sta))]
    status.diff <- apply(combn(sta, 2), 2, diff)

    if(length(unique(state)) > 2){
      w <- sapply(1:length(sta),function(x) sta[x]-
                    mean(abs(frame[which(frame$status!=names(sta)[x]), 2])))
      status.diff <- c(status.diff, w)
      names(status.diff) <- c(apply(combn(names(sta),2), 2, function(x) paste(x[2], x[1], sep = "_")),
                              names(sta))
    } else names(status.diff)<-paste(combn(names(sta), 2)[2:1],collapse="_")


    status.diffS <- matrix(ncol = length(status.diff),
                           nrow = nrep)
    for (i in 1:nrep) {
      s.ran <- sample(frame$status)
      s.frame <- data.frame(s.ran, frame$rate)
      sta <- tapply(abs(frame$rate), s.ran, mean)
      sta <- sta[match(unique(state), names(sta))]
      SD <- apply(combn(sta, 2), 2, diff)
      if(length(unique(state)) > 2){
        w <- sapply(1:length(sta),function(x) sta[x]-
                      mean(abs(frame[which(s.frame$s.ran!=names(sta)[x]), 2])))
        status.diffS[i, ] <- c(SD, w)
      }else status.diffS[i, ] <- SD
    }
    colnames(status.diffS) <- names(status.diff)

    p.status.diff<-sapply(1:length(status.diff),function(i)
      rank(c(status.diff[i],status.diffS[-nrep, i]))[1]/nrep)
    names(p.status.diff) <- names(status.diff)

    if(any(colnames(status.diffS)%in%unique(state)))
      pldata<-status.diffS[,which(colnames(status.diffS)%in%unique(state))] else
      pldata<-status.diffS

    res<-list(state.results=data.frame(rate.difference=status.diff[match(names(p.status.diff),names(status.diff))],p.status.diff),
              plotData=pldata)
  }

  if(!is.null(cov)) {
    if(!is.null(res$plotData)){
      res<-c(res[which(names(res)!="plotData")],rates=list(rates),res[which(names(res)=="plotData")])
    } else res<-c(res,rates=list(rates))
  }

  class(res)<-c("RRphyloList","list")
  attr(res,"hidden")<-"plotData"
  attr(res,"Call")<-funcall

  return(res)
}

