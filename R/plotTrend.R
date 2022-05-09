#'@title Graphical representation of search.trend results
#'@description This function generates customized functions to produce plots of
#'  phenotype versus time and absolute evolutionary rates versus time
#'  regressions for the entire phylogeny and individual clades. Each custom
#'  function takes as the only argument the index or name of the variable (as
#'  in the \code{search.trend} output in \code{$trends.data$phenotypeVStime}) to be plotted.
#'@usage plotTrend(ST)
#'@param ST an object produced by \code{\link{search.trend}}.
#'@return The function returns a list of functions:
#'@return \strong{$plotPhen} returns the plot of rescaled phenotype versus age
#'  regression. The 95\% confidence intervals of slopes produced by regressing
#'  phenotypes simulated under the Brownian motion are plotted as a shaded area.
#'@return \strong{$plotRates} returns the plot of log rescaled rates versus age
#'  regression. The 95\% confidence intervals of slopes produced by regressing
#'  rates simulated under the Brownian motion are plotted as a shaded area.
#'@return \strong{$plotPhenNode} returns plots of rescaled phenotype versus age
#'  regression for individual clades. For each plot, the gray line represents
#'  the slope of phenotype versus age regression for the rest of the tree.
#'@return \strong{$plotRatesNode} returns plots of absolute rates versus age
#'  regression for individual clades. For each plot, the gray line represents
#'  the slope of absolute rates versus age regression for the rest of the tree.
#'@author Silvia Castiglione, Carmela Serio, Pasquale Raia
#'@importFrom graphics points text title polygon pairs plot
#'@export
#'@seealso \href{../doc/search.trend.html}{\code{search.trend} vignette}
#'@references Castiglione, S., Serio, C., Mondanaro, A., Di Febbraro, M.,
#'  Profico, A., Girardi, G., & Raia, P. (2019) Simultaneous detection of
#'  macroevolutionary patterns in phenotypic means and rate of change with and
#'  within phylogenetic trees including extinct species. \emph{PLoS ONE}, 14:
#'  e0210101. https://doi.org/10.1371/journal.pone.0210101
#' @examples
#'  \dontrun{
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#' DataOrnithodirans$massdino->massdino
#' cc<- 2/parallel::detectCores()
#'
#' # Extract Pterosaurs tree and data
#' library(ape)
#' extract.clade(treedino,746)->treeptero
#' massdino[match(treeptero$tip.label,names(massdino))]->massptero
#' massptero[match(treeptero$tip.label,names(massptero))]->massptero
#'
#' RRphylo(tree=treeptero,y=log(massptero))->RRptero
#'
#' search.trend(RR=RRptero, y=log(massptero), nsim=100, node=143, clus=cc,cov=NULL)->ST
#'
#' plotTrends(ST)->plotST
#'
#' plotST$plotPhen(1) # to plot phenotypic trend through time for entire tree
#' plotST$plotPhen("y") # equivalent to the previuous line
#'
#' plotST$plotRates(1) # to plot rates trend through time for entire tree
#'
#' plotST$plotPhenNode("y") # to plot phenotypic trend through time for the clade
#' plotST$plotRatesNode("y") # to plot rates trend through time for the clade
#'
#'    }

plotTrend<-function(ST){
  ST$trend.data$phenotypeVStime->phen.plot
  ST$trend.data$absrateVStime->absrate.plot
  ST$trend.data$rescaledrateVStime->resrate.plot
  ST$ConfInts$phenotype->CIphen
  ST$ConfInts$rescaled_rate->CIresrate
  ST$phenotypic.regression->p.phen
  if(any(grepl(".multi",colnames(phen.plot)))){
    resrate.plot[,-grep(".multi",colnames(phen.plot)),drop=FALSE]->resrate.plot
    absrate.plot[,-grep(".multi",colnames(phen.plot)),drop=FALSE]->absrate.plot
    phen.plot[,-grep(".multi",colnames(phen.plot)),drop=FALSE]->phen.plot
  }
  ncol(phen.plot)-1->ny

  apply(phen.plot[,1:ny,drop=FALSE],2,range01)->phen.plot[,1:ny]
  abs(absrate.plot[,1:ny,drop=FALSE])->absrate.plot[,1:ny]
  phen.plot$age<-max(phen.plot$age)-phen.plot$age

  if(!is.null(absrate.plot$group)){
    phen.plot$group<-absrate.plot$group[match(rownames(phen.plot),rownames(absrate.plot))]
    if("others"%in%unique(absrate.plot$group))
      unique(absrate.plot$group)[-which(unique(absrate.plot$group)=="others")]->groups
    absrate.plot$age<-max(absrate.plot$age)-absrate.plot$age

    car::outlierTest(lm(absrate.plot[,1]~absrate.plot[,2]))->outT
    if(any(outT$p<=0.05))
      max(absrate.plot[-as.numeric(names(outT$p)),1])->maxy else
        max(absrate.plot[,1])->maxy

    cols <- suppressWarnings(RColorBrewer::brewer.pal(length(groups), "Set2"))

    plotPhenNode<-function(variable){
      if(is.character(variable)){
        if(any(is.na(match(variable,colnames(phen.plot))))) stop("variable do not match column names")
        match(variable,colnames(phen.plot))->variable
      }else if(any(variable>ny)) stop("variable are out of y bounds")

      if(length(groups)%in%c(1,2)) {
        mat<-matrix(c(1,length(groups)),ncol=1,nrow=2,byrow=TRUE)

      }else{
        plotind<-1:length(groups)
        if(length(groups)%%2==0) c(2,length(groups)/2)->nn else c(2,(length(groups)+1)/2)->nn
        if(is.integer(mean(nn))) rep(mean(nn),2)->nn
        if((nn[1]*nn[2])>length(groups)) c(plotind,rep(0,(nn[1]*nn[2])-length(groups)))->plotind
        matrix(plotind,nrow=nn[1],ncol=nn[2])->mat

      }
      layout(mat,heights =c(1,1))

      variable->i
      lapply(1:length(groups),function(j){
        phen.plot[which(phen.plot$group==groups[j]),c(i,ny+1)]->node.phen
        phen.plot[which(phen.plot$group!=groups[j]),c(i,ny+1)]->tree.phen

        plot(phen.plot[, c(ny+1, i)],
             xlab = "", ylab = "",cex.axis=0.8,mgp = c(1.2, 0.4, 0),
             xlim=c(max(phen.plot[,ny+1]),min(phen.plot[,ny+1])),
             main= paste("Phenotypic Trend for Variable",colnames(phen.plot)[i]))
        title(xlab = "age", ylab =paste("Clade",gsub("g","",groups[j])),line = 1.5)

        predict(lm(tree.phen[,1]~tree.phen[,2]))->pred.tree

        points(x=tree.phen$age, y=pred.tree,
               lwd=3, col = "gray70",type="l")

        points(node.phen[,2],node.phen[,1],
               pch = 21, col = "black",bg = cols[j],cex=1.4)
        predict(lm(node.phen[,1]~node.phen[,2]))->pred.node
        points(x=node.phen[,2], y=pred.node,
               lwd=5, col = "black",type="l")
        points(x=node.phen[,2], y=pred.node,
               lwd=4, col = cols[j],type="l")
      })
    }

    plotRatesNode<-function(variable){
      if(is.character(variable)){
        if(any(is.na(match(variable,colnames(phen.plot))))) stop("variable do not match column names")
        match(variable,colnames(phen.plot))->variable
      }else if(any(variable>ny)) stop("variable are out of y bounds")

      if(length(groups)%in%c(1,2)) {
        mat<-matrix(c(1,length(groups)),ncol=1,nrow=2,byrow=TRUE)

      }else{
        plotind<-1:length(groups)
        if(length(groups)%%2==0) c(2,length(groups)/2)->nn else c(2,(length(groups)+1)/2)->nn
        if(is.integer(mean(nn))) rep(mean(nn),2)->nn
        if((nn[1]*nn[2])>length(groups)) c(plotind,rep(0,(nn[1]*nn[2])-length(groups)))->plotind
        matrix(plotind,nrow=nn[1],ncol=nn[2])->mat

      }
      layout(mat)

      variable->i
      lapply(1:length(groups),function(j){
        absrate.plot[which(absrate.plot$group==groups[j]),c(i,ny+1)]->node.rates
        absrate.plot[which(absrate.plot$group!=groups[j]),c(i,ny+1)]->tree.rates

        plot(absrate.plot[, c(ny+1, i)],
             xlab = "", ylab = "",cex.axis=0.8,mgp = c(1.2, 0.4, 0),
             ylim=c(min(absrate.plot[,1]),maxy),
             xlim=c(max(absrate.plot[,ny+1]),min(absrate.plot[,ny+1])),
             main= paste("Absolute Rates for Variable",colnames(phen.plot)[i]))
        title(xlab = "age", ylab =paste("Clade",gsub("g","",groups[j])),line = 1.5)

        predict(lm(tree.rates[,1]~tree.rates[,2]))->pred.tree

        points(x=tree.rates[,2], y=pred.tree,
               lwd=3, col = "gray70",type="l")

        points(node.rates[,2],node.rates[,1],
               pch = 21, col = "black",bg = cols[j],cex=1.4)

        predict(lm(node.rates[,1]~node.rates[,2]))->pred.node
        points(x=node.rates[,2], y=pred.node,
               lwd=5, col = "black",type="l")
        points(x=node.rates[,2], y=pred.node,
               lwd=4, col = cols[j],type="l")
      })
    }
  }else plotPhenNode<-NULL


  plotPhen<-function(variable){
    if(is.character(variable)){
      if(any(is.na(match(variable,colnames(phen.plot))))) stop("variable do not match column names")
      match(variable,colnames(phen.plot))->variable
    }else if(any(variable>ny)) stop("variable are out of y bounds")

    variable->i
    ynam.phen<-paste("Rescaled",colnames(phen.plot)[i])
    lm(phen.plot[, i] ~ I(max(phen.plot$age)-phen.plot$age))->plot.reg
    quantile(CIphen[[i]],c(0.025,0.975))->ci.slope
    sapply(ci.slope,function(k) coef(plot.reg)[1]+k*c(max(phen.plot$age),0))->ci.points

    plot(phen.plot[, c(ny+1, i)],
         xlab = "", ylab = "",cex.axis=0.8,mgp = c(1.2, 0.4, 0),
         xlim=c(max(phen.plot[,ny+1]),min(phen.plot[,ny+1])))
    title(xlab = "age", ylab = ynam.phen,main="Phenotypic Trend Test" ,line = 1.5)
    polygon(c(c(0,max(phen.plot$age),c(max(phen.plot$age),0))),
            c(ci.points[,1], rev(ci.points[, 2])), col = rgb(0.5, 0.5,0.5, 0.4), border = NA)
    points(phen.plot[which(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))),ny+1],
           phen.plot[which(suppressWarnings(is.na(as.numeric(rownames(phen.plot))))),i],
           pch = 21, col = "black",
           bg = "red")
    abline(lm(phen.plot[, i] ~ phen.plot$age), lwd = 3,col = "blue")
  }

  plotRates<-function(variable){
    if(is.character(variable)){
      if(any(is.na(match(variable,colnames(phen.plot))))) stop("variable do not match column names")
      match(variable,colnames(phen.plot))->variable
    }else if(any(variable>ny)) stop("variable are out of y bounds")

    resrate.plot$age<-max(na.omit(resrate.plot$age))-resrate.plot$age

    variable->i

    # if(is.null(x1)){
    #   scalrat.age <- max(na.omit(scalrat.data[, ncol(scalrat.data)]))-scalrat.data[, ncol(scalrat.data)]
    #   age <- max(rate.data[, ncol(rate.data)])-rate.data[, ncol(rate.data)]
    #   names(scalrat.age)<-names(age)<-rownames(rate.data)
    # }else{
    #   scalrat.age <- max(scalrat.data[which(rownames(scalrat.data)!=(Ntip(t)+1)), ncol(scalrat.data)])-
    #     scalrat.data[which(rownames(scalrat.data)!=(Ntip(t)+1)), ncol(scalrat.data)]
    #   age <- max(rate.data[which(rownames(rate.data)!=(Ntip(t)+1)), ncol(rate.data)])-
    #     rate.data[which(rownames(rate.data)!=(Ntip(t)+1)), ncol(rate.data)]
    #   names(scalrat.age)<-names(age)<-rownames(rate.data)[which(rownames(rate.data)!=(Ntip(t)+1))]
    # }
    #
    # if(is.null(x1)){
    #   scalrat.bet <- scalrat.data[, i]
    #   bet <- rate.data[, i]
    #   names(bet)<-names(scalrat.bet)<-rownames(rate.data)
    # }else{
    #   bet<-rate.data[-which(rownames(rate.data)==(Ntip(t)+1)),i]
    #   scalrat.bet<-scalrat.data[-which(rownames(scalrat.data)==(Ntip(t)+1)),i]
    #   names(bet)<-names(scalrat.bet)<-rownames(rate.data)[-which(rownames(rate.data)==(Ntip(t)+i))]
    # }

    lm(resrate.plot[,i] ~ I(max(resrate.plot$age)-resrate.plot$age))->plot.reg
    quantile(CIresrate[[i]],c(0.025,0.975))->ci.slope.scal
    sapply(ci.slope.scal,function(k) coef(plot.reg)[1]+k*c(max(resrate.plot$age),0))->ci.points.scal

    ynam.scal<-paste("log rescaled rate for variable",colnames(phen.plot)[i])

    plot(resrate.plot[,i] ~ resrate.plot$age,
         main="Rescaled Rate Trend Test",
         xlab="",ylab="",
         xlim=c(max(na.omit(resrate.plot$age)),min(na.omit(resrate.plot$age))),
         cex.axis=0.8,mgp = c(1.2, 0.4, 0))
    title(xlab = "age", ylab = ynam.scal,line = 1.5,)
    polygon(c(c(0,max(resrate.plot$age),c(max(resrate.plot$age),0))),
            c(ci.points.scal[, 1], rev(ci.points.scal[, 2])),
            col = rgb(0.5, 0.5, 0.5,0.4), border = NA)
    points(resrate.plot[which(suppressWarnings(is.na(as.numeric(rownames(resrate.plot))))),ny+1],
           resrate.plot[which(suppressWarnings(is.na(as.numeric(rownames(resrate.plot))))),i],
           pch = 21, col = "black",bg = "red")
    abline(lm(resrate.plot[,i] ~ resrate.plot$age), lwd = 4, col = "blue")
  }

  list(plotPhen=plotPhen,plotRates=plotRates)->res
  if(!is.null(plotPhenNode)) c(res,list(plotPhenNode=plotPhenNode,plotRatesNode=plotRatesNode))->res
  return(res)
}
