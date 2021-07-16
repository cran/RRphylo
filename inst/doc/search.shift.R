## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),out.width="70%",dpi=220----------
require(ape)
require(phytools)
require(geiger)
require(scales)

set.seed(22)
rtree(100)->Tree
fastBM(Tree)->y

y[tips(Tree,161)]*5->y[tips(Tree,161)]

RRphylo(Tree,y)->RR

RR$tree->tree
RR$rates->rates
round(Ntip(tree)/10)->f
nrep=1000

ST <- subtrees(tree)
len <- array()
for (i in 1:length(ST)) {
  len[i] <- Ntip(ST[[i]])
}
st <- ST[which(len < (Ntip(tree)/2) & len > round(f))]
node <- sapply(st, function(x) getMRCA(tree, x$tip.label))
names(st) <- node

plot(Tree,no.margin = TRUE, show.tip.label = FALSE)
nodelabels(bg="white",frame="none",node=node,col="red")


leaf2N.diff <- array()
p.single <- array()
for (j in 1:length(node)) {
  Cbranch <- getDescendants(tree, node[j])
  Ctips <- tips(tree, node[j])
  Cleaf <- c(Cbranch, Ctips)
  leaf.rates <- rates[match(Cleaf, rownames(rates)),]
  leaf.rates <- na.omit(leaf.rates)
  NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
  leafR <- mean(abs(leaf.rates))
  NCR <- mean(abs(NCrates))
  leaf2N.diff[j] <- leafR - NCR
  NC <- length(rates) - length(leaf.rates)
  C <- length(leaf.rates)
  ran.diffM <- array()
  for (i in 1:nrep) {
    ran.diffM[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates),
                                                              NC))
  }
  p.single[j] <- rank(c(leaf2N.diff[j], ran.diffM[1:(nrep -
                                                       1)]))[1]/nrep
}
names(leaf2N.diff) <- node
names(p.single) <- node

data.frame(difference=leaf2N.diff,p.value=p.single)->tab
require(kableExtra)
knitr::kable(tab,digits=3,align="c", caption="all clades differences") %>%
  kable_styling(full_width = FALSE, position = "float_right")  %>%
  column_spec(1, bold = TRUE)

## ----echo=FALSE,message=FALSE, fig.dim=c(4,4),out.width="70%",dpi=220---------
if (length(p.single[p.single >= 0.975 | p.single <=
                    0.025])==0){
  p.single <- p.single[c(which.max(p.single), which.min(p.single))]
  leaf2N.diff <- leaf2N.diff[match(names(p.single),
                                   names(leaf2N.diff))]
  p.init <- p.single
  l2N.init <- leaf2N.diff[match(names(p.init),
                                names(leaf2N.diff))]
}


if (length(p.single[p.single >= 0.975 | p.single <=
                    0.025]) < 2) {
  p.init <- p.single
  l2N.init <- leaf2N.diff[match(names(p.init),
                                names(leaf2N.diff))]
  
  p.single <- p.single[p.single >= 0.975 | p.single <=
                         0.025]
  leaf2N.diff <- leaf2N.diff[match(names(p.single),
                                   names(leaf2N.diff))]
  
}

if (length(p.single[p.single >= 0.975 | p.single <=
                    0.025]) >= 2)  {
  p.init <- p.single
  l2N.init <- leaf2N.diff[match(names(p.init),
                                names(leaf2N.diff))]
  
  p.single <- p.single[p.single >= 0.975 | p.single <= 0.025]
  
  leaf2N.diff <- leaf2N.diff[match(names(p.single), names(leaf2N.diff))]
  
  ups <- p.single[p.single > 0.975]
  dws <- p.single[p.single < 0.025]
  ups <- ups[na.omit(match(names(leaf2N.diff[order(leaf2N.diff,
                                                   decreasing = FALSE)]), names(ups)))]
  dws <- dws[na.omit(match(names(leaf2N.diff[order(leaf2N.diff,
                                                   decreasing = FALSE)]), names(dws)))]
  if (is.na(mean(dws))) {
    dws = Nnode(tree) * 2
  } else {
    s = 1
    repeat
    {
      d <- which(names(dws) %in% getDescendants(tree, names(dws)[s]))
      if (length(d) > 0) {
        leaf2N.diff[c(match(names(dws[d]),names(leaf2N.diff)),match(names(dws[s]),names(leaf2N.diff)))]->cla
        names(which.max(abs(leaf2N.diff[c(match(names(dws[d]),names(leaf2N.diff)),match(names(dws[s]),names(leaf2N.diff)))])))->IN
        dws[-match(names(cla[which(names(cla)!=IN)]),names(dws))]->dws
        s=1
      } else {
        dws <- dws
        s = s+1
        
      }
      
      if (s > length(dws))  break
      
    }
  }
  
  if (is.na(mean(ups))) {
    ups = Nnode(tree) * 2
  } else {
    
    z = 1
    repeat
    {
      d <- which(names(ups) %in% getDescendants(tree, names(ups)[z]))
      if (length(d) > 0) {
        leaf2N.diff[c(match(names(ups[d]),names(leaf2N.diff)),match(names(ups[z]),names(leaf2N.diff)))]->cla
        names(which.max(abs(leaf2N.diff[c(match(names(ups[d]),names(leaf2N.diff)),match(names(ups[z]),names(leaf2N.diff)))])))->IN
        ups[-match(names(cla[which(names(cla)!=IN)]),names(ups))]->ups
        z=1
      } else {
        ups <- ups
        z = z+1
        
      }
      if (z > length(ups))  break
      
    }
  }
  
  
  # p.init <- p.single
  # l2N.init <- leaf2N.diff[match(names(p.init),
  #                               names(leaf2N.diff))]
  p.single <- p.single[which(names(p.single) %in%
                               names(c(ups, dws)))]
  
  leaf2N.diff <- leaf2N.diff[match(names(p.single),
                                   names(leaf2N.diff))]
}


p.single[order(p.single)]->p.single

data.frame("difference"=leaf2N.diff[match(names(p.single),names(leaf2N.diff))],"p-value"=p.single)->tab2


  plot(tree,no.margin = TRUE, show.tip.label = FALSE)
  xy <- list()
  for (w in 1:length(p.single)) {
    xy[[w]] <- unlist(sapply(get("last_plot.phylo",
                                 envir =ape::.PlotPhyloEnv), function(x) x[as.numeric(names(p.single)[w])]))[c(21,
                                                                                                               22)]
  }


  c(rep("red",length(which(p.single<=0.025))),rep("royalblue",length(which(p.single>=0.975))))->p.col
  symbols(lapply(xy, "[[", 1), lapply(xy, "[[", 2),
          circles = abs(leaf2N.diff[match(names(p.single),names(leaf2N.diff))])^0.5, inches = 0.25,
          add = TRUE, bg = alpha(p.col, 0.5), fg = p.col)

  nodelabels(node = as.numeric(names(p.single)), adj = c(1.5,
                                                         1), text = names(p.single), frame = "none", bg = "white",
             col = "purple")
  

knitr::kable(tab2,digits=3,align="c", caption="single clades differences") %>%
  kable_styling(full_width = FALSE, position = "float_left")  %>%
  column_spec(1, bold = TRUE)

## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),out.width="60%",dpi=220,fig.align="center"----
require(ape)
require(phytools)
require(geiger)
require(scales)

node = 159
Cbranch <- list()
for (i in 1:length(node)) {
  Cbranch[[i]] <- getDescendants(tree, node[i])
}
Cbranch <- unlist(Cbranch)
Cbranch <- Cbranch[-which(Cbranch < Ntip(tree))]
Ctips <- list()
for (i in 1:length(node)) {
  Ctips[[i]] <- tips(tree, node[i])
}
Ctips <- unlist(Ctips)
Ctips <- unique(Ctips)
Cbranch <- unique(Cbranch)
Cleaf <- c(Cbranch, Ctips)
leaf.rates <- rates[match(Cleaf, rownames(rates)),
                    ]
leaf.rates <- na.omit(leaf.rates)
NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
leafR <- mean(abs(leaf.rates))
NCR <- mean(abs(NCrates))
leaf2NC.diff <- leafR - NCR
NC <- length(rates) - length(leaf.rates)
C <- length(leaf.rates)
ran.diffR <- array()
for (i in 1:nrep) {
  ran.diffR[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates),
                                                            NC))
}
p.shift <- rank(c(leaf2NC.diff, ran.diffR[1:(nrep -
                                               1)]))[1]/nrep

par(mar = c(3, 2, 2, 1))
hist(ran.diffR,plot=FALSE)->hi1
hist(ran.diffR, main="",cex.lab=1.5,
     yaxt="n",xaxt="n",ylab=paste("node", node, sep = " "),xlab="",
     xlim = c(1.1 * min(c(leaf2NC.diff,ran.diffR)),1.1 * max(c(leaf2NC.diff,ran.diffR))))
hist(c(leaf2NC.diff,ran.diffR),plot=FALSE)->hi
if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
  c(hi$breaks[seq(1,length(hi$breaks),2)],
    (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
axis(1,at=athi,mgp=c(0.2,0.3,0))
mtext(text="random differences",1,line=1.2)
abline(v = leaf2NC.diff, col = "green", lwd = 3)
text(labels=paste("p-value =",round(p.shift,3)),x=leaf2NC.diff-0.12,y=max(hi1$counts)/2,cex=1.1)
        

## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),out.width='47%',dpi=220----------
require(ape)
require(phytools)
require(geiger)
require(scales)

node = c(159,177)
      Cbranch <- list()
      for (i in 1:length(node)) {
        Cbranch[[i]] <- getDescendants(tree, node[i])
      }
      Cbranch <- unlist(Cbranch)
      Cbranch <- Cbranch[-which(Cbranch < Ntip(tree))]
      Ctips <- list()
      for (i in 1:length(node)) {
        Ctips[[i]] <- tips(tree, node[i])
      }
      Ctips <- unlist(Ctips)
      Ctips <- unique(Ctips)
      Cbranch <- unique(Cbranch)
      Cleaf <- c(Cbranch, Ctips)
      leaf.rates <- rates[match(Cleaf, rownames(rates)),
                          ]
      leaf.rates <- na.omit(leaf.rates)
      NCrates <- rates[-match(names(leaf.rates), rownames(rates))]
      leafR <- mean(abs(leaf.rates))
      NCR <- mean(abs(NCrates))
      leaf2NC.diff <- leafR - NCR
      NC <- length(rates) - length(leaf.rates)
      C <- length(leaf.rates)
      ran.diffR <- array()
      for (i in 1:nrep) {
        ran.diffR[i] <- mean(sample(abs(rates), C)) - mean(sample(abs(rates),
                                                                  NC))
      }
      p.shift <- rank(c(leaf2NC.diff, ran.diffR[1:(nrep -
                                                     1)]))[1]/nrep
      #pdf(file=paste(foldername, "AR results for rate differences.pdf",sep="/"),width=8.3,height=11.7)
        leaf2N.diff <- array()
        p.single <- array()
        ran.diff <- list()
        for (i in 1:length(node)) {
          NOD <- node[-i]
          others <- list()
          mommies <- list()
          for (j in 1:length(NOD)) {
            others[[j]] <- tips(tree, NOD[j])
            mommies[[j]] <- getDescendants(tree, NOD[j])
          }
          others <- unlist(others)
          mommies <- unlist(mommies)
          mommies <- mommies[-which(mommies < Ntip(tree) +
                                      1)]
          otmom <- c(mommies, others)
          Ctips <- tips(tree, node[i])
          Ctips <- unlist(Ctips)
          Cbranch <- getDescendants(tree, node[i])
          Cleaf <- c(Cbranch, Ctips)
          leaf.rates <- rates[match(Cleaf, rownames(rates)),
                              ]
          leaf.rates <- na.omit(leaf.rates)
          NC <- rates[-c(which(rownames(rates) %in% names(leaf.rates)),
                         which(rownames(rates) %in% otmom)), ]
          NR.r <- mean(abs(NC))
          leaf.r <- mean(abs(leaf.rates))
          leaf2N.diff[i] <- leaf.r - NR.r
          NC.l <- length(NC)
          leaf.l <- length(leaf.rates)
          tot.r <- abs(c(NC, leaf.rates))
          RAN.diff <- array()
          for (k in 1:nrep) {
            RAN.diff[k] <- mean(sample(tot.r, leaf.l)) -
              mean(sample(tot.r, NC.l))
            ran.diff[[i]] <- RAN.diff
          }
          p.single[i] <- rank(c(leaf2N.diff[i], RAN.diff[1:(nrep -
                                                              1)]))[1]/nrep
        }
        names(p.single) <- node
        names(leaf2N.diff) <- names(p.single)
        c(0.1,0)->sub
        leaf2N.diff[match(names(p.single),names(leaf2N.diff))]->ldiff

        for (m in 1:length(node)) {
          par(mar = c(3, 2, 2, 1))
          hist(ran.diff[[m]],plot=FALSE)->hi1
          hist(ran.diff[[m]], cex.lab=1.5,
               yaxt="n",ylab="",main=paste("node", node[m], sep = " "),xlab="",xaxt="n",mgp=c(0.2,0.3,0),
               xlim = c(1.1 * min(c(leaf2N.diff[m],ran.diff[[m]])),1.1 * max(c(leaf2N.diff[m],ran.diff[[m]]))))
          hist(c(leaf2N.diff[m],ran.diff[[m]]),plot=FALSE)->hi
          if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
            c(hi$breaks[seq(1,length(hi$breaks),2)],
              (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
          axis(1,at=athi,mgp=c(0.2,0.3,0))
          mtext(text="random differences",1,line=1)
          abline(v = leaf2N.diff[m], col = "blue", lwd = 3)
          text(labels=paste("p-value =",round(p.single[m],3)),x=ldiff[m]-sub[m],y=max(hi1$counts)/2,cex=1.1)
        }

## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),fig.align='center',dpi=220,out.width='47%'----
        hist(ran.diffR,plot=FALSE)->hi1
      p.shift <- rank(c(leaf2NC.diff, ran.diffR[1:(nrep -
                                                     1)]))[1]/nrep
par(mar = c(3, 2, 2, 1))
 hist(ran.diffR, cex.lab=1.5,
             yaxt="n",ylab="",main="All clades together",xlab="",xaxt="n",
             xlim = c(1.1 * min(c(leaf2NC.diff,ran.diffR)),1.1 * max(c(leaf2NC.diff,ran.diffR))))
        hist(c(leaf2NC.diff,ran.diffR),plot=FALSE)->hi
        if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
          c(hi$breaks[seq(1,length(hi$breaks),2)],
            (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
        axis(1,at=athi,mgp=c(0.2,0.3,0))
        mtext(text="random differences",1,line=1)
        abline(v = leaf2NC.diff, col = "green", lwd = 3)
        text(labels=paste("p-value =",round(p.shift,3)),x=leaf2NC.diff-0.07,y=max(hi1$counts)/2,cex=1.1)
          

## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),out.width="60%",dpi=220,fig.align="center"----
require(ape)
require(phytools)
require(geiger)
require(scales)

set.seed(14)
rep("a",100)->categ
categ[sample(1:100,30)]<-"b"
names(categ)<-Tree$tip.label

y[which(categ=="b")]*5->y[which(categ=="b")]

RRphylo(Tree,y)->RR

data.frame(rate=RR$rates[Ntip(Tree):(Ntip(Tree)+Nnode(Tree)),],
           status=categ[match(names(RR$rates[Ntip(Tree):(Ntip(Tree)+Nnode(Tree)),]),names(categ))])->frame

status.diff <- diff(tapply(abs(frame$rate), frame$status,
                           mean))
status.diffS <- array()
for (i in 1:nrep) {
  s.state <- frame$status
  s.ran <- sample(s.state)
  s.frame <- data.frame(s.ran, frame$rate)
  s.frame[, 1] <- as.factor(s.frame[, 1])
  status.diffS[i] <- diff(tapply(abs(s.frame$frame.rate),
                                 s.frame$s.ran, mean))
}
p.status.diff <- rank(c(status.diff, status.diffS[1:(nrep -
                                                       1)]))[1]/nrep

par(mar = c(3, 2, 2, 1))
hist(status.diffS,plot=FALSE)->hi1
hist(status.diffS,xaxt="n",yaxt="n",ylab="",xlab="random differences",
     main=paste("Absolute rate difference \n between",names(p.status.diff),"and",
                unique(frame$status)[which(unique(frame$status)!=names(p.status.diff))],sep=" "),
     cex.lab=1,mgp=c(1.5,0.8,0),
     xlim = c(1.1*min(c(status.diff,status.diffS)), 1.1*max(c(status.diff,status.diffS))))
hist(c(status.diff,status.diffS),plot=FALSE)->hi
if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
  c(hi$breaks[seq(1,length(hi$breaks),2)],
    (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
axis(1,at=athi,mgp=c(0.2,0.5,0))
abline(v = status.diff, lwd = 3, col = "green")
text(labels=paste("p-value =",round(p.status.diff,3)),
     x=2,y=max(hi1$counts)/2,cex=1.1)

## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),out.width="47%",dpi=220----------
#,fig.align="center"}
require(ape)
require(phytools)
require(geiger)
require(scales)

rep("a",100)->categ
categ[sample(1:100,30)]<-"b"
categ[sample(which(categ=="a"),20)]<-"c"
names(categ)<-Tree$tip.label

y[which(categ=="b")]*5->y[which(categ=="b")]


RRphylo(Tree,y)->RR

data.frame(rate=RR$rates[Ntip(Tree):(Ntip(Tree)+Nnode(Tree)),],
           status=as.factor(categ[match(names(RR$rates[Ntip(Tree):(Ntip(Tree)+Nnode(Tree)),]),names(categ))]))->frame

 status.diff <- apply(combn(tapply(abs(frame$rate),
                                      frame$status, mean), 2), 2, diff)
    sta <- tapply(abs(frame$rate), frame$status, mean)
    sta <- sta[match(unique(frame$status), names(sta))]
    w <- array()
    for (x in 1:length(sta)) w[x] <- sta[x] - mean(abs(frame[-which(frame$status ==
                                                                      names(sta)[x]), 1]))
    names(w) <- names(sta)
    status.diff <- c(status.diff, w)
    names(status.diff) <- c(apply(combn(levels(frame$status),
                                        2), 2, function(x) paste(x[2], x[1], sep = "_")),
                            names(sta))
    status.diffS <- matrix(ncol = length(status.diff),
                           nrow = nrep)
    for (i in 1:nrep) {
      s.state <- frame$status
      s.ran <- sample(s.state)
      s.frame <- data.frame(s.ran, frame$rate)
      SD <- apply(combn(tapply(abs(s.frame$frame.rate),
                               s.frame$s.ran, mean), 2), 2, diff)
      sta <- tapply(abs(frame$rate), s.ran, mean)
      sta <- sta[match(unique(frame$status), names(sta))]
      w <- array()
      for (x in 1:length(sta)) w[x] <- sta[x] - mean(abs(frame[-which(s.frame$s.ran ==
                                                                        names(sta)[x]), 1]))
      status.diffS[i, ] <- c(SD, w)
    }
    colnames(status.diffS) <- names(status.diff)
    
    p.status.diff <- array()
    for (i in 1:length(status.diff)) p.status.diff[i] <- rank(c(status.diff[i],
                                                                  status.diffS[1:(nrep - 1), i]))[1]/nrep
    
    idx <- match(unique(frame$status), colnames(status.diffS))
    sub<-c(1,0)
    for (i in 1:2) {
      par(mar = c(3, 2, 2, 1))
      hist(status.diffS[, idx[i]],plot=FALSE)->hi1
      hist(status.diffS[, idx[i]], ylab="",xaxt="n",cex.lab=1.5,
           yaxt="n",main=paste("state", colnames(status.diffS)[idx[i]], sep = " "),
           xlab="random differences",mgp=c(1.5,0.8,0),
           xlim = c(1.1*min(c(status.diff[idx[i]],status.diffS[, idx[i]])),
                    1.1*max(c(status.diff[idx[i]],status.diffS[, idx[i]]))))
      hist(c(status.diff[idx[i]],status.diffS[, idx[i]]),plot=FALSE)->hi
      if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
        c(hi$breaks[seq(1,length(hi$breaks),2)],
          (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
      axis(1,at=athi,mgp=c(0.2,0.5,0))
      abline(v = status.diff[idx[i]], lwd = 3, col = "green")
      text(labels=paste("p-value =",round(p.status.diff[idx[i]],3)),
           x=status.diff[idx[i]]+sub[i],y=max(hi1$counts)/2,cex=1.1)

    }


## ----echo=FALSE,message=FALSE,fig.dim=c(4,4),out.width="47%",dpi=220,fig.align="center"----

i=3
par(mfrow=c(1,1))
par(mar = c(3, 2, 2, 1))
hist(status.diffS[, idx[i]],plot=FALSE)->hi1
hist(status.diffS[, idx[i]], ylab="",xaxt="n",cex.lab=1.5,
     yaxt="n",main=paste("state", colnames(status.diffS)[idx[i]], sep = " "),
     xlab="random differences",mgp=c(1.5,0.8,0),
     xlim = c(1.1*min(c(status.diff[idx[i]],status.diffS[, idx[i]])),
              1.1*max(c(status.diff[idx[i]],status.diffS[, idx[i]]))))
hist(c(status.diff[idx[i]],status.diffS[, idx[i]]),plot=FALSE)->hi
if(length(hi$breaks)%%2==1) (hi$breaks[seq(1,length(hi$breaks),2)])->athi else
  c(hi$breaks[seq(1,length(hi$breaks),2)],
    (hi$breaks[length(hi$breaks)]+abs(diff(hi$breaks[seq(1,length(hi$breaks),2)][1:2]))))->athi
axis(1,at=athi,mgp=c(0.2,0.5,0))
abline(v = status.diff[idx[i]], lwd = 3, col = "green")
text(labels=paste("p-value =",round(p.status.diff[idx[i]],3)),
     x=status.diff[idx[i]]-1.2,y=max(hi1$counts)/2,cex=1.1)


## ----out.width='99%',fig.dim=c(7,6),message=FALSE,dpi=220,warning=FALSE,echo=1:5----
# load the RRphylo example dataset including Ornithodirans tree and data
DataOrnithodirans$treedino->treedino # phylogenetic tree
DataOrnithodirans$massdino->massdino # body mass data
DataOrnithodirans$statedino->statedino # locomotory type data
log(massdino)->lmass

require(geiger)
require(plotrix)
require(RColorBrewer)
statedino->colo
as.character(colo)->colo
names(colo)<-names(statedino)
colo[which(colo=="B")]<-brewer.pal(n = 6, name = "Set1")[6]
colo[which(colo=="BQ")]<-brewer.pal(n = 6, name = "Set1")[5]
colo[which(colo=="F")]<-brewer.pal(n = 6, name = "Set1")[3]
colo[which(colo=="Q")]<-brewer.pal(n = 6, name = "Set1")[1]
ladderize(treedino,FALSE)->treed
#par(mar=c(9,4,1,3))
plot(treed,type="fan",show.tip.label = FALSE,edge.color = "gray40",no.margin = TRUE)
tiplabels(text=rep("",419),bg=colo,frame="circle",cex=.4)
plotinfo<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)
sqrt(plotinfo$xx[416]^2+plotinfo$yy[416]^2)->rad

c(422,516,583,746,689)->clads
c("Ornithischia","Sauropodomorpha","Theropoda","Pterosauria","Avialae")->nams
#c(2,4,3,3,1)->poss
list(c(0.5,0),c(0,1),c(0.5,-0.5),c(0.3,0),c(1,0))->poss
for(k in 1:length(clads)){
  clads[k]->aa
  plotinfo$yy[which(treed$tip.label%in%tips(treed,aa)[c(1,length(tips(treed,aa)))])]->yc
  plotinfo$xx[which(treed$tip.label%in%tips(treed,aa)[c(1,length(tips(treed,aa)))])]->xc
  sqrt(xc[1]^2+yc[1]^2)->r1
  if(yc[1]>=0) acos(xc[1]/r1)->ang1 else 2*pi-acos(xc[1]/r1)->ang1
  sqrt(xc[2]^2+yc[2]^2)->r2
  if(yc[2]>=0) acos(xc[2]/r2)->ang2 else 2*pi-acos(xc[2]/r2)->ang2
  
  if(k<5){
    draw.arc(x=0,y=0,radius=rad+10,angle1 =ang1+0.05,angle2 =ang2,lwd=3,col=brewer.pal(n = 6, name = "Set1")[2]) 
    c(plotinfo$xx[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])],
      plotinfo$yy[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])])->cc
    xcc=sqrt((cc[1]^2*(180+25)^2)/(cc[1]^2+cc[2]^2))
    if(cc[1]>=0) text(x=xcc,y=cc[2]/cc[1]*xcc,labels=nams[k],font=2,adj=poss[[k]]) else text(x=(-xcc),y=cc[2]/cc[1]*-xcc,labels=nams[k],font=2,adj=poss[[k]])
  }else{
    draw.arc(x=0,y=0,radius=rad+20,angle1 =ang1+0.05,angle2 =ang2,lwd=3,col=brewer.pal(n = 6, name = "Set1")[4])
    c(plotinfo$xx[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])],
      plotinfo$yy[which(treed$tip.label==tips(treed,aa)[length(tips(treed,aa))/2])])->cc
    xcc=sqrt((cc[1]^2*(180+35)^2)/(cc[1]^2+cc[2]^2))
    if(cc[1]>=0) text(x=xcc,y=cc[2]/cc[1]*xcc,labels=nams[k],font=2,adj=poss[[k]]) else text(x=(-xcc),y=cc[2]/cc[1]*-xcc,labels=nams[k],font=2,adj=poss[[k]])
  }
}

legend("topleft",legend=c("B","BQ","F","Q"),pt.bg =unique(colo)[c(1,3,4,2)],pch=21,col="black",pt.cex=1.5,bty="n")


## ----message=FALSE,eval=FALSE-------------------------------------------------
#  # check the order of your data: best if data vectors
#  # are sorted in the same order of the species on the phylogeny
#  lmass[match(treedino$tip.label,names(lmass))]->lmass
#  statedino[match(treedino$tip.label,names(statedino))]->statedino
#  
#  # perform RRphylo on the vector of (log) body mass
#  RRphylo(tree=treedino,y=lmass)->RRdinomass
#  
#  # search for clades showing significant shifts in mass specific evolutionary rates
#  # (i.e. using the log body mass itself as a covariate)
#  search.shift(RRdinomass, status.type= "clade",cov=lmass,filename="ss automatic")->SSauto
#  
#  # search for shifts in mass specific evolutionary rates pertaining different locomotory types.
#  search.shift(RRdinomass, status.type= "sparse", state=statedino,cov=lmass, filename="ss state")->SSstate

