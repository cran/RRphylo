## ---- include = FALSE---------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(knitr.kable.NA = '',rmarkdown.html_vignette.check_title = FALSE)

## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,4),fig.align='center',out.width='100%',dpi=220----
require(ape)
require(geiger)
require(phytools)
require(ddpcr)

RRphylo:::range01->range01

sim.bdtree(b = .5, d = 0.2,seed=14)->tree

set.seed(14)
fastBM(tree,a=0,sig2=1)->y

### TREND ###
es=2
yt <- (diag(vcv(tree))^es)/(diag(vcv(tree))) * y
c(0,nodeHeights(tree)[,2])->rootD
hei <- tree$edge[, 2]
hei[which(tree$edge[, 2] < Ntip(tree) + 1)] <- tree$tip.label
names(rootD)<-c(Ntip(tree)+1,hei)

RRphylo(tree,yt)->rr
rr$tip.path->L
rr$node.path->L1
rr$ace[1]->rootV
rr$lambda->lambda
range01(abs(rr$rates))->rts
as.matrix(as.data.frame(rts[match(names(rootD),rownames(rts)),]))->rts
quiet(search.trend(rr,yt,clus=2/parallel::detectCores(),filename=paste(tempdir(),"result", sep = "/"))->st.rates)

BTS<-list()
for(i in 1:100){
  fastBM(tree,sig2=1,a=rootV,bounds=c(min(yt),max(yt)))->yb
  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% 
              t(L)) %*% (as.matrix(yb)-rootV)
  range01(abs(betas))->bts
  as.matrix(as.data.frame(bts[match(names(rootD),rownames(bts)),]))->bts
  bts -> BTS[[i]]
}
do.call(cbind,BTS)->fu

apply(fu,1,function(x) quantile(x,.975))->ful
apply(fu,1,function(x) quantile(x,.025))->mins
predict(loess(ful~rootD))->fk
predict(loess(mins~rootD))->fk.min
names(fk)<-names(ful)
names(fk.min)<-names(ful)

##abline(lm(rts~rootD))
data.frame(rootD,fk,fk.min)->dfk
dfk[order(dfk[,1]),]->dfk

data.frame(rts=rts,col=as.character(rep("red",nrow(rts))))->rts
rts[which(rownames(rts)%in%tree$tip.label),2]<-"green"

par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(NA,ylim=range(rts[,1]),xlim=range(rootD),
     mgp=c(1.8,0.5,0),
     xlab="distance from the root",ylab="rescaled rates",main="Simulation for trend\nin absolute rates")
polygon(c(dfk$rootD,rev(dfk$rootD)),c(dfk$fk,rev(dfk$fk.min)),col = rgb(0.5, 0.5, 0.5,
                                                                        0.3), border = NA)
points(rootD,rts[,1],cex=1.5,bg=rts$col,
       col="black",pch=21)
abline(lm(rts[,1]~rootD),lwd=2,col="#ff00ff")

legend("topleft",legend=c("nodes","tips","brownian range"),fill=c("red","green","#dadad9"),bty="n")

#### DRIFT
ds=1
yd<-range01((y+diag(vcv(tree))*ds))

RRphylo(tree,yd)->rr
rr$lambda->lambda
c(rr$aces,yd)->phen
names(phen)<-c(rownames(rr$aces),names(yd))
rootD[match(names(phen),names(rootD))]->rootP
phen[[1]]->rootV
lm(phen~rootP)->regr
quiet(search.trend(rr,yd,clus=2/parallel::detectCores(),filename=paste(tempdir(),"result", sep = "/"))->st.phen)


phenD<-list()
for(i in 1:100){
  fastBM(tree,sig2=1,a=rootV,bounds=c(min(yd),max(yd)))->yc
  RRphylo:::range01(yc)->yc
  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*% 
              t(L)) %*% (as.matrix(yc) - rootV)
  aceRR <- (L1 %*% betas[1:Nnode(tree), ]) + rootV
  c(aceRR,yc)->phenC
  names(phenC)<-c(rownames(aceRR),names(yc))
  phenC[match(names(rootP),names(phenC))]->phenD[[i]]
}

sapply(phenD,function(x) x[match(names(rootP),names(x))])->fu

apply(fu,1,function(x) quantile(x,.975))->ful
apply(fu,1,function(x) quantile(x,.025))->mins
predict(loess(ful~rootP))->fk
predict(loess(mins~rootP))->fk.min
names(fk)<-names(ful)
names(fk.min)<-names(ful)

data.frame(rootP,fk,fk.min)->dfk
dfk[order(dfk[,1]),]->dfk

data.frame(phen=phen,col=as.character(rep("red",length(phen))))->phen
phen[which(rownames(phen)%in%tree$tip.label),2]<-"green"

plot(NA,ylim=range(phen[,1]),xlim=range(rootP),
     mgp=c(1.8,0.5,0),
     xlab="distance from the root",ylab="rescaled phenotypes",main="Simulation for trend\nin mean phenotypes")
polygon(c(dfk$rootP,rev(dfk$rootP)),c(dfk$fk,rev(dfk$fk.min)),col = rgb(0.5, 0.5, 0.5,
                                                                        0.3), border = NA)
points(rootP,phen[,1],cex=1.5,bg=phen$col,
       col="black",pch=21)
abline(lm(phen[,1]~rootP),lwd=2,col="#ff00ff")



as.data.frame(rbind(c(st.rates[[3]],NA),c(st.phen[[2]][1:3],NA,st.phen[[2]][4])))->res

colnames(res)[4:5]<-c("spread","dev")
rownames(res)<-c("rescaled absolute rate regression","phenotypic regression")

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  search.trend(RR=RR,y=y,filename="st whole")->ST

## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
require(kableExtra)

knitr::kable(res,digits=3,align="c") %>%
kable_styling(full_width = F, position = "center")  


## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(8,4),fig.align='center',out.width='100%',dpi=220----
require(plotrix)

sim.bdtree(b = .5, d = 0.2,seed=14)->T

n=275
n2=198
max(nodeHeights(T))/nodeHeights(T)[,2][which(T$edge[,2]==n)]->dN1
max(nodeHeights(T))/nodeHeights(T)[,2][which(T$edge[,2]==n2)]->dN2

rt=-9.365124
S2=1.04
set.seed(22)
fastBM(T,a=rt,sig2=S2)->yB

### Rates ###
esx1=.01
esx2=.01

yB->y
y[match(tips(T,n),names(y))]->yT
diag(vcv(T)[names(yT),names(yT)])->timeT
yT <- ((timeT^esx1)/timeT)*yT

y[match(tips(T,n2),names(y))]->yT2
diag(vcv(T)[names(yT2),names(yT2)])->timeT2
yT2 <- ((timeT2^esx2)/timeT2)*yT2
y[match(tips(T,n),names(y))]<-yT
y[match(tips(T,n2),names(y))]<-yT2

RRphylo(T,y)->RR
quiet(search.trend(RR,y,node=c(n,n2),clus=2/parallel::detectCores(),filename=paste(tempdir(),"result", sep = "/"))->STrates)
abs(RR$rates[,1])->rats


### Phenotypes ###
set.seed(2293)
fastBM(T,a=rt,sig2=S2)->yB
ds1=.2
ds2=-.2

yB->y
y[match(tips(T,n),names(y))]->yDR
diag(vcv(T)[names(yDR),names(yDR)])->timeDR
coef(lm(yDR~timeDR))[2]->tend
if(tend*ds1>0) ds1->ds1 else ds1*(-1)->ds1
ds1*dN1->dsx1
yD <- (timeDR*dsx1)+ yDR
y[match(tips(T,n),names(y))]<-yD

y[match(tips(T,n2),names(y))]->yDR
diag(vcv(T)[names(yDR),names(yDR)])->timeDR
coef(lm(yDR~timeDR))[2]->tend
# if(tend*ds2>0) ds2->ds2 else 
ds2*(-1)->ds2
ds2*dN2->dsx2
yD <- (timeDR*dsx2)+ yDR
y[match(tips(T,n2),names(y))]<-yD

RRphylo(T,y)->RR
quiet(search.trend(RR,y,node=c(n,n2),clus=2/parallel::detectCores(),filename=paste(tempdir(),"result", sep = "/"))->STphen)
c(RR$aces[,1],y)->phen

c(n,getDescendants(T,n))->desn
T$tip.label[desn[which(desn<=Ntip(T))]]->desn[which(desn<=Ntip(T))]

c(n2,getDescendants(T,n2))->desn2
T$tip.label[desn2[which(desn2<=Ntip(T))]]->desn2[which(desn2<=Ntip(T))]

dist.nodes(T)[Ntip(T)+1,]->ages
T$tip.label[as.numeric(names(ages)[which(as.numeric(names(ages))<=Ntip(T))])]->names(ages)[which(as.numeric(names(ages))<=Ntip(T))]
ages[match(names(phen),names(ages))]->ages


par(mfrow=c(1,2),mar=c(4,3,3,1))
plot(max(ages)-ages[-which(names(ages)%in%c(desn,desn2))],rats[-which(names(rats)%in%c(desn,desn2))],
     ylim=range(rats),xlim=c(max(ages),min(ages)),
     pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="age",ylab="absolute rates",main="Simulation for trend\nin absolute rates")
abline(lm(rats~ages),lwd=3,col="gray20")
points(max(ages)-ages[which(names(ages)%in%desn)],rats[which(names(rats)%in%desn)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="slateblue2",cex=1.5)
ablineclip(lm(rats[which(names(rats)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="black",lwd=4.5)
ablineclip(lm(rats[which(names(rats)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="slateblue2",lwd=3)
points(max(ages)-ages[which(names(ages)%in%desn2)],rats[which(names(rats)%in%desn2)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="#ae9eff",cex=1.5)
ablineclip(lm(rats[which(names(rats)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="black",lwd=4.5)
ablineclip(lm(rats[which(names(rats)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="#ae9eff",lwd=3)
legend("topleft",legend=c("node 275","node 198","entire tree"),fill=c("slateblue2","#ae9eff","gray20"),bty="n",x.intersp = .5)

plot(max(ages)-ages[-which(names(ages)%in%c(desn,desn2))],phen[-which(names(phen)%in%c(desn,desn2))],
     ylim=range(phen),xlim=c(max(ages),min(ages)),
     pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),
     xlab="age",ylab="phenotype",main="Simulation for trend\nin mean phenotypes")
abline(lm(phen~ages),lwd=3,col="gray20")
points(max(ages)-ages[which(names(ages)%in%desn)],phen[which(names(phen)%in%desn)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="aquamarine3",cex=1.5)
ablineclip(lm(phen[which(names(phen)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="black",lwd=4.5)
ablineclip(lm(phen[which(names(phen)%in%desn)]~I(max(ages)-ages[which(names(ages)%in%desn)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn)]),
           col="aquamarine3",lwd=3)

points(max(ages)-ages[which(names(ages)%in%desn2)],phen[which(names(phen)%in%desn2)],
       xlim=c(max(ages),min(ages)),
       pch=21,bg="aquamarine",cex=1.5)
ablineclip(lm(phen[which(names(phen)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="black",lwd=4.5)
ablineclip(lm(phen[which(names(phen)%in%desn2)]~I(max(ages)-ages[which(names(ages)%in%desn2)])),
           x1=max(max(ages)-ages[which(names(ages)%in%desn2)]),
           x2=min(max(ages)-ages[which(names(ages)%in%desn2)]),
           col="aquamarine",lwd=3)
legend("topleft",legend=c("node 275","node 198","entire tree"),fill=c("aquamarine3","aquamarine","gray20"),bty="n",x.intersp = .5)

cbind(data.frame(node=names(STrates[[5]]),do.call(rbind,STrates[[5]])),
      data.frame(node=names(STphen[[4]]),do.call(rbind,STphen[[4]])))->res

## ----eval=FALSE,message=FALSE,warning=FALSE-----------------------------------
#  search.trend(RR=RR,y=y,node=c(275,198),filename="st nodes")->ST

## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
rownames(res)<-NULL
knitr::kable(res,digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
column_spec(5, border_right = TRUE) %>%
add_header_above(c("Trend in absolute rates" = 5, "Trend in phenotypic means" = 5))


## ----echo=FALSE,message=FALSE,warning=FALSE-----------------------------------
knitr::kable(STrates[[6]][[2]],digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
add_header_above(c("Comparison of trends in absolute rates" = 6))

knitr::kable(STphen[[6]][[1]],digits=3,align="c") %>%
kable_styling(full_width = F, position = "center") %>%
add_header_above(c("Comparison of trends in phenotypic means" = 6))


## ----out.width='100%',fig.dim=c(8,7),message=FALSE,dpi=220--------------------
# load the RRphylo example dataset including Cetaceans tree and data
data("DataCetaceans")
DataCetaceans$treecet->treecet # phylogenetic tree
DataCetaceans$masscet->masscet # logged body mass data
DataCetaceans$brainmasscet->brainmasscet # logged brain mass data
DataCetaceans$aceMyst->aceMyst # known phenotypic value for the most recent 
                               # common ancestor of Mysticeti

require(geiger)
par(mar=c(0,0,0,1))
plot(ladderize(treecet,FALSE),show.tip.label = FALSE,edge.color = "gray40")
plotinfo<-get("last_plot.phylo",envir =ape::.PlotPhyloEnv)
nodelabels(text="",node=128,frame="circle",bg="red",cex=0.5)
nodelabels(text="Mystacodon",node=128,frame="n",bg="w",cex=1,adj=c(-0.1,0.5),font=2)
range(plotinfo$yy[which(treecet$tip.label%in%tips(treecet,128))])->yran128
rect(plotinfo$x.lim[2]+0.4,yran128[1],plotinfo$x.lim[2]+0.7,yran128[2],col="red",border="red")
range(plotinfo$yy[which(treecet$tip.label%in%tips(treecet,142))])->yran142
rect(plotinfo$x.lim[2]+0.4,yran142[1],plotinfo$x.lim[2]+0.7,yran142[2],col="blue",border="blue")
mtext(c("Mysticeti","Odontoceti"), side = 4,line=-0.5,at=c(sum(yran128)/2,sum(yran142)/2),
      cex=1.5,adj=0.5,col=c("red","blue"))

## ----eval=FALSE,message=FALSE,dpi=220-----------------------------------------
#  # check the order of your data: best if data vectors
#  # are sorted in the same order of the species on the phylogeny
#  masscet[match(treecet$tip.label,names(masscet))]->masscet
#  
#  ## Example 1: search.trend by setting values at internal nodes
#  # Set the body mass of Mysticetes ancestor (Mystacodon selenensis)
#  # as known value at node and perform RRphylo on the vector of (log) body mass
#  RRphylo(tree=treecet,y=masscet,aces=aceMyst)->RR
#  
#  # Perform search.trend on the RR object and (log) body mass by indicating Mysticeti as focal clade
#  search.trend(RR=RR,y=masscet,node=as.numeric(names(aceMyst)),filename="st mysticetes")
#  
#  
#  ## Example 2: search.trend on multiple regression version of RRphylo
#  # cross-reference the phylogenetic tree and body and brain mass data. Remove from both the tree and
#  # vector of body sizes the species whose brain size is missing
#  drop.tip(treecet,treecet$tip.label[-match(names(brainmasscet),treecet$tip.label)])->treecet.multi
#  masscet[match(treecet.multi$tip.label,names(masscet))]->masscet.multi
#  
#  # check the order of your data: best if
#  # data vectors (i.e. masscet and brainmasscet) are sorted
#  # in the same order of the species on the phylogeny
#  masscet.multi[match(treecet.multi$tip.label,names(masscet.multi))]->masscet.multi
#  brainmasscet[match(treecet.multi$tip.label,names(brainmasscet))]->brainmasscet
#  
#  # perform RRphylo on tree and body mass data
#  RRphylo(tree=treecet.multi,y=masscet.multi)->RRmass.multi
#  
#  # create the predictor vector: retrieve the ancestral character estimates
#  # of body size at internal nodes from the RR object ($aces) and collate them
#  # to the vector of species' body sizes to create
#  c(RRmass.multi$aces[,1],masscet.multi)->x1.mass
#  
#  # perform the multiple regression version of RRphylo by using
#  # the brain size as variable and the body size as predictor
#  RRphylo(treecet.multi,y=brainmasscet,x1=x1.mass)->RRmulti
#  
#  # Perform search.trend on the multiple RR object to inspect the effect of body
#  # size on absolute rates temporal trend only
#  search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,clus=cc,filename="st predictor")
#  
#  # Perform search.trend on the multiple RR object to inspect the effect of body
#  # size on trends in both absolute evolutionary rates and phenotypic values
#  # (by using brain versus body mass regression residuals)
#  search.trend(RR=RRmulti, y=brainmasscet,x1=x1.mass,x1.residuals=TRUE,
#               clus=cc,filename="st residuals")
#  

