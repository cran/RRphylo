## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(RRphylo)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----echo=FALSE,message=FALSE,warning=FALSE,fig.dim=c(6,4),out.width="70%",dpi=220----
require(ape)
require(phytools)
require(geiger)

set.seed(14)
rtree(10)->tree.back

par(mar=c(2,3.5,1,3.5))
plot(tree.back,edge.width=1.2)
title("backbone tree")
axisPhylo()

data.frame(bind=c("a1","a2","a3","a4","a5","a6","a7","t5","a8"),
           reference=c("t6","t10","t9","a2-t5","t10-t2","t2-a3","a1-t10","t7-t10","t6"),
           poly=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE))->dato


require(kableExtra)
knitr::kable(dato,align="c") %>%
  kable_styling(full_width = FALSE, position = "float_right") %>%
add_header_above(c("dato" = 3)) 

## ----echo=1,results = "hide",message=FALSE,fig.dim=c(6,6),out.width="70%",dpi=220,fig.align='center'----
tree.merger(backbone=tree.back,data=dato,plot=FALSE)

tm<-function(backbone,data,source.tree=NULL,tip.ages = NULL, node.ages = NULL,min.branch=NULL,title){
  # require(ape)
  # require(phytools)
  # require(geiger)
  # require(manipulate)

  backbone->tree
  data->dat
  source.tree->tree2
  max(diag(vcv(tree)))->H
  H-diag(vcv(tree))->ages
  if(is.null(min.branch)) min(tree$edge.length)->min.branch

  ### Check data ###
  if(!all(colnames(dat)%in%c("bind","reference","poly"))) {
    if(any(is.na(as.logical(dat[,3])))) stop("Check columns order: it should be 'bind', 'reference', 'poly'")
    warning("Colnames not matching: columns assumed to be ordered as 'bind','reference','poly'",immediate. = TRUE)
    colnames(dat)<-c("bind","reference","poly")
  }

  ### Check for duplicated bind ###
  if(any(duplicated(dat$bind))) stop(paste(paste(dat$bind[duplicated(dat$bind)],collapse = ", "),"names duplicated in supplied tips"))

  if(!is.logical(dat$poly)) as.logical(dat$poly)->dat$poly
  data.frame(dat,bind.type=sapply(strsplit(dat[,1],"-"),length))->dat

  if(all(dat$bind.type==1)&(!is.null(tree2))) tree2<-NULL

  if(any(dat$bind%in%tree2$tip.label)) drop.tip(tree2,dat[which(dat$bind%in%tree2$tip.label),1])->tree2

  dat$bind.tips<-NA
  if(all(dat$bind.type==1)) dat$bind.tips<-dat$bind else{
    dat[which(dat$bind.type==2),]$bind.tips<-lapply(dat$bind[which(dat$bind.type==2)],function(x) tips(tree2,getMRCA(tree2,strsplit(x,"-")[[1]])))
    if(any(is.na(dat$bind.tips))) dat[which(is.na(dat$bind.tips)),]$bind.tips<-dat[which(is.na(dat$bind.tips)),]$bind
  }

  if(!all(unlist(dat[which(dat$bind.type==2),]$bind.tips)%in%tree2$tip.label)){
    stop(paste(paste(unlist(dat[which(dat$bind.type==2),]$bind.tips)[which(
      unlist(dat[which(dat$bind.type==2),]$bind.tips)%in%tree2$tip.label)],
      collapse=", "),"not in source.tree"))
  }


  if(!all(unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,unlist(dat$bind.tips),tree2$tip.label))){
    stop(paste(paste(unlist(strsplit(dat$reference,"-"))[which(!unlist(strsplit(dat$reference,"-"))%in%c(tree$tip.label,unlist(dat$bind.tips),tree2$tip.label))],
                     collapse=","),"missing from the backbone, the source and the tips to be attached"))
  }

  ### bind already on the backbone ###
  if(any(dat$bind%in%tree$tip.label)){
    warning(paste(paste(dat[which(dat$bind%in%tree$tip.label),1],collapse=", "),"removed from the backbone tree"),immediate. = TRUE)
    drop.tip(tree,dat[which(dat$bind%in%tree$tip.label),1])->tree
  }

  ### duplicated reference ###
  table(dat$reference)->tab.ref
  if(any(tab.ref>1)){
    for(j in 1:length(which(tab.ref>1))){
      dat[which(dat$reference==names(which(tab.ref>1)[j])),]->ref.mult
      if(any(isTRUE(ref.mult$poly))) ref.mult[-which(isTRUE(ref.mult$poly)),]->ref.mult
      paste(strsplit(ref.mult$reference[1],"-")[[1]][1],strsplit(ref.mult$bind[1],"-")[[1]][1],sep="-")->ref.mult$reference[-1]
      ref.mult[-1,]$poly<-TRUE
      dat[match(ref.mult$bind,dat$bind),]<-ref.mult
    }
  }

  ### ordering ###
  strsplit(dat$reference,"-")->refs
  dat$ref.tree1<-NA
  dat$ref.tree2<-NA

  lapply(refs,function(x){
    if(length(x[which(!x%in%tree$tip.label)])<1) NA else x[which(!x%in%tree$tip.label)]
  })->ref.tree

  dat[which(is.na(ref.tree)),][order(dat[which(is.na(ref.tree)),]$bind.type),]$ref.tree1<-seq(1,length(which(is.na(ref.tree))))

  if(any(which(!is.na(ref.tree)))){
    dat[which(!is.na(ref.tree)),]->dat.new

    dat.new[,6:7]<-do.call(rbind,ref.tree[-which(is.na(ref.tree))])

    if(any(dat.new[,6]%in%unlist(dat.new$bind.tips)|dat.new[,7]%in%unlist(dat.new$bind.tips))){
      while(nrow(dat.new)>0){
        which(!(dat.new[,6]%in%unlist(dat.new$bind.tips)|dat.new[,7]%in%unlist(dat.new$bind.tips)))->outs
        dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-max(dat$ref.tree1,na.rm=TRUE)+1:length(outs)
        dat.new[-outs,]->dat.new
      }
    }else{
      which(!(dat.new[,6]%in%unlist(dat.new$bind.tips)|dat.new[,7]%in%unlist(dat.new$bind.tips)))->outs
      dat[match(dat.new[outs,1],dat[,1]),]$ref.tree1<-
        max(dat$ref.tree1,na.rm=TRUE)+1:length(which(!dat.new$ref.tree1%in%dat.new[,1]))
    }
  }

  dat[order(dat$ref.tree1),]->dat

  if(!is.null(tree2)){
    dat$MRCAbind<-NA
    sapply(dat[which(dat$bind.type==2),]$bind.tips,function(x) getMRCA(tree2,x))->dat$MRCAbind[which(dat$bind.type==2)]

    if(any(dat$bind.type==2)){
      unlist(sapply(dat[which(dat$bind.type==2),]$MRCAbind,function(x) tree$tip.label[which(tree$tip.label%in%tips(tree2,x))]))->remt
      if(length(remt)>0){
        drop.tip(tree,remt)->tree
        ages[-match(remt,names(ages))]->ages
        warning(paste(paste(remt,collapse=", "),"already on the source tree: removed from the backbone tree"),immediate. = TRUE)
      }
    }
  }

  ### binding ###
  for(k in 1:nrow(dat)){
    if(length(strsplit(dat$reference,"-")[[k]])>1)
      getMRCA(tree,strsplit(dat$reference,"-")[[k]])->where.ref else
        which(tree$tip.label==strsplit(dat$reference,"-")[[k]])->where.ref

    if(where.ref!=(Ntip(tree)+1)) tree$edge.length[which(tree$edge[,2]==where.ref)]->br.len else{
      if(is.null(tree$root.edge)) tree$root.edge<-mean(tree$edge.length)
      tree$root.edge->br.len
    }

    if(dat$bind.type[k]==1){
      if(isTRUE(dat$poly[k])) 0->pos.ref else br.len/2->pos.ref
      bind.tip(tree,dat$bind[k],where.ref,position = pos.ref,edge.length =br.len/2)->tree
    }else{
      extract.clade(tree2,dat$MRCAbind[k])->cla

      if(isTRUE(dat$poly[k])){
        0->pos.ref
        if((max(diag(vcv(cla)))+max(diag(vcv(cla)))/10)>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))
          rescale(cla,"depth",(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]+max(diag(vcv(cla)))/10))->cla
        cla$root.edge<-max(diag(vcv(cla)))/10
      }else {
        if((max(diag(vcv(cla)))+br.len/2)>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),1])){
          rescale(cla,"depth",(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))->cla
          pos.ref<-br.len/2
        }else if((max(diag(vcv(cla)))+br.len/2)>(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2]))
          (max(diag(vcv(cla)))+br.len/2)-(H-nodeHeights(tree)[which(tree$edge[,2]==where.ref),2])->pos.ref else br.len/2->pos.ref
        cla$root.edge<-br.len/2
      }
      bind.tree(tree,cla,where=where.ref,position = pos.ref)->tree
    }

  }

  ### tips calibration ages ###
  if(is.null(tip.ages)){
    rep(0,length(tree$tip.label))->tip.ages
    names(tip.ages)<-tree$tip.label
    tip.ages[match(names(ages),names(tip.ages))]<-ages
  }else{
    if(!all(names(ages)%in%names(tip.ages))) c(ages[which(!names(ages)%in%names(tip.ages))],tip.ages)->tip.ages
    if(!all(tree$tip.label%in%names(tip.ages))){
      rep(0,length(which(!tree$tip.label%in%names(tip.ages))))->tip.add
      names(tip.add)<-tree$tip.label[which(!tree$tip.label%in%names(tip.ages))]
      c(tip.ages,tip.add)->tip.ages
    }
  }

  ### nodes calibration ages ###
  if(!is.null(node.ages)){
    sapply(names(node.ages),function(x){
      getMRCA(tree,unlist(strsplit(x,"-")))
    })->names(node.ages)
  }else{
    node.ages<-c()
  }

  ### age original root ###
  if(!getMRCA(tree,names(ages))%in%names(node.ages)){
    node.ages<-c(node.ages,H)
    names(node.ages)[length(node.ages)]<-getMRCA(tree,names(ages))
  }

  ### time distances inside attached clades ###
  if(any(dat$bind.type==2)){
    unlist(lapply(dat$MRCAbind[which(dat$bind.type==2)],function(x){
      c(x,getDescendants(tree2,x)[which(getDescendants(tree2,x)>Ntip(tree2))])->des
      dist.nodes(tree2)->dn
      max(diag(vcv(tree2)))-dn[which(rownames(dn)==(Ntip(tree2)+1)),match(des,rownames(dn))]->dndes
      names(dndes)<-des
      dndes
    }))->ages.fix

    sapply(names(ages.fix),function(x) getMRCA(tree,tips(tree2,as.numeric(x))))->names(ages.fix)

    if(any(!names(ages.fix)%in%names(node.ages)))
      c(ages.fix[which(!names(ages.fix)%in%names(node.ages))],node.ages)->node.ages
  }

  if(max(diag(vcv(tree)))>H&&(!(Ntip(tree)+1)%in%names(node.ages)))
    warning(paste("Root age not indicated: the tree root arbitrarily set at",round(max(diag(vcv(tree))),2)),immediate.=TRUE)

  scaleTree(tree,node.ages=node.ages,tip.ages =tip.ages,min.branch=min.branch)->tree.final->tree.plot

  if(any(dat$bind.type==2)) lapply(dat$MRCAbind[which(dat$bind.type==2)],function(x)
      c(getMRCA(tree.plot,tips(tree2,x)),getDescendants(tree.plot,getMRCA(tree.plot,tips(tree2,x)))))->cla.plot else cla.plot<-c()

      c(which(tree.plot$tip.label%in%dat$bind[which(dat$bind.type==1)]),unlist(cla.plot))->all.plot
      # colo<-rep(scales::hue_pal()(2)[2],nrow(tree.plot$edge))
      # colo[which(tree.plot$edge[,2]%in%all.plot)]<-scales::hue_pal()(2)[1]
      colo<-rep("gray60",nrow(tree.plot$edge))
      colo[which(tree.plot$edge[,2]%in%all.plot)]<-"red"
      names(colo)<-tree.plot$edge[,2]

      par(mar=c(2,0,1,0))
          plot(tree.plot,edge.color=colo,edge.width=1.5,tip.color=colo[which(tree.plot$edge[,2]<=Ntip(tree.plot))])
title(title)
axisPhylo()
  

  return(tree.final)
}

suppressWarnings(tm(backbone=tree.back,data=dato,title="merged tree"))

## ----echo=7:8,message=FALSE---------------------------------------------------

c(1,2,1.7,1.5,0.8,1.5,0.3,1.2,0.2)->ages.tip
c("a7","a1","t6","a8","a6","a4","a2","a5","a3")->names(ages.tip)
c(2.2,2.9,3.5)->ages.node
c("t2-t1","a1-a8","a7-t10")->names(ages.node)

ages.tip
ages.node

## ----echo=1,results="hide",message=FALSE,fig.dim=c(6,6),out.width="70%",dpi=220,fig.align='center'----
tree.merger(backbone=tree.back,data=dato,tip.ages=ages.tip,node.ages = ages.node,plot=FALSE)

suppressWarnings(tm(backbone=tree.back,data=dato,tip.ages=ages.tip,node.ages = ages.node,title="merged and calibrated tree"))

## ----echo=FALSE,warnings=FALSE,message=FALSE,fig.dim=c(7,3.5),out.width="70%",dpi=220----
set.seed(22)
rtree(8,tip.label = paste("s",1:8,sep=""))->tree.source

par(mfrow=c(1,2),mar=c(2,0,1,0))
plot(tree.back,edge.width=1.8)
title("backbone tree")
axisPhylo()
plot(tree.source,edge.width=1.8)
title("source tree")
axisPhylo()

data.frame(bind=c("a1","s2-s5","s1-s4","s7","a2"),
           reference=c("s3","t10","t3-t9","t3","s2-t7"),
           poly=c(FALSE,FALSE,FALSE,FALSE,FALSE))->dato.clade
knitr::kable(dato.clade,align="c") %>%
  kable_styling(full_width = TRUE, position = "float_right") %>%
add_header_above(c("dato.clade" = 3)) 


## ----echo=1,results = "hide",message=FALSE,fig.dim=c(6,6),out.width="70%",dpi=220,fig.align='center'----
tree.merger(backbone=tree.back,data=dato.clade,source.tree=tree.source,plot=FALSE)

tm(backbone=tree.back,data=dato.clade,source.tree=tree.source,title="merged tree")


## ----echo=c(1:19),message=FALSE,fig.dim=c(6,6),out.width="98%",dpi=220,fig.align='center'----
### load the RRphylo example dataset including Cetaceans tree 
data("DataCetaceans")
DataCetaceans$treecet->treecet # phylogenetic tree

### Select two clades and some species to be removed
tips(treecet,131)->liv.Mysticetes
tips(treecet,193)->Delphininae
c("Aetiocetus_weltoni","Saghacetus_osiris",
  "Zygorhiza_kochii","Ambulocetus_natans",
  "Kentriodon_pernix","Kentriodon_schneideri","Kentriodon_obscurus")->extinct

plot(treecet,show.tip.label = FALSE,no.margin=TRUE)
nodelabels(frame="n",col="blue",font=2,node=c(131,193),text=c("living\nMysticetes","Delphininae"))
tiplabels(frame="circle",bg="red",cex=.3,text=rep("",length(c(liv.Mysticetes,Delphininae,extinct))),
          tip=which(treecet$tip.label%in%c(liv.Mysticetes,Delphininae,extinct)))

### Create the backbone and source trees
drop.tip(treecet,c(liv.Mysticetes[-which(tips(treecet,131)%in%
                                           c("Caperea_marginata","Eubalaena_australis"))],
                   Delphininae[-which(tips(treecet,193)=="Tursiops_aduncus")],extinct))->backtree
drop.tip(treecet,which(!treecet$tip.label%in%
                         c(liv.Mysticetes,Delphininae,extinct)))->sourcetree

### Create the data object
data.frame(bind=c("Balaena_mysticetus-Caperea_marginata",
                  "Aetiocetus_weltoni",
                  "Saghacetus_osiris",
                  "Zygorhiza_kochii",
                  "Ambulocetus_natans",
                  "Kentriodon_pernix",
                  "Kentriodon_schneideri",
                  "Kentriodon_obscurus",
                  "Sousa_chinensis-Delphinus_delphis",
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
                       "Sotalia_fluviatilis",
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

knitr::kable(dato,align="c") %>%
  kable_styling(full_width = TRUE, position = "center") %>%
add_header_above(c("dato" = 3)) 

## ----results="hide",message=FALSE---------------------------------------------
### Merge the backbone and the source trees according to dat without calibrating tip and node ages
tree.merger(backbone = backtree,data=dato,source.tree = sourcetree,plot=FALSE)

### Set tips and nodes calibration ages
c(Aetiocetus_weltoni=28.0,
  Saghacetus_osiris=33.9,
  Zygorhiza_kochii=34.0,
  Ambulocetus_natans=40.4,
  Kentriodon_pernix=15.9,
  Kentriodon_schneideri=11.61,
  Kentriodon_obscurus=13.65)->tipages
c("Ambulocetus_natans-Fucaia_buelli"=52.6,
  "Balaena_mysticetus-Caperea_marginata"=21.5)->nodeages

### Merge the backbone and the source trees and calibrate tips and nodes ages
tree.merger(backbone = backtree,data=dato,source.tree = sourcetree,
            tip.ages=tipages,node.ages=nodeages,plot=FALSE)


## ----echo=c(14:15,32:33), fig.dim=c(6,6), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
library(ape)
library(phytools)
library(geiger)

set.seed(14)
DataFelids$treefel->tree
tree$tip.label<-paste("t",seq(1:Ntip(tree)),sep="")

max(nodeHeights(tree))->H
sample(H-diag(vcv(tree)),8)->sp.ages
sp.ages+tree$edge.length[match(match(names(sp.ages),tree$tip.label),tree$edge[,2])]/2->sp.ages
sp.ages[c(3,7)]<-0

sp.ages
scaleTree(tree,tip.ages=sp.ages)->treeS1

edge.col<-rep("gray60",nrow(tree$edge))
edge.col[which(treeS1$edge[,2]%in%match(names(sp.ages),tree$tip.label))]<-"blue"

par(mfcol=c(2,2),mar=c(0.1,0.1,1,0.1))
plot(tree,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("original",cex.main=1.2)
plot(treeS1,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("species ages rescaled",cex.main=1.2)

set.seed(0)
sample(seq((Ntip(tree)+2),(Nnode(tree)+Ntip(tree))),8)->nods
H-dist.nodes(tree)[(Ntip(tree)+1),nods]->nod.ages

sapply(1:length(nods),function(x) {
H-dist.nodes(tree)[(Ntip(tree)+1),c(getMommy(tree,nods[x])[1],getDescendants(tree,nods[x])[1:2])]->par.ages
nod.ages[x]+((min(abs(nod.ages[x]-par.ages))-0.2)*sample(c(-1,1),1))
})->nod.ages

nod.ages
scaleTree(tree,node.ages=nod.ages)->treeS2
treeS2->treeS1
unlist(lapply(1:length(nods), function(x) c(nods[x],getDescendants(tree,nods[x])[1:2])))->brcol
edge.col<-rep("gray60",nrow(tree$edge))
edge.col[which(treeS1$edge[,2]%in%brcol)]<-"red"

#par(mfrow=c(2,1),mar=c(0.1,0.1,1,0.1))
plot(tree,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("original",cex.main=1.2)
plot(treeS1,edge.color = edge.col,edge.width=1.5,show.tip.label=F)
title("node ages rescaled",cex.main=1.2)

## ----echo=7:8, fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center"----
H-dist.nodes(tree)[(Nnode(tree)+1),91]->sp.ages
names(sp.ages)<-tree$tip.label[1]

H-dist.nodes(tree)[(Nnode(tree)+1),164]->nod.ages
names(nod.ages)<-96

c(sp.ages,nod.ages)
scaleTree(tree,tip.ages = sp.ages,node.ages = nod.ages,min.branch = 1)->treeS

par(mfrow=c(1,2))
plot(tree,edge.color = "gray40",show.tip.label=F,no.margin = TRUE,edge.width=1.5)
plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(plotinfo$xx[1],plotinfo$yy[1],pch=16,col="blue",cex=1.2)
points(plotinfo$xx[91],plotinfo$yy[1],pch=4,col="blue",cex=1.5,lwd=2)
text("target age",x= plotinfo$xx[91]-1,y=plotinfo$yy[1],adj=1,col="blue",cex=0.8)
points(plotinfo$xx[96],plotinfo$yy[96],pch=16,col="red",cex=1.2)
points(plotinfo$xx[164],plotinfo$yy[96],pch=4,col="red",cex=1.5,lwd=2)
text("target age",x= plotinfo$xx[164]-1,y=plotinfo$yy[96],adj=1,col="red",cex=0.8)

edge.col<-rep("gray40",nrow(tree$edge))
edge.col[which(treeS$edge[,2]%in%c(1,getMommy(tree,1)))]<-"blue"
edge.col[which(treeS$edge[,2]%in%c(94,95,96))]<-"red"

plot(treeS,edge.color = edge.col,show.tip.label=F,no.margin = TRUE,edge.width=1.5)
plotinfo <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(plotinfo$xx[1],plotinfo$yy[1],pch=16,col="blue",cex=1.2)
points(plotinfo$xx[96],plotinfo$yy[96],pch=16,col="red",cex=1.2)


## ----echo=c(1:16,26:37,50:56), fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%',fig.align="center"----
# load the RRphylo example dataset including Felids tree
data("DataFelids")
DataFelids$treefel->tree

# get species and nodes ages 
# (meant as distance from the youngest species, that is the Recent in this case)
max(nodeHeights(tree))->H
H-dist.nodes(tree)[(Ntip(tree)+1),(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))]->age.nodes
H-diag(vcv(tree))->age.tips

# apply Pagel's lambda transformation to change node ages only 
geiger::rescale(tree,"lambda",0.8)->tree1

# apply scaleTree to the transformed phylogeny, by setting
# the original ages at nodes as node.ages
scaleTree(tree1,node.ages=age.nodes)->treeS1

par(mfrow=c(1,2),mar=c(1,0.1,1,0.1),mgp=c(3,0.1,0.05))
plot(tree1,edge.color = "black",show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("lambda rescaled",cex.main=1.2)
plot(treeS1,edge.color = "black",show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("scaleTree rescaled",cex.main=1.2)

# change leaf length of 10 sampled species
tree->tree2
set.seed(14)
sample(tree2$tip.label,10)->sam.sp
age.tips[sam.sp]->age.sam
age.sam[which(age.sam>0.1)]<-age.sam[which(age.sam>0.1)]-1.5
age.sam[which(age.sam<0.1)]<-age.sam[which(age.sam<0.1)]+0.2
tree2$edge.length[match(match(sam.sp,tree$tip.label),tree$edge[,2])]<-age.sam

# apply scaleTree to the transformed phylogeny, by setting
# the original ages at sampled tips as tip.ages
scaleTree(tree2,tip.ages=age.tips[sam.sp])->treeS2

edge.col<-rep("black",nrow(tree$edge))
edge.col[which(treeS2$edge[,2]%in%match(sam.sp,tree$tip.label))]<-"red"

par(mfrow=c(1,2),mar=c(1,0.1,1,0.1),mgp=c(3,0.1,0.05))
plot(tree2,edge.color = edge.col,show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("leaves cut",cex.main=1.2)
plot(treeS2,edge.color = edge.col,show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("scaleTree rescaled",cex.main=1.2)

# apply Pagel's kappa transformation to change both species and node ages, 
# including the age at the tree root
geiger::rescale(tree,"kappa",0.5)->tree3

# apply scaleTree to the transformed phylogeny, by setting
# the original ages at nodes as node.ages
scaleTree(tree1,tip.ages = age.tips,node.ages=age.nodes)->treeS3

par(mfrow=c(1,2),mar=c(1,0.1,1,0.1),mgp=c(3,0.1,0.05))
plot(tree3,edge.color = "black",show.tip.label=F)
axisPhylo(tck=-0.02,cex.axis=0.8)
title("kappa rescaled",cex.main=1.2)
plot(treeS3,edge.color = "black",show.tip.label=F)
axis(side=1,at=c(0,4,8,12,16,20,24,28,32),labels=rev(c(0,4,8,12,16,20,24,28,32)),tck=-0.02,cex.axis=0.8)
title("scaleTree rescaled",cex.main=1.2)

## ----echo=FALSE,fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
DataFelids$treefel->tree
max(nodeHeights(tree))->H

par(mfrow=c(1,2),mar=c(1,0.1,1.2,0.1),mgp=c(3,0.1,0.05))
plot(tree,show.tip.label=F)
axis(side=1,at=seq(2,32,5),labels=rev(c(0,5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
abline(v=H-5,col="blue")
title("original",cex.main=1.2)

plot(tree,show.tip.label=F)
axis(side=1,at=seq(2,32,5),labels=rev(c(0,5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
title("original",cex.main=1.2)
points(plotinfo$xx[129],plotinfo$yy[129],pch=16,col="red",cex=1.2)

## ----eval=FALSE---------------------------------------------------------------
#  cutPhylo(tree,age=5,keep.lineage = TRUE)
#  cutPhylo(tree,age=5,keep.lineage = FALSE)

## ----echo=FALSE,fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
par(mfrow=c(1,2),mar=c(1,0.1,1.2,0.1),mgp=c(3,0.1,0.05))
age<-5
{
  distNodes(tree,(Ntip(tree)+1))->dN
  max(nodeHeights(tree))-age->cutT
  dN[,2]->dd
  dd[which(dd>=cutT)]->ddcut
  names(ddcut)->cutter
  
  ### Tips only ###
  tree->tt
  if(all(cutter%in%tree$tip.label)){
    tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]<-
      tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]-(ddcut-cutT)
    #}
  }else{
    ### Tips and nodes ###
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
  
  tt->tt.age
  plot(tt,show.tip.label=FALSE)
  tiplabels(frame="n",col="blue",cex=.6,offset=0.5,
            tip=which(!tt$tip.label%in%tree$tip.label),
            text=tt$tip.label[which(!tt$tip.label%in%tree$tip.label)])
  
  axis(side=1,at=seq(2,27,5),labels=rev(c(5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
  title("cut at 5: keeping lineages",cex.main=1.2)
  
   drop.tip(tt.age,which(!tt.age$tip.label%in%tree$tip.label))->tt.age
   plot(tt.age,show.tip.label=FALSE)
   axis(side=1,at=seq(2,27,5),labels=rev(c(5,10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
   title("cut at 5: removing lineages",cex.main=1.2)
  
}

## ----eval=FALSE---------------------------------------------------------------
#  cutPhylo(tree,node=129,keep.lineage = TRUE)
#  cutPhylo(tree,node=129,keep.lineage = FALSE)

## ----echo=FALSE,fig.dim=c(6,3), message=FALSE, warning=FALSE, dpi=200, out.width='98%'----
node<-129
{
  distNodes(tree,(Ntip(tree)+1))->dN
  dN[match(node,rownames(dN)),2]->cutT
  dN[,2]->dd
  dd[which(dd>=cutT)]->ddcut
  names(ddcut)->cutter
  
  ### Tips only ###
  tree->tt
  if(all(cutter%in%tree$tip.label)){
    tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]<-
      tt$edge.length[match(match(names(ddcut),tt$tip.label),tt$edge[,2])]-(ddcut-cutT)
    #}
  }else{
    ### Tips and nodes ###
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
  
  tt->tt.node
  par(mfrow=c(1,2),mar=c(1,0.1,1.2,0.1),mgp=c(3,0.1,0.05))
  plot(tt,show.tip.label=FALSE)
  tiplabels(frame="n",col="red",cex=.6,offset=0.5,
            tip=which(!tt$tip.label%in%tree$tip.label),
            text=tt$tip.label[which(!tt$tip.label%in%tree$tip.label)])
  
  axis(side=1,at=seq(2,23.5,5),labels=rev(c(10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
  title("cut at node: keeping lineages",cex.main=1.2)
  
   drop.tip(tt.node,which(!tt.node$tip.label%in%tree$tip.label))->tt.node
   plot(tt.node,show.tip.label=FALSE)
   axis(side=1,at=seq(2,23.5,5),labels=rev(c(10,15,20,25,30)),tck=-0.02,cex.axis=0.8)
   title("cut at node: removing lineages",cex.main=1.2)
  
}


## ----warnings=FALSE,message=FALSE,fig.dim=c(6,3),out.width="98%",dpi=220------
 ### load the RRphylo example dataset including Cetaceans tree 
 data("DataCetaceans")
 DataCetaceans$treecet->treecet

 ### Resolve all the polytomies within Cetaceans phylogeny
 fix.poly(treecet,type="resolve")->treecet.fixed
 
 ## Set branch colors
 unlist(sapply(names(which(table(treecet$edge[,1])>2)),function(x) 
   c(x,getDescendants(treecet,as.numeric(x)))))->tocolo
 unlist(sapply(names(which(table(treecet$edge[,1])>2)),function(x) 
   c(getMRCA(treecet.fixed,tips(treecet,x)),
     getDescendants(treecet.fixed,as.numeric(getMRCA(treecet.fixed,tips(treecet,x)))))))->tocolo2
 colo<-rep("gray60",nrow(treecet$edge))
 names(colo)<-treecet$edge[,2]
 colo2<-rep("gray60",nrow(treecet.fixed$edge))
 names(colo2)<-treecet.fixed$edge[,2]
 colo[match(tocolo,names(colo))]<-"red"
 colo2[match(tocolo2,names(colo2))]<-"red"
 
 par(mfrow=c(1,2))
 plot(treecet,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo,edge.width=1.3)
 plot(treecet.fixed,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo2,edge.width=1.3)

 ### Resolve the polytomies pertaining the genus Kentriodon
 fix.poly(treecet,type="resolve",node=221)->treecet.fixed2
 
 ## Set branch colors
 c(221,getDescendants(treecet,as.numeric(221)))->tocolo
 c(getMRCA(treecet.fixed2,tips(treecet,221)),
   getDescendants(treecet.fixed2,as.numeric(getMRCA(treecet.fixed2,tips(treecet,221)))))->tocolo2
 colo<-rep("gray60",nrow(treecet$edge))
 names(colo)<-treecet$edge[,2]
 colo2<-rep("gray60",nrow(treecet.fixed2$edge))
 names(colo2)<-treecet.fixed2$edge[,2]
 colo[match(tocolo,names(colo))]<-"red"
 colo2[match(tocolo2,names(colo2))]<-"red"
 
 
 par(mfrow=c(1,2))
 plot(treecet,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo,edge.width=1.3)
 plot(treecet.fixed2,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo2,edge.width=1.3)

 ### Collapse Delphinidae into a polytomous clade
 fix.poly(treecet,type="collapse",node=179)->treecet.collapsed
 
 # Set branch colors
 c(179,getDescendants(treecet,as.numeric(179)))->tocolo
 c(getMRCA(treecet.collapsed,tips(treecet,179)),
   getDescendants(treecet.collapsed,as.numeric(getMRCA(treecet.collapsed,tips(treecet,179)))))->tocolo2
 colo<-rep("gray60",nrow(treecet$edge))
 names(colo)<-treecet$edge[,2]
 colo2<-rep("gray60",nrow(treecet.collapsed$edge))
 names(colo2)<-treecet.collapsed$edge[,2]
 colo[match(tocolo,names(colo))]<-"red"
 colo2[match(tocolo2,names(colo2))]<-"red"
 
 par(mfrow=c(1,2))
 plot(treecet,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo,edge.width=1.3)
 plot(treecet.collapsed,no.margin=TRUE,show.tip.label=FALSE,edge.color = colo2,edge.width=1.3)


