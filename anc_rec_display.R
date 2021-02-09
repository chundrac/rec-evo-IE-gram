require(rstan)
require(phytools)
require(phangorn)
require(expm)
require(tikzDevice)
require(RColorBrewer)

if (!dir.exists('graphics/')) {
  dir.create('graphics/')
}

#feature indices {1..65}
feat.ind = as.integer(commandArgs(trailingOnly=TRUE)[1])

diacl.data <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)

IE.trees <- read.newick('IE.newick')

diacl.data <- droplevels(diacl.data[IE.trees[[1]]$tip.label,]) #exclude non-IE

feat <- colnames(diacl.data)[feat.ind]
feat.old <- feat
print(feat)

F <- length(levels(as.factor(diacl.data[,feat.ind])))

load(paste('stanfits/',feat,'.Rdata',sep=''))

diacl.data <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1, check.names = FALSE)

diacl.data <- droplevels(diacl.data[IE.trees[[1]]$tip.label,]) #exclude non-IE

feat <- colnames(diacl.data)[feat.ind]

anc.states <- function(t,s,f) {
  tree <- IE.trees[[s]]
  tree <- reorder.phylo(tree,'pruningwise')
  parent <- tree$edge[,1]
  child <- tree$edge[,2]
  Q <- extract(fit.full)$Q[t,,]
  brlens <- tree$edge.length/1000
  B <- length(brlens)
  P <- list() #each branch's transition probability matrix
  #get likelihoods	
  curr.states <- likelihood[tree$tip.label,]
  lambda <- matrix(data=0,nrow=N,ncol=F)
  for (n in 1:T) {
    for (f in 1:F) {
      lambda[n,f] <- log(curr.states[n,f])
    }
  }
  for (b in 1:B) {
    P[[b]] <- expm(Q*brlens[b])
    for (f in 1:F) {
      lambda[parent[b],f] <- lambda[parent[b],f] + log(P[[b]][f,]%*%exp(lambda[child[b],]))
    }
  }
  pi = expm(1000*Q)[1,]
  lambda[parent[B],] <- lambda[parent[B],] + log(pi)
  #normalize likelihoods
  phi <- prop.table(exp(lambda),margin=1)
  #draw state and root and move down the tree
  states.sim <- matrix(data=0,nrow=N,ncol=F)
  states.sim[parent[B],] <- t(rmultinom(1,1,phi[parent[B],]))
  for (b in B:1) {
    states.sim[child[b],] <- t(rmultinom(1,1,(states.sim[parent[b],]%*%P[[b]]) * phi[child[b],]))
  }
  return(states.sim)
  #    return(states.sim[(T+1):N,])
}

#set.seed(1)
#states.list <- list()
#for (i in 1:100000) {
#	print(i)
#    j = sample(1:dim(extract(fit.full)$Q)[1],1)
#	k = sample(1:length(IE.trees),1)
#    states.list[[i]] <- anc.states(j,k,f)
#	states.agg <- prop.table(Reduce('+',states.list),margin=1)
#	colnames(states.agg) <- levels(as.factor(diacl.data[,feat.ind]))
#	saveRDS(states.agg,file=paste('graphics/',feat.old,'.rds',sep=''))
#}




##states.agg <- data.frame(prop.table(Reduce('+',states.list),margin=1))
##states.agg <- data.frame(Reduce('+',states.list))
#states.agg <- prop.table(Reduce('+',states.list),margin=1)
##rownames(states.agg) <- tree$node.label
#colnames(states.agg) <- levels(as.factor(diacl.data[,feat.ind]))
##states.agg$node <- c(1:N)

#saveRDS(states.agg,file=paste('graphics/',feat,'.rds',sep=''))
states.agg=readRDS(file=paste('graphics/',feat.old,'.rds',sep=''))

node.colors <- brewer.pal(n = 8, name = "Set2")
#tikz('state_tree')
cairo_pdf(paste('graphics/',feat.old,'for_display.pdf',sep=''),height=15,width=10)
mcc.tree <- maxCladeCred(IE.trees)


boxlabel<-function(x,y,text,cex=.6,bg="transparent",offset=0){
  w<-strwidth(text)*cex*1.1
  h<-strheight(text)*cex*1.4
  os<-offset*strwidth("W")*cex
  rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
  text(x,y,text,cex=cex,pos=4,offset=offset,font=1,direction="downwards")
}

#par(fg="transparent")
#plotTree(mcc.tree,lwd = .5,fsize=.4,direction="downwards")
plotTree(mcc.tree,lwd = .5,fsize=.4)
#plot(mcc.tree,edge.width = .1,fsize=.4)

#add.scale.bar()
#tiplabels(pie=states.agg[1:T,],cex=.5,piecol=node.colors)
#nodelabels(pie=states.agg[(T+1):N,],cex=.5,piecol=node.colors)
non.advent<-c((T+1):N)[sapply(c((T+1):N),function(x) {return(!any(mcc.tree$edge.length[mcc.tree$edge[,1]==x]==1))})]

#tiplik<-as.matrix(likelihood[mcc.tree$tip.label,]/rowSums(likelihood[mcc.tree$tip.label,]))
#tiplabels(pie=tiplik,node=mcc.tree$tip.label,cex=.5,piecol=node.colors)


mcc.feat <- as.factor(diacl.data[mcc.tree$tip.label,feat.ind])
print(length(levels(mcc.feat)))
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
N<-Ntip(mcc.tree)
par(fg="black")
for(i in 1:Ntip(mcc.tree)) {
  #if (sum(likelihood[mcc.tree$tip.label[i],])==1) {
  #if !(is.na(diacl.data[mcc.tree$tip.label[i],f])) {
  #print(node.colors[as.numeric(mcc.feat[i])])
  #print(as.numeric(mcc.feat[i]) %in% c(1:length(levels(mcc.feat))))
  #if (as.numeric(mcc.feat[i]) %in% c(1:length(levels(mcc.feat)))) {
  #  boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg=node.colors[as.numeric(mcc.feat[i])])
  #}
  #else {
  #	boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg='white')
  #}
  if (is.na(mcc.feat[i])) {
    boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg='white')
  }
  else {
    boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg=node.colors[as.numeric(mcc.feat[i])])
  }
  #}
  #else {
  #    boxlabel(pp$xx[i],pp$yy[i],tree$tip.label[i],bg='gray')
  #}
}

nodelabels(pie=states.agg[non.advent,],node=non.advent,cex=.6,col=NULL,piecol=node.colors)

nodelabels(text=mcc.tree$node.label[c(1:4,21,55,68,69,39)],frame='none',node=T+c(1:4,21,55,68,69,39),adj=c(0,-2))

#odelabels(text=mcc.tree$node.label[c(1:4,21,55,68,69,39)],frame='none',node=T+c(1:4,21,55,68,69,39),adj=c(1,1))

#cladelabels(text=mcc.tree$node.label[non.advent-T],node=non.advent)

#nodelabels(text=mcc.tree$node.label[non.advent],frame='none',pie=states.agg,cex=.4,col=NULL,piecol=node.colors)

#nodelabels(text=mcc.tree$node.label[non.advent-T],frame='none',thermo=states.agg[non.advent,],node=non.advent,cex=.7,col=NULL,piecol=node.colors)

#nodelabels(text=mcc.tree$node.label[non.advent-T],frame='none',pie=states.agg[non.advent,],node=non.advent,cex=.7,col=NULL,piecol=node.colors,label.offset=1)

feat.title <- paste(unlist(strsplit(unlist(strsplit(feat,';'))[1],"|",fixed=T))[1],unlist(strsplit(unlist(strsplit(feat,';'))[1],"|",fixed=T))[2],sep=' — ')
print(feat.title)
#CAREFUL!!
legend('bottomleft',legend=levels(as.factor(diacl.data[,feat.ind])),title=feat.title,fill=node.colors,bty='o',cex=1.5)
#title(gsub(feat,'|','; '))
#title(feat,cex=.25)
#legend('topleft',legend=levels(states$V2),fill=node.colors,bty='n',cex=1)
dev.off()

