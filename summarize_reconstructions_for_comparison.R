require(rstan)
require(phytools)

diacl.data1 <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)
IE.trees <- read.newick('IE.newick')
diacl.data1 <- droplevels(diacl.data1[IE.trees[[1]]$tip.label,])

results <- NULL
for (feat.ind in c(1:65)[-45]) {
	feat1 <- colnames(diacl.data1)[feat.ind]
	load(paste('stanfits/',feat1,'.Rdata',sep=''))
	probs.0 <- colMeans(extract(fit.full)$z)
	load(paste('stanfits_dirprior/',feat1,'.Rdata',sep=''))
	probs.dir <- colMeans(extract(fit.full)$z)
	load(paste('stanfits_unifprior/',feat1,'.Rdata',sep=''))
	probs.unif <- colMeans(extract(fit.full)$z)
	values <- levels(as.factor(diacl.data1[,feat.ind]))
	results <- rbind(results,cbind(rep(feat1,length(probs.0)),values,probs.0,probs.dir,probs.unif))
}

write.table(file='reconstructions_compare.tsv',results,row.names=F,sep='\t',quote=F)