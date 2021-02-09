require(rstan)
require(phytools)
require(expm)

diacl.data1 <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)
diacl.data2 <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1, check.names = FALSE)
IE.trees <- read.newick('IE.newick')
diacl.data1 <- droplevels(diacl.data1[IE.trees[[1]]$tip.label,])
diacl.data2 <- droplevels(diacl.data2[IE.trees[[1]]$tip.label,])


results <- NULL
rate.pairs <- NULL
all.entry.rates <- NULL
all.exit.rates <- NULL
for (feat.ind in c(1:65)[-45]) {
	feat1 <- colnames(diacl.data1)[feat.ind]
	feat2 <- colnames(diacl.data2)[feat.ind]
	load(paste('stanfits/',feat1,'.Rdata',sep=''))
	probs <- colMeans(extract(fit.full)$z)
	rates <- colMeans(extract(fit.full)$Q)
	stat.probs <- expm(rates*10000)[1,]
	values <- levels(as.factor(diacl.data1[,feat.ind]))
	results <- rbind(results,cbind(rep(feat2,length(probs)),values,probs))
	for (i in 1:length(values)) {
		for (j in 1:length(values)) {
			if (i != j) {
				rate.pairs <- rbind(rate.pairs,c(values[i],values[j],rates[i,j]))
			}
		}
	}
	entry.rates <- c()
	exit.rates <- c()
	for (i in 1:length(values)) {
		i.rates <- rates[-i,i]
		i.probs <- stat.probs[-i]
		entry.rates <- c(entry.rates, sum(i.probs*i.rates)/sum(i.probs))
		exit.rates <- c(exit.rates, -rates[i,i])
	}
	all.entry.rates <- rbind(all.entry.rates, cbind(rep(feat2,length(values)),values,entry.rates))
	all.exit.rates <- rbind(all.exit.rates, cbind(rep(feat2,length(values)),values,exit.rates))
}

write.table(file='reconstructions.tsv',results,row.names=F,sep='\t',quote=F)
write.table(file='rates.tsv',rate.pairs,row.names=F,sep='\t',quote=F)
write.table(file='entry_rates.tsv',all.entry.rates,row.names=F,sep='\t',quote=F)
write.table(file='exit_rates.tsv',all.exit.rates,row.names=F,sep='\t',quote=F)