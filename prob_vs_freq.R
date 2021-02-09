require(phytools)
require(ggplot2)

diacl.data <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)

IE.trees <- read.newick('IE.newick')

diacl.data <- droplevels(diacl.data[IE.trees[[1]]$tip.label,]) #exclude non-IE

rec <- read.csv('reconstructions.tsv',header=T,sep='\t')

#probs <- NULL

#for (i in c(1:65)[-45]) {
#	#probs <- c(probs,unname(prop.table(xtabs(~na.omit(diacl.data[,i])))))
#	probs <- rbind(probs,cbind(names(prop.table(xtabs(~na.omit(diacl.data[,i])))),unname(prop.table(xtabs(~na.omit(diacl.data[,i]))))))
#}

#rownames(probs) <- probs[,1]
#probs <- probs[,2]

#probs <- probs[rec$value,]

ancestral <- c(
  'Classical_Greek',
  'Old_Norse',
  'Sanskrit',
  'Latin',
  'Middle_Greek',
  'Old_Dutch',
  'Old_English',
  'Old_Frisian',
  'Old_High_German',
  'Old_Saxon',
  'Old_Irish',
  'Old_Russian',
  'Old_Swedish',
  'Middle_Dutch',
  'Middle_English',
  'Middle_High_German',
  'Middle_Low_German',
  'Old_Persian',
  'Middle_Breton',
  'Middle_Welsh',
  'Middle_Irish',
  'Middle_Persian',
  'Old_Italian',
  'Old_French',
  'Old_Portuguese',
  'Old_Provençal',
  'Old_Spanish'
)

names.p <- c()
probs <- c()

for (i in c(1:65)[-45]) {
	#probs <- c(probs,unname(prop.table(xtabs(~na.omit(diacl.data[,i])))))
	probs <- c(probs,unname(prop.table(xtabs(~na.omit(diacl.data[,i])))))
	names.p <- c(names.p,names(prop.table(xtabs(~na.omit(diacl.data[,i])))))
}

names(probs) <- names.p

rec$freq <- probs[as.character(rec$values)]

rec$maxprob <- 1/xtabs( ~ rec$X)[rec$X]

#plot(probs ~ freq, rec[rec$probs > rec$maxprob,])

pdf('prob_v_freq.pdf')

#ggplot(rec[rec$probs > rec$maxprob,]) + geom_point(aes(x=freq,y=probs)) + labs(x='Relative frequency',y='Reconstruction probability') + theme(text = element_text(size=25))

ggplot() + geom_point(data = rec[rec$probs <= rec$maxprob,], aes(x=freq,y=probs), col='gray') + 
  geom_point(data = rec[rec$probs > rec$maxprob,], aes(x=freq,y=probs), col='black') + 
  geom_text(data = rec[rec$freq > .75 & rec$probs < .3,], aes(x=freq,y=probs,label=values), nudge_x = .1, col='gray25', size=5) + 
  geom_text(data = rec[rec$freq < .3 & rec$probs > .75,], aes(x=freq,y=probs,label=values), nudge_x = -.1, col='gray25', size=5) + 
  labs(x='Relative frequency',y='Reconstruction probability') + theme(text = element_text(size=25))

dev.off()

with(rec[rec$probs > rec$maxprob,],cor.test(freq,probs,method='spearman'))

#entry.exit[log(entry.exit$Gain) <= median(log(entry.exit$Gain)) & log(entry.exit$Loss) <= median(log(entry.exit$Loss)),]$fullfeat



#entry.exit[log(entry.exit$Gain) <= median(log(entry.exit$Gain)) & log(entry.exit$Loss) <= median(log(entry.exit$Loss)),c('fullfeat','Probability')]



names.p <- c()
probs <- c()

for (i in c(1:65)[-45]) {
  #probs <- c(probs,unname(prop.table(xtabs(~na.omit(diacl.data[,i])))))
  probs <- c(probs,unname(prop.table(xtabs(~na.omit(diacl.data[ancestral,i])))))
  names.p <- c(names.p,names(prop.table(xtabs(~na.omit(diacl.data[ancestral,i])))))
}

names(probs) <- names.p

rec$freq <- probs[as.character(rec$values)]

rec$maxprob <- 1/xtabs( ~ rec$X)[rec$X]

#plot(probs ~ freq, rec[rec$probs > rec$maxprob,])

pdf('prob_v_freq_ancestral.pdf')

#ggplot(rec[rec$probs > rec$maxprob,]) + geom_point(aes(x=freq,y=probs)) + labs(x='Relative frequency',y='Reconstruction probability') + theme(text = element_text(size=25))

ggplot() + geom_point(data = rec[rec$probs <= rec$maxprob,], aes(x=freq,y=probs), col='gray') + geom_point(data = rec[rec$probs > rec$maxprob,], aes(x=freq,y=probs), col='black') + labs(x='Relative frequency',y='Reconstruction probability') + theme(text = element_text(size=25))


dev.off()

with(rec[rec$probs > rec$maxprob,],cor.test(freq,probs,method='spearman'))

#entry.exit[log(entry.exit$Gain) <= median(log(entry.exit$Gain)) & log(entry.exit$Loss) <= median(log(entry.exit$Loss)),]$fullfeat



#entry.exit[log(entry.exit$Gain) <= median(log(entry.exit$Gain)) & log(entry.exit$Loss) <= median(log(entry.exit$Loss)),c('fullfeat','Probability')]