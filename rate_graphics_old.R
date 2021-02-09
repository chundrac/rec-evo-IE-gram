require(ggplot2)
require(gridExtra)

entry <- read.csv('entry_rates.tsv',sep='\t')
exit <- read.csv('exit_rates.tsv',sep='\t')
reconstruction <- read.csv('reconstructions.tsv',sep='\t')

entry.exit <- merge(entry,exit,by='values')
entry.exit <- entry.exit[,c(2,1,3,5)]
entry.exit$grid <- unname(sapply(as.character(entry.exit[,1]), function(x) {return(unlist(strsplit(x,"|",fixed=T))[1])}))
entry.exit[entry.exit$grid=='ALIGNMENT',]$grid <- 'Alignment'
entry.exit[entry.exit$grid=='TENSE',]$grid <- 'Tense'

colnames(entry.exit) <- c('fullfeat','values','Gain','Loss','Type')

pdf('rate_gain_loss.pdf',width=12,height=8)
ggplot(data=entry.exit) + geom_point(aes(x=log(1/Gain),y=log(1/Loss),shape=Type,color=Type)) + labs(x='average time before gain (log 1000 years)',y='average time before loss (log 1000 years)') + scale_shape_manual(values=c(15:19))
dev.off()


entry.exit <- merge(entry.exit,reconstruction,by='values')

#pdf('rate_prob.pdf',width=12,height=8)

#par(mfrow=c(1,2))
g1 <- ggplot(data=entry.exit) + geom_point(aes(x=log(1/Gain),y=probs)) + labs(x='average time before gain (log 1000 years)',y='probability of reconstruction')# + scale_shape_manual(values=c(15:19)) + theme(legend.title = element_blank())
#dev.off()


#pdf('rate_loss_prob.pdf',width=12,height=8)
g2 <- ggplot(data=entry.exit) + geom_point(aes(x=log(1/Loss),y=probs)) + labs(x='average time before loss (log 1000 years)',y='probability of reconstruction')# + scale_shape_manual(values=c(15:19))
#dev.off()

pdf('rate_prob.pdf',width=12,height=8)
grid.arrange(
  g1,
  g2,
  nrow = 1)
  dev.off()