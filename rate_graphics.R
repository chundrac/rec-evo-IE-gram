require(ggplot2)
require(gridExtra)
require(ggrepel)

entry <- read.csv('entry_rates.tsv',sep='\t')
exit <- read.csv('exit_rates.tsv',sep='\t')
reconstruction <- read.csv('reconstructions.tsv',sep='\t')

entry.exit <- merge(entry,exit,by='values')
entry.exit <- entry.exit[,c(2,1,3,5)]
entry.exit$grid <- unname(sapply(as.character(entry.exit[,1]), function(x) {return(unlist(strsplit(x,"|",fixed=T))[1])}))
entry.exit[entry.exit$grid=='ALIGNMENT',]$grid <- 'Alignment'
entry.exit[entry.exit$grid=='TENSE',]$grid <- 'Tense'

#colnames(entry.exit) <- c('fullfeat','values','Gain','Loss','Type')

entry.exit <- merge(entry.exit,reconstruction,by='values')

colnames(entry.exit) <- c('fullfeat','values','Gain','Loss','Type','X','Probability')

#pdf('rate_gain_loss.pdf',width=12,height=8)
#ggplot(data=entry.exit) + geom_point(aes(x=log(1/Gain),y=log(1/Loss),alpha=Probability,shape=Type,color=Type)) + labs(x='average time before gain (log 1000 years)',y='average time before loss (log 1000 years)') + scale_shape_manual(values=c(15:19))
#dev.off()

pdf('rate_gain_loss.pdf',width=12,height=8)
ggplot(data=entry.exit) + 
  geom_vline(xintercept = median(log(entry.exit$Gain)),#color='gray',
             linetype="dashed") + 
  geom_hline(yintercept = median(log(entry.exit$Loss)),#color='gray',
             linetype="dashed") +
  geom_point(aes(x=log(Gain),y=log(Loss),size=Probability,shape=Type,fill=Type),alpha=.5) + 
  geom_point(aes(x=log(Gain),y=log(Loss),size=Probability,shape=Type),fill=NA,color='black') + 
  scale_shape_manual(values=c(21:25)) + 
  theme(text = element_text(size=15)) + 
  labs(x='log gain rate (per 1000 years)',y='log loss rate (per 1000 years)')
  # + 
  #geom_text_repel(aes(x=log(Gain),y=log(Loss),label=fullfeat))
dev.off()
