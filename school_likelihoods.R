require(ggplot2)

#kl <- read.csv('kldivs.tsv',header=F,sep='\t')
#wide <- reshape(kl,idvar='V2',timevar='V1',direction='wide')

liks <- read.csv('likelihoods.tsv',header=F,sep='\t')

wide <- reshape(liks,idvar='V2',timevar='V1',direction='wide')
wide <- na.omit(wide)

colnames(liks) <- c('Model','Feature','Likelihood')

#pairwise.wilcox.test(c(wide$V3.Canonical,wide$V3.Isolating,wide$V3.Active-Stative),c(rep('c',46),rep('i',46),rep('a',46)),paired=T,'bonferroni')
pairwise.wilcox.test(c(wide[,2],wide[,3],wide[,4]),c(rep('c',46),rep('i',46),rep('a',46)),paired=T,'bonferroni')


pdf('school_likelihoods.pdf')
ggplot(data=liks,aes(x=Model,y=Likelihood)) + geom_violin() + #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange") #
				 geom_boxplot(width=0.1)

dev.off()