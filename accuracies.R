require(ggplot2)
require(phytools)

filelist <- paste('accuracies',dir('accuracies'),sep='/')

acc <- NULL
for (file in filelist) {
  acc <- rbind(acc,read.csv(file,sep='\t',header=F))
}

#node.list <- read.delim('nodes_ordered.txt',header=F)

#acc$node <- node.list[acc$V2,]

colnames(acc) <- c('feature','node','accuracy')

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

acc$ancestral <- ifelse(acc$node %in% ancestral, 'ancestral', 'non-ancestral')

diacl.data <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)

IE.trees <- read.newick('IE.newick')

diacl.data <- droplevels(diacl.data[IE.trees[[1]]$tip.label,]) #exclude non-IE

keep <- c()
for (i in 1:nrow(acc)) {
  if (!is.na(diacl.data[acc[i,]$node,acc[i,]$feature])) {
    keep <- c(keep,i)
  }
}

acc <- acc[keep,]

#colnames(acc) <- c('feature','node.num','accuracy','node','ancestral')

print(median(acc$accuracy))
#print(mean(acc$accuracy))
print(median(acc[acc$ancestral=='ancestral',]$accuracy))
print(median(acc[acc$ancestral=='non-ancestral',]$accuracy))

wilcox.test(acc[acc$ancestral=='ancestral',]$accuracy,acc[acc$ancestral=='non-ancestral',]$accuracy)

require(lme4)
require(MASS)

dropterm(lmer( log(accuracy) ~ ancestral + (1|feature) + (1|node), acc),test='Chisq')

pdf('LOO_accuracy.pdf')
ggplot(data=acc,aes(y=accuracy,x=ancestral))+geom_violin() + stat_summary(fun.y=median, geom="point", size=2) + labs(x = "language type", y = "leave-one-out accuracy")
dev.off()

#ggplot(data=acc,aes(x=accuracy)) + #geom_histogram(data=subset(acc,ancestral=='ancestral'),fill='red',alpha=.2,bins=100) + #geom_histogram(data=subset(acc,ancestral=='non-ancestral'),fill='blue',alpha=.2,bins=100)