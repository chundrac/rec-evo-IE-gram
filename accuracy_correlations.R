require(phytools)

filelist <- paste('accuracies',dir('accuracies'),sep='/')

acc <- NULL
for (file in filelist) {
  acc <- rbind(acc,read.csv(file,sep='\t',header=F))
}

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

diacl.data <- read.csv('diacl_qualitative_wide.tsv',sep='\t',row.names=1)

IE.trees <- read.newick('IE.newick')

diacl.data <- droplevels(diacl.data[IE.trees[[1]]$tip.label,]) #exclude non-IE

tree <- IE.trees[[1]]

ancestral.nodes <- NULL
non.ancestral.nodes <- NULL
for (i in 1:nrow(acc)) {
	feat <- acc[i,1]
	lang <- acc[i,2]
	p_acc <- acc[i,3]
	true.feat <- diacl.data[as.character(lang),as.character(feat)]
	if (!is.na(true.feat)) {
		if (lang %in% ancestral) {
			langnode <- paste(lang,'advent',sep='_')
			nodes <- getDescendants(tree,length(tree$tip.label)+which(tree$node.label==langnode))
			desc <- nodes[nodes<=length(tree$tip.label)]
			desc <- tree$tip.label[desc]
			desc <- desc[which(desc!=lang)]
			prop <- prop.table(xtabs( ~ diacl.data[desc,feat]))[feat]
			ancestral.nodes <- rbind(ancestral.nodes,c(feat,lang,p_acc,prop))
		}
		else {
			prop <- prop.table(xtabs( ~ diacl.data[-which(rownames(diacl.data)==lang),feat]))[feat]
			non.ancestral.nodes <- rbind(non.ancestral.nodes,c(feat,lang,p_acc,prop))
		}
	}
}


cor.test(ancestral.nodes[,3],ancestral.nodes[,4],method='spearman')

cor.test(non.ancestral.nodes[,3],non.ancestral.nodes[,4],method='spearman')