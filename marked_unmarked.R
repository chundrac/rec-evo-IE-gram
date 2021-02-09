unmarked = "Pronoun, Present progressive: Nominative-accusative
Pronoun, Present progressive: No marking
Pronoun, Simple past: Nominative-accusative
Pronoun, Simple past: No marking
Pronoun, Simple past: Ergative
Pronoun, Simple past: Tripartite
Agglutination for case (pronouns)
Agglutination for number (pronouns)
No agglutination for case (pronouns)
No agglutination for number (pronouns)
Difference A and O (pronouns)
No peripheral cases (pronouns)
Difference, O and Dative (pronouns)
Vocative (pronouns)
No Vocative (pronouns)
More than 7 pronominal cases (pronouns)
Fewer than 7 pronominal cases (pronouns)
Synthetic Present progressive
No synthetic Present progressive
Present progressive by auxiliary
No Present progressive by auxiliary
Pronoun, Present progressive: Nominative-accusative
Pronoun, Present progressive: No marking
Noun, Present progressive: Nominative-accusative
Noun, Present progressive: No marking
Noun, Present progressive: Tripartite
Verb, Present progressive: Nominative-Accusative
Verb, Present progressive: No marking
Verb, Present progressive: Tripartite
Present progressive: Syncretic A Agreement
Present progressive: No A Agreement
Present progressive: Gender A Agreement
Present progressive: Full A Agreement
Present progressive: Full and Gender A Agreement
Present progressive: Syncretic Dative Agreement
Present progressive: No Dative Agreement
Present progressive: Full Dative Agreement
Present progressive: Syncretic O Agreement
Present progressive: No O Agreement
Reflexive with Agent
Reflexive not with Agent
Case difference A and O
Difference A and O (pronouns)
No difference A and O (pronouns)
Masculine/feminine distinction
No masculine/feminine distinction"

marked = "Noun, Present progressive: Nominative-accusative
Noun, Present progressive: No marking
Noun, Simple past: Nominative-accusative
Noun, Simple past: No marking
Noun, Simple past: Ergative
Noun, Simple past: Tripartite
Agglutination for case
Agglutination for number
No agglutination for case
No agglutination for number
Case difference A and O
Peripheral cases
No peripheral cases
Vocative
No Vocative
More than 7 cases
Not more than 7 cases
Synthetic Future
No synthetic Future
Future by auxiliary
No Future by auxiliary
Pronoun, Simple past: Nominative-accusative
Pronoun, Simple past: No marking
Noun, Simple past: Nominative-accusative
Noun, Simple past: No marking
Noun, Simple past: Tripartite
Verb, Simple past: Nominative-accusative
Verb, Simple past: No marking
Verb, Simple past: Tripartite
Simple past: Syncretic A Agreement
Simple past: No A Agreement
Simple past: Gender A Agreement
Simple past: Full A Agreement
Simple past: Full and Gender A Agreement
Simple past: Syncretic Dative Agreement
Simple past: No Dative Agreement
Simple past: Full Dative Agreement
Simple past: Syncretic O Agreement
Simple past: No O Agreement
Reflexive with Object
Reflexive not with Object
Genitive and dative
Difference, O and Dative (pronouns)
No difference O and Dative (pronouns)
Neuter gender
No neuter gender"

unmarked <- read.delim(textConnection(unmarked),header=F)
marked <- read.delim(textConnection(marked),header=F)

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

#levels(unmarked$V1) <- levels(entry.exit$fullfeat)
#levels(marked$V1) <- levels(entry.exit$fullfeat)

#unmarked$V1 <- factor(unmarked$V1,levels=entry.exit$fullfeat)
#marked$V1 <- factor(marked$V1,levels=entry.exit$fullfeat)

unmarked$V1 <- as.character(unmarked$V1)
marked$V1 <- as.character(marked$V1)
entry.exit$fullfeat <- as.character(entry.exit$fullfeat)


marked.full <- NULL
for (i in 1:nrow(marked)) {
  marked.full <- rbind(marked.full,entry.exit[entry.exit$fullfeat==marked[i,],][1,]) }

unmarked.full <- NULL
for (i in 1:nrow(unmarked)) {
  unmarked.full <- rbind(unmarked.full,entry.exit[entry.exit$fullfeat==unmarked[i,],][1,]) }

#wilcox.test(marked.full$Gain,unmarked.full$Gain,paired=T)
#wilcox.test(marked.full$Loss,unmarked.full$Loss,paired=T)

#for (type in levels(as.factor(marked.full$Type))) {
#  print(wilcox.test(marked.full[marked.full$Type==type,]$Loss,unmarked.full[unmarked.full$Type==type,]$Loss,paired=T))
#  print(wilcox.test(marked.full[marked.full$Type==type,]$Gain,unmarked.full[unmarked.full$Type==type,]$Gain,paired=T))
#  }

require(ggplot2)
#full.df <- data.frame('Unmarked'=unmarked.full$Loss,'Marked'=marked.full$Loss)
full.df <- data.frame('Hierarchical.level'=c(rep('Unmarked',nrow(unmarked.full)),rep('Marked',nrow(unmarked.full))),'Loss.rate'=c(unmarked.full$Loss,marked.full$Loss))

pdf('marked_unmarked_loss.pdf')
ggplot(data=full.df,aes(x=Hierarchical.level,y=Loss.rate)) + geom_violin() + stat_summary(fun.y=median, geom="point", size=2) + labs(x='Hierarchical level',y='Loss rate')
dev.off()

