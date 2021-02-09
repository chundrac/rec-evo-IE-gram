rec <- read.csv('reconstructions_compare.tsv',header=T,sep='\t')

print('uniform')


unif <- NULL
for (level in levels(rec$X)) {
  subframe <- rec[rec$X==level,]
  if (!(which(subframe$probs.0==max(subframe$probs.0)) 
        %in% 
        which(subframe$probs.unif==max(subframe$probs.unif)))) {
    #print(level)
    #print(subframe[,c('values','probs.0','probs.unif')])
    f1 <- subframe[which(subframe$probs.0==max(subframe$probs.0))[1],]$values
    p1 <- subframe[which(subframe$probs.0==max(subframe$probs.0))[1],]$probs.0
    f2 <- subframe[which(subframe$probs.unif==max(subframe$probs.unif))[1],]$values
    p2 <- subframe[which(subframe$probs.unif==max(subframe$probs.unif))[1],]$probs.unif
    unif <- rbind(unif,c(as.character(f1),p1,as.character(f2),p2))
  }
}

colnames(unif) <- c('sf','sp','uf','up')
  
dir <- NULL
for (level in levels(rec$X)) {
  subframe <- rec[rec$X==level,]
  if (!(which(subframe$probs.0==max(subframe$probs.0)) 
        %in% 
        which(subframe$probs.dir==max(subframe$probs.dir)))) {
    #print(level)
    #print(subframe[,c('values','probs.0','probs.unif')])
    f1 <- subframe[which(subframe$probs.0==max(subframe$probs.0))[1],]$values
    p1 <- subframe[which(subframe$probs.0==max(subframe$probs.0))[1],]$probs.0
    f2 <- subframe[which(subframe$probs.dir==max(subframe$probs.dir))[1],]$values
    p2 <- subframe[which(subframe$probs.dir==max(subframe$probs.dir))[1],]$probs.dir
    dir <- rbind(dir,c(as.character(f1),p1,as.character(f2),p2))
  }
}

colnames(dir) <- c('sf','sp','df','dp')


merged <- merge(unif,dir,by='sf',all=T)


merged <- merged[,-5]

colnames(merged)[2] <- 'sp'

write.table(file='prior_disagreement.csv',merged,quote=F,row.names = F,sep='\t')