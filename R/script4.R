args = commandArgs( trailingOnly = TRUE )
MER_SCORE=args[1]

bed <- read.table(file='Final.bed', sep='\t')
load(file='Final.RData')

score <- read.delim(MER_SCORE ,header=FALSE,as.is=TRUE)
sum(score[,1]=='')
# 0
# score <- score[!score[,1]=='',c(1,2,3,7)]
score <- score[,c(1,2,3,7)]
id <- paste0(score[,1],"-",score[,2])
sum(duplicated(id))
score <- score[!duplicated(id),]
id <- paste0(score[,1],"-",score[,2])
rownames(score) <- id
ID <- paste0(bed[,1],"-",bed[,2])
Score <- score[ID,4]
Final <- cbind(Score,Final)
sum(is.na(Score))
#0
# View(Final[is.na(Score),]) #strange contigs (ie: chr1_gl000192_random)
toRemove <- which(as.numeric(Final[,'Score'])<0.3)
Final[toRemove,"fiftymer"] <- "reject"
table(Final$fiftymer)
# . reject 

bed <- Final[,c("Chr","Start","End")]
write.table(bed,file='Finalmask.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
save(Final,file='Final.RData')