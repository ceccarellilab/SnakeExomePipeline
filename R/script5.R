args = commandArgs( trailingOnly = TRUE )
MER_SCORE=args[1]
FLAGS_FILE=args[2]
MASK_SCORE=args[3]

bed <- read.table(file='Finalmask.bed', sep='\t')
load(file='Final.RData')

score <- read.delim(MASK_SCORE,header=FALSE,as.is=TRUE)
id <- paste0(score[,1],":",score[,2],"-",score[,3])
table(duplicated(id))
score <- score[!duplicated(id),]
id <- paste0(score[,1],":",score[,2],"-",score[,3])
rownames(score) <- id
ID <- paste0(bed[,1],":",bed[,2],"-",bed[,3])
table(ID%in%id)
Score <- score[ID,7]
Final <- data.frame(MaskedFilter=Score,Final,stringsAsFactors=FALSE)
Final[Final[,"MaskedFilter"]=="pilot","MaskedFilter"] <- "pass"
Final[Final[,"MaskedFilter"]==".","MaskedFilter"] <- "fail"
table(Final$MaskedFilter)


# FLAGS genes annotation (https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0064-y#Sec29)

Final <- data.frame(FLAGS=".",FLAGSscore=0,Final,stringsAsFactors = FALSE)
FLAGS <- read.delim(FLAGS_FILE,sep="\t",as.is=TRUE,header = FALSE)
# quantile(FLAGS[,2])
# # 0%  25%  50%  75% 100% 
# # 1   15   28   49 2659 
# quantile(FLAGS[,2],0.99) #195
# View(FLAGS[FLAGS[,2]>195,])
tmp <- FLAGS$V1
# sum(duplicated(tmp))
# tmp[duplicated(tmp)]
# # "01-Mar" "02-Mar"
# View(FLAGS[order(tmp),])
FLAGS <- FLAGS[!duplicated(tmp),]
rownames(FLAGS) <- FLAGS$V1
genes <- sapply(Final[,"Gene.refGene"],function(x) unlist(strsplit(x,","))[1])
Final$FLAGSscore <- FLAGS[genes,2]
table(is.na(Final$FLAGSscore))
Final$FLAGSscore[is.na(Final$FLAGSscore)] <- 0
# View(Final[order(Final[,1],decreasing=TRUE),])
sum(Final[,"FLAGSscore"]>195) #2606
Final[Final[,"FLAGSscore"]>195,"FLAGS"] <- "reject"
table(Final$FLAGS)
# . reject 
# 25277   2606
save(Final,file="Final.tmp.rda")

bed <- Final[,c("Chr","Start","End")]
write.table(bed,file="SomaticVariants.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
save(Final,file='Final.RData')