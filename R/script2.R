args = commandArgs( trailingOnly = TRUE )
BASE_DIR=args[1]


fileDir = paste0(BASE_DIR,"/Analysis/")
fpfilterDir = paste0(BASE_DIR, "/fpfilter/")

load(paste0(fileDir,"vcf.rda"))

setwd(fpfilterDir)
failfiles <- dir(path=".",pattern = "fpfilter.fail")

print(failfiles)

FPfail <- NULL
for(i in 1:length(failfiles)){
  print(i)
  patient <- unlist(strsplit(failfiles[i],".",fixed=TRUE))[1]
  tmp <- readLines(failfiles[i])
  tmp <- sapply(tmp,function(x){
    Tmp <- unlist(strsplit(x,"\t"))
    length(Tmp) <- 9
    Tmp
  })
  tmp <- t(tmp)
  failID <- paste0(patient,"-",tmp[,1],":",tmp[,2],"-",tmp[,4])
  tmp <- cbind(failID,tmp[,9])
  rownames(tmp) <- failID
  FPfail <- rbind(FPfail,tmp)
}

passfiles <- dir(path=".",pattern = "fpfilter.pass")
FPpass <- NULL
for(i in 1:length(passfiles)){
  print(i)
  patient <- unlist(strsplit(passfiles[i],".",fixed=TRUE))[1]
  tmp <- read.delim(passfiles[i],header=FALSE,as.is=TRUE)
  tmp <- paste0(patient,"-",tmp[,1],":",tmp[,2],"-",tmp[,4])
  FPpass <- c(FPpass,tmp)
}

length(FPpass) #8967
nrow(FPfail) #18886

length(FPpass)+nrow(FPfail) #27853

all <- c(FPpass,FPfail[,1])
table(vcf$fpFid%in%all)
#View(vcf[!vcf$fpFid%in%all,])

vcf <- data.frame(vcf,FPfilter=".",FPfail=".",stringsAsFactors=FALSE)
vcf[vcf$fpFid%in%FPpass,"FPfilter"] <- "PASS"
vcf[vcf$fpFid%in%rownames(FPfail),"FPfilter"] <- "FAIL"
vcf[vcf$FPfilter=="FAIL","FPfail"] <- FPfail[vcf[vcf$FPfilter=="FAIL","fpFid"],2]

vcf[is.na(vcf$FPfail),"FPfail"] <- ""
vcf[vcf$FPfail=="","FPfilter"] <- "no_readcounts" # alt absent in the bam readcount file
vcf[vcf$FPfail=="","FPfail"] <- "."
vcf[vcf$FPfilter==".","FPfilter"] <- "no_readcounts"# alt = 0 in the bam readcount file
table(vcf$FPfilter)
# FAIL no_readcounts          PASS 
# 15577          3350          8967  
rm(FPpass,FPfail,all)

save(vcf,file = paste(fileDir,"vcf.rda",sep=""))

setwd(fileDir)
ann.in <- vcf[,c(1,2,2,3,4)]
colnames(ann.in) <- c("Chr","Start","End","Ref","Alt")

cond <- nchar(ann.in[,"Ref"])>1&nchar(ann.in[,"Alt"])==1
ann.in[cond,"Start"] <- as.numeric(ann.in[cond,"Start"])+1
ann.in[cond,"Ref"] <- substring(ann.in[cond,"Ref"],2,nchar(ann.in[cond,"Ref"]))
ann.in[cond,"Alt"] <- "-"
ann.in[cond,"End"] <- as.numeric(ann.in[cond,"Start"])+nchar(ann.in[cond,"Ref"])-1

cond <- nchar(ann.in[,"Alt"])>1&nchar(ann.in[,"Ref"])==1
ann.in[cond,"Ref"] <- "-"
ann.in[cond,"Alt"] <- substring(ann.in[cond,"Alt"],2,nchar(ann.in[cond,"Alt"]))

cond <- nchar(ann.in[,"Ref"])>1&nchar(ann.in[,"Alt"])>1
ann.in[cond,"End"] <- as.numeric(ann.in[cond,"Start"])+nchar(ann.in[cond,"Ref"])-1

write.table(ann.in,file=paste(fileDir,"ann.in",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
