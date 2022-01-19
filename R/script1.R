args = commandArgs( trailingOnly = TRUE )
BASE_DIR=args[1]
WORK_DIR=args[2]
VCF_DIR=args[3]
VCF_FILE=args[4]
NUM_FORKS=args[5]

fileDir = paste0(BASE_DIR,"/Analysis/")
fpfilterDir = paste0(BASE_DIR, "/fpfilter/")

print(fileDir)
print(WORK_DIR)

library(vcfR)
library(doParallel)
setwd("/storage/gluster/vol1/data/PUBLIC/SCAMBIO/ABT414_WES_Analysis/")
dir.create(fileDir)
setwd(WORK_DIR)

# VCF-------------------
filenames <- dir(VCF_DIR,VCF_FILE)

print(filenames)
#vcfDir <- "ensemble/vcf/" MY
#filenames <- dir(vcfDir,"ensemble.snpEff.vcf.gz")

#filenames <- filenames[!grepl(".tbi",filenames)]

vcfRead <- function(filename,VCF_DIR){
  tumorID <- sub(paste("-",VCF_FILE,sep=""),"",filename) #"-ensemble.snpEff.vcf.gz"
  tmp <- vcfR::read.vcfR(paste0(VCF_DIR,filename), verbose=FALSE)
  tmpVcf <- tmp@fix
  tmpGt <- tmp@gt
  if(length(colnames(tmpGt))==3){
    if(!grepl("_N\\b",colnames(tmpGt)[3])){ #CHECK THE NAME OF TUMOR AND NORMAL SAMPLES IN VCF!
      tmpGt <- tmpGt[,c(1,3,2)]
    }
  }
  info <- tmpVcf[,"INFO"]
  func <- sapply(info,function(x){
    temp <- unlist(strsplit(x,","))
    paste0(unique(na.omit(unlist(sapply(temp,function(x) unlist(strsplit(x,"|",fixed=TRUE))[2])))),collapse=",")
  })
  func <- unname(func)
  callers <- sapply(info,function(x){
    temp <- unlist(strsplit(x,";"))
    sub("CALLERS=","",temp[grepl("CALLERS=",temp)])
  })
  callers <- unname(callers)
  id <- paste0(tumorID,"-",tmpVcf[,1],":",tmpVcf[,2])
  if(length(colnames(tmpGt))==3){
    tmp <- data.frame(tmpVcf[,c(1,2,4,5,7)],tmpGt,tumorID,info,func,callers,id,stringsAsFactors=FALSE)
  } else {
    tmp <- data.frame(tmpVcf[,c(1,2,4,5,7)],tmpGt,normal=NA,tumorID,info,func,callers,id,stringsAsFactors=FALSE)
  }
  colnames(tmp) <- c("CHROM","POS","REF","ALT","FILTER","FORMAT","Tumor","Normal","Tumor_ID","snpEff","snpEff_func","caller","id")
  rownames(tmp) <- NULL
  return(tmp)
}

vcf <- lapply(filenames,function(x) vcfRead(x,VCF_DIR))
vcf <- do.call(rbind,vcf)
nrow(vcf) #306237


#CORRECT VCF FOR MULTIPLE ALT

table(is.na(vcf$ALT))
# vcf <- vcf[!is.na(vcf$ALT),]
table(grepl(",",vcf$ALT))
tmp <- sapply(vcf$ALT,function(x) unlist(strsplit(x,",",fixed=TRUE))[1])
vcf$ALT <- tmp





# SELECT CODING/SPLICING

tmp <- unique(vcf$snpEff_func)
tmp <- unique(unlist(sapply(tmp,function(x) unlist(strsplit(x,",")))))
tmp <- unique(unlist(sapply(tmp,function(x) unlist(strsplit(x,"&")))))
tmp <- tmp[!grepl("ENS",tmp)]
tmp <- tmp[!grepl("=",tmp)]
sort(tmp)
sel <- c("5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_truncation", #NEW
         "3_prime_UTR_variant", "5_prime_UTR_variant", # Added
         "bidirectional_gene_fusion","conservative_inframe_deletion","conservative_inframe_insertion",
         "disruptive_inframe_deletion","disruptive_inframe_insertion",
         "exon_loss_variant","initiator_codon_variant","gene_fusion",
         "frameshift_variant","missense_variant","splice_acceptor_variant",
         "splice_donor_variant","splice_region_variant","start_lost","stop_gained",
         "stop_lost","stop_retained_variant")
setdiff(tmp,sel)
setdiff(sel,tmp)
sel <- sapply(sel,function(x) grep(x,vcf$snpEff_func))
sel <- unique(unlist(sel))
sel <- sort(sel)
sum(is.na(sel))
nrow(vcf) #227168
rownames(vcf) <- NULL
length(sel) #36938
vcf <- vcf[sel,]
sort(table(vcf$Tumor_ID))

callerN <- sapply(vcf$caller,function(x) length(unlist(strsplit(x,",",fixed=TRUE))))
callerN <- unname(callerN)
vcf <- data.frame(vcf,callerN,stringsAsFactors = FALSE)
table(vcf$callerN)
#   1     2     3     4     5     6 
# 29690  3753   704   328   438  2025  
rownames(vcf) <- NULL
#----------------------------------------------

save(vcf,file = paste(fileDir,"vcf.rda",sep=""))




# COUNTS and ALLELE FREQUENCY -----------------

#countsVcf <- function(vcf,i){
countsVcf <- function(i){
  Tmp <- rep(NA,4)
  names(Tmp) <- c("t_ref_count","t_alt_count","n_ref_count","n_alt_count")
  vcfLine <- as.character(vcf[i,])
  names(vcfLine) <- colnames(vcf)
  caller <- unlist(strsplit(vcfLine["caller"],",",fixed=TRUE))[1]
  if(caller=="varscan"){
    tmp <- unlist(strsplit(vcfLine["Tumor"],":"))
    names(tmp) <- unlist(strsplit(vcfLine["FORMAT"],":"))
    Tmp["t_alt_count"] <- as.numeric(tmp["AD"])
    Tmp["t_ref_count"] <- as.numeric(tmp["RD"])
    tmp <- unlist(strsplit(vcfLine["Normal"],":"))
    names(tmp) <- unlist(strsplit(vcfLine["FORMAT"],":"))
    Tmp["n_alt_count"] <- as.numeric(tmp["AD"])
    Tmp["n_ref_count"] <- as.numeric(tmp["RD"])
  } else {
    tmp <- unlist(strsplit(vcfLine["Tumor"],":"))
    names(tmp) <- unlist(strsplit(vcfLine["FORMAT"],":"))
    if(tmp["AD"]!=".") Tmp[c("t_ref_count","t_alt_count")] <- as.numeric(unlist(strsplit(tmp["AD"],","))[1:2])
    if(!is.na(vcfLine["Normal"])) {
      tmp <- unlist(strsplit(vcfLine["Normal"],":"))
      names(tmp) <- unlist(strsplit(vcfLine["FORMAT"],":"))[1:length(tmp)]
      if(tmp["AD"]!=".") Tmp[c("n_ref_count","n_alt_count")] <- as.numeric(unlist(strsplit(tmp["AD"],","))[1:2])
    }
  }
  return(Tmp)
}

#cl <- makePSOCKcluster(NUM_FORKS,outfile="")  #60
#registerDoParallel(cl)
#tmp <- foreach(i = 1:nrow(vcf)) %dopar% {
#  tmp <- countsVcf(vcf,i)
#  return(tmp)
#}
#stopCluster(cl)


if(Sys.info()["sysname"]=="Windows"){

cl <- parallel::makeCluster(getOption("cl.cores", NUM_FORKS))

tmp <- parallel::parLapply(cl, 1:nrow(vcf), countsVcf)

parallel::stopCluster(cl)
}else{
tmp <- parallel::mclapply(1:nrow(vcf), countsVcf, mc.cores = NUM_FORKS)
}


tmp <- do.call(rbind,tmp)
dim(tmp) #36938     4

counts <- cbind(tmp[,1:2],NA,NA,tmp[,3:4],NA,NA)
colnames(counts) <- c("t_ref_count","t_alt_count","t_depth","t_vaf","n_ref_count","n_alt_count","n_depth","n_vaf")
counts[,"t_depth"] <- counts[,"t_alt_count"] + counts[,"t_ref_count"]
counts[,"t_vaf"] <- round(counts[,"t_alt_count"]/counts[,"t_depth"],2)
counts[,"n_depth"] <- counts[,"n_alt_count"] + counts[,"n_ref_count"]
counts[,"n_vaf"] <- round(counts[,"n_alt_count"]/counts[,"n_depth"],2)
rownames(counts) <- NULL

#View(counts[sample(1:nrow(counts),100),])

vcf <- data.frame(vcf, counts, stringsAsFactors = FALSE)

table(cond <- is.na(vcf$n_depth)&!is.na(vcf$Normal))
# FALSE  TRUE 
# 27894  9044 
# table(cond <- is.na(vcf$n_depth))

tmp <- vcf[cond,]

# table(tmp$caller)
# freebayes 
# 9044

vcf <- vcf[!cond,]

rownames(vcf) <- NULL
nrow(vcf) #27894
#rm(counts)

#----------------------------------------------
# save(vcf,file = "Analysis/vcf.rda")




############## FP FILTER #########
# FP FILTER
#system("mkdir Analysis/fpfilter")
#setwd("Analysis/fpfilter/")
dir.create(fpfilterDir)
setwd(fpfilterDir)

print(fpfilterDir)

vcf <- data.frame(vcf,fpFid="",stringsAsFactors=FALSE)
patients <- unique(vcf$Tumor_ID)
for(i in 1:length(patients)){
  print(patients[i])
  name <- patients[i]
  var <- vcf[vcf$Tumor_ID==name,1:4]
  var[nchar(var$ALT)>1,"ALT"] <- paste0("+",substring(var[nchar(var$ALT)>1,"ALT"],2))
  var[nchar(var$REF)>1,"POS"] <- as.numeric(var[nchar(var$REF)>1,"POS"])+1
  var[nchar(var$REF)>1,"ALT"] <- paste0("-",substring(var[nchar(var$REF)>1,"REF"],2))
  var[nchar(var$REF)>1,"REF"] <- substring(var[nchar(var$REF)>1,"REF"],2,2)
  vcf[vcf$Tumor_ID==name,"fpFid"] <- paste0(name,"-",var[,1],":",var[,2],"-",var[,4])
  write.table(var,file=paste0(name,".var"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  # n <- ceiling(nrow(var)/60)
  # system(paste0("split -l ",n," --additional-suffix .var ",paste0(name,".var")))
  loc <- paste(var$CHROM,var$POS,var$POS)
  write.table(loc,file=paste0(name,".loc"),quote=FALSE,row.names=FALSE,col.names=FALSE)
  # system(paste0("split -l ",n," --additional-suffix .loc ",paste0(name,".loc")))
}
