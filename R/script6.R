args = commandArgs( trailingOnly = TRUE )
CBIO_DB=args[1]
CBIO_SCORE=args[2]
NUM_FORKS=args[3]
CENSUS_FILE=args[4]
NATURE_FILE=args[5]
ONCOKB_FILE=args[6]
SMG_FILE=args[7]
library(parallel)
load(file='Final.RData')

selVar_cBio <- read.delim(CBIO_SCORE,header=FALSE,as.is=TRUE)
load(CBIO_DB,verbose=TRUE)
ID <- paste0(selVar_cBio[,1],":",selVar_cBio[,2],"-",selVar_cBio[,3])
id <- paste0(Final[,"Chr"],":",Final[,"Start"],"-",Final[,"End"])

cBio <- function(x){
  cBioID <- as.numeric(sub("cBio_","",unique(selVar_cBio[ID==x,7])))
  a <- paste0(cBioID,collapse="|")
  b <- paste0(cBioMutDb[cBioID,"AAchange"],collapse="|")
  c <- paste0(cBioMutDb[cBioID,"KnownVar"],collapse="|")
  d <- paste0(cBioMutDb[cBioID,"StudyID"],collapse="|")
  e <- sum(as.numeric(cBioMutDb[cBioID,"SamplesNum"]))
  f <- paste0(cBioMutDb[cBioID,"Freq"],collapse="|")
  e <- c(a,b,c,d,e,f)
  return(e)
}

#cl <- makePSOCKcluster(NUM_FORKS,outfile="") #20
#registerDoParallel(cl)
##clusterExport(cl, c("cBioMutDb","id"))
#cBioAnn <- foreach(i = 1:length(id)) %dopar% {
#  tmp <- cBio(id[i])
#  return(tmp)
#}
#stopCluster(cl)


cBioFunc <- function(i){
  return(cBio(id[i]))
}

if(Sys.info()["sysname"]=="Windows"){
  
  cl <- parallel::makeCluster(getOption("cl.cores", NUM_FORKS))
  
  cBioAnn <- parallel::parLapply(cl, 1:length(id), cBioFunc)
  
  parallel::stopCluster(cl)
}else{
  cBioAnn <- parallel::mclapply(1:length(id), cBioFunc, mc.cores = NUM_FORKS)
}



cBioAnn <- do.call(rbind,cBioAnn)
nrow(cBioAnn) #27883
colnames(cBioAnn) <- c("cBioId","cBioAAchange","cBioKnownVar","cBioStudyId","cBioSamplesNum","cBioFreq")
Final <- data.frame(Final,cBioAnn,stringsAsFactors = FALSE)


#ONCO ANNOTATION


census <- read.delim(CENSUS_FILE,as.is=TRUE)
nature <- read.delim(NATURE_FILE,as.is=TRUE)
oncoKB <- read.delim(ONCOKB_FILE,as.is=TRUE)
oncoKB <- oncoKB[oncoKB$Oncogenicity=="Oncogenic",]

onco <- function(i){
  x <- Final[i,]
  gene <- unlist(strsplit(as.character(x["Gene.refGene"]),";"))[1]
  #census
  if(gene%in%census$Gene.Symbol){
    a <- census[census$Gene.Symbol%in%gene,"Role.in.Cancer"]
    b <- paste0(census[census$Gene.Symbol%in%gene,c("Tumour.Types.Germline.","Cancer.Syndrome")],collapse="|")
  } else {
    a <- "."
    b <- "."
  }
  #nature
  if(gene%in%nature$Gene..Symbol){
    c <- paste0(nature[nature$Gene..Symbol%in%gene,"Mechanism.of.action.of.CPG.mutations"],collapse="|")
    d <- nature[nature$Gene..Symbol%in%gene,"Cancer.syndrome.s."]
  } else {
    c <- "."
    d <- "."
  }
  #oncoKB
  if(gene%in%oncoKB$Gene) {
    e <- paste0(oncoKB[oncoKB$Gene%in%gene,"Mutation.Effect"],collapse="|")
    f <- paste0(oncoKB[oncoKB$Gene%in%gene,"Alteration"],collapse="|")
  } else {
    e <- "."
    f <- "."
  }
  g <- c(a,b,c,d,e,f)
  return(g)
}


oncoFunc <- function(i){
  return(onco(i))
}

if(Sys.info()["sysname"]=="Windows"){
  
  cl <- parallel::makeCluster(getOption("cl.cores", NUM_FORKS))
  
  Onco <- parallel::parLapply(cl, 1:nrow(Final), oncoFunc)
  
  parallel::stopCluster(cl)
}else{
  Onco <- parallel::mclapply(1:nrow(Final), oncoFunc, mc.cores = NUM_FORKS)
}

#cl <- makePSOCKcluster(NUM_FORKS,outfile="") #30
#registerDoParallel(cl)
#Onco <- foreach(i = 1:nrow(Final)) %do% {
#  tmp <- onco(i)
#  return(tmp)
#}
#stopCluster(cl)

Onco <- do.call(rbind,Onco)
nrow(Onco) #27883
colnames(Onco) <- c("census","censusGermline","nature","natureGermline","oncoKB","oncoKBmut")

Final <- data.frame(Final,Onco ,stringsAsFactors = FALSE)
#save(Final,file="Final.tmp.rda")

# Panglioma recurrent mutations


smg <- read.csv(SMG_FILE,as.is=TRUE)
genes <- sapply(Final[,"Gene.refGene"],function(x) unlist(strsplit(x,","))[1])
smgOcc <- NULL
for(i in 1:nrow(Final)){
  if(genes[i]%in%smg[,1]){
    smgOcc <- c(smgOcc,smg[smg[,1]==genes[i],2])
  }
  else smgOcc <- c(smgOcc,0)
}
Final <- data.frame(Final,smgOcc,stringsAsFactors = FALSE)


# Damaging prediction counts
# pred <- Final[,colnames(Final)[grep("_pred",colnames(Final))]]
# pred <- pred[,c("SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")]
unique(Final$SIFT_pred) # "D" "." "T"
unique(Final$Polyphen2_HDIV_pred) #"P" "." "B" "D"
Final$Polyphen2_HDIV_pred[Final$Polyphen2_HDIV_pred=="P"] <- "D"
Final$Polyphen2_HDIV_pred[Final$Polyphen2_HDIV_pred=="B"] <- "T"
unique(Final$MutationTaster_pred) #"D" "." "N" "A" "P"
Final$MutationTaster_pred[Final$MutationTaster_pred=="A"] <- "D"
Final$MutationTaster_pred[Final$MutationTaster_pred=="N"] <- "T"
Final$MutationTaster_pred[Final$MutationTaster_pred=="P"] <- "T"
unique(Final$PROVEAN_pred) #"N" "." "D"
Final$PROVEAN_pred[Final$PROVEAN_pred=="N"] <- "T"
tmp <- apply(Final[,c("SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")],1,function(x){
  d <- sum(x%in%c("D"))
  t <- sum(x%in%c("T"))
  u <- sum(x%in%c("."))
  c(d,t,u)
})
tmp <- t(tmp)
colnames(tmp) <- c("Damaging","Tolerated","Unpredicted")

Final <- data.frame(Pathogenicity=".",MissenseDamaging=tmp[,1],Final,stringsAsFactors = FALSE)

table(Final$Func.refGene)

Final[Final$Func.refGene=="nonsynonymous SNV"&as.numeric(Final$MissenseDamaging)>=2,"Pathogenicity"] <- "pathogenic missense"
Final[Final$Func.refGene=="nonsynonymous SNV"&as.numeric(Final$MissenseDamaging)<2,"Pathogenicity"] <- "tolerated missense"

Final[Final$Pathogenicity=="."&as.numeric(Final$MissenseDamaging)>=2,"Pathogenicity"] <- "pathogenic"

Final[Final$Func.refGene%in%c("frameshift deletion","frameshift insertion","frameshift substitution"),"Pathogenicity"] <- "frameshift"
Final[Final$Func.refGene%in%c("frameshift deletion","frameshift insertion","frameshift substitution"),
      c("MissenseDamaging","SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")] <- "."

Final[Final$Func.refGene%in%c("splicing","exonic;splicing"),"Pathogenicity"] <- "splicing"
Final[Final$Func.refGene%in%c("splicing","exonic;splicing"),
      c("MissenseDamaging","SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")] <- "."

Final[Final$Func.refGene%in%c("stopgain"),"Pathogenicity"] <- "truncating"
Final[Final$Func.refGene%in%c("stopgain"),
      c("MissenseDamaging","SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")] <- "."

Final[Final$Func.refGene%in%c("stoploss"),"Pathogenicity"] <- "stoploss"
Final[Final$Func.refGene%in%c("stoploss"),
      c("MissenseDamaging","SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")] <- "."

Final[Final$Func.refGene%in%c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution"),"Pathogenicity"] <- "inframe"
Final[Final$Func.refGene%in%c("nonframeshift deletion","nonframeshift insertion","nonframeshift substitution"),
      c("MissenseDamaging","SIFT_pred","Polyphen2_HDIV_pred","MutationTaster_pred","PROVEAN_pred")] <- "."

table(Final$Pathogenicity)


### DOWNSTREAM SELECTION

# 1000Genomes and dbNSFP filter AF ≤ 0.0004
Final <- data.frame(PopFreq=".",Final,stringsAsFactors = FALSE)
toReject <- as.numeric(Final[,"ExAC_ALL"])>0.0004
toReject[is.na(toReject)] <- FALSE
Final[toReject,"PopFreq"] <- "reject"
toReject <- as.numeric(Final[,"X1000g2015aug_all"])>0.0004
toReject[is.na(toReject)] <- FALSE
Final[toReject,"PopFreq"] <- "reject"
table(Final$PopFreq)
# . reject 
# 15887  11996  

# Low confidence n_alt_count > 1 | t_depth < 15 | t_alt_count <= 3 | T_vaf < 0.05
Final <- data.frame(LowConf=".",Final,stringsAsFactors = FALSE)

sum(sel <- is.na(Final$n_vaf))
unique(Final[sel,"Patient_ID"])
# sum(sel <- grepl("NA",Final$n_vaf))
# unique(Final[sel,"caller"])

lowConf <- function(i){
  if(!is.na(Final$n_vaf[i])){
    tmp <- Final$n_alt_count[i]
    # tmp <- as.numeric(unlist(strsplit(Final$n_alt_count[i],"|",fixed=TRUE)))
    # if(length(tmp)>1){
    #   toRej1 <- sum(tmp>1)>1
    # } else {
    toRej1 <- as.numeric(tmp)>1
    # } 
  } else toRej1 <- FALSE
  tmp <- Final$t_depth[i]
  # tmp <- as.numeric(unlist(strsplit(Final$t_depth[i],"|",fixed=TRUE)))
  # if(length(tmp)>1) {
  #   toRej2 <- sum(tmp<15)>1
  # } else {
  toRej2 <- tmp<15
  # }
  tmp <- Final$t_alt_count[i]
  # tmp <- as.numeric(unlist(strsplit(Final$t_alt_count[i],"|",fixed=TRUE)))
  # if(length(tmp)>1) {
  #   toRej3 <- sum(tmp<=3)>1
  # } else {
  toRej3 <- tmp<=3
  # }
  tmp <- Final$t_vaf[i]
  # tmp <- as.numeric(unlist(strsplit(Final$t_vaf[i],"|",fixed=TRUE)))
  # if(length(tmp)>1) {
  #   toRej4 <- sum(tmp<0.05)>1
  # } else {
  toRej4 <- tmp<=0.05
  # }
  toRej <- toRej1|toRej2|toRej3|toRej4
  return(toRej)
}

# res <- NULL
# for(i in 1:nrow(Final)){
#   print(i)
#   tmp <- lowConf(i)
#   res <- c(res,tmp)
# }


# cl <- makePSOCKcluster(NUM_FORKS,outfile="") #30
#registerDoParallel(cl)
# res <- foreach(i = 1:nrow(Final)) %dopar% {
#  tmp <- lowConf(i)
#  return(tmp)
#}
#stopCluster(cl)


lowConfFunc <- function(i){
  return(lowConf(i))
}

if(Sys.info()["sysname"]=="Windows"){
  
  cl <- parallel::makeCluster(getOption("cl.cores", NUM_FORKS))
  
  res <- parallel::parLapply(cl, 1:nrow(Final), lowConfFunc)
  
  parallel::stopCluster(cl)
}else{
  res <- parallel::mclapply(1:nrow(Final), lowConfFunc, mc.cores = NUM_FORKS)
}

res <- unlist(res)

table(res)
# FALSE   TRUE 
# 12495 15388 
# table(is.na(res))
Final[res,"LowConf"] <- "reject"
table(Final$LowConf)
# . reject 
#  12495  15388 


# TABLE FINALIZATION

# tmp <- paste0(colnames(Final),collapse = "','")
# 
# tmp <- c('Description','LowConf','PopFreq','Pathogenicity','MissenseDamaging','FLAGSscore','FLAGS',
#          'MaskedFilter','Score','fiftymer','SNP','blackGenes','effectFilt','effectSnpEff','effectAnnovar','Patient_ID',
#          'Type','Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','Interpro_domain',
#          'avsnp150','snp138NonFlagged','ExAC_ALL','X1000g2015aug_all','esp6500siv2_all','Kaviar_AF','HRC_AF','cosmic80','CLINSIG',
#          'CLNDBN','CLNACC','CLNDSDB','CLNDSDBID','SIFT_pred','Polyphen2_HDIV_pred','MutationTaster_pred','PROVEAN_pred',
#          'dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','dann','Eigen','gerp..gt2','CADD','CADD_Phred',
#          'snpEff','snpEff_func','caller','t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','t_vaf','n_vaf',
#          'FPfilter','FPfail','cBio','cBioId','cBioAAchange','cBioKnownVar','cBioStudyId','cBioSamplesNum','cBioFreq',
#          'census','censusGermline','nature','natureGermline','oncoKB','oncoKBmut')

tmp <- colnames(Final)

sel <- c("Patient_ID",#'Tumor_ID',"Normal_ID", #SAMPLE
         'caller','callerN',"duplicated",'Type','Chr','Start','End','Ref','Alt',
         #"HGVSc","HGVSp", #VARIANT
         'Gene.refGene','Description','Func.refGene','GeneDetail.refGene','snpEff','snpEff_func','Interpro_domain','effectFilt','effectSnpEff','effectAnnovar', #GENE ANNOTATION
         'SNP','PopFreq','avsnp150','snp138NonFlagged','ExAC_ALL','X1000g2015aug_all','esp6500siv2_all','Kaviar_AF','HRC_AF', #POPULATION FREQUENICIES
         'FPfilter','FPfail','MaskedFilter','Score','fiftymer','FLAGSscore','FLAGS','blackGenes', #FALSE POSITIVE FILTERS
         #"VN_occurrences","VN_frequency","VN_fullycalled_count","VN_fullycalled_frequency","vn",
         "PON_filter","PON_filter_pass", #FALSE POSITIVE FILTERS
         'LowConf','t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','t_vaf','n_vaf', #READ COUNTS
         'Pathogenicity','MissenseDamaging','SIFT_pred','Polyphen2_HDIV_pred','MutationTaster_pred','PROVEAN_pred', #PATHOGENICITY
         'dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','dann','Eigen','gerp..gt2','CADD','CADD_Phred', #WGS PATHOGENICITY
         'cosmic80','CLINSIG','CLNDBN','CLNACC','CLNDSDB','CLNDSDBID', #CLINICAL ANNOTATION
         'census','censusGermline','nature','natureGermline','oncoKB','oncoKBmut',"smgOcc", #ONCO ANNOTATION
         'cBioId','cBioAAchange','cBioKnownVar','cBioStudyId','cBioSamplesNum','cBioFreq' #cBio ANNOTATION
)
setdiff(tmp,sel)
setdiff(sel,tmp)
Final <- Final[,sel]

Final <- Final[order(as.numeric(Final$Start)),]
Final <- Final[order(Final$Gene.refGene),]
#Final <- Final[order(Final$Tumor_ID),]
rownames(Final) <- NULL
# tmp <- sapply(Final$Patient_ID, function(x) unlist(strsplit(x,"_"))[1])
# Final <- Final[order(tmp),]

Final[Final==""] <- "."

# HGVS from SNPEFF

snpEff <- sapply(Final$snpEff,function(x){
  temp <- unlist(strsplit(x,"ANN"))[2]
  # temp <- unlist(strsplit(temp,","))[1]
  temp <- unlist(strsplit(temp,","))
  if(sum(grepl("p.",temp,fixed=TRUE))>0){
    temp <- temp[grep("p.",temp,fixed=TRUE)][1]
  } else temp <- temp[1]
  info <- unlist(strsplit(temp,"|",fixed=TRUE))
  HGVSc <- info[grep("c.",info,fixed=TRUE)][1]
  HGVSp <- info[grep("p.",info,fixed=TRUE)][1]
  return(c(HGVSc,HGVSp))
})
dim(snpEff)
snpEff <- t(snpEff)
colnames(snpEff) <- c("HGVSc","HGVSp")
rownames(snpEff) <- NULL
Final <- data.frame(snpEff,Final,stringsAsFactors = FALSE)
# Final[,c("HGVSc","HGVSp")] <- snpEff[,c("HGVSc","HGVSp")]



Final <- data.frame(selected=".",Final,stringsAsFactors = FALSE)
Final$selected[grepl("path",Final$CLINSIG,ignore.case = TRUE)] <- "yes"
Final$selected[Final$Gene.refGene%in%c("NF1","SUZ12","EED")] <- "yes"
table(Final$selected) # 36 
# View(Final[Final$selected=="yes",])
reject <- Final$effectFilt=="reject" |  Final$SNP=="reject" | Final$fiftymer=="reject" | Final$FLAGS=="reject" | Final$blackGenes=="reject"
# View(Final[reject&Final$selected=="yes",])
table(reject&Final$selected=="yes") #3
Final$selected[reject&Final$selected=="."] <- "no"
table(Final$selected)
# .    no   yes 
# 8591 19261    36 

FinalS <- Final[Final$selected!="no",]
rownames(FinalS) <- NULL

table(FP_1 <- FinalS$FPfilter=="FAIL" & FinalS$MaskedFilter=="fail") #472
table(FP_tmp <- FinalS$PopFreq=="reject" | FinalS$FPfilter=="FAIL" | FinalS$MaskedFilter=="fail" | FinalS$LowConf=="reject") #5707
table(FP_2 <- FP_tmp & (FinalS$PON_filter_pass==FALSE | FinalS$callerN<3)) #5329
table(FP_3 <- FinalS$effectAnnovar=="noncoding") #2414
table(FP_4 <- FinalS$FPfilter=="no_readcounts"&FinalS$callerN<3&(nchar(FinalS$Ref)>3|nchar(FinalS$Alt)>3)) #155
table(FP_5 <- FinalS$callerN<3 & FinalS$PON_filter_pass==FALSE) #3565 
#table(FP_6 <-  (FinalS$PON_filter > -10 | (FinalS$avsnp150!="." & FinalS$cosmic80=="."))) #5442  

# tmp <- FinalS[!(FP_1 | FP_2 | FP_3 | FP_4 | FP_5) & FP_6,]
# View(tmp)

table(FinalS$toCheck <- (FP_1 | FP_2 | FP_3 | FP_4 | FP_5 ) & grepl("path",FinalS$CLINSIG,ignore.case = TRUE)) #11 
#table(FinalS$toCheck <- (FP_1 | FP_2 | FP_3 | FP_4 | FP_5 | FP_6) & grepl("path",FinalS$CLINSIG,ignore.case = TRUE)) #11 

#View(FinalS[FinalS$toCheck,])
sort(unique(FinalS$CLNDBN[FinalS$toCheck]))
sort(unique(FinalS$Gene.refGene[FinalS$toCheck]))

table(toReject1 <- FinalS$toCheck & !FinalS$n_alt_count%in%c("0",".")) #1
table(toReject2 <- FinalS$toCheck  & FinalS$PopFreq=="reject") #6
table(toKeep <- FinalS$toCheck & !(toReject1 | toReject2)) #5
table(toKeep <- toKeep | (FinalS$toCheck & FinalS$Gene.refGene%in%c("TP53","IDH1"))) #5

#View(FinalS[toKeep,])

FP <- (FP_1 | FP_2 | FP_3 | FP_4 | FP_5) & !toKeep
#FP <- (FP_1 | FP_2 | FP_3 | FP_4 | FP_5 | FP_6) & !toKeep

table(FP)
# FALSE  TRUE 
#  1109  6459 

FinalF <- FinalS[!FP,]
FinalS <- data.frame(FP=".",FinalS,stringsAsFactors = FALSE)
FinalS$FP[FP] <- "reject"
table(FinalS$FP)
#    . reject 
#  1109   6459 
# FinalF <- FinalS[FinalS$FP==".",]
nrow(FinalF) #1104
rownames(FinalF) <- NULL


# HOTSPOST LIKE

if(nrow(FinalF)>0){
  FinalF <- data.frame(Hotspot=".",FinalF,stringsAsFactors = FALSE)
  mut <- paste0(FinalF$Chr,",",FinalF$Start)
  Mut <- table(mut)
  Mut <- sort(Mut,decreasing=TRUE)
  Mut <- cbind(Mut,1)
  for(i in 1:nrow(Mut)){
    print(i)
    tmp <- unique(FinalF[mut==rownames(Mut)[i],"Patient_ID"])
    # tmp <- unique(sapply(tmp,function(x) unlist(strsplit(x,"_"))[1]))
    Mut[i,2] <- length(tmp)
  }
  
  for(i in 1:nrow(FinalF)){
    FinalF[i,"Hotspot"] <- Mut[mut[i],2]
  }
  table(FinalF$Hotspot)
  
  # 1  10  12  13  14  15  16  17  18  19   2  20   3   4   5   6   7   8 
  # 42  10  12  13  28  45  48  51  18  38  22 700   9  12   5  12   7  32 
  
  #View(FinalF[FinalF$Hotspot!=1,])
  table(FinalF$Gene.refGene[FinalF$Hotspot!=1])



### REDUNDANT MUTATIONS

  rownames(FinalF) <- NULL
  id <- paste0(FinalF$Tumor_ID,"|",FinalF$Gene.refGene)
  table(red <- duplicated(id))
  # FALSE  TRUE 
  #  112   992 
  red <- id%in%id[red]
  FinalF <- data.frame(red,FinalF,stringsAsFactors = FALSE)
  #View(FinalF[FinalF$red==TRUE,c("Patient_ID","caller","callerN","Chr","Start","End","Ref","Alt","Gene.refGene")])
  toRemove <- FinalF$red & FinalF$callerN==1
  #View(FinalF[toRemove,c("Patient_ID","caller","callerN","Chr","Start","End","Ref","Alt","Gene.refGene")])
  FinalF <- FinalF[!toRemove,]
  rownames(FinalF) <- NULL
  nrow(FinalF) #1078
  
  # Mutation recurrence (correct to count only one mutation for patient)
  genes <- sapply(FinalF[,"Gene.refGene"],function(x) unlist(strsplit(x,","))[1])
  MutGenes <- table(genes)
  MutGenes <- sort(MutGenes,decreasing=TRUE)
  MutGenes <- cbind(MutGenes,1)
  for(i in 1:nrow(MutGenes)){
    print(rownames(MutGenes)[i])
    tmp <- FinalF[genes==rownames(MutGenes)[i],"Patient_ID"]
    tmp <- unique(tmp)
    MutGenes[i,2] <- length(tmp)
  }
  FinalF <- data.frame(FinalF,GeneRec=1,GenePzRec=1,stringsAsFactors = FALSE)
  for(i in 1:nrow(FinalF)){
    FinalF[i,"GeneRec"] <- MutGenes[genes[i],1]
    FinalF[i,"GenePzRec"] <- MutGenes[genes[i],2]
  }
  table(FinalF$GeneRec)
  table(FinalF$GenePzRec)
  #View(FinalF[FinalF$GenePzRec!=1,])
  
  mutCounts <- sort(table(FinalF$Tumor_ID))
  
  tmp <- colnames(FinalF)
  sel <- c("Patient_ID", #"Tumor_ID", "Normal_ID",
           "caller",  "callerN",  "Type", "Chr",  "Start",  "End", "Ref", "Alt",  
           "Gene.refGene",  "Description",  "Func.refGene", "GeneDetail.refGene","snpEff_func", "HGVSc",  "HGVSp",  "Interpro_domain",  
           "Pathogenicity", "MissenseDamaging", "SIFT_pred",  "Polyphen2_HDIV_pred",  "MutationTaster_pred",  "PROVEAN_pred", 
           "census",  "censusGermline", "nature", "natureGermline", "oncoKB", "oncoKBmut",  "smgOcc", 
           "FPfilter",  "FPfail", "MaskedFilter", "PON_filter", "PON_filter_pass", "Hotspot",
           # "VN_occurrences","vn",
           "avsnp150",  "ExAC_ALL", "X1000g2015aug_all", "esp6500siv2_all", "Kaviar_AF",  "HRC_AF", 
           "cosmic80",  "CLINSIG",  "CLNDBN", "CLNACC", "CLNDSDB",  "CLNDSDBID",
           "t_depth", "t_ref_count",  "t_alt_count", "n_depth", "n_ref_count",  "n_alt_count",  "t_vaf",  "n_vaf",
           "cBioId",  "cBioAAchange", "cBioKnownVar", "cBioStudyId",  "cBioSamplesNum", "cBioFreq",
           "GeneRec", "GenePzRec")
  setdiff(sel, tmp)
  setdiff(tmp, sel)  
  FinalF <- FinalF[,sel]
}

save(Final,FinalS,FinalF,file="Final_no_PoN_filter.rda")
FinalS <- FinalS[,!colnames(FinalS)%in%c("snpEff")]
write.table(FinalS,file = "SomaticVariants.txt",sep="\t",row.names = FALSE)
write.table(FinalF,file = "SomaticVariants_FPfilt.txt",sep="\t",row.names = FALSE)
write.table(FinalF$Gene.refGene,file = "SomaticVariants_FPfilt_genes.txt",sep="\t",row.names = FALSE)