args = commandArgs( trailingOnly = TRUE )
BASE_DIR=args[1]
PONFILTER_FILE=args[2]
PON_DIR=args[3]


httr::set_config(httr::config(ssl_verifypeer = FALSE))

fileDir = paste0(BASE_DIR,"/Analysis/")
load(paste0(fileDir,"vcf.rda"))

ann.out <- read.delim(paste0(fileDir,"ann.out.hg19_multianno.txt"),as.is=TRUE)
paste0(colnames(ann.out),collapse="','")
tmp <- c('Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene', #refGene
         'avsnp150', #avsnp150
         'snp138NonFlagged', #snp138NonFlagged
         'ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS', #exac03
         'X1000g2015aug_all', #1000g2015aug_all
         'esp6500siv2_all', #esp6500siv2_all
         'Kaviar_AF','Kaviar_AC','Kaviar_AN', #kaviar_20150923
         'HRC_AF','HRC_AC','HRC_AN','HRC_non1000G_AF','HRC_non1000G_AC','HRC_non1000G_AN', #hrcr1
         'cosmic80', #cosmic80
         'CLINSIG','CLNDBN','CLNACC','CLNDSDB','CLNDSDBID', #clinvar_20170905
         'SIFT_score','SIFT_converted_rankscore','SIFT_pred', #dbnsfp33a
         'Polyphen2_HDIV_score','Polyphen2_HDIV_rankscore','Polyphen2_HDIV_pred', #dbnsfp33a
         'Polyphen2_HVAR_score','Polyphen2_HVAR_rankscore','Polyphen2_HVAR_pred', #dbnsfp33a
         'LRT_score','LRT_converted_rankscore','LRT_pred', #dbnsfp33a
         'MutationTaster_score','MutationTaster_converted_rankscore','MutationTaster_pred', #dbnsfp33a
         'MutationAssessor_score','MutationAssessor_score_rankscore','MutationAssessor_pred', #dbnsfp33a
         'FATHMM_score','FATHMM_converted_rankscore','FATHMM_pred', #dbnsfp33a
         'PROVEAN_score','PROVEAN_converted_rankscore','PROVEAN_pred', #dbnsfp33a
         'VEST3_score','VEST3_rankscore','MetaSVM_score','MetaSVM_rankscore','MetaSVM_pred', #dbnsfp33a
         'MetaLR_score','MetaLR_rankscore','MetaLR_pred', #dbnsfp33a
         'M.CAP_score','M.CAP_rankscore','M.CAP_pred', #dbnsfp33a
         'CADD_raw','CADD_raw_rankscore','CADD_phred', #dbnsfp33a
         'DANN_score','DANN_rankscore', #dbnsfp33a
         'fathmm.MKL_coding_score','fathmm.MKL_coding_rankscore','fathmm.MKL_coding_pred', #dbnsfp33a
         'Eigen_coding_or_noncoding','Eigen.raw','Eigen.PC.raw', #dbnsfp33a
         'GenoCanyon_score','GenoCanyon_score_rankscore', #dbnsfp33a
         'integrated_fitCons_score','integrated_fitCons_score_rankscore','integrated_confidence_value', #dbnsfp33a
         'GERP.._RS','GERP.._RS_rankscore', #dbnsfp33a
         'phyloP100way_vertebrate','phyloP100way_vertebrate_rankscore','phyloP20way_mammalian','phyloP20way_mammalian_rankscore', #dbnsfp33a
         'phastCons100way_vertebrate','phastCons100way_vertebrate_rankscore','phastCons20way_mammalian','phastCons20way_mammalian_rankscore', #dbnsfp33a
         'SiPhy_29way_logOdds','SiPhy_29way_logOdds_rankscore', #dbnsfp33a
         'Interpro_domain', #dbnsfp33a
         'GTEx_V6_gene','GTEx_V6_tissue', #dbnsfp33a
         'dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE', #dbscsnv11 SPLICE SITE
         'dann', #dann
         'Eigen', #eigen
         'gerp..gt2', #gerp++gt2
         'CADD','CADD_Phred') #cadd
sel <- c('Chr','Start','End','Ref','Alt',
         'Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene', 'Interpro_domain',
         'avsnp150','snp138NonFlagged',
         'ExAC_ALL','X1000g2015aug_all','esp6500siv2_all','Kaviar_AF','HRC_AF',
         'cosmic80','CLINSIG','CLNDBN','CLNACC','CLNDSDB','CLNDSDBID',
         'SIFT_pred','Polyphen2_HDIV_pred','MutationTaster_pred','PROVEAN_pred',
         'dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE',
         'dann','Eigen','gerp..gt2','CADD','CADD_Phred')

ann.out <- ann.out[,sel]
ann.out[ann.out$ExonicFunc.refGene!=".","Func.refGene"] <- ann.out[ann.out$ExonicFunc.refGene!=".","ExonicFunc.refGene"]
ann.out[ann.out$AAChange.refGene!=".","GeneDetail.refGene"] <- ann.out[ann.out$AAChange.refGene!=".","AAChange.refGene"]
ann.out <- ann.out[,!colnames(ann.out)%in%c("ExonicFunc.refGene","AAChange.refGene")]


colnames(vcf)
sel <- c("snpEff","snpEff_func","caller","callerN","t_depth","t_ref_count", "t_alt_count",
         "n_depth","n_ref_count", "n_alt_count", "t_vaf","n_vaf","FPfilter","FPfail"#,
         #"VN_occurrences","VN_frequency","VN_fullycalled_count","VN_fullycalled_frequency","vn"
)
setdiff(colnames(vcf),sel)
setdiff(sel,colnames(vcf))

# tumor <- vcf$Tumor_ID
# info <- NULL
# for(i in 1:length(tumor)){
#   print(i)
#   tmp <- samples[samples$TUMOR_WES_ID==tumor[i],c("Patient_ID", "TUMOR_WES_ID", "NORMAL_WES_ID")]
#   info <- rbind(info,tmp)
# }
# colnames(info) <- c("Patient_ID","Tumor_ID","Normal_ID")

final <- data.frame(Patient_ID = vcf[, "Tumor_ID"], Type= "SNV", ann.out, vcf[, sel], stringsAsFactors = FALSE)
final[final$Ref=="-","Type"] <- "INS"
final[final$Alt=="-","Type"] <- "DEL"

id <- paste0(final$Patient_ID,".",final$Chr,".",final$Start,".",final$End,".",final$Alt)
sum(duplicated(id)) #11
final <- data.frame(final,duplicated=duplicated(id),stringsAsFactors = FALSE)
final <- final[!duplicated(id),]
Final <- final
rownames(Final) <- NULL
# save(Final,file="Final.tmp.rda")

if(sum(grepl("chr",Final$Chr))==0) Final$Chr <- paste0("chr", Final$Chr) 

# PoNs FILTER (reference_PoNs_final_summed_tokens.hist)----------
#source(PONFILTER_FILE)

ponFilter <- function(maf, ncores=60, threshold=-2.5, pondir){
  #maf=Final
  #ncores=1 
  #threshold=-2.5
  #pondir=PON_DIR
  
  maf$PON_filter <- rep(100,nrow(maf))
  maf$PON_filter_pass <- rep(TRUE,nrow(maf))
  chrs <- paste0("chr", c(1:22,"X"))
  #chrs <- c(1:22,"X","Y")
  x <- c(0.001, 0.003, 0.01, 0.03, 0.2, 1)
  for(chr in chrs){
    
    f_name<- paste0(pondir,"/pon_",chr,".rds")
    print(paste0("reading pon file ",f_name))
    pon <- readRDS(f_name)
   
    loc <- maf$Chr==chr
    
    sel_mut <- maf[loc,]
    require(doMC)
    registerDoMC(ncores)
    ans <- foreach(i=1:nrow(sel_mut)) %dopar% {
    #for(i in 1:nrow(sel_mut)){
      n_alt <- sel_mut$t_alt_count[i]
      n_ref <- sel_mut$t_ref_count[i]
      d <- pbeta(x, shape1=n_alt+1, shape2=n_ref+1)
      d <- c(0,d)
      f <- d[-1] - d[-length(d)]
      #f[length(f)+1]<-f[length(f)]
      position <- sel_mut$Start[i]
      p1 <- 8*(position-1)+1
      p2 <- p1+7
      h <- pon[(p1+2):p2]
      s <- f%*% h
      s <- log10(s)
      #print(s)
      #maf[loc,]$PON_filter[i]<- s
      if(is.na(s)){
        print(paste0(i, ": NA occurred f: ",f, " h:", h))
      }
      return(s)
    }
    ss <- unlist(lapply(ans, function(x) x[1,1]))
    maf$PON_filter[loc] <- ss
  }
  maf$PON_filter_pass <- maf$PON_filter < threshold  
  return(maf)
}

Final <- ponFilter(maf=Final, ncores=20, threshold=-2.5, 
                   pondir=PON_DIR)
table(Final$PON_filter_pass)
# FALSE  TRUE 
# 19207  8676


library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host="uswest.ensembl.org")
genes <- unique(Final$Gene.refGene)
# tmp <- listAttributes(ensembl)
# tmp[grep("description",tmp$name),]
# tmp <- listFilters(ensembl)
# tmp[grep("external_gene_name",tmp$name),]
info <- getBM(attributes = c("external_gene_name","description"),filters = c("external_gene_name"),values = genes, mart = ensembl)
Final <- data.frame(Description=".",Final,stringsAsFactors = FALSE)
info <- info[!duplicated(info$external_gene_name),]
rownames(info) <- info$external_gene_name
Final$Description <- info[Final$Gene.refGene,"description"]
table(is.na(Final$Description))


sel <- grep(";",Final$Gene.refGene)
genes <- Final$Gene.refGene[sel]
genes <- sapply(genes,function(x) unlist(strsplit(x,";")))
genes <- unique(unlist(genes))
info <- getBM(attributes = c("external_gene_name","description"),filters = c("external_gene_name"),values = genes, mart = ensembl)
description <- sapply(Final$Gene.refGene[sel],function(x){
  tmp <- unlist(strsplit(x,";"))
  if(tmp[1]%in%info$external_gene_name){
    a <- info[info$external_gene_name==tmp[1],"description"][1]
  } else {
    a <- "."
  }
  if(tmp[2]%in%info$external_gene_name){
    b <- info[info$external_gene_name==tmp[2],"description"][1]
  } else {
    b <- "."
  }
  paste0(a," ; ",b)
})

Final$Description[sel] <- description
table(is.na(Final$Description))


Final <- data.frame(effectSnpEff=".",effectAnnovar=".",Final,stringsAsFactors = FALSE)

#snpEff
tmp <- unique(Final$snpEff_func)
tmp <- unique(unlist(sapply(tmp,function(x) unlist(strsplit(x,",")))))
tmp <- unique(unlist(sapply(tmp,function(x) unlist(strsplit(x,"&")))))
sort(tmp)
coding <- c("conservative_inframe_deletion","conservative_inframe_insertion",
            "disruptive_inframe_deletion","disruptive_inframe_insertion",
            "frameshift_variant","initiator_codon_variant","missense_variant",
            "start_lost","stop_gained","stop_lost")
sort(setdiff(tmp,coding))
coding <- sapply(coding,function(x) grep(paste0("\\b",x,"\\b"),Final[,"snpEff_func"]))
coding <- sort(unique(unlist(coding)))
Final[coding,"effectSnpEff"] <- "coding"
splicing <- c("splice_acceptor_variant","splice_donor_variant","splice_region_variant")
splicing <- sapply(splicing,function(x) grep(paste0("\\b",x,"\\b"),Final[,"snpEff_func"]))
splicing <- sort(unique(unlist(splicing)))
Final[rownames(Final)%in%splicing & Final$effectSnpEff==".","effectSnpEff"] <- "splicing"
Final$effectSnpEff[Final$effectSnpEff=="."] <- "noncoding"
table(Final$effectSnpEff)
# coding noncoding  splicing 
# 13201      8874      5808 
#Annovar
tmp <- unique(Final[,"Func.refGene"])
tmp <- sort(tmp)
paste0(tmp,collapse="','")
tmp <- c('downstream','frameshift deletion','frameshift insertion','intergenic','intronic',
         'ncRNA_exonic','ncRNA_exonic;splicing','ncRNA_intronic','ncRNA_splicing',
         'nonframeshift deletion','nonframeshift insertion','nonsynonymous SNV','splicing',
         'stopgain','stoploss','synonymous SNV','unknown','upstream','upstream;downstream',
         'UTR3','UTR5','UTR5;UTR3')
Final[Final$Func.refGene%in%c('frameshift deletion','frameshift insertion','frameshift substitution',
                              'nonframeshift deletion','nonframeshift insertion','nonframeshift substitution',
                              'nonsynonymous SNV','stopgain','stoploss'),"effectAnnovar"] <- "coding" 
Final[Final$Func.refGene%in%c('splicing','ncRNA_exonic;splicing'),"effectAnnovar"] <- 'splicing' 
Final[Final$effectAnnovar==".","effectAnnovar"] <- "noncoding"
table(Final$effectAnnovar)
# 
# coding noncoding  splicing 
# 11149     16450       284 

Final <- data.frame(effectFilt="reject",Final,stringsAsFactors = FALSE)
Final[Final$effectSnpEff%in%c("coding","splicing")|Final$effectAnnovar%in%c("coding","splicing"),"effectFilt"] <- "keep"
table(Final$effectFilt)
# keep reject 
# 19054   8829


##Annotation of gene paralogs with high homology (http://massgenomics.org/2013/06/ngs-false-positives.html)

Final <- data.frame(blackGenes=".",Final,stringsAsFactors = FALSE)
artifacts <- c("LOC","ENS","FAM","GOL","PRA","NBP","POT","DEF","MUC","KRT","WAS","ANK","TRI","FRG","HLA",paste0("OR",1:9))
# artifacts <- c("LOC","ENS","FAM","GOL","PRA","NBP","DEF","MUC","KRT","WAS","ANK","FRG","HLA",paste0("OR",1:9)) #KEEP POT AND TRIM
genesSuff <- sapply(Final[,"Gene.refGene"],function(x) substring(unlist(strsplit(x,","))[1],1,3))
toRemove <- genesSuff%in%artifacts
Final[toRemove,"blackGenes"] <- "reject"
genesSuff <- sapply(Final[,"Gene.refGene"],function(x) substring(unlist(strsplit(x,","))[1],1,4))
toRemove <- genesSuff%in%c("PLIN","CELA","SRA1")
Final[toRemove,"blackGenes"] <- "reject"
artifacts <- c("ATXN1","PBRM1","ZNF814","MSH3","TTN","USH2A")
genes <- sapply(Final[,"Gene.refGene"],function(x) unlist(strsplit(x,","))[1])
toRemove <- genes%in%artifacts
Final[toRemove,"blackGenes"] <- "reject"
table(Final$blackGenes)
# . reject 
# 22994   4889

#ExAC, 1000Genomes, esp6500siv2_all, Kaviar_AF, HRC_AF allele frequency < 5% annotation

Final <- data.frame(SNP=".",Final,stringsAsFactors = FALSE)
toReject <- as.numeric(Final[,"ExAC_ALL"])>=0.05
toReject[is.na(toReject)] <- FALSE
table(toReject)
# FALSE  TRUE 
# 21888  5995  
Final[toReject,"SNP"] <- "reject"
toReject <- as.numeric(Final[,"X1000g2015aug_all"])>=0.05
toReject[is.na(toReject)] <- FALSE
table(toReject)
# FALSE   TRUE 
# 25356  2527 
Final[toReject,"SNP"] <- "reject"
toReject <- as.numeric(Final[,"esp6500siv2_all"])>=0.05
toReject[is.na(toReject)] <- FALSE
table(toReject)
# FALSE   TRUE 
# 26081  1802
Final[toReject,"SNP"] <- "reject"
toReject <- as.numeric(Final[,"Kaviar_AF"])>=0.05
toReject[is.na(toReject)] <- FALSE
table(toReject)
# FALSE   TRUE 
# 26529  1354 
Final[toReject,"SNP"] <- "reject"
toReject <- as.numeric(Final[,"HRC_AF"])>=0.05
toReject[is.na(toReject)] <- FALSE
table(toReject)
# FALSE   TRUE 
# 27154   729 
Final[toReject,"SNP"] <- "reject"
table(Final$SNP)
# . reject 
#  20729   7154



#snp138NonFlagged annotation

sum(toRemove <- Final[,"snp138NonFlagged"]!=".") #8935
Final[toRemove,"SNP"] <- "reject"
table(Final$SNP)
# . reject 
# 15986  11897

#50-mer alignability annotation

Final <- data.frame(fiftymer=".",Final,stringsAsFactors = FALSE)
bed <- Final[,c("Chr","Start","End")]
write.table(bed,file='Final.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
save(Final,file='Final.RData')