export PATH=$PATH:/storage/qnap_vol1/SHARED/NGSTOOLS/BCBIO/bcbio/anaconda/bin
DATADIR=/home/adabbo/gdc/SnakeFromBAM/SnakeExomePipeline/Data
FASTA=/home/adabbo/gdc/reference/Homo_sapiens_assembly19.fasta
JAVAPATH=/storage/qnap_vol1/bcbio/anaconda/pkgs/java-jdk-8.0.92-1/bin/java
BCBIOPATH=/storage/qnap_vol1/SHARED/NGSTOOLS/BCBIO/bcbio/anaconda/bin/bcbio-variation-recall


for file in $DATADIR/*-filtTNscope.snpEff.vcf.gz
do
	#echo $file
	name=${file##*/}
	name=${name%%-filtTNscope.snpEff.vcf.gz}
	echo $name
	echo $JAVAPATH -Xms512m -Xmx1g -jar $BCBIOPATH ensemble -c 80 --numpass 1 --names TNhaplotyper,TNscope,TNsnv --nofiltered $DATADIR/ensemble/${name}-ensemble.snpEff.vcf.gz $FASTA $DATADIR/${name}-TNhaplotyper.snpEff.vcf.gz $DATADIR/${name}-filtTNscope.snpEff.vcf.gz $DATADIR/${name}-TNsnv.snpEff.vcf.gz
	$JAVAPATH -Xms512m -Xmx1g -jar $BCBIOPATH ensemble -c 1 --numpass 1 --names TNhaplotyper,TNscope,TNsnv --nofiltered $DATADIR/ensemble/${name}-ensemble.snpEff.vcf.gz $FASTA $DATADIR/${name}-TNhaplotyper.snpEff.vcf.gz $DATADIR/${name}-filtTNscope.snpEff.vcf.gz $DATADIR/${name}-TNsnv.snpEff.vcf.gz
done
