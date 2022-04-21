fpfilterDir=$1
fastaPath=$2
workPath=$3 
bamreadcount=$4
fileloc=$5
filevar=$6
fpfilter=$7

cd $fpfilterDir

#for file in *.loc
#do
echo $fileloc
name=${fileloc%.loc}
echo $name
echo "$name bam-readcount"
bam=`ls $workPath/${name}T_recal.bam`
echo $bam
echo '\n'
$bamreadcount -f $fastaPath -w 1 -l $fileloc $bam  > ${name}.readcount
#done 
echo "done bamreadcount"

#for file in *.var
#do
echo $filevar
perl $fpfilter ${filevar} ${filevar%.var}.readcount --output-basename ${filevar%.var}.fpfilter > ${filevar%.var}.log

echo "done fpfilter"
#done  