fpfilterDir=$1
filevar=$5

cd $fpfilterDir

#for file in *.var
#do
	echo $filevar
	perl $fpfilter ${filevar} ${filevar%.var}.readcount --output-basename ${filevar%.var}.fpfilter > ${filevar%.var}.log &
#done  