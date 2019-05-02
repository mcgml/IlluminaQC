# temp sample sheet location for testing
CURRENTDIRECTORY=$(pwd)
INPUTFOLDER=${CURRENTDIRECTORY%/*}/

# clear sample name list
>"$INPUTFOLDER"/samplenames.list

# obtain list of samples from sample sheet
for line in $(sed "1,/Sample_ID/d" "$INPUTFOLDER""SampleSheet.csv" | tr -d " ")
    do
        # obtain sample name and patient name		
        samplename=$(printf "$line" | cut -d, -f1 | sed 's/[^a-zA-Z0-9]+/-/g')
		echo $samplename
		echo $samplename >> "$INPUTFOLDER"/samplenames.list
		#mkdir $INPUTFOLDER/$samplename
		#touch $INPUTFOLDER/$samplename".variables"
		#touch $INPUTFOLDER/$samplename"_1_.fastq"
		#touch $INPUTFOLDER/$samplename"_2_.fastq"
		#touch $INPUTFOLDER/$samplename"_3_.fastq"
		#touch $INPUTFOLDER/$samplename"_4_.fastq"

    done
	
# compare samples list to generated variables files
if [ $(find $INPUTFOLDER -maxdepth 1 -mindepth 1 -name *".variables" | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort $INPUTFOLDER/samplenames.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ] 
	then
		echo "expected number of variables files created"
	else
		echo "number of variables files not as expected"
		#break
	fi
	
# compare samples list to fastqs
# create file of generated fastqs


#compare to samples list