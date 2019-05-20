# temp sample sheet location for testing
CURRENTDIRECTORY=$(pwd)
INPUTFOLDER=${CURRENTDIRECTORY%/*}/

# clear sample name list
>samplenames.list

# obtain list of samples from sample sheet
for line in $(sed "1,/Sample_ID/d" "$INPUTFOLDER""SampleSheet.csv" | tr -d " ")
do
    # obtain sample name and patient name		
    samplename=$(printf "$line" | cut -d, -f1 | sed 's/[^a-zA-Z0-9]+/-/g')
	echo $samplename >> samplenames.list
    done

	
# compare samples list to generated variables files by count 
if ! [ $(find . -maxdepth 1 -mindepth 1 -name '*.variables' | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort samplenames.list | wc -l | sed 's/^[[:space:]]*//g') ] 
then
	echo "Number of variables files not as expected"
	exit 1
else
	echo "Expected number of variables files created"
fi
	
# compare samples list to fastqs generated
for sn in $( cat samplenames.list )
do
	if ! [[ $(find . -name $sn'*fastq*') ]]
	then
		echo "Missing fastqs for sample " $sn
		exit 1
	else
		echo "Fastqs generated for sample " $sn
	fi
done 


# compare directory number in data/results/ to list of samples from sample sheet
if ! [[ $(ls ../data/results/ | wc -l) -eq $(sort samplenames.list | wc -l | sed 's/^[[:space:]]*//g') ]]
then
	echo "Number of sample directories created not as expected"
	exit 1
else
	echo "Expected number of sample directories created"
fi