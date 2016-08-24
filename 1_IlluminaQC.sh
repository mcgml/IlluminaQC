#!/bin/bash -e
#PBS -l walltime=04:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Quality control for Illumina sequencing data. Not for use with other instruments.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/results/"$seqId" && cd /data/results/"$seqId" && qsub -v seqId="$seqId",sourceDir="/data/archive/miseq/$seqId" /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-"$version"/1_IlluminaQC.sh
version="dev"

#TODO highest unmatched index seq -- check for sample contamination
#TODO update interop

#log with Trello
/share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
/data/diagnostics/scripts/TrelloAPI.js \
"$seqId" "Starting QC"

### Set up ###
passedSeqId="$seqId"
passedSourceDir="$sourceDir"
. /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-"$version"/variables

#get SAV metrics & check %Q30 passed QC
/share/apps/interop-distros/interop/build/bin/bin/imaging_table "$sourceDir" | grep -vP "#|Lane|^$" | \
awk -F, '{ density[$1]+=$6; pf[$1]+=$10; q30[$1]+=$15; n++ } END { print "Lane\tClusterDensity\tPctPassingFilter\tPctGtQ30"; for(i in density) print i"\t"density[i]/n"\t"pf[i]/n"\t"q30[i]/n; }' | \
tee "$seqId"_sav.txt | awk '{ if (NR > 1 && $4 < 80) { print "Run generated insufficient Q30 data"; exit -1 } }'

#convert bcls to FASTQ
/usr/local/bin/bcl2fastq -l WARNING -R "$sourceDir" -o .

#link SampleSheet & runParameters.xml
ln -s "$passedSourceDir"/SampleSheet.csv
ln -s "$passedSourceDir"/?unParameters.xml

#Make variable files
/share/apps/jre-distros/jre1.8.0_101/bin/java -jar /data/diagnostics/apps/MakeVariableFiles-2.0.0/MakeVariableFiles-2.0.0.jar \
SampleSheet.csv \
?unParameters.xml

#move fastq & variable files into project folders
for variableFile in $(ls *.variables); do
	
	#load variables into scope
	. "$variableFile"

	#make project folders
	mkdir -p "$panel"
	mkdir "$panel"/"$sampleId"
	mv "$variableFile" "$panel"/"$sampleId"

	#move FASTQs into sample folder
    for fastqPair in $(ls "$sampleId"_*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
        
        #parse fastq filenames
        laneId=$(echo "$fastqPair" | cut -d_ -f3)
        read1Fastq=$(ls "$fastqPair"_R1*fastq.gz)
        read2Fastq=$(ls "$fastqPair"_R2*fastq.gz)

        mv "$read1Fastq" "$panel"/"$sampleId"
		mv "$read2Fastq" "$panel"/"$sampleId"

    done

	#cp pipeline scripts
	cp /data/diagnostics/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/*sh "$panel"/"$sampleId"

    #record job workdirs
    echo $PWD/"$panel"/"$sampleId" >> workdirs.list

done

#move unindexed data into sample folder 
mkdir Undetermined
mv  Undetermined_*.fastq.gz Undetermined

### QC ###
cd Undetermined

for fastqPair in $(ls Undetermined_S0_*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2*fastq.gz)
    
    #trim adapters and remove short reads
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    -o $(echo "$read1Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') \
    -p $(echo "$read2Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') \
    --minimum-length 30 \
    "$read1Fastq" \
    "$read2Fastq"
    
    #Align reads to reference genome, retain only proper pairs, sort by coordinate and convert to BAM
    /share/apps/bwa-distros/bwa-0.7.15/bwa mem \
    -M \
    -R '@RG\tID:'"$passedSeqId"_"$laneId"'_PhiX\tSM:PhiX\tPL:ILLUMINA\tLB:'"$passedSeqId"'_PhiX\tPU:'"$passedSeqId"_"$laneId" \
    -t 12 \
    /data/db/phix/mappers/bwa/genome.fa \
    $(echo "$read1Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') $(echo "$read2Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') | \
    /share/apps/samtools-distros/samtools-1.3.1/samtools view -h -f2 | \
    /share/apps/samtools-distros/samtools-1.3.1/samtools sort -l0 -o "$passedSeqId"_PhiX_"$laneId"_sorted.bam

done

#merge mulitple lanes
/share/apps/samtools-distros/samtools-1.3.1/samtools merge \
-@ 12 \
-u "$passedSeqId"_PhiX_all_sorted.bam \
"$passedSeqId"_PhiX_*_sorted.bam

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$passedSeqId"_PhiX_all_sorted.bam \
OUTPUT="$passedSeqId"_PhiX_rmdup.bam \
METRICS_FILE="$passedSeqId"_PhiX_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#check PhiX has been loaded
if [ $(/share/apps/samtools-distros/samtools-1.3.1/samtools view -c -F 0x400 "$passedSeqId"_PhiX_rmdup.bam) < 50000 ]; then
    echo "Insufficient PhiX reads for error modelling"
    exit -1
fi

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$passedSeqId"_PhiX_rmdup.bam \
-o "$passedSeqId"_PhiX_realign.intervals \
-nt 12 \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /data/db/phix/Illumina/1.1/genome.fa \
-targetIntervals "$passedSeqId"_PhiX_realign.intervals \
-I "$passedSeqId"_PhiX_rmdup.bam \
-o "$passedSeqId"_PhiX_realigned.bam \
-compress 0 \
-dt NONE

#Get quality vs error rate data
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$passedSeqId"_PhiX_realigned.bam \
-knownSites /data/db/phix/phix.vcf \
-o "$passedSeqId"_PhiX_BaseRecalibrator.txt \
-nct 12 \
-dt NONE

#Calculate pearson correlation with 95% confidence
/share/apps/R-distros/R-3.3.1/bin/Rscript \
/data/diagnositcs/pipelines/IlluminaQC/IlluminaQC-"$version"/bqsrAnalysis.R \
-r "$passedSeqId"_PhiX_BaseRecalibrator.txt | \
tee "$passedSeqId"_bqsr_pearson.txt | \
awk '{if ($1 < 0.95) { print "Poor correlation between reported and emperical Q scores"; exit -1; }}'

### Clean up ###
rm -r tmp
rm *.fastq "$passedSeqId"_PhiX_*_sorted.bam "$passedSeqId"_PhiX_*_sorted.bai
rm "$passedSeqId"_PhiX_MarkDuplicatesMetrics.txt
rm "$passedSeqId"_PhiX_realign.intervals

#log with Trello
/share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
/data/diagnostics/scripts/TrelloAPI.js \
"$seqId" "Passed QC. Starting analysis"

#launch analyses
for i in $(sort workdirs.list | uniq); do
    cd "$i"
    qsub 1_*.sh
done