#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Quality control for Illumina sequencing data. Not for use with other instruments.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/results/"$seqId" && cd /data/results/"$seqId" && qsub -v seqId="160725_M02641_0122_000000000-AU4LF",sourceDir="/data/archive/miseq/160725_M02641_0122_000000000-AU4LF" /data/diagnositcs/pipelines/QC/QC-1.0.0/QC.sh
version="dev"

#TODO get metrics from bcl2fastq output: clusterDensity, clusterDensityPassingFilter, pctPassingFilter, pctGtQ30, highest unmatched index seq
#TODO print metrics
#TODO qsub if run passes QC then launch analysis
#TODO report to trello

### Set up ###
passedSeqId="$seqId"
passedSourceDir="$sourceDir"

#convert bcls to FASTQ
bcl2fastq -l WARNING -R "$sourceDir" -o .

#link SampleSheet & runParameters.xml
ln -s "$passedSourceDir"/SampleSheet.csv
ln -s "$passedSourceDir"/?unParameters.xml

#Make variable files
/share/apps/jre-distros/jre1.8.0_71/bin/java -jar /share/apps/MakeVariableFiles-2.0.0/MakeVariableFiles.jar \
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
	for fastqPair in fastqPairs; do
		read1Fastq=$(echo "$fastqPair" | cut -d, -f1)
		read2Fastq=$(echo "$fastqPair" | cut -d, -f2)

		mv "$read1Fastq" "$panel"/"$sampleId"
		mv "$read2Fastq" "$panel"/"$sampleId"
	done

	#cp pipeline scripts
	cp /data/diagnostics/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/*sh "$panel"/"$sampleId"

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
    -R '@RG\tID:'"$passedSeqId"_"$laneId"_"PhiX"'\tSM:'"PhiX"'\tPL:ILLUMINA\tLB:'"$passedSeqId_PhiX" \
    -t 6 \
    /data/db/phix/mappers/bwa/genome.fa \
    $(echo "$read1Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') $(echo "$read2Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') | \
    /share/apps/samtools-distros/samtools-1.3.1/samtools view -h -f2 | \
    /share/apps/samtools-distros/samtools-1.3.1/samtools sort -@4 -l0 -o "$passedSeqId"_PhiX_"$laneId"_sorted.bam

done

#merge mulitple lanes
if [ $(ls "$passedSeqId"_PhiX_*_sorted.bam | wc -l | sed 's/^[[:space:]]*//g') -gt 1 ]; then
    /share/apps/samtools-distros/samtools-1.3.1/samtools merge -u "$passedSeqId"_PhiX_all_sorted.bam "$passedSeqId"_PhiX_*_sorted.bam
else
    mv "$passedSeqId"_PhiX_*_sorted.bam "$passedSeqId"_PhiX_all_sorted.bam
fi

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$passedSeqId"_PhiX_all_sorted.bam \
OUTPUT="$passedSeqId"_PhiX_rmdup.bam \
METRICS_FILE="$seqId"_PhiX_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$passedSeqId"_PhiX_rmdup.bam \
-o "$passedSeqId"_PhiX_realign.intervals \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /data/db/phix/Illumina/1.1/genome.fa \
-targetIntervals "$passedSeqId"_PhiX_realign.intervals \
-I "$passedSeqId"_PhiX_rmdup.bam \
-o "$passedSeqId"_PhiX_realigned.bam \
-compress 0 \
-dt NONE

#Get quality vs error rate data
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$passedSeqId"_PhiX_realigned.bam \
-knownSites /data/db/phix/phix.vcf \
-o "$passedSeqId"_PhiX_BaseRecalibrator.txt \
-dt NONE

#Calculate pearson correlation
pearson=$(/share/apps/R-distros/R-3.3.1/bin/Rscript bqsrAnalysis.R -r "$passedSeqId"_PhiX_BaseRecalibrator.txt)

### Clean up ###

#log run complete
#/share/apps/node-distros/node-v0.12.7-linux-x64/bin/node \
#/data/diagnostics/scripts/TrelloAPI.js \
#"$seqId" "$worksheetId"