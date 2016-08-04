#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l ncpus=1
#PBS -l mem=2gb
#PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: QC Script
#Author: Matthew Lyon
#Status: Development
#Mode: BY_LANE

#load sample variables
. *.variables

if [[ $(gunzip -t "$R1Filename") -ne 0 ]] || [[ $(gunzip -t "$R2Filename") -ne 0 ]] ;then
	echo "FASTQ file(s) are corrupt, cannot proceed."
	exit -1
fi

#Get cluster density (K/mm2)
clusterDensity=$(cut -d, -f23 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')

#Get density passing filter (K/mm2)
clusterDensityPassingFilter=$(cut -d, -f24 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')

#Get clusters passing filter
pctPassingFilter=$(cut -d, -f44 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')

#Get % >Q30
pctGtQ30=$(cut -d, -f43 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')

#trim read 1 adapter from R1
/share/apps/cutadapt-distros/cutadapt-1.8/bin/cutadapt \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-e 0.2 \
-o "$SampleID"_R1_trimmed.fastq \
"$R1Filename"

#trim read 2 adapter from R2
/share/apps/cutadapt-distros/cutadapt-1.8/bin/cutadapt \
-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-e 0.2 \
-o "$SampleID"_R2_trimmed.fastq \
"$R2Filename"

#remove short reads
/usr/java/jdk1.7.0_51/bin/java -Xmx2g -jar /share/apps/trimmomatic/trimmomatic-0.32.jar PE \
"$SampleID"_R1_trimmed.fastq \
"$SampleID"_R2_trimmed.fastq \
"$SampleID"_R1_paired.fastq \
"$SampleID"_R1_unpaired.fastq \
"$SampleID"_R2_paired.fastq \
"$SampleID"_R2_unpaired.fastq \
MINLEN:36

#run fastqc
fastqc "$SampleID"_R1_paired.fastq
fastqc "$SampleID"_R2_paired.fastq

#locally align reads to reference genome
/share/apps/bwa-distros/bwa-0.7.10/bwa mem \
-M \
-R '@RG\tID:'"$RunID"'\tSM:'"$RunID"'\tPL:'"$Platform"'\tLB:'"$ExperimentName" \
/data/db/phix/mappers/bwa/genome.fa \
"$SampleID"_R1_paired.fastq "$SampleID"_R2_paired.fastq \
> "$RunID"_"$SampleID".sam

#Sort reads and convert to BAM
/usr/java/jdk1.7.0_51/bin/java -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-1.129/picard.jar SortSam \
INPUT="$RunID"_"$SampleID".sam \
OUTPUT="$RunID"_"$SampleID"_sorted.bam \
SORT_ORDER=coordinate \
COMPRESSION_LEVEL=0

#flagstat
samtools flagstat "$RunID"_"$SampleID"_sorted.bam

#Mark duplicate reads
/usr/java/jdk1.7.0_51/bin/java -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-1.129/picard.jar MarkDuplicates \
INPUT="$RunID"_"$SampleID"_sorted.bam \
OUTPUT="$RunID"_"$SampleID"_rmdup.bam \
CREATE_INDEX=true \
METRICS_FILE="$RunID"_"$SampleID"_dupMetrics.txt \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/usr/java/jdk1.7.0_51/bin/java -Xmx2g -jar /share/apps/GATK-distros/GATK_3.3.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$RunID"_"$SampleID"_rmdup.bam \
-o "$RunID"_"$SampleID".intervals \
-dt NONE

#Realign around indels
/usr/java/jdk1.7.0_51/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.3.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /data/db/phix/Illumina/1.1/genome.fa \
-targetIntervals "$RunID"_"$SampleID".intervals \
-I "$RunID"_"$SampleID"_rmdup.bam \
-o "$RunID"_"$SampleID"_realigned.bam \
-dt NONE

#calculate error rate
/usr/java/jdk1.7.0_51/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.3.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$RunID"_"$SampleID"_realigned.bam \
-knownSites phix.vcf \
-o recal_data.table

#clean up
rm "$SampleID"_R1_paired_fastqc.zip
rm "$SampleID"_R2_paired_fastqc.zip
rm "$SampleID"_R1_trimmed.fastq
rm "$SampleID"_R2_trimmed.fastq
rm "$SampleID"_R1_unpaired.fastq
rm "$SampleID"_R2_unpaired.fastq
rm "$RunID"_"$SampleID".sam
rm "$RunID"_"$SampleID"_sorted.bam
rm "$RunID"_"$SampleID"_rmdup.bam
rm "$RunID"_"$SampleID"_rmdup.bai
