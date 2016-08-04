#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l ncpus=12
#PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
#cd $PBS_O_WORKDIR

#Description: Quality control for Illumina sequencing data
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_CHIP
version="dev"
seqId="160728_M02641_0123_000000000-ARBU9"

#TODO support multiple lanes
#TODO support HiSeq
#TODO launch sample analysis on success
#TODO move all archived NGS runs into a single folder
#TODO replace SAV metrics with Stats folder from bcl2fastq
#TODO mkdir run folder, cd and qsub this script

#convert Bcls to FASTQ
bcl2fastq -l WARNING -R /data/archive/miseq/"$seqId" -o .

#Get metrics from SAV
clusterDensity=$(cut -d, -f23 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')
clusterDensityPassingFilter=$(cut -d, -f24 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')
pctPassingFilter=$(cut -d, -f44 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')
pctGtQ30=$(cut -d, -f43 /data/archive/metrics/"$seqId"_SAV.txt | awk '{ if (NR>1) total += $1 } END { print total/(NR-1)}')

#move fastq files into folders
for i in $(ls *fastq.gz); do
	 sampleId=$(echo "$i" | cut -d_ -f1);

	 mkdir -p "$sampleId"
	 mv "$i" "$sampleId"
done

#start QC
cd Undetermined

#trim adapters and remove short reads
/share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-o "$seqId"_Undetermined_R1_trimmed.fastq \
-p "$seqId"_Undetermined_R2_trimmed.fastq \
--minimum-length 30 \
Undetermined_S0_L001_R1_001.fastq.gz \
Undetermined_S0_L001_R2_001.fastq.gz

#Align reads to reference genome, sort by coordinate and convert to BAM
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-R '@RG\tID:'"$seqId"'_PhiX\tSM:PhiX\tPL:ILLUMINA\tLB:'"$seqId"'_PhiX' \
/data/db/phix/mappers/bwa/genome.fa \
"$seqId"_Undetermined_R1_trimmed.fastq "$seqId"_Undetermined_R2_trimmed.fastq | \
/share/apps/samtools-distros/samtools-1.3.1/samtools view -h -f2 | \
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -l0 -o "$seqId"_Undetermined_sorted.bam

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$seqId"_Undetermined_sorted.bam \
OUTPUT="$seqId"_Undetermined_rmdup.bam \
METRICS_FILE="$seqId"_Undetermined_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$seqId"_Undetermined_rmdup.bam \
-o "$seqId"_Undetermined_realign.intervals \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /data/db/phix/Illumina/1.1/genome.fa \
-targetIntervals "$seqId"_Undetermined_realign.intervals \
-I "$seqId"_Undetermined_rmdup.bam \
-o "$seqId"_Undetermined_realigned.bam \
-compress 0 \
-dt NONE

#Get quality vs error rate data
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$seqId"_Undetermined_realigned.bam \
-knownSites /data/db/phix/phix.vcf \
-o "$seqId"_Undetermined_BaseRecalibrator.txt

#TODO calculate pearson correlation


#Calculate error-rate by cycle
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T ErrorRatePerCycle \
-R /data/db/phix/Illumina/1.1/genome.fa \
-I "$seqId"_Undetermined_realigned.bam \
-knownSites /data/db/phix/phix.vcf \
-o "$seqId"_Undetermined_ErrorRatePerCycle.txt

#clean up
rm -r tmp
rm "$seqId"_Undetermined_sorted.bam
rm "$seqId"_Undetermined_rmdup.bam "$seqId"_Undetermined_rmdup.bai