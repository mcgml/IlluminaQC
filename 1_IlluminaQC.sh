#!/bin/bash
set -euo pipefail
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Quality control for Illumina sequencing data. Not for use with other instruments.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/archive/ubam/"$seqId" && cd /data/archive/ubam/"$seqId" && qsub -v seqId="$seqId",sourceDir="/data/archive/miseq/$seqId" /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-"$version"/1_IlluminaQC.sh
version="dev"

#TODO highest unmatched index seq -- check for sample contamination
#TODO update interop command

phoneTrello() {
    /share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2" #seqId & message
}

### Preparation ###
passedSeqId="$seqId"
passedSourceDir="$sourceDir"

#log with Trello
phoneTrello "$passedSeqId" "Starting QC"

#convert BCLs to FASTQ
/usr/local/bin/bcl2fastq -l WARNING -R "$passedSourceDir" -o .

#copy files to keep to long-term storage
cp "$passedSourceDir"/SampleSheet.csv .
cp "$passedSourceDir"/?unParameters.xml RunParameters.xml
cp "$passedSourceDir"/RunInfo.xml .
cp -R "$passedSourceDir"/InterOp .

#Make variable files
/share/apps/jre-distros/jre1.8.0_101/bin/java -jar /data/diagnostics/apps/MakeVariableFiles/MakeVariableFiles-2.1.0.jar \
SampleSheet.csv \
RunParameters.xml

#make results folder
mkdir /data/results/"$passedSeqId"

#move fastq & variable files into project folders
for variableFile in $(ls *.variables); do

{
    #load variables into scope
	. "$variableFile"

    #make sample folder
    mkdir "$sampleId"
    mv "$variableFile" "$sampleId"

	#convert FASTQ to uBAM with RGIDs
    for fastqPair in $(ls "$sampleId"_*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
        
        #parse fastq filenames
        laneId=$(echo "$fastqPair" | cut -d_ -f3)
        read1Fastq=$(ls "$fastqPair"_R1*fastq.gz)
        read2Fastq=$(ls "$fastqPair"_R2*fastq.gz)
        
        #convert fastq to ubam
        /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar FastqToSam \
        F1="$read1Fastq" \
        F2="$read2Fastq" \
        O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
        READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
        SAMPLE_NAME="$sampleId" \
        LIBRARY_NAME="$worklistId"_"$panel"_"$sampleId" \
        PLATFORM_UNIT="$seqId"_"$laneId" \
        PLATFORM="ILLUMINA"

    done

    #merge lane ubams
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MergeSamFiles \
    I="$seqId"_"$sampleId"_*_unaligned.bam \
    O="$sampleId"/"$seqId"_"$sampleId"_unaligned.bam \
    SORT_ORDER=queryname \
    USE_THREADING=true

    #make project folders
    if [[ -z "$panel" || -z "$pipelineVerison" || -z "$pipelineName" ]]; then
        continue;
    fi

    mkdir -p /data/results/"$seqId"/"$panel"
    mkdir /data/results/"$seqId"/"$panel"/"$sampleId"
    ln -s "$sampleId"/"$seqId"_"$sampleId"_unaligned.bam /data/results/"$seqId"/"$panel"/"$sampleId"
    ln -s "$variableFile" /data/results/"$seqId"/"$panel"/"$sampleId"
    cp /data/diagnostics/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVerison"/*sh /data/results/"$seqId"/"$panel"/"$sampleId"
    echo /data/results/"$seqId"/"$panel"/"$sampleId" >> ../workdirs.list

}

done

### Clean up ###
rm *.fastq.gz *unaligned.bam

### QC ###

#get SAV metrics & check %Q30 passed QC
/share/apps/interop-distros/interop-1.0.11/build/bin/usr/local/bin/imaging_table "$passedSourceDir" | grep -vP "#|Lane|^$" | \
awk -F, '{ density[$1]+=$6; pf[$1]+=$10; q30[$1]+=$16; n[$1]++ } END { print "Lane\tClusterDensity\tPctPassingFilter\tPctGtQ30"; for(i in density) print i"\t"density[i]/n[i]"\t"pf[i]/n[i]"\t"q30[i]/n[i]; }' | \
tee "$passedSeqId"_sav.txt | awk '{ if (NR > 1 && $4 < 80) { print "Run generated insufficient Q30 data"; phoneTrello "$passedSeqId" "Failed QC. Insufficient data"; exit -1 } }'

### Launch analysis ###
phoneTrello "$passedSeqId" "Passed QC. Launching analysis..."
for i in $(sort workdirs.list | uniq); do
    bash -c cd "$i" && qsub 1_*.sh
done