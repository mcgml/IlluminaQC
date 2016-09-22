#!/bin/bash
set -euo pipefail
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Quality control for paired-end Illumina sequencing data. Not for use with other instruments or run configurations.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/archive/ubam/"$seqId" && cd /data/archive/ubam/"$seqId" && qsub -v sourceDir=/data/archive/miseq/"$seqId" /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-"$version"/1_IlluminaQC.sh
version="dev"

phoneTrello() {
    /share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2" #seqId & message
}

countQCFlagFails() {
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

### Preparation ###

#log with Trello
phoneTrello $(basename "$sourceDir") "Starting QC"

#make output folder and change dir
mkdir Data

#convert BCLs to FASTQ
/usr/local/bin/bcl2fastq -l WARNING -R "$sourceDir" -o .
rm Undetermined_S0_*.fastq.gz

#copy files to keep to long-term storage
cp "$sourceDir"/SampleSheet.csv .
cp "$sourceDir"/?unParameters.xml RunParameters.xml
cp "$sourceDir"/RunInfo.xml .
cp -R "$sourceDir"/InterOp .

#Make variable files
/share/apps/jre-distros/jre1.8.0_101/bin/java -jar /data/diagnostics/apps/MakeVariableFiles/MakeVariableFiles-2.1.0.jar \
SampleSheet.csv \
RunParameters.xml

#move fastq & variable files into project folders
for variableFile in $(ls *.variables); do

    #reset variables
    unset sampleId fastqPair laneId read1Fastq read2Fastq seqId worklistId pipelineVersion pipelineName panel
    failed=false

    #load variables into local scope
	. "$variableFile"

    #make sample folder
    mkdir Data/"$sampleId"
    mv "$variableFile" Data/"$sampleId"

	#convert FASTQ to uBAM with RGIDs
    for fastqPair in $(ls "$sampleId"_*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
        
        #parse fastq filenames
        laneId=$(echo "$fastqPair" | cut -d_ -f3)
        read1Fastq=$(ls "$fastqPair"_R1*fastq.gz)
        read2Fastq=$(ls "$fastqPair"_R2*fastq.gz)
        
        #convert fastq to ubam
        /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar FastqToSam \
        F1="$read1Fastq" \
        F2="$read2Fastq" \
        O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
        READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
        SAMPLE_NAME="$sampleId" \
        LIBRARY_NAME="$worklistId"_"$panel"_"$sampleId" \
        PLATFORM_UNIT="$seqId"_"$laneId" \
        PLATFORM="ILLUMINA" \
        SORT_ORDER=queryname \
        MAX_RECORDS_IN_RAM=5000000 \
        TMP_DIR=/state/partition1/tmpdir

        #fastqc
        /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract -o Data/"$sampleId" "$read1Fastq"
        /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract -o Data/"$sampleId" "$read2Fastq"

        #check FASTQ output
        if [ $(countQCFlagFails Data/"$sampleId"/"$(echo $read1Fastq | sed 's/\.fastq\.gz/_fastqc/g')/summary.txt") -gt 0 ] || [ $(countQCFlagFails Data/"$sampleId"/"$(echo $read2Fastq | sed 's/\.fastq\.gz/_fastqc/g')/summary.txt") -gt 0 ]; then
            failed=true
        fi

    done
    
    #merge lane ubams
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MergeSamFiles \
    $(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
    SORT_ORDER=queryname \
    USE_THREADING=true \
    MAX_RECORDS_IN_RAM=5000000 \
    TMP_DIR=/state/partition1/tmpdir \
    CREATE_MD5_FILE=true \
    O=Data/"$sampleId"/"$seqId"_"$sampleId"_unaligned.bam

    #clean up
    rm "$sampleId"_*.fastq.gz "$seqId"_"$sampleId"_*_unaligned.bam Data/*/*_fastqc.html Data/*/*_fastqc.zip

    #skip failed samples
    if [ "$failed" = true ] ; then
        phoneTrello "$seqId" "$sampleId has failed QC"
        continue;
    fi
    
    #analysis
    if [[ ! -z ${pipelineVersion-} && ! -z ${pipelineName-} && ! -z ${panel-} ]]; then

        #make project folders
        mkdir -p /data/results/"$seqId"
        mkdir -p /data/results/"$seqId"/"$panel"

        #make sample folder & link uBam
        mkdir /data/results/"$seqId"/"$panel"/"$sampleId"
        ln -s $PWD/Data/"$sampleId"/"$seqId"_"$sampleId"_unaligned.bam /data/results/"$seqId"/"$panel"/"$sampleId"
        ln -s $PWD/Data/"$variableFile" /data/results/"$seqId"/"$panel"/"$sampleId"

        #copy scripts
        cp /data/diagnostics/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/*sh /data/results/"$seqId"/"$panel"/"$sampleId"

        #create worklist
        echo /data/results/"$seqId"/"$panel"/"$sampleId" >> workdirs.list

    fi

done

### QC ###

#get SAV metrics & check %Q30 passed QC
/share/apps/interop-distros/interop-1.0.11/build/bin/usr/local/bin/imaging_table "$sourceDir" | grep -vP "#|Lane|^$" | \
awk -F, '{ density[$1]+=$6; pf[$1]+=$10; q30[$1]+=$15; n[$1]++ } END { print "Lane\tClusterDensity\tPctPassingFilter\tPctGtQ30"; for(i in density) print i"\t"density[i]/n[i]"\t"pf[i]/n[i]"\t"q30[i]/n[i]; }' | \
tee $(basename "$sourceDir")_sav.txt | awk '{ if (NR > 1 && $4 < 80) { print "Run generated insufficient Q30 data"; phoneTrello $(basename "$sourceDir") "Failed QC. Insufficient data"; exit -1 } }'

### Launch analysis ###
phoneTrello $(basename "$sourceDir") "Passed QC. Launching analysis..."
for i in $(sort workdirs.list | uniq); do
    bash -c cd "$i" && qsub 1_*.sh
done