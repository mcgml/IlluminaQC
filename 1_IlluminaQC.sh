#!/bin/bash
set -euo pipefail
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: FASTQ generation and quality control for paired-end Illumina sequencing data. Not for use with other instruments or run configurations.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/archive/fastq/"$seqId" && cd /data/archive/fastq/"$seqId" && qsub -v sourceDir=/data/archive/miseq/"$seqId" /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-"$version"/1_IlluminaQC.sh
version="dev"

phoneTrello() {
    /share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2" #seqId & message
}

### Preparation ###

#log with Trello
phoneTrello $(basename "$sourceDir") "Starting QC"

#convert BCLs to FASTQ
/usr/local/bin/bcl2fastq -l WARNING -R "$sourceDir" -o .
rm Undetermined_S0_*.fastq.gz

#copy files to keep to long-term storage
mkdir Data
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
    unset sampleId seqId worklistId pipelineVersion pipelineName panel

    #load variables into local scope
	. "$variableFile"

    #make sample folder
    mkdir Data/"$sampleId"
    mv "$variableFile" Data/"$sampleId"
    mv "$sampleId"_S*.fastq.gz Data/"$sampleId"
    
    #launch analysis
    if [[ ! -z ${pipelineVersion-} && ! -z ${pipelineName-} && ! -z ${panel-} ]]; then

        #make project folders
        mkdir -p /data/results/"$seqId"
        mkdir -p /data/results/"$seqId"/"$panel"

        #make sample folder
        mkdir /data/results/"$seqId"/"$panel"/"$sampleId"

        #soft link files
        ln -s $PWD/Data/"$variableFile" /data/results/"$seqId"/"$panel"/"$sampleId"
        for i in $(ls Data/"$sampleId"/"$sampleId"_S*.fastq.gz); do
            ln -s $PWD/Data/"$sampleId"/"$i" /data/results/"$seqId"/"$panel"/"$sampleId"
        done

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
tee $(basename "$sourceDir")_sav.txt | awk '{ if (NR > 1 && $4 < 80) { print "Run generated insufficient Q30 data"; phoneTrello $(basename "$sourceDir") "Failed QC. Insufficient data"; print "Low % Q30" > /data/results/"$seqId"/QC_FAIL } }'
ln -s $(basename "$sourceDir")_sav.txt /data/results/"$seqId"

### Launch analysis ###
phoneTrello $(basename "$sourceDir") "Launching analysis..."
for i in $(sort workdirs.list | uniq); do
    bash -c cd "$i" && qsub 1_*.sh
done