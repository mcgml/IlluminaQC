#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: FASTQ generation and quality control for paired-end Illumina sequencing data. Not for use with other instruments or run configurations.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/archive/fastq/"$seqId" && cd /data/archive/fastq/"$seqId" && qsub -v sourceDir=/data/archive/miseq/"$seqId" /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-"$version"/1_IlluminaQC.sh
version="1.0.1"

### Preparation ###

#get SAV metrics & check %Q30 passed QC
/share/apps/interop-distros/interop-1.0.11/build/bin/usr/local/bin/imaging_table "$sourceDir" | grep -vP "#|Lane|^$" | \
awk -F, '{ density[$1]+=$6; pf[$1]+=$10; q30[$1]+=$15; n[$1]++ } END { print "Lane\tClusterDensity\tPctPassingFilter\tPctGtQ30"; for(i in density) print i"\t"density[i]/n[i]"\t"pf[i]/n[i]"\t"q30[i]/n[i]; }' | \
tee $(basename "$sourceDir")_sav.txt | awk '{ if (NR > 1 && $4 < 80) { print "Lane $1 generated low Q30% ($4%)" } }'

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
        ln -s $PWD/Data/"$sampleId"/"$variableFile" /data/results/"$seqId"/"$panel"/"$sampleId"
        for i in $(ls Data/"$sampleId"/"$sampleId"_S*.fastq.gz); do
            ln -s $PWD/"$i" /data/results/"$seqId"/"$panel"/"$sampleId"
        done

        #copy scripts
        cp /data/diagnostics/pipelines/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/*sh /data/results/"$seqId"/"$panel"/"$sampleId"

        #queue pipeline
        bash -c "cd /data/results/$seqId/$panel/$sampleId && qsub 1_*.sh"

    fi

done