#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: FASTQ generation and quality control for paired-end Illumina sequencing data. Not for use with other instruments or run configurations.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_RUN
#Usage: mkdir /data/archive/fastq/<seqId> && cd /data/archive/fastq/<seqId> && qsub -v sourceDir=/data/archive/miseq/<seqId> /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-<version>/1_IlluminaQC.sh
version="1.0.4"

### Preparation ###

#collect interop data
summary=$(/share/apps/interop-distros/InterOp-1.0.25-Linux-GNU-4.8.2/bin/summary --level=3 --csv=1 "$sourceDir")

#extract fields
yieldGb=$(echo "$summary" | grep ^Total | cut -d, -f2)
q30Pct=$(echo "$summary" | grep ^Total | cut -d, -f7)
avgDensity=$(echo "$summary" | grep -A999 "^Level" | grep ^[[:space:]]*[0-9] | awk -F',| ' '{print $1"\t"$4}' | sort | uniq | awk -F'\t' '{total += $2; count++} END {print total/count}')
avgPf=$(echo "$summary" | grep -A999 "^Level" |grep ^[[:space:]]*[0-9] | awk -F',| ' '{print $1"\t"$7}' | sort | uniq | awk -F'\t' '{total += $2; count++} END {print total/count}')
totalReads=$(echo "$summary" | grep -A999 "^Level" | grep ^[[:space:]]*[0-9] | awk -F',| ' '{print $1"\t"$19}' | sort | uniq | awk -F'\t' '{total += $2} END {print total}')

#print metrics (headers)
if [ ! -e /data/temp/Metrics.txt ]; then
    echo -e "Run\tTotalGb\tQ30\tAvgDensity\tAvgPF\tTotalMReads" > /data/temp/Metrics.txt
fi

#print metrics (values)
echo -e "$run\t$yieldGb\t$q30Pct\t$avgDensity\t$avgPf\t$totalReads" >> /data/temp/Metrics.txt

#update run log
#TODO

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