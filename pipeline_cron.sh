#!/bin/sh
set -euo pipefail

#Description: shell script to launch bioinformatics analysis pipelines. Run as root cron job without .sh extension
#Author:Matt Lyon
#Date: 28/02/17
version="1.0.5"

function processJobs {
    echo "checking for jobs in $1 ..."

    for path in $(find "$1" -maxdepth 2 -mindepth 2 -type f -name "RTAComplete.txt" -exec dirname '{}' \;); do

        #extract run info from path
        instrumentType=$(basename $(dirname "$path"))
        run=$(basename "$path")

        #log
	    echo "path: $path"
        echo "run: $run"
        echo "instrumentType: $instrumentType"

        #move run to archive
        mv "$path" /data/archive/"$instrumentType"

        #change access permissions
        chown -R transfer /data/archive/"$instrumentType"/"$run"
        chgrp -R transfer /data/archive/"$instrumentType"/"$run"
        chmod -R 755 /data/archive/"$instrumentType"/"$run"

        #launch IlluminaQC for demultiplexing and QC
        ssh transfer@10.59.210.245 "mkdir /data/archive/fastq/$run && cd /data/archive/fastq/$run && qsub -v sourceDir=/data/archive/$instrumentType/$run /data/diagnostics/pipelines/IlluminaQC/IlluminaQC-$version/1_IlluminaQC.sh"

    done

}

processJobs "/data/raw/miseq"
processJobs "/data/raw/hiseq"