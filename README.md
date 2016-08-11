# IlluminaQC

Quality control workflow for Illumina sequencing data


Launch with

mkdir /data/results/"$seqId" && cd /data/results/"$seqId" && qsub -v seqId="$seqId",sourceDir="/data/archive/miseq/$seqId" /data/diagnositcs/pipelines/IlluminaQC/IlluminaQC-"$version"/IlluminaQC.sh
