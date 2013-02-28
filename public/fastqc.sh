#!/bin/sh
#PARAMETER INPUT
#PARAMETER OUTPUT_DIR
#PARAMETER THREAD


INPUT=/srv/GT/analysis/masaomi/sushi/work/SushiFabric/public/20120905-I8_970_20081205_R1_sample.fastq.gz
OUTPUT_DIR=/srv/GT/analysis/masaomi/sushi/work/fastqc_out
RESULT_LINK=file:///srv/GT/analysis/masaomi/sushi/work/fastqc_out/ 
THREAD=4

mkdir -p $OUTPUT_DIR
fastqc -t $THREAD -o $OUTPUT_DIR $INPUT

echo RESULT_LINK=$RESULT_LINK
