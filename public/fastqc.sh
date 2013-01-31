#!/bin/sh
#PARAMETER INPUT
#PARAMETER OUTPUT_DIR

INPUT=/srv/GT/analysis/masaomi/sushi/work/SushiFabric/20120905-I8_970_20081205_R1_sample.fastq.gz
OUTPUT_DIR=/srv/GT/analysis/masaomi/sushi/work/fastqc_out
mkdir -p $OUTPUT_DIR
fastqc -t 4 -o $OUTPUT_DIR $INPUT

