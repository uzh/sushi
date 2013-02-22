#!/bin/sh
#PARAMETER INPUT
#PARAMETER BOWTIE_INDEX
#PARAMETER OUTPUT_DIR
#PARAMETER THREAD 

INPUT='/srv/GT/analysis/masaomi/sushi/work/SushiFabric/public/20120905-I8_970_20081205_R1_sample.fastq.gz /srv/GT/analysis/masaomi/sushi/work/SushiFabric/public/20120905-I8_970_20081205_R2_sample.fastq.gz'
OUTPUT_DIR=/srv/GT/analysis/masaomi/sushi/work/tophat_out
BOWTIE_INDEX=/srv/GT/analysis/masaomi/sushi/work/bowtie_index/flowering_genes
THREAD=1
RESULT_LINK=file:///srv/GT/analysis/masaomi/sushi/work/tophat_out

echo tophat -p $THREAD -o $OUTPUT_DIR $BOWTIE_INDEX $INPUT
tophat -p $THREAD -o $OUTPUT_DIR $BOWTIE_INDEX $INPUT
echo RESULT_LINK=$RESULT_LINK

