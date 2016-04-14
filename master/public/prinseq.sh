#!/bin/sh
#PARAMETER INPUT
#PARAMETER WORK_DIR
#PARAMETER OUTPUT_DIR


WORK_DIR=/srv/GT/analysis/masaomi/sushi/work
INPUT='public/20120905-I8_970_20081205_R1_sample.fastq.gz public/20120905-I8_970_20081205_R2_sample.fastq.gz'
OUTPUT_DIR=prinseq_out
RESULT_LINK=file:///$WORK_DIR/prinseq_out/ 

cd $WORK_DIR
echo cd $WORK_DIR
mkdir $OUTPUT_DIR
echo mkdir $OUTPUT_DIR
for input in $INPUT
do
  output=`basename $input | sed 's/\.fastq\.gz/_filtered/'`
  echo "gunzip -c $input | prinseq-lite.pl -fastq stdin -out_good $OUTPUT_DIR/$output -out_bad null -min_qual_score 20"
  gunzip -c $input | prinseq-lite.pl -fastq stdin -out_good $OUTPUT_DIR/$output -out_bad null -min_qual_score 20
done

echo RESULT_LINK=$RESULT_LINK
for input in $INPUT
do
  output=`basename $input | sed 's/\.fastq\.gz/_filtered.fastq/'`
  echo RESULT_PATH=$WORK_DIR/$OUTPUT_DIR/$output
done
