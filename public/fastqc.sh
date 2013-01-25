#!/bin/sh
#PARAMETER TEST_PARAMETER1
#PARAMETER TEST_PARAMETER2

TEST_PARAMETER1=123
TEST_PARAMETER2=test_parameter

INPUT=/srv/GT/analysis/masaomi/sushi/20120921-hir_hay1_800_R1.fastq.gz
OUTPUT=/srv/GT/analysis/masaomi/sushi/SushiFabric/public/fastqc_out.zip
EXDIR=/srv/GT/analysis/masaomi/sushi/SushiFabric/public
fastqc_wrapper $INPUT $OUTPUT
unzip -o $OUTPUT -d $EXDIR
