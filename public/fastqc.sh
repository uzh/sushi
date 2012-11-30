#!/bin/sh

INPUT=/srv/GT/analysis/masaomi/sushi/20120921-hir_hay1_800_R1.fastq.gz
OUTPUT=/srv/GT/analysis/masaomi/sushi/SushiFabric/public/fastqc_out.zip
EXDIR=/srv/GT/analysis/masaomi/sushi/SushiFabric/public
fastqc_wrapper $INPUT $OUTPUT
unzip -o $OUTPUT -d $EXDIR
