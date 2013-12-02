#!/bin/bash

PICKRELL=/srv/gs1/projects/snyder/ccenik/PICKRELL_RNASEQ
cd $PICKRELL
for dir in $(find . -maxdepth 1 -mindepth 1 -name 'GM*' -type d -printf '%f\n')
do
cd $PICKRELL/$dir/
pwd
zcat *  > Merged_Reads.fastq
qsub  ~/SCRIPTS/BASH_SCRIPTS/alignment_strategy_pickrell_qsub.sh
done